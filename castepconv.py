#!/usr/bin/env python

"""
CASTEP convergence automation tool
by Simone Sturniolo

Copyright 2013-2018 Science and Technology Facilities Council
This software is distributed under the terms of the GNU General Public
License (GNU GPL)
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# Main script

import os
import sys
import time
import glob
import json
import pickle
import shutil
import numpy as np
from copy import deepcopy
import subprocess as sp
from collections import OrderedDict, namedtuple
from ase.calculators.castep import CastepOption, Castep
from ase.io.castep import (read_castep_cell, read_param, write_castep_cell,
                           write_param)

from cconv import utils
from cconv.input import parse_convfile
from cconv.output import gp_plot, agr_plot, write_dat, write_report

__version__ = "2.0.1"

__intromsg__ = """
CASTEPconv v. {version}
by Simone Sturniolo
Copyright 2014-2018 Science and Technology Facilities Council

=======

""".format(version=__version__)

# Default folder paths
_in_dir = '_cconv_in'
_out_dir = '_cconv_out'
_serial_dir = 'serial'
_serial_path = os.path.join(_in_dir, 'serial')
_pspot_dir = 'pspot'
_pspot_path = os.path.join(_in_dir, _pspot_dir)

_json_file = '_cconv.json'


def make_ranges(params, kpnbase):

    # Create ranges for the three main scanning variables
    ranges = OrderedDict()

    cutvals = np.arange(params['cutmin'], params['cutmax'] + params['cutstep'],
                        params['cutstep'])
    ranges['cut'] = {
        'values': cutvals.tolist(),
        'labels': list(map(str, cutvals)),
    }

    kpnvals = np.arange(params['kpnmin'], params['kpnmax'] + params['kpnstep'],
                        params['kpnstep'])
    kpnvals = (kpnbase[None, :]*kpnvals[:, None]).astype(int)
    ranges['kpn'] = {
        'values': kpnvals.tolist(),
        'labels': list(map(lambda kpn: 'x'.join(map(str, kpn)), kpnvals))
    }

    if params['fgmmode'] is not None:
        fgmvals = np.arange(params['fgmmin'],
                            params['fgmmax'] + params['fgmstep'],
                            params['fgmstep'])
    else:
        fgmvals = np.array([None])

    ranges['fgm'] = {
        'values': fgmvals.tolist(),
        'labels': list(map(str, fgmvals))
    }

    return ranges


def find_pspots(cell, basename, seedpath='.'):

    pspot_exts = ['.usp', '.uspcc', '.recpot']

    pspot_block = cell.calc.cell.species_pot.value

    if pspot_block is None:
        return

    try:
        elems = set(cell.get_array('castep_custom_species'))
    except KeyError:
        elems = set(cell.get_chemical_symbols())

    pspots = {}

    pspot_lines = pspot_block.split('\n')
    pspot_block = ''  # New block

    to_copy = []  # Potential files to actually copy

    for l in pspot_lines:
        l_s = l.split()
        if len(l_s) != 2:
            # Not a two-word pspot line (e.g. a library)
            pspot_block += l + '\n'
            continue
        # Only extract it if it's a file that can be found; otherwise leave
        # everything as it is
        el, pspsrc = l_s
        if el not in elems:
            # Not a known element
            pspot_block += l + '\n'
            continue
        pspsrc = os.path.join(seedpath, pspsrc)
        if not os.path.isfile(pspsrc):
            # Not a file (e.g. a string)
            ext = os.path.splitext(pspsrc)[1]
            if ext in pspot_exts:
                raise IOError('Pseudopotential file '
                              '{0} not found'.format(pspsrc))
            pspot_block += l + '\n'
            continue

        # If we're here, it's an actual pseudopotential present in the
        # right folder

        pspfile = os.path.split(pspsrc)[1]
        pspdest = os.path.join(basename + _pspot_path, pspfile)
        pspreldest = os.path.join('..', _pspot_dir, pspfile)
        to_copy.append((pspsrc, pspdest))

        pspot_block += '{0} {1}\n'.format(el, pspreldest)

    if len(to_copy) > 0:
        try:
            os.mkdir(basename + _pspot_path)
        except OSError:
            pass  # Directory already exists

        for s, d in to_copy:
            shutil.copy2(s, d)

    cell.calc.cell.species_pot.value = pspot_block


def add_castepopts(a, c8p=True):
    # These adjustments work for convenience
    # since we're not using a castep_keywords.json
    sgen = CastepOption('symmetry_generate', 'B', 'defined', True)
    a.calc.cell._options['symmetry_generate'] = sgen
    kpmp = CastepOption('kpoints_mp_grid', 'B', 'integer vector', [1, 1, 1])
    a.calc.cell._options['kpoints_mp_grid'] = kpmp
    kpoff = CastepOption('kpoints_mp_offset', 'B', 'real vector')
    a.calc.cell._options['kpoints_mp_offset'] = kpoff

    if c8p:
        # This is only valid for Castep 8+
        wcheck = CastepOption('write_checkpoint', 'B', 'string', 'MINIMAL')
        a.calc.param._options['write_checkpoint'] = wcheck


def compile_cmd_line(jobname, cmd_line):

    # Check for redirections

    stdin_file = None
    stdout_file = None

    cmd_line = cmd_line.replace('<seedname>', jobname)

    if ('<' in cmd_line):
        # Take the last of the filenames given after a < but before a >
        stdin_file = cmd_line.split('<')[-1].split('>')[0].split()[-1]

    if ('>' in cmd_line):
        # Take the last of the filenames given after a > but before a <
        stdout_file = cmd_line.split('>')[-1].split('<')[0].split()[-1]

    cmd_line = cmd_line.split('<')[0].split('>')[0].split()

    return cmd_line, stdin_file, stdout_file


WorktreeFolder = namedtuple('WorktreeFolder',
                            ['name', 'dir', 'seed', 'cell',
                             'param', 'castep', 'values', 'labels'])

# Check states
C_READY = 0         # Input ready, no output
C_COMPLETE = 1      # Job completed successfully
C_ERROR = 2         # Some kind of error has occurred


class Worktree(object):

    # A Worktree is a representation of the directory structure holding
    # input files for calculations

    def __init__(self, basename, convpars, convranges):

        self._has_fgm = convranges['fgm']['values'][0] is not None
        jobstr = '{base}_cut_{cut}_kpn_{kpn}' + ('_fgm_{fgm}'
                                                 if self._has_fgm else '')

        self._reuse = convpars['rcalc']
        self._runmode = convpars['rmode']
        self._sreuse = convpars['sruse'] and self._runmode == 'serial'
        self._maxjobs = convpars['maxjobs']
        self._sscript = convpars['subs']
        self._gamma = convpars['gamma']

        # The 'base' values for each structure (first of each range)
        self._basevals = OrderedDict()
        self._baselabels = OrderedDict()

        for key, crange in convranges.items():
            self._basevals[key] = crange['values'][0]
            self._baselabels[key] = crange['labels'][0]

        self._worktree = OrderedDict()
        self._ranges = OrderedDict()  # Store the jobs for each range

        for key, crange in convranges.items():

            self._ranges[key] = []

            for v_i, v in enumerate(crange['values']):
                if v is None:
                    continue

                jobvals = OrderedDict(self._basevals)
                joblabels = OrderedDict(self._baselabels)

                jobvals[key] = v
                joblabels[key] = crange['labels'][v_i]

                # Name?
                jobname = jobstr.format(base=basename, **joblabels)

                self._ranges[key].append(jobname)

                if key != 'cut' and v_i == 0:
                    # Don't repeat first structure pointlessly
                    continue

                # Now write this stuff
                if self._runmode == 'serial':
                    jobdir = basename + _serial_path
                else:
                    jobdir = os.path.join(basename + _in_dir, jobname)

                jobseed = os.path.join(jobdir, jobname)

                self._worktree[jobname] = WorktreeFolder(
                    jobname,
                    jobdir,
                    jobseed,
                    jobseed + '.cell',
                    jobseed + '.param',
                    jobseed + '.castep',
                    jobvals,
                    joblabels)

    @property
    def tree(self):
        return self._worktree

    def write(self, a):

        def set_vals(a, cut, kpn, fgm):
            a.calc.param.cut_off_energy = cut
            a.calc.cell.kpoints_mp_grid = kpn
            if self._gamma:
                # Compute the necessary offset for the gamma point
                offset = [0 if k%2 == 1 else 1.0/(2*k) for k in kpn]
                a.calc.cell.kpoints_mp_offset = offset
            if fgm is not None:
                a.calc.param.fine_gmax = fgm

        set_vals(a, **self._basevals)

        prevjob = None

        for name, job in self._worktree.items():

            set_vals(a, **job.values)

            try:
                os.mkdir(job.dir)
            except OSError:
                pass

            if self._reuse and os.path.isfile(job.cell):
                # If files exist and we shouldn't overwrite, skip
                print('Reusing files for {0}'.format(job.seed))
                continue

            if self._sscript:
                script = open(self._sscript).read()
                script = script.replace('<seedname>', name)
                with open(job.seed, 'w') as outf:
                    outf.write(script)

            if self._sreuse:
                a.calc.param.write_checkpoint = 'ALL'
                if prevjob is not None:
                    a.calc.param.reuse = prevjob.name + '.check'

            print('Writing files for {0}'.format(job.seed))
            write_castep_cell(open(job.cell, 'w'), a)
            write_param(job.param, a.calc.param, force_write=True)

            prevjob = job

    def check(self, names=None):
        """Status code:
        0 - still only running
        1 - complete
        2 - error
        """

        # Our token to find the end of the file
        endstr = 'Total time          ='

        if names is None:
            names = self._worktree.keys()

        # Check if the jobs in the worktree are complete or not
        complete = OrderedDict()
        for name in names:

            complete[name] = C_READY  # Default

            job = self._worktree[name]
            errf = glob.glob(job.seed + '.*.err')

            if len(errf) > 0:
                complete[name] = C_ERROR
                continue

            end_i = -1
            try:
                print(job.castep)
                castlines = open(job.castep).readlines()
                has_end = any(map(lambda l: endstr in l, castlines[-10:]))
            except IOError:
                continue  # No castep file at all, it's C_READY

            complete[name] = C_COMPLETE if has_end else C_READY

        return complete

    def run(self, castep_command='castep <seedname>', wait=True):

        # First, which ones are to run?
        to_run = list(self._worktree.keys())
        if self._reuse:
            jobstate = self.check()
            to_run = [jn for jn, finished in jobstate.items()
                      if finished != C_COMPLETE]
        running = []

        # How many at once?
        maxjobs = 1
        if self._runmode == 'parallel':
            maxjobs = self._maxjobs if self._maxjobs > 0 else len(to_run)

        # Now start running
        while len(to_run) > 0 or len(running) > 0:
            if len(running) < maxjobs and len(to_run) > 0:
                # Push another!
                name = to_run.pop(0)
                job = self._worktree[name]

                # Remove CASTEP and any error files
                to_rm = glob.glob(job.seed + '.*.err')
                to_rm += [job.castep]

                for f in to_rm:
                    try:
                        os.remove(f)
                    except OSError:
                        pass

                cmd_line, stdin, stdout = compile_cmd_line(
                    name, castep_command)

                if (stdin is not None):
                    stdin = open(stdin)
                if (stdout is not None):
                    stdout = open(stdout, 'w')

                print('Running job {0}...'.format(name))
                sp.Popen(cmd_line, stdin=stdin, stdout=stdout, cwd=job.dir)

                if wait:
                    running.append(name)  # If not we're just launching all
            elif wait:
                # Wait for one to finish...
                jobstate = self.check(running)
                # Which ones are finished?
                complete = [(jn, finished) for jn, finished
                            in jobstate.items()
                            if finished != C_READY]
                for jn, res in complete:
                    print('Job {0} completed {1}.'.format(jn,
                                                          'successfully'
                                                          if res == 1 else
                                                          'with an error'))
                    running.remove(jn)
                time.sleep(1)

    def read_data(self):

        # Read the output from all completed jobs
        jobstate = self.check()

        jobdata = OrderedDict([('E', {}), ('F', {}), ('S', {})])

        def get_vals(c):
            nrg = c._energy_total
            frc = c._forces
            strs = c._stress

            return nrg, frc, strs

        for name, job in self._worktree.items():
            if jobstate[name] != C_COMPLETE:
                # Job isn't finished
                utils.warn('Results for {0} missing, skipping.'.format(name))
                continue

            ccalc = Castep(keyword_tolerance=3)

            ccalc.read(job.castep)
            nrg, frc, strs = get_vals(ccalc)

            jobdata['E'][name] = nrg
            jobdata['F'][name] = max(np.linalg.norm(frc, axis=1))
            jobdata['S'][name] = np.linalg.norm(strs)

        # Organise it by ranges
        wtree = self._worktree
        data_curves = OrderedDict()
        for X, jobrange in self._ranges.items():
            # Get an effective jobrange
            jobcomplete = [j for j in jobrange if jobstate[j] == C_COMPLETE]
            data_curves[X] = {'values': np.array([wtree[j].values[X]
                                                  for j in jobcomplete]),
                              'labels': [wtree[j].labels[X]
                                         for j in jobcomplete],
                              'Ys': OrderedDict()
                              }
            for Y, data in jobdata.items():
                data_curves[X]['Ys'][Y] = np.array([data[j]
                                                    for j in jobcomplete])

        return data_curves


"""
The full workflow operates in the following manner:

- read in command line arguments
- read in input files 
- generate input structures
- run calculations (or let the user do it)
- process output

The folder structure generated works as follows:

-|
 |- cconv_in  --|
 |              | serial (only if in serial mode)
 |              | cut_X_kpn_Y (if in parallel mode, for each X and Y)
 |              | ppots (optional)
 |
 |- cconv_out --|
 |              | .dat files
 |              | .gp and .agr files
 |              | report.txt
 |- .json       | <seedname>_cconv.json file, stores the created and used
 |              | ranges from last calculation
"""


def main(seedname, cmdline_task):

    seedpath, basename = os.path.split(seedname)

    utils.check_pyversion()

    print(__intromsg__)

    print('Reading {0}.conv'.format(seedname))

    try:
        with open('{0}.conv'.format(seedname)) as f:
            convpars = parse_convfile(f.read())
    except IOError:
        utils.warn('.conv file not found - using default parameters')
        convpars = parse_convfile()

    task = cmdline_task if cmdline_task is not None else convpars['ctsk']

    # Now open the base cell and param files
    cname = '{0}.cell'.format(seedname)
    print('Reading ' + cname)

    # Necessary because of ASE's annoying messages...
    cfile = read_castep_cell(open(cname),
                             calculator_args={'keyword_tolerance': 3})

    pname = '{0}.param'.format(seedname)
    print('Reading ' + pname)
    try:
        read_param(pname, calc=cfile.calc)
    except FileNotFoundError:
        print('File {0} not found, skipping'.format(pname))

    print('')

    # Now go for clearing
    if task == 'clear':

        to_del = [basename + f
                  for f in [_in_dir, _out_dir, _json_file]]

        print('The following files and folders will be removed:\n\t' +
              '\n\t'.join(to_del))
        ans = utils.safe_input('Continue (y/N)? ')

        if ans.lower() == 'y':
            for f in to_del:
                if not os.path.exists(f):
                    continue
                try:
                    os.remove(f)
                except OSError:
                    shutil.rmtree(f)

        return

    # Strip the files
    cfile.calc.param.task = 'SinglePoint'
    cfile.calc.param.calculate_stress = convpars['cnvstr']

    # Clean up all references to kpoints
    kclean = ['kpoints_mp_grid', 'kpoint_mp_grid', 'kpoints_mp_spacing',
              'kpoint_mp_spacing', 'kpoints_list', 'kpoint_list']
    # These, clean up only in the presence of the relevant option
    if convpars['gamma']:
        kclean += ['kpoints_mp_offset', 'kpoint_mp_offset']

    for k in kclean:
        cfile.calc.cell.__setattr__(k, None)

    # Get the kpoint basis
    invcell = cfile.get_reciprocal_cell()
    kpnbase = np.linalg.norm(invcell, axis=1)
    kpnbase = kpnbase/min(kpnbase)

    # Convergence ranges?
    convranges = make_ranges(convpars, kpnbase)

    # Ask for confirmation
    if task in ('input', 'inputrun', 'all'):
        try:
            os.mkdir(basename + _in_dir)
        except OSError:
            utils.warn('Input directory existing - some files could be '
                       'overwritten.')
            ans = utils.safe_input('Continue (y/N)? ')
            if ans.lower() != 'y':
                return

    if task in ('output', 'all'):
        try:
            os.mkdir(basename + _out_dir)
        except OSError:
            utils.warn('Output directory existing - some files could be '
                       'overwritten.')
            ans = utils.safe_input('Continue (y/N)? ')
            if ans.lower() != 'y':
                return

    # The Worktree object is useful for a number of things in all tasks
    wtree = Worktree(basename, convpars, convranges)

    ### PHASE 1: Input ###
    if task in ('input', 'inputrun', 'all'):
        # Now look for pseudopotentials
        find_pspots(cfile, basename, seedpath)

        add_castepopts(cfile, convpars['c8plus'])

        # If required, rattle the atoms
        if convpars['displ'] != 0:
            cfile.rattle(abs(convpars['displ']))
            cfile.calc.cell.symmetry_generate = None
            cfile.calc.cell.symmetry_ops = None

        wtree.write(cfile)

    ### PHASE 2: Running ###
    if task in ('inputrun', 'all'):

        # Not waiting only makes sense for inputrun
        wait = convpars['jwait'] or (task == 'all')
        wtree.run(convpars['rcmd'], wait)

    ### PHASE 3: Output processing ###
    if task in ('output', 'all'):

        data_curves = wtree.read_data()

        print('Writing output to ' + basename + _out_dir)

        write_dat(basename, data_curves, basename + _out_dir)

        write_report(basename, data_curves, convpars['nrgtol'],
                     convpars['fortol'], convpars['strtol'],
                     basename + _out_dir)

        if convpars['outp'] == 'gnuplot':
            gp_plot(basename, data_curves, basename + _out_dir)
        elif convpars['outp'] == 'grace':
            agr_plot(basename, data_curves, basename + _out_dir)


if __name__ == '__main__':

    seedname, cmdline_task = utils.parse_cmd_args()
    main(seedname, cmdline_task)
