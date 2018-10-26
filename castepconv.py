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

"""Main script"""

import os
import sys
import json
import shutil
import numpy as np
from copy import deepcopy
from collections import OrderedDict
from ase.calculators.castep import CastepOption
from ase.io.castep import (read_castep_cell, read_param, write_castep_cell,
                           write_param)

from cconv import utils
from cconv.io import parse_convfile


__version__ = "2.0"

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
        'name': 'Cut Off Energy (eV)',
        'values': cutvals.tolist(),
        'labels': list(map(str, cutvals)),
    }

    kpnvals = np.arange(params['kpnmin'], params['kpnmax'] + params['kpnstep'],
                        params['kpnstep'])
    kpnvals = (kpnbase[None, :]*kpnvals[:, None]).astype(int)
    ranges['kpn'] = {
        'name': 'K-points grid',
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
        'name': 'Fine GMax (eV)',
        'values': fgmvals.tolist(),
        'labels': list(map(str, fgmvals))
    }

    return ranges


def find_pspots(cell, basename, seedpath='.'):

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


def add_castepopts(cell):
    # These adjustments work for convenience
    # since we're not using a castep_keywords.json
    sgen = CastepOption('symmetry_generate', 'B', 'defined', True)
    cell.calc.cell._options['symmetry_generate'] = sgen
    kpmp = CastepOption('kpoints_mp_grid', 'B', 'integer vector', [1, 1, 1])
    cell.calc.cell._options['kpoints_mp_grid'] = kpmp


def create_worktree(cell, basename, convpars, convranges):

    # Create filename/folder string
    has_fgm = convranges['fgm']['values'][0] is not None
    jobstr = 'cut_{cut}_kpn_{kpn}' + ('_fgm_{fgm}' if has_fgm else '')

    def set_vals(cell, cut, kpn, fgm):
        cell.calc.param.cut_off_energy = cut
        cell.calc.cell.kpoints_mp_grid = kpn
        if fgm is not None:
            cell.calc.param.fine_gmax = fgm

    basevals = OrderedDict()
    baselabels = OrderedDict()
    for key, crange in convranges.items():
        basevals[key] = crange['values'][0]
        baselabels[key] = crange['labels'][0]

    for key, crange in convranges.items():
        # Start by setting the base values
        set_vals(cell, **basevals)

        for v_i, v in enumerate(crange['values']):

            if v is None:
                continue

            jobvals = OrderedDict(basevals)
            joblabels = OrderedDict(baselabels)

            jobvals[key] = v
            joblabels[key] = crange['labels'][v_i]

            set_vals(cell, **jobvals)

            # Name?
            jobname = jobstr.format(**joblabels)

            # Now write this stuff
            if convpars['rmode'] == 'serial':
                jobdir = basename + _serial_path
            else:
                jobdir = os.path.join(basename + _in_dir, jobname)

            try:
                os.mkdir(jobdir)
            except OSError:
                pass

            write_castep_cell(open(os.path.join(jobdir, jobname + '.cell'),
                                   'w'),
                              cell)


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
    if (task in ('input', 'inputrun', 'all')):

        cname = '{0}.cell'.format(seedname)
        print('Reading ' + cname)

        # Necessary because of ASE's annoying messages...
        cfile = read_castep_cell(open(cname),
                                 calculator_args={'keyword_tolerance': 3})

        pname = '{0}.param'.format(seedname)
        print('Reading ' + pname)
        read_param(pname, calc=cfile.calc)

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
    for k in ('kpoints_mp_grid', 'kpoint_mp_grid', 'kpoints_mp_spacing',
              'kpoint_mp_spacing', 'kpoints_list', 'kpoint_list'):
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

    # Now look for pseudopotentials
    find_pspots(cfile, basename, seedpath)

    add_castepopts(cfile)

    # If required, rattle the atoms
    if convpars['displ'] != 0:
        cfile.rattle(abs(convpars['displ']))
        cfile.calc.cell.symmetry_generate = None
        cfile.calc.cell.symmetry_ops = None

    create_worktree(cfile, basename, convpars, convranges)

    # Store the ranges for future reference
    json.dump(convranges, open(seedname + '_convtab.json', 'w'))


if __name__ == '__main__':

    seedname, cmdline_task = utils.parse_cmd_args()
    main(seedname, cmdline_task)
