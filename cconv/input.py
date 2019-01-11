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

"""Input/Output that needs internal handling"""

import os
import numpy as np
from cconv import utils


class ConvError(Exception):
    pass


def parsebool(v):
    return (str(v).lower().strip() in ('true', 't'))


def parsestr(v):
    return str(v).lower().strip()


def parsepath(v):
    return str(v).strip()


class ConvPar(object):
    """A definition of a convergence parameter"""

    def __init__(self, name, shortname, parser, default=None,
                 validrange=None,
                 validoptions=None):
        """
        There are no checks, but it's important:

        validrange to be used only for numbers, as [min, max]

        validoptions to be used for strings, and they all have to be stripped
        and lowercase

        """
        self.name = name
        self.sname = shortname
        self.parser = parser
        self.default = default
        self.range = validrange
        self.options = validoptions

    def parse(self, v):

        v = self.parser(v)

        if self.range is not None:
            if ((self.range[0] is not None and v < self.range[0]) or
                    (self.range[1] is not None and v > self.range[1])):
                raise ValueError(('Value of parameter {0} outside of '
                                  'admissible range').format(self.name))

        if self.options is not None:
            if (v not in self.options):
                raise ValueError(
                    'Invalid value for parameter {0}'.format(self.name))

        return v

    def __repr__(self):
        return '{name} [{sname}, {parser}, {default}]'.format(**self.__dict__)


# List of parameters: entries are [short name, type, default value]
_conv_parameters = {
    # String parameters
    'convergence_task': ConvPar('convergence_task', 'ctsk', parsestr, 'input',
                                validoptions=['input', 'inputrun', 'output',
                                              'all', 'clear']),
    'running_mode': ConvPar('running_mode', 'rmode', parsestr, 'serial',
                            validoptions=['serial', 'parallel']),
    'output_type': ConvPar('output_type', 'outp', parsestr, 'gnuplot',
                           validoptions=['gnuplot', 'grace']),
    'fine_gmax_mode': ConvPar('fine_gmax_mode', 'fgmmode', parsestr,
                              validoptions=['min', 'max']),
    # The following are case sensitive, so we don't use parsestr
    'running_command': ConvPar('running_command', 'rcmd', parsepath,
                               'castep <seedname> -dryrun'),
    'dryrun_command': ConvPar('dryrun_command', 'drcmd', parsepath),
    'submission_script': ConvPar('submission_script', 'subs', parsepath),
    # Float parameters
    #   eV
    'cutoff_min': ConvPar('cutoff_min', 'cutmin', float, 400.0,
                          validrange=[0, None]),
    'cutoff_max': ConvPar('cutoff_max', 'cutmax', float, 800.0,
                          validrange=[0, None]),
    'cutoff_step': ConvPar('cutoff_step', 'cutstep', float, 100.0,
                           validrange=[0, None]),
    'fgmoff_min': ConvPar('fgmoff_min', 'fgmmin', float,
                          validrange=[0, None]),
    'fgmoff_max': ConvPar('fgmoff_max', 'fgmmax', float,
                          validrange=[0, None]),
    'fgmoff_step': ConvPar('fgmoff_step', 'fgmstep', float, 100.0,
                           validrange=[0, None]),
    #   Ang
    'displace_atoms': ConvPar('displace_atoms', 'displ', float, 0.0),
    #   eV/Atom
    'final_energy_delta': ConvPar('final_energy_delta', 'nrgtol', float, 1e-5,
                                  validrange=[0, None]),
    #   eV/Ang
    'forces_delta': ConvPar('forces_delta', 'fortol', float, 5e-2,
                            validrange=[0, None]),
    #   GPa
    'stresses_delta': ConvPar('stresses_delta', 'strtol', float, 0.1,
                              validrange=[0, None]),
    # Integer parameters
    'kpoint_n_min': ConvPar('kpoint_n_min', 'kpnmin', int, 1,
                            validrange=[0, None]),
    'kpoint_n_max': ConvPar('kpoint_n_max', 'kpnmax', int, 4,
                            validrange=[0, None]),
    'kpoint_n_step': ConvPar('kpoint_n_step', 'kpnstep', int, 1,
                             validrange=[0, None]),
    'max_parallel_jobs': ConvPar('max_parallel_jobs', 'maxjobs', int, 0,
                                 validrange=[-1, None]),
    # Boolean parameters
    'job_wait': ConvPar('job_wait', 'jwait', parsebool, True),
    'converge_stress': ConvPar('converge_stress', 'cnvstr', parsebool, False),
    'reuse_calcs': ConvPar('reuse_calcs', 'rcalc', parsebool, False),
    'serial_reuse': ConvPar('serial_reuse', 'sruse', parsebool, True),
    'castep_8plus': ConvPar('castep_8plus', 'c8plus', parsebool, True),
    'include_gamma': ConvPar('include_gamma', 'gamma', parsebool, False)
}


def param_check(params):
    """Check the values of the parameters for internal contradictions"""

    if params['cutmax'] < params['cutmin']:
        raise ConvError('Invalid parameter - must be cutmax > cutmin')

    if params['kpnmax'] < params['kpnmin']:
        raise ConvError('Invalid parameter - must be kpnmax > kpnmin')

    # Fine grid max checks
    if params['fgmmode'] is not None:
        gscale = 1.75   # This is a CASTEP default
        lbound = gscale**2*(params['cutmin'] if params['fgmmode'] == 'min'
                            else params['cutmax'])
        params['fgmmin'] = (lbound if params['fgmmin'] is None
                            else params['fgmmin'])
        params['fgmmax'] = (params['fgmmin']+3*params['fgmstep']
                            if params['fgmmax'] is None
                            else params['fgmmax'])
        if any([params['fgmmax'] <= lbound, params['fgmmin'] < lbound,
                params['fgmmin'] > params['fgmmax']]):
            raise ConvError('Invalid parameter - must be fgmmax > fgmmin >= '
                            '{0}*cutoff_{1}'.format(gscale**2,
                                                    params['fgmmode']))

    if '<seedname>' not in params['rcmd']:
        utils.warn('Running command does not contain a <seedname> tag.'
                   ' This is likely erroneous and needs checking.')

    if params['subs'] is not None:
        try:
            with open(params['subs']) as f:
                if '<seedname>' not in f.read():
                    utils.warn('Submission script does not contain a'
                               '<seedname> tag. This is likely erroneous and'
                               ' needs checking.')
        except IOError:
            raise ConvError('Submission script file does not exist')


def parse_convfile(cfile=''):
    """Parse .conv file. Accepts file object or contents. If no argument is
    passed, will return default values"""

    if hasattr(cfile, 'readlines'):
        clines = cfile.readlines()
    else:
        # It should be the contents already
        clines = cfile.split('\n')

    params = {}
    for cp in _conv_parameters.values():
        params[cp.sname] = cp.default

    for i, l in enumerate(clines):

        l = l.split('#', 1)[0].strip()

        if l == '':
            continue

        token, value = l.split(':', 1)
        token = parsestr(token)

        pardef = _conv_parameters.get(token, None)

        if pardef is None:
            raise ConvError('Bad token in .conv file at line {0}'.format(i+1))

        value = pardef.parse(value)

        params[pardef.sname] = value

    param_check(params)

    return params
