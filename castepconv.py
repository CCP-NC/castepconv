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
from ase.io.castep import read_castep_cell, read_param

from cconv import utils
from cconv.io import parse_convfile


__version__ = "2.0"

__intromsg__ = """
CASTEPconv v. {version}
by Simone Sturniolo
Copyright 2014-2018 Science and Technology Facilities Council

=======

""".format(version=__version__)


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
        fgmvals = np.array([])

    ranges['fgm'] = {
        'name': 'Fine GMax (eV)',
        'values': fgmvals.tolist(),
        'labels': list(map(str, fgmvals))
    }

    return ranges


def main():

    utils.check_pyversion()

    print(__intromsg__)

    seedname, cmdline_task = utils.parse_cmd_args()

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

        to_del = ['input', 'output']

        print('The following files and folders will be removed:\n\t' +
              '\n\t'.join(to_del))
        ans = utils.safe_input('Continue (y/N)?')

        if ans == 'y':
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

    # Store the ranges for future reference
    json.dump(convranges, open(seedname + '_convtab.json', 'w'))


if __name__ == '__main__':
    main()
