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

"""Plotting functions"""

import os
import numpy as np
from cconv import utils

## Default properties ##

_x_types = {
    'cut': {
        'xlabel': 'Cut Off Energy (eV)',
    },
    'kpn': {
        'xlabel': 'K-Points Grid'
    },
    'fgm': {
        'xlabel': 'Fine Gmax (eV)'
    }
}

_y_types = {
    'E': {
        'lc': 1,
        'col': 2,
        'legend': 'Final energy'
    },
    'F': {
        'lc': 2,
        'col': 3,
        'legend': 'Max force'
    },
    'S': {
        'lc': 3,
        'col': 4,
        'legend': 'Stress (norm)'
    }
}

## Templates ##

_gp_template = """
set xlabel '{xlabel}'
set ylabel 'Deviation'
set xtics nomirror
set ytics nomirror
{xticsline}
plot {plot}, 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
"""

_gp_plot_template = ('"{file}" using 1:((${col}-({ref}))*{scale})'
                     ' with linespoints pt 7 lc {lc} ti "{legend}"')


def find_scale(data_curves):
    # Find common scales for the residuals of the given data
    scales = {}

    for xtype, xdata in data_curves.items():

        scales[xtype] = {}

        if len(xdata['values']) == 0:
            continue

        # Maximum deviations from final value
        Ms = {}

        for ytype, yvals in xdata['Ys'].items():
            Ms[ytype] = max(abs(yvals-yvals[-1]))

        for ytype, M in Ms.items():
            if M > 0:
                scales[xtype][ytype] = 10**np.floor(np.log10(Ms['E']/M))
            else:
                scales[xtype][ytype] = 1

    return scales


def gp_plot(seedname, data_curves, cwd='.'):

    scales = find_scale(data_curves)

    for xtype, xdata in data_curves.items():

        if len(xdata['values']) == 0:
            # No data
            continue

        # Do we have the corresponding .dat file?
        dataf = '{0}_{1}_conv.dat'.format(seedname, xtype)

        if not os.path.isfile(os.path.join(cwd, dataf)):
            utils.warn('Data file {0} not found, skipping.'.format(dataf))
            continue

        # Plot lines
        plines = []
        for ytype, yvals in xdata['Ys'].items():

            args = dict(_y_types[ytype])

            args['file'] = dataf
            args['ref'] = yvals[-1]
            args['scale'] = scales[xtype][ytype]

            plines.append(_gp_plot_template.format(**args))

        args = dict(_x_types[xtype])
        if xtype == 'kpn':
            xvals = np.prod(xdata['values'], axis=1)
            args['xticsline'] = ('set xtics (' +
                                 ', '.join(['"{1}" {0}'.format(v, l)
                                            for v, l in zip(xvals,
                                                            xdata['labels'])]
                                           ) + ') rotate by 45 right')
        else:
            args['xticsline'] = ''

        args['plot'] = ', '.join(plines)

        gp_file = _gp_template.format(**args)

        with open(os.path.join(cwd, '{0}_{1}_conv.gp'.format(seedname,
                                                             xtype)),
                  'w') as f:
            f.write(gp_file)

    return
