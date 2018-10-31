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
        'title': 'Convergence vs. plane wave cut off'
    },
    'kpn': {
        'xlabel': 'K-Points Grid',
        'title': 'Convergence vs. k-points'
    },
    'fgm': {
        'xlabel': 'Fine Gmax (eV)',
        'title': 'Convergence vs. fine grid cut off'
    }
}

_y_types = {
    'E': {
        'set': 0,
        'lc': 1,
        'col': 2,
        'legend': 'Final energy (eV)'
    },
    'F': {
        'set': 1,
        'lc': 2,
        'col': 3,
        'legend': 'Max force ({scale} ev/Ang)'
    },
    'S': {
        'set': 2,
        'lc': 3,
        'col': 4,
        'legend': 'Stress (norm, {scale} GPa)'
    }
}

## Templates ##

_gp_template = """
set xlabel '{xlabel}'
set ylabel 'Deviation'
set title '{title}'
set xtics nomirror
set ytics nomirror
{xticsblock}
plot {plot}, 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
"""

_gp_plot_template = ('"{file}" using 1:((${col}-({ref}))*{scale})'
                     ' with linespoints pt 7 lc {lc} ti "{legend}"')

_agr_template = """
@version 50123
@title "{title}"
@g0 on
@g0 hidden false
@with g0
@    world {world}
@    view 0.150000, 0.150000, 1.150000, 0.850000
@    xaxis label "{xlabel}"
{xticsblock}
@    xaxis offset 0.0, 1.0
{plot}
{data}
"""

_agr_plot_template = """
@    yaxis label "Final energy (eV/atom)"
@    yaxis tick major {ytics}
@    yaxis offset 0.0, 1.0
@    s{set} hidden false
@    s{set} on
@    s{set} legend "{legend}"
@    s{set} line color {lc}
@    s{set} symbol 1
@    s{set} symbol color {lc}
@    s{set} symbol fill color {lc}
@    s{set} symbol fill pattern 1
"""

_agr_data_template = """
@target G0.S{set}
@type xy
{tabdata}
&
"""


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
            lsc = '10^{{{0}}}'.format(-int(np.log10(scales[xtype][ytype])))
            args['legend'] = args['legend'].format(scale=lsc)

            if (max(yvals)-min(yvals)) == 0:
                continue

            plines.append(_gp_plot_template.format(**args))

        args = dict(_x_types[xtype])
        if xtype == 'kpn':
            xvals = np.prod(xdata['values'], axis=1)
            args['xticsblock'] = ('set xtics (' +
                                  ', '.join(['"{1}" {0}'.format(v, l)
                                             for v, l in zip(xvals,
                                                             xdata['labels'])]
                                            ) + ') rotate by 45 right')
        else:
            args['xticsblock'] = ''

        args['plot'] = ', '.join(plines)

        gp_file = _gp_template.format(**args)

        with open(os.path.join(cwd, '{0}_{1}_conv.gp'.format(seedname,
                                                             xtype)),
                  'w') as f:
            f.write(gp_file)

    return


def agr_plot(seedname, data_curves, cwd='.'):

    # data needs to be a dict formatted as: {'cut':  {'range': [...], 'nrg':
    # [...], 'for': [...] etc.}}

    scales = find_scale(data_curves)

    for xtype, xdata in data_curves.items():

        xvals = xdata['values']
        if xtype == 'kpn':
            xvals = np.prod(xvals, axis=1)

        if len(xvals) == 0:
            # No data
            continue

        xrng = (np.amin(xvals),
                np.amax(xvals))
        world = [xrng[0], 0, xrng[1], 0]

        # Plot lines and data blocks
        plines = []
        dblocks = []
        for ytype, yvals in xdata['Ys'].items():

            args = dict(_y_types[ytype])

            # Scaled data
            yerr = (yvals - yvals[-1])*scales[xtype][ytype]
            yrng = (np.clip(np.amin(yerr), None, 0),
                    np.clip(np.amax(yerr), 0, None))
            yspan = np.diff(yrng)[0]

            if yspan == 0:
                continue

            world[1] = min(world[1], yrng[0])
            world[3] = max(world[3], yrng[1])

            lsc = '1E{0}'.format(-int(np.log10(scales[xtype][ytype])))
            args['legend'] = args['legend'].format(scale=lsc)

            if yspan > 0:
                args['ytics'] = 0.5*10**(np.floor(np.log10(yspan)))
            else:
                args['ytics'] = 1

            plines.append(_agr_plot_template.format(**args))

            args['tabdata'] = '\n'.join(['{0}\t{1}'.format(x, y)
                                         for x, y in zip(xvals, yerr)])

            dblocks.append(_agr_data_template.format(**args))

        # Set the tics specially in case of k-points
        if xtype == 'kpn':
            xticsblock = ('@    xaxis  tick spec type both\n' +
                          '@    xaxis  tick spec {0}\n' +
                          '@    xaxis ticklabel angle 315\n'
                          ).format(len(xvals))

            for i, (v, l) in enumerate(zip(xvals, xdata['labels'])):
                xticsblock += '@    xaxis tick major {0},{1}\n'.format(i, v)
                xticsblock += '@    xaxis ticklabel {0},"{1}"\n'.format(i, l)
        else:
            xticsblock = '@    xaxis tick major {0}\n'.format(xvals[1] -
                                                              xvals[0])

        args = dict(_x_types[xtype])
        args['world'] = ','.join(map(str, world))
        args['xticsblock'] = xticsblock
        args['plot'] = '\n'.join(plines)
        args['data'] = '\n'.join(dblocks)

        agr_file = _agr_template.format(**args)

        with open(os.path.join(cwd, '{0}_{1}_conv.agr'.format(seedname,
                                                              xtype)),
                  'w') as f:
            f.write(agr_file)
