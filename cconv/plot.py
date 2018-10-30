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

import numpy as np

def find_scale(data):
    # Find common scales for the residuals of the given data
    scales = {}

    for xtype, datarange in data.items():

        if len(datarange['values']) == 0:
            continue

        # Maximum deviations from final value
        Ms = {}

        for ytype, yvals in datarange['Ys'].items():
            Ms[ytype] = max(abs(yvals-yvals[-1]))

        scales = {}
        for ytype, M in Ms.items():
            if M > 0:
                scales[ytype] = 10**np.floor(np.log10(Ms['E']/M))
            else:
                scales[ytype] = 1

    return scales


def gp_plot(seedname, data, cwd='.'):

    scales = find_scale(data)

    for xtype, datarange in data.items():

        pass

    return

    plot_details = {'nrg': {
        'col': '3',
        'lc': '1',
    },
        'for': {
        'col': '4',
        'lc': '2',
    },
        'str': {
        'col': '5',
        'lc': '3',
    }
    }

    for x, xlabel in list(_x_types.items()):

        # Check if source .dat file exists
        file_seed = '{seedname}_{x}_conv'.format(seedname=seedname, x=x)
        if not os.path.isfile(file_seed + '.dat'):
            continue
        out_file = open(file_seed + '.gp', 'w')

        out_file.write('set xlabel "{xlabel}"\n'.format(xlabel=xlabel))
        out_file.write('set ylabel "Deviation"\n')
        out_file.write('set xtics nomirror\n')

        # Set the tics specially in case of k-points
        if x == 'kpn':
            out_file.write('set xtics (' +
                           ', '.join(['"{0}" {1}'.format(rs, r)
                                      for r, rs in zip(data[x]['range'],
                                                       data[x]['rangestr'])]) +
                           ') rotate by 45 right\n')

        out_file.write('set ytics nomirror\n')
        out_file.write('plot ')
        for y, legend in list(_y_types.items()):

            if y == 'str' and not cnvstr:
                continue

            scale = scales[x][y]
            lsc = '10^{{{0}}}'.format(-int(log10(scale)))
            out_file.write(('"{file_seed}.dat" using '
                            '1:((${col}-({ref}))*{scale})'
                            ' with linespoints pt 7 lc {lc} '
                            'ti "{legend}",').format(file_seed=file_seed,
                                                     scale=scale,
                                                     legend=legend.format(
                                                         scale=lsc),
                                                     ref=data[x][y][-1],
                                                     **plot_details[y])
                           )
        out_file.write(' 0 lt 0 lc 0 notitle\n')

        out_file.write('pause -1 "Hit return to continue"\n')
        out_file.close()
