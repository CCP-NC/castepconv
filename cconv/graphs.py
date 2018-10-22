# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
from math import log10, floor

_x_types = {'cut': 'Cutoff (eV)', 'kpn': 'k-points', 'fgm': 'Fine Gmax (eV)'}
_y_types = {'nrg': 'Final energy (eV/Atom)',
            'for': 'Maximum force ({scale} ev/Ang)',
            'str': 'Maximum stress ({scale} GPa)'}

# Templates for printing
_gp_template = """
set xlabel "{xlabel}"\n
set ylabel "Deviation"
set xtics nomirror
set xtics ({xtics}) {xticspos}
set ytics nomirror
plot {plotlines} 0 lt 0 lc 0 notitle
pause -1 "Hit return to continue"
"""

_gp_plot_template = ('"{file_seed}.dat" using 1:((${col}-({ref}))*{scale})'
                     ' with linespoints pt 7 lc {lc} ti "{legend}",')


def find_scale(data):
    # Find common scales for the residuals of the given data
    scales = {}

    for x in _x_types:

        if data[x]['range'][0] is None:
            continue

        # Maximum deviations from final value
        Ms = {y: max([abs(d-data[x][y][-1]) for d in data[x][y]])
              if y in data[x] else 0
              for y in _y_types}
        # Proposed scaling factors
        scales[x] = {y: 10**floor(log10(Ms['nrg']/m)) if m > 0 else 1
                     for y, m in Ms.items()}

    return scales


def gp_graph(seedname, data, cnvstr=False):

    # data needs to be a dict formatted as: {'cut':  {'range': [...], 'nrg':
    # [...], 'for': [...] etc.}}

    scales = find_scale(data)

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

    for x, xlabel in _x_types.items():

        # Check if source .dat file exists
        file_seed = '{seedname}_{x}_conv'.format(seedname=seedname, x=x)
        if not os.path.isfile(file_seed + '.dat'):
            continue
        out_file = open(file_seed + '.gp', 'w')

        # X Tics
        xtics = ', '.join(['"{0}" {1}'.format(rs, r)
                           for r, rs in zip(data[x]['range'],
                                            data[x]['rangestr'])])
        xticspos = 'rotate by 45 right' if x == 'kpn' else ''

        # Individual line plots
        plotlines = ''

        for y, legend in _y_types.items():

            if y == 'str' and not cnvstr:
                continue

            scale = scales[x][y]
            lsc = '10^{{{0}}}'.format(-int(log10(scale)))

            plotlines += _gp_plot_template.format(file_seed=file_seed,
                                                  scale=scale,
                                                  legend=legend.format(
                                                      scale=lsc),
                                                  ref=data[x][y][-1],
                                                  **plot_details[y])

        out_file.write(_gp_template.format(xlabel=xlabel,
                                           xtics=xtics,
                                           xticspos=xticspos,
                                           plotlines=plotlines))

        out_file.close()


def agr_graph(seedname, data, cnvstr=False):

    # data needs to be a dict formatted as: {'cut':  {'range': [...], 'nrg':
    # [...], 'for': [...] etc.}}

    scales = find_scale(data)

    plot_details = {'nrg': {
        'set': '0',
        'lc': '1',
    },
        'for': {
        'set': '1',
        'lc': '2',
    },
        'str': {
        'set': '2',
        'lc': '3',
    }
    }

    for x, xname in _x_types.items():

        if x == 'fgm' and len(data[x]['nrg']) == 0:
            # No fine grid convergence was performed
            continue

        file_seed = '{seedname}_{x}_conv'.format(seedname=seedname, x=x)
        out_file = open(file_seed + '.agr', 'w')

        x1 = data[x]['range']
        x1rng = (min(x1), max(x1))

        # Set up the graphics

        out_file.write('@version 50123\n')
        out_file.write(('@title "{seedname} - Convergence '
                        'vs {x}"\n').format(seedname=seedname, x=xname))
        out_file.write('@g0 on\n@g0 hidden false\n@with g0\n')

        # Find the optimal size for the world
        y1rng = [0, 0]
        for y in _y_types:
            if y == 'str' and not cnvstr:
                continue
            y1 = data[x][y]
            s = scales[x][y]
            y1rng[0] = min(y1rng[0], min([(p-y1[-1])*s for p in y1]))
            y1rng[1] = max(y1rng[1], max([(p-y1[-1])*s for p in y1]))

        y1span = y1rng[1]-y1rng[0]
        y1rng = [y1rng[0]-0.1*y1span, y1rng[1]+0.1*y1span]

        out_file.write('@world ' + str(x1rng[0]) + ',' + str(
            y1rng[0]) + ',' + str(x1rng[1]) + ',' + str(y1rng[1]) + '\n')

        out_file.write(
            '@    view 0.150000, 0.150000, 1.150000, 0.850000\n')
        out_file.write('@    xaxis label "' + _x_types[x] + '"\n')

        # Set the tics specially in case of k-points
        if x == 'kpn':
            out_file.write("""@    xaxis  tick spec type both
@    xaxis  tick spec {0}
@    xaxis ticklabel angle 315\n""".format(len(data[x]['range'])))
            for i, (r, rs) in enumerate(zip(data[x]['range'],
                                            data[x]['rangestr'])):
                out_file.write('@    xaxis tick major {0},{1}\n'.format(i, r))
                out_file.write('@    xaxis ticklabel {0},"{1}"\n'.format(i,
                                                                         rs))
        else:
            out_file.write('@    xaxis tick major ' +
                           str(data[x]['step']) + '\n')
        out_file.write('@    xaxis offset 0.0, 1.0\n')

        for y, legend in _y_types.items():
            if y == 'str' and not cnvstr:
                continue

            # Cutoff vs energy and forces

            y1 = data[x][y]
            set_header = '@    s{set}'.format(**plot_details[y])
            lsc = legend.format(scale='10\\S{0}\\N'.format(
                -int(log10(scales[x][y]))))

            out_file.write("""
@    yaxis label "Final energy (eV/atom)"
@    yaxis tick major {ytick}
@    yaxis offset 0.0, 1.0
{set_header} hidden false
{set_header} on
{set_header} legend "{legend}"
{set_header} line color {lc}
{set_header} symbol 1
{set_header} symbol color {lc}
{set_header} symbol fill color {lc}
{set_header} symbol fill pattern 1
""".format(set_header=set_header, legend=lsc,
                ytick=0.25*10**(floor(log10((y1rng[1]-y1rng[0])))),
                **plot_details[y]))

        # Input the actual data

        for y in _y_types:

            if y == 'str' and not cnvstr:
                continue

            out_file.write('@target G0.S{0}\n@type xy\n'.format(
                plot_details[y]['set']))
            for i, c in enumerate(x1):
                out_file.write('{x}\t{y}\n'.format(
                    x=c, y=(data[x][y][i]-data[x][y][-1])*scales[x][y]))
            out_file.write('&\n')

        out_file.close()
