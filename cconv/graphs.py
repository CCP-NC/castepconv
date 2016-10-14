import os

_x_types = {'cut': 'Cutoff (eV)', 'kpn': 'k-points', 'fgm': 'Fine Gmax (eV)'}
_y_types = {'for': 'Maximum force (ev/Ang)', 'str': 'Maximum stress (GPa)'}

def gp_graph(seedname, cnvstr=False):

    for x in _x_types:
        for y in _y_types:
            if y == 'str' and not cnvstr:
                continue

            # Check if source .dat file exists
            if not os.path.isfile(os.path.join(seedname + '_' + x  + '_conv.dat')):
                continue
            out_file = open(seedname + '_' + x  + ('_str' if y == 'str' else '') + '_conv.gp', 'w')

            out_file.write('set xlabel "' + _x_types[x] + '"\n')
            out_file.write('set ylabel "Final energy (eV/atom)"\n')
            out_file.write('set y2label "' + _y_types[y] + '"\n')
            out_file.write('set ytics nomirror\n')
            out_file.write('set y2tics\n')
            out_file.write('plot "' + seedname + '_' + x + '_conv.dat" using 1:2 with linespoints pt 7 lc 1 ti "Final energy (eV)",')
            out_file.write('"" using 1:' + {'for': '3', 'str': '4'}[y] +
            ' with linespoints pt 7 lc ' + {'for': '2', 'str': '3'}[y] + ' axes x1y2 ti "' + _y_types[y] + '"\n')
            out_file.write('pause -1 "Hit return to continue"\n')

            out_file.close()

def agr_graph(seedname, data, cnvstr=False):

    # data needs to be a dict formatted as: {'cut':  {'range': [...], 'nrg': [...], 'for': [...] etc.}}


    for x in _x_types:

        if x == 'fgm' and len(data[x]['nrg']) == 0:
            # No fine grid convergence was performed
            continue

        for y in _y_types:
            if y == 'str' and not cnvstr:
                continue

            # Cutoff vs energy and forces

            out_file = open(seedname + '_' + x + ('_str' if y == 'str' else '') + '_conv.agr', 'w')

            x1 = data[x]['range']
            y1 = data[x]['nrg']
            y2 = data[x][y]

            x1rng = (min(x1), max(x1))
            y1rng = (min(y1)-0.1*(max(y1)-min(y1)), max(y1)+0.1*(max(y1)-min(y1)))
            y2rng = (min(y2)-0.1*(max(y2)-min(y2)), max(y2)+0.1*(max(y2)-min(y2)))

            # Set up the graphics

            out_file.write('@version 50123\n')
            out_file.write('@title "' + seedname + ' - Energy and ' + {'for': 'forces', 'str': 'stresses'}[y] + 'vs ' +
            {'cut': 'cutoff', 'kpn': 'k-points', 'fgm': 'fine Gmax'}[x] + '"\n')

            out_file.write('@g0 on\n@g0 hidden false\n@with g0\n')
            out_file.write('@world ' + str(x1rng[0]) + ',' + str(y1rng[0]) + ',' + str(x1rng[1]) + ',' + str(y1rng[1]) + '\n')
            out_file.write('@    view 0.150000, 0.150000, 1.150000, 0.850000\n')
            out_file.write('@    xaxis label "' + _x_types[x] + '"\n')
            out_file.write('@    xaxis tick major ' + str(data[x]['step']) + '\n')
            out_file.write('@    xaxis offset 0.0, 1.0\n')
            out_file.write('@    yaxis label "Final energy (eV/atom)"\n')
            out_file.write('@    yaxis tick major ' + str((y1rng[1]-y1rng[0])/8.0) + '\n')
            out_file.write('@    yaxis offset 0.0, 1.0\n')
            out_file.write('@    s0 hidden false\n@    s0 on\n')
            out_file.write('@    s0 legend "Final energy"\n')
            out_file.write('@    s0 line color 1\n')
            out_file.write('@    s0 symbol 1\n')
            out_file.write('@    s0 symbol size 0.7\n')
            out_file.write('@    s0 symbol color 1\n')
            out_file.write('@    s0 symbol fill color 1\n')
            out_file.write('@    s0 symbol fill pattern 1\n')

            out_file.write('@g1 on\n@g1 hidden false\n@with g1\n')
            out_file.write('@world ' + str(x1rng[0]) + ',' + str(y2rng[0]) + ',' + str(x1rng[1]) + ',' + str(y2rng[1]) + '\n')
            out_file.write('@    view 0.150000, 0.150000, 1.150000, 0.850000\n')
            out_file.write('@    xaxis label ""\n')
            out_file.write('@    xaxis tick off\n')
            out_file.write('@    xaxis tick major ' + str(data[x]['step']) + '\n')
            out_file.write('@    yaxis label "' + _y_types[y] + '"\n')
            out_file.write('@    yaxis tick major ' + str((y2rng[1]-y2rng[0])/8.0) + '\n')
            out_file.write('@    yaxis ticklabel format exponential\n')
            out_file.write('@    yaxis ticklabel prec 1\n')
            out_file.write('@    yaxis offset 1.0, 0.0\n')
            out_file.write('@    yaxis label place opposite\n')
            out_file.write('@    yaxis tick place opposite\n')
            out_file.write('@    yaxis ticklabel place opposite\n')
            out_file.write('@    s0 hidden false\n@    s0 on\n')
            out_file.write('@    s0 legend "' + _y_types[y] + '"\n')
            out_file.write('@    s0 line color 2\n')
            out_file.write('@    s0 symbol 1\n')
            out_file.write('@    s0 symbol size 0.7\n')
            out_file.write('@    s0 symbol color 2\n')
            out_file.write('@    s0 symbol fill color 2\n')
            out_file.write('@    s0 symbol fill pattern 1\n')

            # Input the actual data

            out_file.write('@target G0.S0\n@type xy\n')
            for i, c in enumerate(x1):
                out_file.write(str(c) + '\t' + str(y1[i]) + '\n')
            out_file.write('&\n')

            out_file.write('@target G1.S0\n@type xy\n')
            for i, c in enumerate(x1):
                out_file.write(str(c) + '\t' + str(y2[i]) + '\n')
            out_file.write('&\n')

            out_file.close()

