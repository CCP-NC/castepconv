"""
A specialised version of the freeform file class for castep CELL files
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import re
import math
import random
from collections import namedtuple, OrderedDict

from cconv.units import castep_length_units
from cconv.io.freeform import IOFreeformFile, Keyword

# A definition for an entry in the pseudopotential block
PPotEntry = namedtuple('pspot_entry',
                       ['element', 'entry', 'is_file'])


class CellError(Exception):
    pass


class IOCellFile(IOFreeformFile):

    def __init__(self, fname):

        # This being a cellfile, it is by default tolerant but has certain
        # keywords

        cell_kw = [
            Keyword('positions_frac', 'B'),
            Keyword('positions_abs', 'B'),
            Keyword('lattice_cart', 'B'),
            Keyword('lattice_abc', 'B'),
            Keyword('species_pot', 'B'),
            Keyword('kpoints_mp_grid', 'W')
        ]

        IOFreeformFile.__init__(self, fname=fname, keywords=cell_kw)

        # First, parse the lattice

        hkl_data = []

        u = castep_length_units['ang']  # Default
        kbase = [1, 1, 1]

        if self.freeform_present("lattice_cart"):
            abc_block = self.freeform_block("lattice_cart")
            for l in abc_block:
                l_split = l.lower().strip().split()
                if l_split[0] in castep_length_units:
                    u = castep_length_units[l_split[0]]
                else:
                    try:
                        hkl_data.append([u*float(x) for x in l_split])
                    except ValueError:
                        raise CellError(
                            'Bad formatting in .cell file LATTICE_CART block')
                    if (len(hkl_data) == 3):
                        abc = tuple([math.sqrt(sum([x**2.0 for x in hkl]))
                                     for hkl in hkl_data])
                    if (len(hkl_data) > 3):
                        raise CellError(
                            'Bad formatting in .cell file LATTICE_CART block')
        elif self.freeform_present("lattice_abc"):
            if abc is not None:
                raise CellError('Duplicated LATTICE_* block in .cell file')
            abc_block = self.freeform_block("lattice_abc")
            for l in abc_block:
                l_split = l.lower().strip().split()
                if l_split[0] in castep_length_units:
                    u = castep_length_units[l_split[0]]
                else:
                    try:
                        hkl_data.append([u*float(x) for x in l_split])
                    except ValueError:
                        raise CellError(
                            'Bad formatting in .cell file LATTICE_ABC block')
                    if (len(hkl_data) > 2):
                        raise CellError(
                            'Bad formatting in .cell file LATTICE_ABC block')
            abc = tuple(hkl_data[0])
        else:
            raise CellError('No lattice parameters found in .cell file')

        # Save lattice parameters
        self.abc = abc

        if len(hkl_data) == 2:
            for i in range(0, 3):
                kbase[i] = abs(hkl_data[0][i-1]*hkl_data[0]
                               [i-2]*math.sin(hkl_data[1][i]))
        elif len(hkl_data) == 3:
            for i in range(0, 3):
                kbase[i] = math.sqrt(sum([(hkl_data[i-2][j-2] *
                                           hkl_data[i-1][j-1] -
                                           hkl_data[i-2][j-1] *
                                           hkl_data[i-1][j-2]
                                           )**2.0
                                          for j in range(0, 3)]))

        self.kbase = tuple([k/min(kbase) for k in kbase])

        # Read pseudopotentials

        self.ppots = []
        if self.freeform_present("species_pot"):
            ppot = []
            pspotline_regex = re.compile('[0-9\\.]+\\|[0-9\\.]+'
                                         '\\|[0-9\\.]+\\|')
            pot_block = self.freeform_block("species_pot")
            try:
                for l in pot_block:
                    l_split = l.strip().split()
                    # What type is it?
                    if len(l_split) == 1:
                        # It's a single library
                        seed, ext = os.path.splitext(l_split[0])
                        if ext == '':
                            # Just a library, leave unchanged
                            ppot.append(PPotEntry(None, seed, False))
                        elif ext.lower() == '.otfglib':
                            # Library file
                            ppot.append(PPotEntry(None, l_split[0],
                                                  True))
                        else:
                            raise RuntimeError()
                    elif len(l_split) == 2:
                        # File, library or string?
                        seed, ext = os.path.splitext(l_split[1])
                        if (len(pspotline_regex.findall(l_split[1])) > 1 or
                                ext == ''):
                            # String or library
                            ppot.append(
                                PPotEntry(l_split[0], l_split[1], False))
                        elif ext.lower() in ('.usp', '.recpot', '.ucf'):
                            ppot.append(PPotEntry(l_split[0],
                                                  l_split[1],
                                                  True))
                        else:
                            raise RuntimeError()

            except RuntimeError:
                raise CellError(
                    'Bad formatting in .cell file SPECIES_POT block')

            self.ppots = ppot

    def displace_cell(self, d=0.05):

        u = castep_length_units['ang']  # Default

        # Are positions absolute or fractionary?
        pos_is_abs = self.freeform_present('positions_abs')
        if pos_is_abs:
            pos_block = self.freeform_block('positions_abs')
        elif self.freeform_present('positions_frac'):
            pos_block = self.freeform_block('positions_frac')
        else:
            raise CellError('No positions found in file')

        pos_block_displ = []

        for i, l in enumerate(pos_block):
            l_split = l.strip().split()
            if pos_is_abs and l_split[0].lower() in castep_length_units:
                u = length_units[l_split[0].lower()]
                pos_block_displ.append('ang\n')
            else:
                if len(l_split) != 4:
                    raise CellError(
                        'Bad formatting in .cell file POSITION_* block')
                try:
                    xyz = [float(x) for x in l_split[1:]]
                    displ = [(random.random()-0.5)*2*d for i in range(3)]
                    if pos_is_abs:
                        l = '{0} {1} {2} {3}'.format(
                            *[l_split[0]] + [u*x+displ[j]
                                             for j, x in enumerate(xyz)])
                    else:
                        l = '{0} {1} {2} {3}'.format(
                            *[l_split[0]] + [x+displ[j]/self.abc[j]
                                             for j, x in enumerate(xyz)])
                except ValueError:
                    raise CellError(
                        'Bad formatting in .cell file POSITION_* block')
                pos_block_displ.append(l)

        if pos_is_abs:
            self.freeform_block('positions_abs', pos_block_displ)
        else:
            self.freeform_block('positions_frac', pos_block_displ)

        return self

    def set_kpoint_grid(self, small_k=1):

        kpn_grid = tuple([int(small_k * e)
                          for e in self.kbase])

        # Clean all pre-existing kpoints keywords
        for k in ("kpoints_mp_grid", "kpoint_mp_grid", "kpoints_mp_spacing",
                  "kpoint_mp_spacing", "kpoints_list", "kpoint_list"):
            self.freeform_remove(k)

        self.freeform_integer_vector('kpoints_mp_grid', kpn_grid)
