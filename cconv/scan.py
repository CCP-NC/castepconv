"""
A Scan class, representing a series of structures with different parameters
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from cconv.io.cell import IOCellFile


class CastepScan(object):
    """CastepScan class - represents a series of structures corresponding to
    various settings (cutoff, k-points etc.) to explore.
    """

    def __init__(self, base_cell, base_param, base_cut=400, base_kpn=1,
                 use_stress=False):

        if not isinstance(base_cell, IOCellFile):
            raise ValueError('base_cell is not IOCellFile')

        self._cell = base_cell.copy()
        self._param = base_param.copy()

        # These are necessary for everything
        self._param.freeform_string('task', 'SinglePoint')
        if use_stress:
            self._param.freeform_boolean('calculate_stress', True)
        # Set the base cutoff
        self._param.freeform_physical('cut_off_energy', 'E', base_cut, 'ev')
        self._cell.set_kpoint_grid(base_kpn)

        self._ranges = {}

    def set_cut_range(self, cut_range=[]):
        self._ranges['cut'] = {
            'longname': 'Cut Off Energy (eV)'
            'values': cut_range,
            'labels': [str(c) for c in cut_range],
        }

        files = []

        for i, cut in enumerate(cut_range):
            cf = self._cell.copy()
            pf = self._param.copy()

            pf.freeform_physical('cut_off_energy', 'E', cut, 'ev')
            files.push((cf, pf))

        self._ranges['cut']['files'] = files

    def set_kpn_range(self, kpn_range=[]):
        self._ranges['kpn'] = {
            'longname': 'K-points Grid'
            'values': kpn_range
        }

        files = []

        for i, kpn in enumerate(kpn_range):
            cf = self._cell.copy()
            pf = self._param.copy()

            cf.set_kpoint_grid(kpn)

            files.push((cf, pf))

        labels = ['x'.join([str(k) for k in
                            c.freeform_integer_vector('kpoints_mp_grid')])
                  for c, p in files],

        self._ranges['kpn']['labels'] = labels
        self._ranges['kpn']['files'] = files
