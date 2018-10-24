"""Unit testing routines"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import unittest
import tempfile

# To import castepconv, set the parent directory in the PYTHONPATH
sys.path = [os.path.abspath('../')] + sys.path
import cconv

# Test files directory
testdir = os.path.join(os.path.split(__file__)[0], 'files')


class CConvTests(unittest.TestCase):

    def test_keyword(self):
        from cconv.io.freeform import (IOKeywordError, Keyword)

        with self.assertRaises(IOKeywordError):
            kw = Keyword('dummy', 'none')

    def test_freeform(self):
        from cconv.io.freeform import (IOKeywordError, IOFreeformError,
                                       Keyword, IOFreeformFile)

        # Write a test temporary file
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write(b"""
                keyw: this
                keyn: 333
            """)
            tmp.flush()

            ioff = IOFreeformFile(tmp.name)

            self.assertEqual(ioff.freeform_string('KEYW'), 'this')
            self.assertEqual(ioff.freeform_integer('KEYN'), 333)

            with self.assertRaises(IOFreeformError):
                ioff.freeform_integer('KEYW')

            # Now try a keyworded approach
            ioff = IOFreeformFile(tmp.name, keywords=[
                Keyword('keyw', 'S:B'),
                Keyword('keyn', 'I:B')])

            # And this should fail
            with self.assertRaises(IOFreeformError):
                ioff = IOFreeformFile(tmp.name, keywords=[
                    Keyword('keyw', 'S:B')], tolerant=False)

        # Now test an empty, tolerant file
        ioff = IOFreeformFile()
        ioff.freeform_integer_vector('Test', value=[3, 3, 3])
        self.assertEqual(ioff.keyvals['TEST'], '3 3 3')
        ioff = IOFreeformFile()
        ioff.freeform_real_vector('Test', value=[2.0, 3.0, 4.5])
        self.assertListEqual(ioff.freeform_real_vector('TEST'), [2, 3, 4.5])

    def test_cell(self):

        from cconv.io.cell import IOCellFile

        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write(b"""
                %block lattice_cart
                    2   0   0
                    0   2   0
                    0   0   2
                %endblock lattice_cart

                %block positions_frac
                H 0 0 0
                H 0.5 0.5 0.5
                %endblock positions_frac

                kpoints_mp_grid 3 1 2

                %block species_pot
                Si  Si_00.usp
                %endblock species_pot
            """)
            tmp.flush()

            iocf = IOCellFile(tmp.name)

            self.assertListEqual(iocf
                                 .freeform_integer_vector('kpoints_mp_grid'),
                                 [3, 1, 2])

            iocf.displace_cell()

            displ = sum([float(x)**2 for x in
                         iocf.freeform_block('positions_frac')[0].split()[1:]])
            self.assertTrue(displ**0.5 <= 0.15)

            iocf.set_kpoint_grid(2)

            self.assertListEqual(iocf
                                 .freeform_integer_vector('kpoints_mp_grid'),
                                 [2, 2, 2])

            # Check copying
            iocf_copy = iocf.copy()
            iocf_copy.ppots[0] = 'N/A'
            # Verify that they're pointing at different locations
            self.assertNotEqual(iocf_copy.ppots[0], iocf.ppots[0])

    def test_scan(self):
        from cconv.scan import CastepScan
        from cconv.io import IOFreeformFile, IOCellFile

        cf = IOCellFile(os.path.join(testdir, 'Si.cell'))
        pf = IOFreeformFile(os.path.join(testdir, 'Si.param'))

        sc = CastepScan(cf, pf)

        sc.set_cut_range([300, 400, 500])
        sc.set_kpn_range([1, 2, 4])

        cut_params = [p for c, p in sc._ranges['cut']['files']]
        self.assertAlmostEqual(cut_params[2]
                               .freeform_physical('cut_off_energy',
                                                  'E', unit='ev'),
                               500)
        kpn_cells = [c for c, p in sc._ranges['kpn']['files']]
        self.assertListEqual(kpn_cells[2]
                             .freeform_integer_vector('kpoints_mp_grid'),
                             [4, 4, 4])

    def test_utils(self):

        from cconv.utils import floatrange

        frng = floatrange(300, 400, 50)
        self.assertEqual(len(frng), 3)
        self.assertListEqual(frng, [300, 350, 400])


if __name__ == "__main__":
    unittest.main()
