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

    def test_utils(self):

        from cconv.utils import floatrange

        frng = floatrange(300, 400, 50)
        self.assertEqual(len(frng), 3)
        self.assertListEqual(frng, [300, 350, 400])


if __name__ == "__main__":
    unittest.main()
