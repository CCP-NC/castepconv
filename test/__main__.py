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
        from cconv.io_freeform import (IOKeywordError, Keyword)

        with self.assertRaises(IOKeywordError):
            kw = Keyword('dummy', 'none')

    def test_freeform(self):
        from cconv.io_freeform import (IOKeywordError, IOFreeformError,
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
                    Keyword('keyw', 'S:B')])

    def test_utils(self):

        from cconv.utils import floatrange

        frng = floatrange(300, 400, 50)
        self.assertEqual(len(frng), 3)
        self.assertListEqual(frng, [300, 350, 400])


if __name__ == "__main__":
    unittest.main()
