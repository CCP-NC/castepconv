import os
import sys
import unittest
import tempfile

# To import castepconv, set the parent directory in the PYTHONPATH
sys.path = [os.path.abspath('../')] + sys.path
import cconv


class CConvTests(unittest.TestCase):

    def test_keyword(self):
        from cconv.io_freeform import io_keyword_error, keyword

        with self.assertRaises(io_keyword_error):
            kw = keyword(1, 'S')

        with self.assertRaises(io_keyword_error):
            kw = keyword('dummy', 'none')

        with self.assertRaises(io_keyword_error):
            kw = keyword('dummy', 'S', 123)

    def test_freeform(self):
        from cconv.io_freeform import (io_freeform_error, io_freeform_file,
                                       keyword)

        # Write a test temporary file
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write("""
                keyw:   this
                keyn: 333
            """)
            tmp.flush()

            ioff = io_freeform_file(tmp.name)

            self.assertEqual(ioff.freeform_string('KEYW'), 'this')
            self.assertEqual(ioff.freeform_integer('KEYN'), 333)

            with self.assertRaises(io_freeform_error):
                ioff.freeform_integer('KEYW')

            # Now try a keyworded approach
            ioff = io_freeform_file(tmp.name, keywords=[
                keyword('keyw', 'S:B'),
                keyword('keyn', 'I:B')])

            # And this should fail
            with self.assertRaises(io_freeform_error):
                ioff = io_freeform_file(tmp.name, keywords=[
                    keyword('keyw', 'S:B')])


if __name__ == "__main__":
    unittest.main()
