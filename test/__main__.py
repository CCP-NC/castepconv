import os
import sys
import unittest

# To import castepconv, set the parent directory in the PYTHONPATH
sys.path = [os.path.abspath('../')] + sys.path
import cconv


class CConvTests(unittest.TestCase):

    def test_one(self):
        from cconv.io_freeform import (io_keyword_error, io_freeform_error,
                                       keyword)

        with self.assertRaises(io_keyword_error):
            kw = keyword(1, 'S')

        with self.assertRaises(io_keyword_error):
            kw = keyword('dummy', 'none')

        with self.assertRaises(io_keyword_error):
            kw = keyword('dummy', 'S', 123)

                


if __name__ == "__main__":
    unittest.main()
