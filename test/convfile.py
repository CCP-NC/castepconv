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


class ConvfileTests(unittest.TestCase):

    def test_parsebool(self):

        from cconv.input import parsebool

        self.assertTrue(parsebool('TRUE'))
        self.assertTrue(parsebool('T'))
        self.assertFalse(parsebool('F'))

    def test_parsestr(self):

        from cconv.input import parsestr

        self.assertEqual(parsestr(' this AND tHaT  '),
                         'this and that')

    def test_convpar(self):

        from cconv.input import ConvPar, parsestr

        cp = ConvPar('number', 'n', float, 1.0, validrange=(0, 2))

        with self.assertRaises(ValueError):
            cp.parse('3')

        with self.assertRaises(ValueError):
            cp.parse('-1')

        cp = ConvPar('string', 's', parsestr, 'ok', validoptions=['ok', 'no'])

        self.assertEqual(cp.parse('OK'), 'ok')

        with self.assertRaises(ValueError):
            cp.parse('test')

    def test_parseconv(self):

        from cconv.input import parse_convfile, ConvError

        # Create a mock convfile
        mock_file = """
        # Comment
        convergence_task: ALL
        cutoff_max: 800.0
        """

        parsed = parse_convfile(mock_file)

        self.assertEqual(parsed['ctsk'], 'all')
        self.assertEqual(parsed['cutmax'], 800.0)

        # Now one that should fail
        mock_file = """
        # Comment
        invented_keyword: NO
        """

        with self.assertRaises(ConvError):
            parsed = parse_convfile(mock_file)

        mock_file = """
        # Comment
        cutoff_max: yes
        """

        with self.assertRaises(ValueError):
            parsed = parse_convfile(mock_file)

        mock_file = """
        # Comment
        cutoff_max: 400
        cutoff_min: 500
        """

        with self.assertRaises(ConvError):
            parsed = parse_convfile(mock_file)

        mock_file = """
        # Comment
        cutoff_max: 800
        fine_gmax_mode: max
        fgmoff_min: 2000
        """

        with self.assertRaises(ConvError):
            parsed = parse_convfile(mock_file)
