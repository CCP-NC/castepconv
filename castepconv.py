#!/usr/bin/env python

"""
CASTEP convergence automation tool
by Simone Sturniolo

Copyright 2013-2018 Science and Technology Facilities Council
This software is distributed under the terms of the GNU General Public
License (GNU GPL)
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""Main script"""

import sys
from cconv import utils

__version__ = "2.0"

__intromsg__ = """
CASTEPconv v. {version}
by Simone Sturniolo
Copyright 2014 Science and Technology Facilities Council

=======

""".format(version=__version__)

if __name__ == '__main__':

    utils.check_pyversion()

    print(__intromsg__)

    seedname, cmdline_task = utils.parse_cmd_args()

    print('Reading {0}.conv'.format(seedname))

    