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

"""Utility functions"""

import sys
import traceback


def check_pyversion():
    """Check Python version"""
    if (sys.version_info[0] < 2 or (sys.version_info[0] == 2 and
                                    sys.version_info[1] < 6)):
        raise RuntimeError("Python version 2.6 or higher required "
                           "to run the script")


def safe_input(msg):
    try:
        return raw_input(msg)
    except NameError:
        return input(msg)


def parse_cmd_args():
    """Parse command line arguments"""

    # Try importing argparse - if the Python version is too old, use optparse
    try:
        import argparse as ap
        parser = ap.ArgumentParser()
        parser.add_argument("seedname", type=str, default=None,
                            help=("Seedname of the convergence job to run - " +
                                  "must be the name of the .cell file before" +
                                  " the extension."))
        parser.add_argument('-t', '--task', type=str, default=None,
                            dest="task", choices=['c', 'i', 'ir', 'o', 'a'],
                            help=("Task to run - can be c (clear), i (input),"
                                  " ir (inputrun), o (output) or a (all). "
                                  "Overrides the .conv file value"),
                            )
        args = parser.parse_args()
        seed = args.seedname
        task = args.task
    except ImportError:
        import optparse as op
        parser = op.OptionParser()
        parser.add_option('-t', '--task', default=None,
                          dest="task", choices=['c', 'i', 'ir', 'o', 'a'],
                          help=("Task to run - can be c (clear), i (input),"
                                " ir (inputrun), o (output) or a (all)."
                                " Overrides the .conv file value"),
                          )
        (options, args) = parser.parse_args()
        if len(args) != 1:
            raise RuntimeError('Error: too few arguments')
        else:
            seed = args[0]
        task = options.task

    # Interpret the task
    task = {'i': 'input', 'c': 'clear', 'ir': 'inputrun', 'o': 'output',
            'a': 'all'}[task]

    return seed, task       

def warn(msg):
    """Print a highlighted warning"""
    print('\033[93mWARNING\033[0m : ' + msg)
