"""Utility functions"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import math


def floatrange(vmin, vmax, vstep):
    # A range operation with floating points numbers
    n = int(math.ceil((vmax-vmin)/vstep))+1
    return [vmin + i * vstep for i in range(0, n)]
