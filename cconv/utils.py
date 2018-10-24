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


_cut2k = 0.512316724383056

def cut_to_k(E): return _cut2k*E**0.5

def k_to_cut(k): return (k/_cut2k)**2