# -*- coding: utf-8 -*-

# Here are stored the unit convertion functions for the various quantities, classified by dimension

import math

# Physical constants

__m_e__    = 9.10938291E-31    # Electron mass, kg
__hbar__   = 1.054571726E-34   # Reduced Planck constant, J*s
__eV__     = 1.602176565E-19   # electronVolt, J
__Ang__    = 1E-10             # Angstrom, m

def uconv(a):    # A decorating function to create conversion functions
    return (lambda x: x*a)

def cut_to_k(E): return math.sqrt(2*__m_e__*E*__eV__)/__hbar__*__Ang__
def k_to_cut(k): return __hbar__**2*(k/__Ang__)**2/(2.0*__m_e__)/__eV__

cconv_def_units = {
    'adim': '',
    'L': 'ang',
    'E': 'eV',
    'F': 'eV/ang',
    'S': 'GPa',
    'Sh': 'ppm',
    'G': '1/ang',
    }

cconv_units = {
    'adim': {    # Adimensional, no units
            '': uconv(1.0),
        },
    'L':    {    # Length units, default is angstrom
            'ang': uconv(1.0), 
            'm': uconv(1.0e10), 
            'cm': uconv(1.0e8), 
            'nm': uconv(10.0), 
            'bohr': uconv(0.529), 
            'a0': uconv(0.529),
        },
    'E':    {    # Energy units, default is electronVolt
            'eV': uconv(1.0),
            'J': uconv(6.24150934e18),
            'Ry': uconv(13.6057),
        },
    'F':    {    # Force units, default is electronVolt/Angstrom
            'eV/ang': uconv(1.0),
            'N':      uconv(6.24150934e8),
        },
    'S':    {    # Stress units, default is gigaPascal
            'GPa': uconv(1.0),
            'Pa':  uconv(1e9),
        },
    'Sh':    {    # Shielding units, default is ppm
            'ppm': uconv(1.0),
        },
    'G':    {    # A reciprocal G vector, default is 1/Angstrom
            '1/ang': uconv(1.0),
            'eV':    cut_to_k,
        }
    }