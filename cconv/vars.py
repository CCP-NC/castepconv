# This file holds all the possible convergence variables (on X and Y axes) and the methods necessary to have them work

import math
import types
from units import cconv_units, cconv_def_units

def wrap_range_calc(rfunc):

    def inner(self):
        if (self.vars['max'] < self.vars['min'] or self.vars['min'] <= 0 or self.vars['step'] <= 0):
            raise CconvVarError('Incompatible boundary values defined for ' + self.name)
        
        self.n = int(math.ceil(self.vars['max']-self.vars['min'])/self.vars['step'])+1
        self.range = rfunc(self)

    return inner

@wrap_range_calc
def range_float(self):
    return [self.vars['min'] + i * self.vars['step'] for i in range(0, self.n)]

@wrap_range_calc
def range_kpn(self):
    try:
        return [tuple([int((self.vars['min'] + i * self.vars['step']) * e) for e in self.kpn_base]) for i in range(0, self.n)]
    except AttributeError:
        return [(0, 0, 0)]

class CconvVarError(Exception):
    pass

class cconv_var(object):
    
    def __init__(self, symbol, name, Name, vtype, dim, active):

        args = [symbol, name, Name, vtype, dim, active]
        typemask = [str, str, str, type, str, bool] # Check that proper types were used

        if (any([type(args[i]) != t for i,t in enumerate(typemask)])):
            raise TypeError()
        if vtype not in (int, float):
            raise ValueError()
        if dim not in cconv_units:
            raise ValueError()

        super(cconv_var, self).__init__()
        self.sym  = symbol
        self.name = name
        self.Name = Name
        self.vtype = vtype
        self.dim  = dim
        self.unit = cconv_def_units[dim]
        self.active = active
        self.fullName = self.Name + (' (' + self.unit + ')' if (dim != 'adim') else '')

class cconv_xvar(cconv_var):
    
    def __init__(self, symbol, name, Name, vtype, dim, active):
        super(cconv_xvar, self).__init__(symbol, name, Name, vtype, dim, active)

        self.vars = {
            'min': None,
            'max': None,
            'step': None,
        }
        
        self.cross_xvars = {}    # The values of the OTHER xvars when converging on this one, if specified


class cconv_yvar(cconv_var):
    
    def __init__(self, symbol, name, Name, vtype, dim, active):
        super(cconv_yvar, self).__init__(symbol, name, Name, vtype, dim, active)

        self.vars = {
            'delta': None
        }

xvar_list = {
        'cut': cconv_xvar('cut', 'cutoff', 'Cut-off', float, 'E', True),
        'kpn': cconv_xvar('kpn', 'kpoint', 'k-points', int, 'adim', True),
        'fgm': cconv_xvar('fgm', 'fine_gmax', 'Fine Gmax', float, 'G', False),
    }

# Set the default values
xvar_list['cut'].vars = {'min': 400.0, 'max': 800.0, 'step': 100.0}
xvar_list['cut'].range_calc = types.MethodType(range_float, xvar_list['cut'])
xvar_list['kpn'].vars = {'min': 1, 'max': 4, 'step': 1}
xvar_list['kpn'].range_calc = types.MethodType(range_kpn, xvar_list['kpn'])
xvar_list['fgm'].vars = {'min': None, 'max':  None, 'step': 100.0}
xvar_list['fgm'].range_calc = types.MethodType(range_float, xvar_list['fgm'])

yvar_list = {
        'nrg': cconv_yvar('nrg', 'energy', 'Final energy', float, 'E', True),
        'for': cconv_yvar('for', 'force', 'Total force', float, 'F', True),
        'str': cconv_yvar('str', 'stress', 'Total stress', float, 'S', False),
    }

yvar_list['nrg'].vars['delta'] = 0.1000E-04    # eV/atom
yvar_list['for'].vars['delta'] = 0.5000E-01    # eV/angstrom
yvar_list['str'].vars['delta'] = 0.1           # GPa

if __name__ == "__main__":

    xvar_list['cut'].range_calc()
    print xvar_list['cut'].range
    xvar_list['kpn'].kpn_base = (1, 1.4, 2.0)
    xvar_list['kpn'].range_calc()
    print xvar_list['kpn'].range
