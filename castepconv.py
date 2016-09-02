#!/usr/bin/env python

# CASTEP convergence automation tool
# by Simone Sturniolo
#
# Copyright 2013 Science and Technology Facilities Council
# This software is distributed under the terms of the GNU General Public License (GNU GPL)

import sys, time, math, os, shutil, glob, re, copy
import subprocess as sp

from cconv.graphs import gp_graph, agr_graph
from cconv.io_freeform import io_freeform_file, io_freeform_error

__vers_number__ = "0.9.6"

# Try importing argparse - if the Python version is too old, use optparse

try:
    import argparse as ap
    has_ap = True
except ImportError:
    has_ap = False
    import optparse as op

class ConvError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class JobError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class CellError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class CastepError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class OptError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class PotError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

# Util: ANSI colored WARNING text

__WARNING__ = '\033[93mWARNING\033[0m'

# Length units allowed by CASTEP. Internally used unit is always Angstroms

length_units = {'ang':1.0, 'm':1.0e10, 'cm':1.0e8, 'nm':10.0, 'bohr':0.529, 'a0':0.529}

# Parameters from CONV files - names and values

str_par_names = {
"convergence_task"      : "ctsk",
"running_mode"          : "rmode",
"output_type"           : "outp",
"running_command"       : "rcmd",
"submission_script"     : "subs",
"fine_gmax_mode"        : "fgmmode",
}

str_par_vals = {
"ctsk"    : "input",                         # Can be INPUT, INPUTRUN, OUTPUT or ALL
"rmode"   : "serial",                        # Can be PARALLEL or SERIAL
"outp"    : "gnuplot",                       # Can be GNUPLOT or XMGRACE
"rcmd"    : "castep <seedname> -dryrun",
"subs"    : None,
"fgmmode" : None,                            # Can be MIN or MAX
}

float_par_names = {
"cutoff_min"        : "cutmin",
"cutoff_max"        : "cutmax",
"cutoff_step"       : "cutstep",
"fine_gmax_min"     : "fgmmin",
"fine_gmax_max"     : "fgmmax",
"fine_gmax_step"    : "fgmstep",
"displace_atoms"    : "displ",
"final_energy_delta": "nrgtol",
"forces_delta"      : "fortol",
"stresses_delta"    : "strtol"
}

float_par_vals = {
"cutmin" : 400.0,           # eV
"cutmax" : 800.0,           # eV
"cutstep": 100.0,           # eV
"fgmmin" : None,            # eV
"fgmmax" : None,            # eV
"fgmstep": 100.0,           # eV
"displ"  : 0.0,             # Ang
"nrgtol" : 0.1000E-04,      # eV/atom
"fortol" : 0.5000E-01,      # eV/Ang
"strtol" : 0.1              # GPa
}

int_par_names = {
"kpoint_n_min"      : "kpnmin",
"kpoint_n_max"      : "kpnmax",
"kpoint_n_step"     : "kpnstep",
"max_parallel_jobs" : "maxjobs"
}

int_par_vals = {
"kpnmin"  : 1,
"kpnmax"  : 4,
"kpnstep" : 1,
"maxjobs" : 0
}

bool_par_names = {
"converge_stress"   : "cnvstr",
"reuse_calcs"       : "rcalc",
"serial_reuse"      : "sruse"
}

bool_par_vals = {
"cnvstr" : False,
"rcalc"  : False,
"sruse"  : True
}

# Physical constants

__m_e__    = 9.10938291E-31    # Electron mass, kg
__hbar__   = 1.054571726E-34   # Reduced Planck constant, J*s
__eV__     = 1.602176565E-19   # electronVolt, J
__Ang__    = 1E-10             # Angstrom, m

# CASTEP constants

__gscale__ = 1.75              # Default value for grid_scale

# CELL and PARAM files

cellfile_lines = None
paramfile_lines = None

cutrange = None
kpnrange = None
fgmrange = None
allrange = None

abc_len = None
kpn_base = (1, 1, 1)
pseudo_pots = None
has_fix_occ = False
sscript = None

ovwrite_files = False

__CASTEP_HEADER__       = "+-------------------------------------------------+"
__CASTEP_TIME__         = "Total time          ="
__CASTEP_ATOMN__        = "Total number of ions in cell = "
__CASTEP_CUTOFF__       = "plane wave basis set cut-off                   :"
__CASTEP_KPOINTS__      = "MP grid size for SCF calculation is"
__CASTEP_FINEGMAX__     = "size of   fine   gmax                          :"
__CASTEP_ENERGY__       = "Final energy, E             ="
__CASTEP_ENERGY_FIX__   = "Final energy ="
__CASTEP_FORCES__       = "***************** Symmetrised Forces *****************"
__CASTEP_FORCES_ALT__   = "*********************** Forces ***********************"
__CASTEP_FORCES_END__   = "*                                                    *"
__CASTEP_STRESSES__     = "*********** Symmetrised Stress Tensor ***********"
__CASTEP_STRESSES_ALT__ = "***************** Stress Tensor *****************"
__CASTEP_STRESSES_END__ = "*                                               *"

# Just a useful snippet to return a kpoint grid from a 3-ple

def kgrid(t, sep=' '): return str(t[0]) + sep + str(t[1]) + sep + str(t[2])

# Another convenient snippet - this time to build a 'jobname'

def jname(s, c, k, f=None): return s + "_cut_" + str(c) + "_kpn_" + str(min(k)) + ("_fgm_" + str(f) if f is not None else "")

# Yet another snippet - just a decimal rounding utility

def round_digits(x, n=0): return int(x*10**n+0.5)/10.0**n

# Useful conversion utilities - go from a cutoff energy to a wave vector and back

def cut_to_k(E): return math.sqrt(2*__m_e__*E*__eV__)/__hbar__*__Ang__
def k_to_cut(k): return __hbar__**2*(k/__Ang__)**2/(2.0*__m_e__)/__eV__

# Another utility - find the LAST (rather than the first) occurrence that contains an element in a list

def rindex_cont(my_list, my_el):

    for i in range(len(my_list), 0, -1):
        if my_el in my_list[i-1]:
            return i-1
    return None

# Parse .conv file for convergence calculation options

def parse_convfile(cfile):

    global par_names
    global par_vals

    cfile.seek(0)
    clines = cfile.readlines()

    for i, l in enumerate(clines):
        # Skip comments
        if ('#' in l):
            l = l[:l.index('#')]

        #Parse options
        cline = l.strip()
        if len(cline) == 0:
            #Skip empty lines
            continue
        else:
            cline = cline.split(':', 1)
        if (len(cline) < 2):
            raise ConvError("Bad formatting in .conv file at line " + str(i))
        par_name = cline[0].strip().lower()
        if (par_name in str_par_names):
            if str_par_names[par_name] not in ('rcmd', 'subs'): # If it IS rcmd or subs, we don't need to alter it
                cline[1] = cline[1].lower()
            else:
                # A basic sanity check
                if '<seedname>' not in cline[1]:
                    print(__WARNING__ + ": {0} does not contain a <seedname> tag.\n"
                          "This is likely erroneous. Please check.\n".format(par_name))
            str_par_vals[str_par_names[par_name]] = cline[1].strip()
        elif (par_name in float_par_names):
            float_par_vals[float_par_names[par_name]] = float(cline[1])
        elif (par_name in int_par_names):
            int_par_vals[int_par_names[par_name]] = int(cline[1])
        elif (par_name in bool_par_names):
            bool_par_vals[bool_par_names[par_name]] = (cline[1].lower().strip() == "true")
        else:
            raise ConvError("Unrecognized option in .conv file at line " + str(i))


    fgmax_validate()

# Parse .conv_tab file to generate cutoff and k points ranges

def parse_conv_tab_file(cfile):

    cutrange = []
    kpnrange = []
    fgmrange = []

    cfile.seek(0)
    clines = cfile.readlines()

    if len(clines) == 3:
        try:
            cutline = clines[0].split(':')[1].split('eV')
            kpnline = clines[1].split(':')[1].split('|')
            fgmline = clines[2].split(':')[1].split('eV')
            if (cutline is not None and kpnline is not None):
                cutrange = [float(cut.strip()) for cut in cutline[:-1]]
                kpnrange = [tuple([int(k) for k in kpn.strip().split()]) for kpn in kpnline[:-1]]
            if (fgmline[0].strip() != 'N/A'):
                fgmrange = [float(fgm.strip()) for fgm in fgmline[:-1]]
            else:
                fgmrange = [None]
        except IndexError:
            return [], [], []
    else:
        return [], [], []

    return cutrange, kpnrange, fgmrange

# Strip .cell file from unnecessary lines to get only the ones we care for (i.e. remove all reference to kpoints, we're going to put those in ourselves)
# Also read cell parameters and construct the proper kpn_base

def strip_cellfile(fname):

    stripped = io_freeform_file(fname)

    # Parse/remove all undesired blocks

    abc   = None
    u     = 1.0
    kbase = [1.0, 1.0, 1.0]
    ppot  = None

    hkl_data = []

    # Remove everything having to do with kpoints
    
    for k in ("kpoints_mp_grid", "kpoint_mp_grid", "kpoints_mp_spacing", "kpoint_mp_spacing", "kpoints_list", "kpoint_list"):
        stripped.freeform_remove(k)

    # Read lattice parameters
    if stripped.freeform_present("lattice_cart"):
        abc_block = stripped.freeform_block("lattice_cart")
        for l in abc_block:
            l_split = l.lower().strip().split()
            if l_split[0] in length_units:
                u = length_units[l_split[0]]
            else:
                try:
                    hkl_data.append([u*float(x) for x in l_split])
                except ValueError:
                    raise CellError('Bad formatting in .cell file LATTICE_CART block')
                if (len(hkl_data) == 3):
                    abc = tuple([math.sqrt(sum([x**2.0 for x in hkl])) for hkl in hkl_data])
                if (len(hkl_data) > 3):
                    raise CellError('Bad formatting in .cell file LATTICE_CART block')
    
    if stripped.freeform_present("lattice_abc"):
        if abc is not None:
            raise CellError('Duplicated LATTICE_* block in .cell file')
        abc_block = stripped.freeform_block("lattice_abc")
        for l in abc_block:
            l_split = l.lower().strip().split()
            if l_split[0] in length_units:
                u = length_units[l_split[0]]
            else:
                try:
                    hkl_data.append([u*float(x) for x in l_split])
                except ValueError:
                    raise CellError('Bad formatting in .cell file LATTICE_ABC block')
                if (len(hkl_data) == 1):
                    abc = tuple(hkl_data[0])
                if (len(hkl_data) > 2):
                    raise CellError('Bad formatting in .cell file LATTICE_ABC block')

    if abc is None:
        raise CellError('.cell file does not contain a LATTICE_* block')
    
    # Read pseudopotentials
    if stripped.freeform_present("species_pot"):
        ppot = []
        pot_block = stripped.freeform_block("species_pot")
        for l in pot_block:
            l_split = l.strip().split()
            try:
                ppot.append([l_split[0], l_split[1]])
            except IndexError:
                raise CellError('Bad formatting in .cell file SPECIES_POT block')

    if len(hkl_data) == 2:
        for i in range(0, 3):
            kbase[i] = abs(hkl_data[0][i-1]*hkl_data[0][i-2]*math.sin(hkl_data[1][i]))
    elif len(hkl_data) == 3:
        for i in range(0, 3):
            kbase[i] = math.sqrt(sum([(hkl_data[i-2][j-2]*hkl_data[i-1][j-1] - hkl_data[i-2][j-1]*hkl_data[i-1][j-2])**2.0 for j in range(0, 3)]))

    kbase = tuple([k/min(kbase) for k in kbase])

    return stripped, abc, kbase, ppot

# Strip .param file from unnecessary lines to get only the ones we care for (i.e. remove all reference to task and cutoff)

def strip_paramfile(fname):
    
    global bool_par_vals, has_fix_occ

    stripped = io_freeform_file(fname)

    # Is occupancy fixed?
    if stripped.freeform_present('fix_occupancy'):
        has_fix_occ = stripped.freeform_boolean('fix_occupancy')

    # Remove keywords that need to be stripped
    for k in ("cut_off_energy", "cutoff_energy"):
        stripped.freeform_remove(k)

    if bool_par_vals['cnvstr']:
        stripped.freeform_remove('calculate_stress')
    if bool_par_vals['sruse']:
        stripped.freeform_remove('continuation')
        stripped.freeform_remove('reuse')
    if str_par_vals['fgmmode'] is not None:
        stripped.freeform_remove('fine_gmax')

    return stripped

# Displace atoms by a small amount in cellfile if needed to have nonzero forces

def displace_cell_atoms(cfile, abc, d):

    u = None

    pos_is_abs = cfile.freeform_present('positions_abs')    # Are positions absolute or fractionary?
    if pos_is_abs:
        pos_block = cfile.freeform_block('positions_abs')
    else:
        pos_block = cfile.freeform_block('positions_frac')

    pos_block_displ = []

    for i, l in enumerate(pos_block):
        l_split = l.strip().split()
        if pos_is_abs and l_split[0].lower() in length_units:
            u = length_units[l_split[0].lower()]
            pos_block_displ.append('ang\n')
        else:
            if len(l_split) != 4:
                raise CellError('Bad formatting in .cell file POSITION_* block')
            try:
                xyz = [float(x) for x in l_split[1:]]
                # The displacement is in alternated directions, since otherwise you'd still get an equilibrium structure
                # Won't work for single atom unit cells. This needs fixing. It might also break if there are symmetries in the system.
                fac = (i%2)*2-1
                if pos_is_abs:
                    #l = l_split[0] + '\t' + str(u*xyz[0]+fac*d) + '\t' + str(u*xyz[1]+fac*d) + '\t' + str(u*xyz[2]+fac*d) + '\n'
                    l = '{0} {1} {2} {3}'.format(*[l_split[0]] + [u*x+fac*d for x in xyz])
                else:
                    #l = l_split[0] + '\t' + str(xyz[0]+fac*d/abc[0]) + '\t' + str(xyz[1]+fac*d/abc[1]) + '\t' + str(xyz[2]+fac*d/abc[2]) + '\n'
                    l = '{0} {1} {2} {3}'.format(*[l_split[0]] + [x+fac*d/abc[j] for j, x in enumerate(xyz)])
            except ValueError:
                raise CellError('Bad formatting in .cell file POSITION_* block')
            pos_block_displ.append(l)
    
    if pos_is_abs:
        cfile.freeform_block('positions_abs', pos_block_displ)
    else:
        cfile.freeform_block('positions_frac', pos_block_displ)

    return cfile

# Waits for a running job until it's finished

def jobfinish_wait(foldname, jobname):

    is_finished = False
    while not is_finished:
        try:
            is_finished = jobfinish_check(foldname, jobname)
        except JobError:
            raise
        time.sleep(1)

# Checks whether a job is over or not

def jobfinish_check(foldname, jobname):

    if not os.path.isfile(os.path.join(foldname, jobname + '.castep')):
        return False
    # In case of errors, abort
    if os.path.isfile(os.path.join(foldname, jobname + '.0001.err')):
        raise JobError("Job " + jobname + " failed")
    cast_file = open(os.path.join(foldname, jobname + '.castep'), 'r')
    cast_lines = cast_file.readlines()
    cast_lines.reverse()
    cast_file.close()
    if cast_lines is None:
        return False
    for l in cast_lines:
        if __CASTEP_HEADER__ in l:
            return False
        elif __CASTEP_TIME__ in l:
            return True

# Read forces and find the maximum in .castep file

def parse_forces(cfile):

    global __CASTEP_FORCES_END__

    max_for = 0.0

    for l in cfile[6:]:

        if __CASTEP_FORCES_END__ in l:
            return max_for

        try:
            cur_for = math.sqrt(sum([float(x)**2.0 for x in l.split()[3:6]]))
        except ValueError:
            pass
        if cur_for > max_for:
            max_for = cur_for

    raise CastepError("Corrupted forces block")

# Read stresses and find the maximum in .castep file

def parse_stresses(cfile):

    global __CASTEP_STRESSES_END__

    max_stress = 0.0

    for l in cfile[6:]:

        if __CASTEP_STRESSES_END__ in l:
            return max_stress

        try:
            cur_stress = math.sqrt(sum([float(x)**2.0 for x in l.split()[2:5]]))
        except ValueError:
            pass

        if cur_stress > max_stress:
            max_stress = cur_stress

    raise CastepError("Corrupted stresses block")

# Parse a .castep file as a whole

def parse_castep_file(cfile, filepath):

    global bool_par_vals, has_fix_occ
    global __CASTEP_ATOMN__, __CASTEP_CUTOFF__, __CASTEP_ENERGY_FIX__, __CASTEP_ENERGY__, __CASTEP_FINEGMAX__, __CASTEP_FORCES_ALT__
    global __CASTEP_FORCES_END__, __CASTEP_FORCES__, __CASTEP_HEADER__, __CASTEP_KPOINTS__, __CASTEP_STRESSES_ALT__
    global __CASTEP_STRESSES_ALT__, __CASTEP_STRESSES_END__, __CASTEP_STRESSES__

    atom_n = None
    i_nrg = None
    i_for = None
    i_str = None
    cut_check = None
    kpn_check = None
    fgm_check = None

    calc_str = bool_par_vals["cnvstr"]

    start_l = rindex_cont(castepfile, __CASTEP_HEADER__)

    for i, l in enumerate(castepfile[start_l:]):

        try:
            if __CASTEP_ATOMN__ in l and atom_n is None:
                atom_n = int(l.split()[7])
            elif __CASTEP_CUTOFF__ in l:
                cut_check = float(l.split()[6])
            elif __CASTEP_KPOINTS__ in l:
                kpn_check = tuple([int(x) for x in l.split()[7:]])
            elif __CASTEP_FINEGMAX__ in l:
                fgm_check = float(l.split()[5])
            elif __CASTEP_ENERGY__ in l and not has_fix_occ:
                i_nrg = float(l.split()[4])
            elif __CASTEP_ENERGY_FIX__ in l and has_fix_occ:
                i_nrg = float(l.split()[3])
            elif __CASTEP_FORCES__ in l or __CASTEP_FORCES_ALT__ in l:
                i_for = parse_forces(castepfile[start_l+i:])
            elif calc_str and (__CASTEP_STRESSES__ in l or __CASTEP_STRESSES_ALT__ in l):
                i_str = parse_stresses(castepfile[start_l+i:])
        except ValueError:
            raise CastepError("Corrupted " + filepath + " file detected")
        except CastepError as CE:
            raise CE

        if atom_n is not None and i_nrg is not None and i_for is not None and cut_check is not None and kpn_check is not None and fgm_check is not None \
        and ((calc_str and i_str is not None) or not calc_str):
            break

    return {
        'atom_n': atom_n,
        'nrg':    i_nrg,
        'for':    i_for,
        'str':    i_str,
        'cut':    cut_check,
        'kpn':    kpn_check,
        'fgm':    fgm_check,
        }

# Find reasonable estimates for convergence values

def conv_estimates(seedname, data, cnvstr=False):

    global float_par_vals

    ordered_x = ['cut', 'kpn', 'fgm']    # Just for the sake of printing them in the right order
    ordered_y = ['nrg', 'for', 'str']    # Same as above
    data_x = {'cut': ['cutoff', 'eV'], 'kpn': ['k-point grid', 'points'], 'fgm': ['fine Gmax', 'eV']}
    data_y = {'nrg': ['total energy', 'eV per atom'], 'for': ['maximum force', 'eV/Ang'], 'str': ['maximum stress', 'GPa']}

    opt_vals = {}
    
    out_string = ""
    
    
    out_string += "Convergence results:\n"

    for x in ordered_x:

        x_name = data_x[x][0]
        x_unit = data_x[x][1]

        # Skip if there is no X convergence performed
        if None in data[x]['range']:
            continue

        opt_vals[x] = {}
        
        for y in ordered_y:

            if y == 'str' and not cnvstr:
                continue
                
            opt_vals[x][y] = None

            dataset = data[x][y]
            y_name = data_y[y][0]
            y_unit = data_y[y][1]
            tol = float_par_vals[y+"tol"]

            # Cases where we can't say anything meaningful
            if len(dataset) < 2:
                # 1. If there's only one point
                out_string += "Impossible to give a convergence estimate with a single " + x_name + " point\n"
            elif tol <= 0.0:
                # 2. If the tolerance value is unphysical
                out_string +=  "Impossible to give a convergence estimate with a null or negative value for " + y_name + " tolerance\n"
            elif all([abs(dpoint) < tol for dpoint in dataset]):
                # 3. The values are so small (in absolute sense) that there is nothing to converge (e.g. forces with an equilibrated system)
                out_string +=  "Impossible to give a convergence estimate with " + y_name + " values lower than the tolerance \n"
            else:

                delta = dataset[0]

                for i, val in enumerate(dataset[1:]):

                    delta = abs(val - delta)
                    if delta < tol:
                        out_string += "Based on converging " + y_name + " to " + str(tol) + " " + y_unit + ", minimum " + x_name + " suggested is " + str(data[x]['rangestr'][i]) + " " + x_unit
                        if (x == 'fgm'):
                            out_string +=  " (corresponding to a value of " + str(round_digits(cut_to_k(data[x]['range'][i]), 4)) + " 1/Ang)\n"
                        else:
                            out_string += "\n"
                        
                        opt_vals[x][y] = data[x]['range'][i]
                        break

                    delta = val

                    if i == len(dataset[1:])-1:
                        out_string +=  "Unable to converge " + x_name + " with " + y_name + " within given range.  Try increasing maximum " + x_name + " for your search\n"

    out_string += "\n"

    # Generate a global estimate

    max_vals = {}
    max_str  = {}

    for x in ordered_x:

        if x not in opt_vals:
            continue

        max_vals[x] = -1
        max_str[x] = ''

        for y in opt_vals[x]:

            if opt_vals[x][y] > max_vals[x]:
                max_vals[x] = opt_vals[x][y]
                max_str[x] = data[x]['rangestr'][data[x]['range'].index(opt_vals[x][y])]
            elif opt_vals[x][y] is None:
                max_str[x] = "N/A"
                break
        
        out_string += "Overall suggested " + data_x[x][0] + " is " + max_str[x] + " " + data_x[x][1] + "\n"
    
    sys.stdout.write(out_string)

    # Also save as report file

    optfile = open(seedname + "_report.txt", 'w')
    optfile.write(out_string)
    optfile.close()

    return opt_vals


# Parse command line arguments

def parse_cmd_args():

    global has_ap

    if has_ap:
        parser = ap.ArgumentParser()
        parser.add_argument("seedname", type=str, default=None, help="Seedname of the convergence job to run - must be the name of the .cell file before the extension.")
        parser.add_argument('-t', '--task', type=str, default=None, help="Task to run - can be c (clear), i (input), ir (inputrun), o (output) or a (all). Overrides the .conv file value", \
        dest="task", choices=['c', 'i', 'ir', 'o', 'a'])
        args = parser.parse_args()
        return args.seedname, args.task
    else:
        parser = op.OptionParser()
        parser.add_option('-t', '--task', default=None, help="Task to run - can be c (clear), i (input), ir (inputrun), o (output) or a (all). Overrides the .conv file value", \
        dest="task", choices=['c', 'i', 'ir', 'o', 'a'])
        (options, args) = parser.parse_args()
        if len(args) != 1:
            seedname = None
        else:
            seedname = args[0]
        return seedname, options.task

# Write a folder with given cutoff and kpoints

def create_conv_folder(foldname, jobname, cut, kpn, fgm, prev_jobname=None):

    global ovwrite_files, stripped_cell, stripped_param, bool_par_vals, str_par_vals, pseudo_pots, sscript

    if not os.path.exists(foldname):
        print "Creating folder " + foldname
        os.makedirs(foldname)
    elif not ovwrite_files and len(os.listdir(foldname)) > 0 and prev_jobname is None:
        to_del = raw_input(__WARNING__ + ": folder " + foldname + " already exists. \n"
                           "Some files might be overwritten. Continue (y/N/y-all)?")
        if to_del.lower() == 'y-all':
            ovwrite_files = True
        elif to_del.lower() != 'y':
            sys.exit("Aborting")

    print "Creating files for job " + jobname
    icell = open(os.path.join(foldname, jobname + '.cell'), 'w')
    iparam = open(os.path.join(foldname, jobname + '.param'), 'w')

    if sscript is not None:

        print "Copying submission script"

        iscript = open(os.path.join(foldname, jobname), 'w')
        for l in sscript:
            mod_l = l.replace('<seedname>', jobname)
            iscript.write(mod_l)
        iscript.flush()
        iscript.close()

    icell_fform = copy.copy(stripped_cell)
    # Write pseudopotentials block
    if pseudo_pots is not None:
        icell_fform.freeform_block("species_pot", ['{0}\t{1}'.format(*p) for p in pseudo_pots])
    # Write kpoint grid
    icell_fform.freeform_integer_vector("kpoints_mp_grid", kpn)
    
    icell.write(icell_fform.freeform_print())
    icell.close()

    # Write param file

    iparam_fform = copy.copy(stripped_param)
    iparam_fform.freeform_string('task', 'SinglePoint')
    iparam_fform.freeform_physical('cut_off_energy', 'E', cut, 'ev')
    if prev_jobname is not None and bool_par_vals["sruse"]:
        iparam_fform.freeform_string('reuse', prev_jobname + '.check')
    if bool_par_vals["cnvstr"]:
        iparam_fform.freeform_boolean('calculate_stress', True)
    if fgm is not None:
        iparam_fform.freeform_real('fine_gmax', cut_to_k(fgm))
    iparam.write(iparam_fform.freeform_print())
    iparam.close()

# Run a job, eventually wait for its completion, etc.

def job_run(foldname, jobname, running_jobs=None):

    global bool_par_vals, int_par_vals

    if bool_par_vals["rcalc"]:
        try:
            if jobfinish_check(foldname, jobname):
                return
        except JobError:
            pass

    os.chdir(foldname)

    if os.path.isfile(jobname + ".castep") or os.path.isfile(jobname + ".check") or os.path.isfile(jobname + ".0001.err"):
        print "Removing output files from previous jobs for " + jobname
        if os.path.isfile(jobname + ".castep"):
            os.remove(jobname + ".castep")
        if os.path.isfile(jobname + ".check"):
            os.remove(jobname + ".check")
        if os.path.isfile(jobname + ".0001.err"):
            os.remove(jobname + ".0001.err")
    
    cmd_line, stdin_file, stdout_file = compile_cmd_line(jobname)

    if (stdin_file is not None):
        if os.path.isfile(stdin_file):
            stdin_file = open(stdin_file, 'r')
        else:
            sys.exit("ERROR - STDIN redirected from inexistent file " + stdin_file)
    if (stdout_file is not None):
        stdout_file = open(stdout_file, 'w')

    try:
        print "Running job " + jobname
        sp.Popen(cmd_line, stdin=stdin_file, stdout=stdout_file)
    except OSError:
        sys.exit("ERROR - Command:\n>\t" + cmd_line[0] + "\ndoes not exist on this system")
    os.chdir("..")
    if (running_jobs is not None):
        running_jobs.append([foldname, jobname])

    if(stdin_file is not None):
        stdin_file.close()
    if(stdout_file is not None):
        stdout_file.close()

# Wait for multiple jobs to finish, respecting the maxjobs roof. If wait_all is True, wait for all of them to complete

def multijob_wait(running_jobs=None, wait_all=False):

    global int_par_vals

    if (running_jobs is not None):
        # PARALLEL operation case
        if int_par_vals["maxjobs"] > 0 or wait_all:
            while (not wait_all and len(running_jobs) >= int_par_vals["maxjobs"]) or (wait_all and (len(running_jobs) > 0)):
                for job_ind in range(len(running_jobs)-1, -1, -1):
                    try:
                        if jobfinish_check(running_jobs[job_ind][0], running_jobs[job_ind][1]):
                            del running_jobs[job_ind]
                    except JobError as JE:
                        del running_jobs[job_ind]
                        print __WARNING__ + " - " + str(JE)
    else:
        # SERIAL operation case
        try:
            jobfinish_wait(foldname, jobname)
        except JobError as JE:
            print __WARNING__ + " - " + str(JE)

# Find and if needed copy the pseudo potential files

def find_pseudopots(seedname, pseudo_pots):

    global ovwrite_files

    # First: check whether there are pseudo pots at all

    if (pseudo_pots is None) or (len(pseudo_pots) == 0):
        return

    pp_folder = seedname + "_pspot"

    # Second: check whether a system path for pseudopotentials is indicated

    try:
        default_path = os.environ['PSPOT_DIR']
    except KeyError:
        default_path = None

    # Third: check whether the pseudopotentials do indeed exist and update their paths

    pspotline_regex = re.compile("[0-9\.]+\|[0-9\.]+\|[0-9\.]+\|")

    for i in range(0, len(pseudo_pots)):

        if default_path is None or not os.path.isfile(os.path.join(default_path, pseudo_pots[i][1])):

            # Exclude that pseudopotentials are just pseudopotential LINES

            if len(pspotline_regex.findall(pseudo_pots[i][1])) > 0:
                continue

            if os.path.isfile(pseudo_pots[i][1]):
                if not os.path.exists(pp_folder):
                    os.makedirs(pp_folder)
                elif not ovwrite_files:
                    to_del = raw_input(__WARNING__ + ": folder " + pp_folder + " already exists.\n"
                                       "Some files might be overwritten. Continue (y/N/y-all)?")
                    if to_del.lower() == 'y-all':
                        ovwrite_files = True
                    elif to_del.lower() != 'y':
                        sys.exit("Aborting")
                shutil.copy2(pseudo_pots[i][1], os.path.join(pp_folder, pseudo_pots[i][1]))
                pseudo_pots[i][1] = os.path.join('..', pp_folder, pseudo_pots[i][1])
            else:
                raise PotError("Pseudo potential file " + pseudo_pots[i][1] + " could not be found")

# Fill in the fine gmax parameters after parsing the .conv file

def fgmax_validate():

    global str_par_vals, float_par_vals

    fgmmode = str_par_vals['fgmmode']

    if fgmmode is None:
        pass
    elif fgmmode == 'min_cutoff':
        if float_par_vals['fgmmin'] is None:
            float_par_vals['fgmmin'] = __gscale__**2*float_par_vals['cutmin']
        elif float_par_vals['fgmmin'] < __gscale__**2*float_par_vals['cutmin']:
            raise ConvError("fine_gmax_min must be greater or equal than " + str(__gscale__**2) + "*cutoff_min for fine_gmax_mode = MIN")
        if float_par_vals['fgmmax'] is None:
            float_par_vals['fgmmax'] = float_par_vals['fgmmin']+3.0*float_par_vals['fgmstep']
        elif float_par_vals['fgmmax'] <= __gscale__**2*float_par_vals['cutmin']:
            raise ConvError("fine_gmax_max must be greater than " + str(__gscale__**2) + "*cutoff_min for fine_gmax_mode = MIN")
        elif float_par_vals['fgmmax'] <= float_par_vals['fgmmin']:
            raise ConvError("Invalid fine Gmax range defined in .conv file")
    elif fgmmode ==  'max_cutoff':
        if float_par_vals['fgmmin'] is None:
            float_par_vals['fgmmin'] = __gscale__**2*float_par_vals['cutmax']
        elif float_par_vals['fgmmin'] < __gscale__**2*float_par_vals['cutmax']:
            raise ConvError("fine_gmax_min must be greater or equal than " + str(__gscale__**2) + "*cutoff_max for fine_gmax_mode = MAX")
        if float_par_vals['fgmmax'] is None:
            float_par_vals['fgmmax'] = float_par_vals['fgmmin']+3.0*float_par_vals['fgmstep']
        elif float_par_vals['fgmmax'] <= __gscale__**2*float_par_vals['cutmax']:
            raise ConvError("fine_gmax_max must be greater than " + str(__gscale__**2) + "*cutoff_max for fine_gmax_mode = MAX")
        elif float_par_vals['fgmmax'] <= float_par_vals['fgmmin']:
            raise ConvError("Invalid fine Gmax range defined in .conv file")
    elif fgmmode == 'none':
        str_par_vals['fgmmode'] = None
    else:
        raise ConvError("Invalid value for fine_gmax_mode parameter in .conv file")

# Interpret the command line pipelining and such in running_command

def compile_cmd_line(jname):

    global str_par_vals

    cmd_line = str_par_vals["rcmd"]

    # Check for redirections

    stdin_file = None
    stdout_file = None

    cmd_line = cmd_line.replace('<seedname>', jname)

    if ('<' in cmd_line):
        # Take the last of the filenames given after a < but before a >
        stdin_file = cmd_line.split('<')[-1].split('>')[0].split()[-1]

    if ('>' in cmd_line):
        # Take the last of the filenames given after a > but before a <
        stdout_file = cmd_line.split('>')[-1].split('<')[0].split()[-1]

    cmd_line = cmd_line.split('<')[0].split('>')[0].split()

    return cmd_line, stdin_file, stdout_file

###### -- MAIN PROGRAM -- ######

if (sys.version_info[0] < 2 or (sys.version_info[0] == 2 and sys.version_info[1] < 6)):
    sys.exit("ERROR - Python version 2.6 or higher required to run the script")

print "CASTEPconv v. " + __vers_number__ + "\n"
print "by Simone Sturniolo"
print "Copyright 2014 Science and Technology Facilities Council\n"

seedname, cmdline_task = parse_cmd_args()

if seedname is None:
    sys.exit("ERROR - <seedname> is a required argument")

# PHASE 0 - Check for existence of all required files and read the necessary information

print "Reading " + seedname + ".conv"

try:
    job_convfile = open(seedname + ".conv", 'r')
except IOError:
    print __WARNING__ + ": - Convergence parameter file for job " + seedname + " not found. Using default parameters"
    job_convfile = None

if job_convfile is not None:
    try:
        parse_convfile(job_convfile)
    except ConvError as CE:
        sys.exit("ERROR - " + str(CE))

#    CELL and PARAM files need to be opened only if the task is not OUTPUT

if (str_par_vals['ctsk'] != "output"):

    print "Reading " + seedname + ".cell"

    try:
        assert os.path.isfile(seedname + ".cell")
    except AssertionError:
        sys.exit("ERROR - .cell file for job " + seedname + " not found")

    print "Reading " + seedname + ".param"

else:
    cellfile_lines = None

# This is needed even for output tasks since the program needs to be aware whether fix_occupancy is used

try:
    stripped_param = strip_paramfile(seedname + ".param")
except io_freeform_error:
    print(__WARNING__ + " - .param file for job " + seedname + " not found")
    stripped_param = io_freeform_file()

# Override task in .conv file with command line options

if cmdline_task is not None:
    str_par_vals['ctsk'] = {'c': 'clear', 'i': 'input', 'ir': 'inputrun', 'o': 'output', 'a': 'all'}[cmdline_task]

try:
    assert(str_par_vals['ctsk'] in ("clear", "input", "inputrun", "output", "all"))
except AssertionError as e:
    sys.exit("ERROR - Invalid convergence_task parameter")

# PHASE 0.5 - CLEAR
# Clear all files and folders from previous jobs

if (str_par_vals['ctsk'] in ("clear")):

    # Create a list of files and folders to delete
    
    to_del_fold = set(glob.glob(seedname + "_cut_*_kpn_*_fgm_*"))
    to_del_fold = to_del_fold.union(set(glob.glob(seedname + "_cut_*_kpn_*")))
    
    if os.path.exists(seedname + "_conv"):
        to_del_fold.add(seedname + "_conv")
    if os.path.exists(seedname + "_pspot"):
        to_del_fold.add(seedname + "_pspot")

    # Create a list of files to delete
    to_del_suffixes = [".conv_tab", "_report.txt",
        "_cut_conv.dat", "_kpn_conv.dat", "_fgm_conv.dat",
        "_cut_conv.gp",  "_kpn_conv.gp", "_fgm_conv.gp",   "_cut_str_conv.gp",  "_kpn_str_conv.gp",  "_fgm_str_conv.gp",
        "_cut_conv.agr", "_kpn_conv.agr", "_fgm_conv.agr", "_cut_str_conv.agr", "_kpn_str_conv.agr", "_fgm_str_conv.agr"]   # A list of all possible files to delete

    to_del_files = []

    for suff in to_del_suffixes:
        if os.path.isfile(seedname + suff):
            to_del_files += [seedname + suff]

    if (len(to_del_fold) > 0):
        print "The following folders will be deleted:"
        for f in to_del_fold:
            sys.stdout.write(f + ' ')
        print ""

    if (len(to_del_files) > 0):
        print "The following files will be deleted:"
        for f in to_del_files:
            sys.stdout.write(f + ' ')
        print ""

    if ((len(to_del_fold) + len(to_del_files)) == 0):
        print "No folders or files to delete found"
    else:
        to_del = raw_input("Continue (y/N)?")
    
        if to_del.lower() == 'y':
    
            for f in to_del_fold:
                shutil.rmtree(f)
            for f in to_del_files:
                os.remove(f)

# PHASE 1 - INPUT
# Produce a batch of folders and files with the various configurations

if (str_par_vals['ctsk'] in ("input", "inputrun", "all")):

    # First create a "stripped" version of cell and param files, to use as a template for the new files to generate

    try:
        stripped_cell, abc_len, kpn_base, pseudo_pots = strip_cellfile(seedname + ".cell")
    except CellError as CE:
        sys.exit("ERROR - " + str(CE))

    # Check the positions and existence of pseudopotentials

    try:
        find_pseudopots(seedname, pseudo_pots)
    except PotError as PE:
        sys.exit("ERROR - " + str(PE))

    # Check the existence of any eventual submission scripts

    if str_par_vals["subs"] is not None:

        if not os.path.isfile(str_par_vals["subs"]):
            print __WARNING__ + ": submission script " + str_par_vals["subs"] + " not found, skipping"
        else:
            sscript = open(str_par_vals["subs"], 'r').readlines()

    # Apply displacements to .cell atoms if needed to have non-zero forces

    if float_par_vals["displ"] > 0.0:
        print "Displacing atoms in .cell file of " + str(float_par_vals["displ"]) + " Angstroms"
        try:
            stripped_cell = displace_cell_atoms(stripped_cell, abc_len, float_par_vals["displ"])
        except CellError as CE:
            sys.exit("ERROR - " + str(CE))

    # If reuse of previous calculations is required, try to open the old .conv_tab file

    if bool_par_vals["rcalc"] and os.path.isfile(seedname + ".conv_tab"):
        print "Reusing results from previous convergence calculations"
        old_cutrange, old_kpnrange, old_fgmrange = parse_conv_tab_file(open(seedname + ".conv_tab", 'r'))
    else:
        old_cutrange = []
        old_kpnrange = []
        old_fgmrange = []

    # Open a .conv_tab file to keep track of the created files and folders. Will be read if output is done as a separate operation

    if os.path.isfile(seedname + ".conv_tab") and not ovwrite_files:
        to_del = raw_input(__WARNING__ + ": " + seedname + ".conv_tab already exists. This file will be overwritten. Continue (y/N/y-all)?")
        if to_del.lower() == 'y-all':
            ovwrite_files = True
        elif to_del.lower() != 'y':
            sys.exit("Aborting")

    conv_tab_file = open(seedname + ".conv_tab", 'w')

    if (float_par_vals["cutmin"] <= 0.0 or float_par_vals["cutstep"] <= 0.0 or float_par_vals["cutmax"] < float_par_vals["cutmin"]):
        sys.exit("ERROR - Invalid cutoff range defined in .conv file")
    if (int_par_vals["kpnmin"] <= 0 or int_par_vals["kpnstep"] <= 0 or int_par_vals["kpnmax"] < int_par_vals["kpnmin"]):
        sys.exit("ERROR - Invalid k-points range defined in .conv file")

    cut_n = int(math.ceil(float_par_vals["cutmax"]-float_par_vals["cutmin"])/float_par_vals["cutstep"])+1
    cutrange = [float_par_vals["cutmin"] + i * float_par_vals["cutstep"] for i in range(0, cut_n)]

    kpn_n = int(math.ceil(int_par_vals["kpnmax"]-int_par_vals["kpnmin"])/int_par_vals["kpnstep"])+1
    kpnrange = [tuple([int((int_par_vals["kpnmin"] + i * int_par_vals["kpnstep"]) * e) for e in kpn_base]) for i in range(0, kpn_n)]

    if (str_par_vals["fgmmode"] is not None):
        fgm_n = int(math.ceil(float_par_vals["fgmmax"]-float_par_vals["fgmmin"])/float_par_vals["fgmstep"])+1
        fgmrange = [float_par_vals["fgmmin"] + i * float_par_vals["fgmstep"] for i in range(0, fgm_n)]
    else:
        fgmrange = [None]
    
    # Build an array with the simulation conditions to employ condensed into one

    allrange = []
    fgm_cutoff = cutrange[-1] if str_par_vals['fgmmode'] == 'max_cutoff' else cutrange[0]
    conv_tab_file.write("cutoff:\t")
    for cut in cutrange:
        conv_tab_file.write(str(cut) + " eV\t")
        allrange.append({'cut': cut, 'kpn': kpnrange[0], 'fgm': fgmrange[0], 'scan': 'cut'})
    conv_tab_file.write("\nkpoint_n:\t" + kgrid(kpnrange[0]) + "\t|\t")
    for kpn in kpnrange[1:]:
        conv_tab_file.write(kgrid(kpn) + "\t|\t")
        allrange.append({'cut': cutrange[0], 'kpn': kpn, 'fgm': fgmrange[0], 'scan': 'kpn'})
    conv_tab_file.write("\nfine_gmax:\t" + (str(fgmrange[0]) + " eV\t" if fgmrange[0] is not None else "N/A"))
    for fgm in fgmrange[1:]:
        conv_tab_file.write(str(fgm) + " eV\t")
        allrange.append({'cut': fgm_cutoff, 'kpn': kpnrange[0], 'fgm': fgm, 'scan': 'fgm'})
    conv_tab_file.close()

    foldname = None
    jobname = None
    prev_jobname = None

    try:
        assert(str_par_vals["rmode"] in ("serial", "parallel"))
    except AssertionError as e:
        sys.exit("ERROR - Invalid value for the running_mode parameter")

    if str_par_vals["rmode"] == "parallel":

        print "Creating folders for parallel convergence run"

    elif str_par_vals["rmode"] == "serial":

        foldname = seedname + "_conv"

        print "Creating folder " + foldname + " for serial convergence run"

        # Check whether the folder exists from previous jobs or you're just creating it now...

        if not os.path.exists(foldname):
            os.makedirs(foldname)
        elif not ovwrite_files:
            to_del = raw_input(__WARNING__ + ": folder " + foldname + " already exists. "
                               "\nSome files might be overwritten. Continue (y/N/y-all)?")
            if to_del.lower() == 'y-all':
                ovwrite_files = True
            elif to_del.lower() != 'y':
                sys.exit("Aborting")

    for pars in allrange:

        cut = pars['cut']
        kpn = pars['kpn']
        fgm = pars['fgm']

        jobname = jname(seedname, cut, kpn, fgm)

        if str_par_vals["rmode"] == "parallel":
            foldname = jobname

        # If we're reusing data, skip recreating the folder only if we already have some results

        if bool_par_vals["rcalc"]:
            try:
                if jobfinish_check(foldname, jobname):
                    continue
            except JobError:
                pass

        create_conv_folder(foldname, jobname, cut, kpn, fgm, prev_jobname)

        if str_par_vals["rmode"] == "serial":
            prev_jobname = jobname


# PHASE 2 - EXECUTION
# Run through the files and obtain the results

if (str_par_vals["ctsk"] in ("all", "inputrun")):

    if str_par_vals["rmode"] == "parallel":

        print "Running parallel convergence jobs"

        running_jobs = []

    elif str_par_vals["rmode"] == "serial":

        print "Running serial convergence jobs"

        running_jobs = None

        foldname = seedname + "_conv"

    else:
        sys.exit("ERROR - Invalid value for the running_mode parameter")

    for i, pars in enumerate(allrange):

        cut = pars['cut']
        kpn = pars['kpn']
        fgm = pars['fgm']

        jobname = jname(seedname, cut, kpn, fgm)
        if str_par_vals["rmode"] == "parallel":
            foldname = jobname

        job_run(foldname, jobname, running_jobs)
        # Wait for everyone to finish
        multijob_wait(running_jobs, wait_all=(i==(len(allrange)-1)))

    if (str_par_vals["ctsk"] == "all"):
        print "\n -- All jobs finished. Proceeding to analyze output --\n"

# PHASE 3 - OUTPUT PARSING
# Here is where we parse the .CASTEP files and draw our conclusions

if (str_par_vals["ctsk"] in ("all", "output")):

    # If this is an output job, try opening a .conv_tab file and rebuild cutrange and kpnrange from that - otherwise just use the parameters from the .conv file

    if (cutrange is None or kpnrange is None or fgmrange is None):
        if os.path.isfile(seedname + ".conv_tab"):
            cutrange, kpnrange, fgmrange = parse_conv_tab_file(open(seedname + ".conv_tab", 'r'))
        else:
            cut_n = int(math.ceil(float_par_vals["cutmax"]-float_par_vals["cutmin"])/float_par_vals["cutstep"])+1
            cutrange = [float_par_vals["cutmin"] + i * float_par_vals["cutstep"] for i in range(0, cut_n)]
    
            kpn_n = int(math.ceil(int_par_vals["kpnmax"]-int_par_vals["kpnmin"])/int_par_vals["kpnstep"])+1
            kpnrange = [tuple([int((int_par_vals["kpnmin"] + i * int_par_vals["kpnstep"]) * e) for e in kpn_base]) for i in range(0, kpn_n)]
    
            if (str_par_vals["fgmmode"] is not None):
                fgm_n = int(math.ceil(float_par_vals["fgmmax"]-float_par_vals["fgmmin"])/float_par_vals["fgmstep"])+1
                fgmrange = [float_par_vals["fgmmin"] + i * float_par_vals["fgmstep"] for i in range(0, fgm_n)]
            else:
                fgmrange = [None]

    if (len(cutrange) + len(kpnrange) + len(fgmrange) == 0):
        sys.exit("ERROR - Corrupted .conv_tab file. Output task can not be carried out")
    
    calc_str = bool_par_vals["cnvstr"]
    conv_fgm = (str_par_vals["fgmmode"] is not None)

    if (allrange is None):
        fgm_cutoff = cutrange[-1] if str_par_vals['fgmmode'] == 'max_cutoff' else cutrange[0]
        allrange = []
        for cut in cutrange:
            allrange.append({'cut': cut, 'kpn': kpnrange[0], 'fgm': fgmrange[0], 'scan': 'cut'})
        for kpn in kpnrange[1:]:
            allrange.append({'cut': cutrange[0], 'kpn': kpn, 'fgm': fgmrange[0], 'scan': 'kpn'})
        for fgm in fgmrange[1:]:
            allrange.append({'cut': fgm_cutoff, 'kpn': kpnrange[0], 'fgm': fgm, 'scan': 'fgm'})

    cutnrg = []
    cutfor = []
    cutstr = []
    kpnnrg = []
    kpnfor = []
    kpnstr = []
    fgmnrg = []
    fgmfor = []
    fgmstr = []

    # Try opening the various .castep files and collect energy and forces

    for pars in allrange:

        cut = pars['cut']
        kpn = pars['kpn']
        fgm = pars['fgm']

        scan = pars['scan']

        jobname = jname(seedname, cut, kpn, fgm)

        if (str_par_vals["rmode"] == "serial"):
            foldname = seedname + "_conv"
        elif (str_par_vals["rmode"] == "parallel"):
            foldname = jobname
        else:
            sys.exit("ERROR - Invalid value for the running_mode parameter")

        filepath = os.path.join(foldname, jobname + '.castep')

        if not os.path.isfile(filepath):
            sys.exit("ERROR - File " + filepath + " not found. Please check the .conv file and that all calculations have actually finished")

        castepfile = open(filepath, 'r').readlines()

        try:
            castep_data = parse_castep_file(castepfile, filepath)
        except CastepError as CE:
            sys.exit("ERROR - " + str(CE))

        atom_n = castep_data['atom_n']

        if castep_data['cut'] is None or castep_data['kpn'] is None or castep_data['cut'] != cut or castep_data['kpn'] != kpn:
            sys.exit("ERROR - Simulation parameters in " + filepath + " do not correspond to expected values")

        if conv_fgm and (castep_data['fgm'] is None or castep_data['fgm'] != round_digits(cut_to_k(fgm), 4)):
            if scan in ('fgm', 'kpn'):
                sys.exit("ERROR - Simulation parameters in " + filepath + " do not correspond to expected values")
            else:
                print __WARNING__ + " - fine_Gmax value used in " + filepath + " is greater than the set value of fine_gmax_min due to the higher cutoff"

        if castep_data['nrg'] is None or castep_data['for'] is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing energy or forces)")

        if calc_str and castep_data['str'] is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing stresses)")

        if castep_data['atom_n'] is None:
            sys.exit("ERROR - Corrupted " + filepath + " file")

        if (scan == 'cut'):
            cutnrg.append(castep_data['nrg'])
            cutfor.append(castep_data['for'])
            if calc_str:
                cutstr.append(castep_data['str'])
        elif (scan == 'kpn'):
            if (len(kpnnrg) == 0):
                kpnnrg.append(cutnrg[0])
                kpnfor.append(cutfor[0])
                if calc_str:
                    kpnstr.append(cutstr[0])
            kpnnrg.append(castep_data['nrg'])
            kpnfor.append(castep_data['for'])
            if calc_str:
                kpnstr.append(castep_data['str'])
        elif (scan == 'fgm'):
            # This could be useful if there are NO kpn points
            if (len(kpnnrg) == 0):
                kpnnrg.append(cutnrg[0])
                kpnfor.append(cutfor[0])
                if calc_str:
                    kpnstr.append(cutstr[0])
            if (len(fgmnrg) == 0):
                fgmnrg.append(cutnrg[0])
                fgmfor.append(cutfor[0])
                if calc_str:
                    fgmstr.append(cutstr[0])
            fgmnrg.append(castep_data['nrg'])
            fgmfor.append(castep_data['for'])
            if calc_str:
                fgmstr.append(castep_data['str'])

    # Save the results in data files

    cutfile = open(seedname + "_cut_conv.dat", 'w')
    cutfile.write("Cutoff (eV)\tEnergy (eV)\tForces (eV/A)" + ("\tTotal stresses (GPa)" if calc_str else "") + "\n")

    for i in range(0, len(cutrange)):
        cutfile.write(str(cutrange[i]) + '\t\t' + str(cutnrg[i]) + '\t\t' + str(cutfor[i]) + ('\t\t' + str(cutstr[i]) if calc_str else '') + '\n')
    
    cutfile.close()

    kpnfile = open(seedname + "_kpn_conv.dat", 'w')
    kpnfile.write("kpoints (tot)\tEnergy (eV)\tForces (eV/A)" + ("\tTotal stresses (GPa)" if calc_str else "") + "\n")

    for i in range(0, len(kpnrange)):
        kpnfile.write(str(kpnrange[i][0]*kpnrange[i][1]*kpnrange[i][2]) + '\t\t' + str(kpnnrg[i]) + '\t\t' + str(kpnfor[i]) + ('\t\t' + str(kpnstr[i]) if calc_str else '') + '\n')

    kpnfile.close()

    if (conv_fgm):
        fgmfile = open(seedname + "_fgm_conv.dat", 'w')
        fgmfile.write("Fine Gmax (eV)\tEnergy (eV)\t\tMax force (eV/A)" + ("\tTotal stress (GPa)" if calc_str else "") + "\tFine Gmax (1/A)\n")
    
        for i in range(0, len(fgmrange)):
            fgmfile.write(str(fgmrange[i]) + '\t\t' + str(fgmnrg[i]) + '\t\t' + str(fgmfor[i]) + ('\t\t' + str(fgmstr[i]) if calc_str else '') + '\t\t' + str(round_digits(cut_to_k(fgmrange[i]), 4)) + '\n')
        
        fgmfile.close()
    
    # Make an estimate of the best values

    conv_data = {}
    conv_data['cut'] = {'range': cutrange,
        'rangestr': [str(c) for c in cutrange],
        'step': float_par_vals['cutstep'],
        'nrg': [x/atom_n for x in cutnrg],
        'for': cutfor}
    conv_data['kpn'] = {'range': [k[0]*k[1]*k[2] for k in kpnrange],
        'rangestr': [kgrid(k, 'x') for k in kpnrange],
        'step': int_par_vals['kpnstep']*3,
        'nrg': [x/atom_n for x in kpnnrg],
        'for': kpnfor}
    conv_data['fgm'] = {'range': fgmrange,
        'rangestr': [str(f) for f in fgmrange],
        'step': float_par_vals['fgmstep'],
        'nrg': [x/atom_n for x in fgmnrg],
        'for': fgmfor}

    if bool_par_vals['cnvstr']:
        conv_data['cut']['str'] = cutstr
        conv_data['kpn']['str'] = kpnstr
        conv_data['fgm']['str'] = fgmstr

    opt_estimates = conv_estimates(seedname, conv_data, calc_str)
    
    # Finally, let's do the plotting

    if str_par_vals["outp"] == "gnuplot":
        gp_graph(seedname, calc_str)

    elif str_par_vals["outp"] in ("xmgrace", "grace"):
        agr_graph(seedname, conv_data, calc_str)
