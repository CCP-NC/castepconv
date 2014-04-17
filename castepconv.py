#!/usr/bin/env python

# CASTEP convergence automation tool
# by Simone Sturniolo
#
# Copyright 2013 Science and Technology Facilities Council
# This software is distributed under the terms of the GNU General Public License (GNU GPL)

import sys, time, math, os, shutil, glob, re
import subprocess as sp

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

# Length units allowed by CASTEP. Internally used unit is always Angstroms 

length_units = {'ang':1.0, 'm':1.0e10, 'cm':1.0e8, 'nm':10.0, 'bohr':0.529, 'a0':0.529}

# Parameters from CONV files - names and values

str_par_names = {
"convergence_task"  : "ctsk",
"running_mode"      : "rmode",
"output_type"       :  "outp",
"running_command"   : "rcmd",
}

str_par_vals = {
"ctsk"   : "input",                         # Can be INPUT, INPUTRUN, OUTPUT or ALL
"rmode"  : "parallel",                      # Can be PARALLEL or SERIAL
"outp"   : "gnuplot",                       # Can be GNUPLOT or XMGRACE
"rcmd"   : "castep <seedname> -dryrun",
}

float_par_names = {
"cutoff_min"        : "cutmin",
"cutoff_max"        : "cutmax",
"cutoff_step"       : "cutstep",
"displace_atoms"    : "displ",
"final_energy_delta": "nrgtol",
"forces_delta"      : "fortol",
"stresses_delta"    : "strtol"
}

float_par_vals = { 
"cutmin" : 400.0,           # eV
"cutmax" : 800.0,           # eV
"cutstep": 100.0,           # eV
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
"rcalc"  : True,
"sruse"  : True
}

# CELL and PARAM files

cellfile_lines = None
paramfile_lines = None

cutrange = None
kpnrange = None

abc_len = None
kpn_base = (1, 1, 1)
pseudo_pots = None
has_fix_occ = False

ovwrite_files = False

__CASTEP_HEADER__       = "+-------------------------------------------------+"
__CASTEP_TIME__         = "Total time          ="
__CASTEP_ATOMN__        = "Total number of ions in cell = "
__CASTEP_CUTOFF__       = "plane wave basis set cut-off                   :"
__CASTEP_KPOINTS__      = "MP grid size for SCF calculation is"
__CASTEP_ENERGY__       = "Final energy, E             ="
__CASTEP_ENERGY_FIX__   = "Final energy ="                                                      
__CASTEP_FORCES__       = "***************** Symmetrised Forces *****************"
__CASTEP_FORCES_ALT__   = "*********************** Forces ***********************"
__CASTEP_FORCES_END__   = "*                                                    *"
__CASTEP_STRESSES__     = "*********** Symmetrised Stress Tensor ***********"
__CASTEP_STRESSES_ALT__ = "***************** Stress Tensor *****************"
__CASTEP_STRESSES_END__ = "*                                               *"

# Just a useful snippet to return a kpoint grid from a 3-ple

def kgrid(t): return str(t[0]) + '\t' + str(t[1]) + '\t' + str(t[2])

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
        #Skip comments
        if (l[0] == '#'):
            continue
        
        #Parse options
        cline = l.strip()
        if len(cline) == 0:
            #Skip empty lines
            continue
        else:
            cline = cline.split(':')
        if (len(cline) < 2):
            raise ConvError("Bad formatting in .conv file at line " + str(i))
        par_name = cline[0].strip().lower()
        if (par_name in str_par_names):
            # A condition added to take into account weird command line instructions
            if par_name == 'rcmd':
                cline[1] = ':'.join(cline[1:])
            else:
                cline[1] = cline[1].lower()
            str_par_vals[str_par_names[par_name]] = cline[1].strip()
        elif (par_name in float_par_names):
            float_par_vals[float_par_names[par_name]] = float(cline[1])
        elif (par_name in int_par_names):
            int_par_vals[int_par_names[par_name]] = int(cline[1])
        elif (par_name in bool_par_names):
            bool_par_vals[bool_par_names[par_name]] = (cline[1].lower().strip() == "true")
        else:
            raise ConvError("Unrecognized option in .conv file at line " + str(i))

# Parse .conv_tab file to generate cutoff and k points ranges

def parse_conv_tab_file(cfile):
    
    cutrange = []
    kpnrange = []
    
    cfile.seek(0)
    clines = cfile.readlines()
    
    if len(clines) == 2:
        try:
            cutline = clines[0].split(':')[1].split('eV')
            kpnline = clines[1].split(':')[1].split('|')
            if (cutline is not None and kpnline is not None):
                cutrange = [float(cut.strip()) for cut in cutline[:-1]]
                kpnrange = [tuple([int(k) for k in kpn.strip().split()]) for kpn in kpnline[:-1]]
        except IndexError:
            return [], []
    else:
        return [], []
    
    return cutrange, kpnrange

# Strip .cell file from unnecessary lines to get only the ones we care for (i.e. remove all reference to kpoints, we're going to put those in ourselves)
# Also read cell parameters and construct the proper kpn_base

def strip_cellfile(clines):
    
    stripped = []
    
    to_strip = False
    to_read_abc = False
    to_read_cart = False
    to_read_pot = False
    
    abc = None
    u = 1.0
    kbase = (1, 1, 1)
    ppot = None
    
    for l in clines:
        
        if not l[0] == '#':
            
            l_low = l.lower()
            
            if to_read_abc:
                l_split = l_low.strip().split()
                if l_split[0] in length_units:
                    u = length_units[l_split[0]]
                else: 
                    try:
                        abc = [u*float(x) for x in l_split]
                    except ValueError:
                        raise CellError('Bad formatting in .cell file LATTICE_ABC block')
                    to_read_abc = False
            if to_read_cart:
                l_split = l_low.strip().split()
                if l_split[0] in length_units:
                    u = length_units[l_split[0]]
                else:
                    try:
                        for i in range(0, 3):
                            if abc[i] < 0:
                                abc[i] = math.sqrt(sum([(u*float(x))**2.0 for x in l_split]))
                                break
                    except ValueError:
                        raise CellError('Bad formatting in .cell file LATTICE_CART block')
                    to_read_cart = (min(abc) < 0.0)
            if to_read_pot:
                l_split = l.strip().split()
                if "%endblock" in l_low:
                    if not "species_pot" in l_low:
                        raise CellError('Bad formatting in .cell file SPECIES_POT block')
                    else:
                        to_read_pot = False
                else:
                    try:
                        ppot.append([l_split[0], l_split[1]])
                    except ValueError:
                        raise CellError('Bad formatting in .cell file SPECIES_POT block')
                continue
            
            if "kpoints_mp_grid" in l_low:
                continue
            if "kpoints_mp_spacing" in l_low:
                continue
            if "%block" in l_low:
                if "kpoints_list" in l_low:
                    to_strip = True
                    continue
                elif "lattice_abc" in l_low:
                    if abc is not None:
                        raise CellError('Duplicated LATTICE_* block in .cell file')
                    to_read_abc = True
                    abc = [-1.0, -1.0, -1.0]
                elif "lattice_cart" in l_low:
                    if abc is not None:
                        raise CellError('Duplicated LATTICE_* block in .cell file')
                    to_read_cart = True
                    abc = [-1.0, -1.0, -1.0]
                elif "species_pot" in l_low:
                    if ppot is not None:
                        raise CellError('Duplicated SPECIES_POT block in .cell file')                        
                    to_read_pot = True
                    ppot = []
                    continue
            if to_strip:
                if "%endblock" in l_low:
                    to_strip = False
                continue
        
        stripped.append(l)
        
    if abc is None:
        raise CellError('.cell file does not contain a LATTICE_* block')
    
    max_e = max(abc)
    kbase = tuple([int(max_e/x) for x in abc])
    
    return stripped, tuple(abc), kbase, ppot

# Strip .param file from unnecessary lines to get only the ones we care for (i.e. remove all reference to task and cutoff)

def strip_paramfile(plines):
    
    global bool_par_vals, has_fix_occ
    
    stripped = []
    
    for l in plines:
        
        l_split = l.strip().lower().split()
        
        # Skip empty lines
        
        if len(l_split) == 0:
            continue
        if "fix_occupancy" in l_split[0]:
            if len(l_split) > 1 and l_split[1] == "true":
                has_fix_occ = True
        if "task" in l_split[0]:
            continue
        if "cut_off_energy" in l_split[0]:
            continue
        if "calculate_stress" in l_split[0] and bool_par_vals["cnvstr"]:
            continue
        if "reuse" in l_split[0] and bool_par_vals["sruse"]:
            continue
        
        stripped.append(l)
    
    return stripped

# Displace atoms by a small amount in cellfile if needed to have nonzero forces

def displace_cell_atoms(clines, abc, d):
    
    to_displ_abs = False
    to_displ_frac = False
    u = None
    
    clines_displ = []
    
    for i, l in enumerate(clines):
        
        if not l[0] == '#':
            
            l_low = l.lower()
            
            if to_displ_abs or to_displ_frac:
                l_split = l.strip().split()
                if "%endblock" in l_low:
                    to_displ_abs = False
                    to_displ_frac = False
                elif to_displ_abs and l_split[0].lower() in length_units:
                    u = length_units[l_split[0].lower()]
                    l = "ang\n"
                else:
                    if len(l_split) != 4:
                        raise CellError('Bad formatting in .cell file POSITION_* block')
                    try:
                        xyz = [float(x) for x in l_split[1:]]
                        # The displacement is in alternated directions, since otherwise you'd still get an equilibrium structure
                        # Won't work for single atom unit cells. This needs fixing. It might also break if there are symmetries in the system.
                        fac = (i%2)*2-1
                        if to_displ_abs:
                            l = l_split[0] + '\t' + str(u*xyz[0]+fac*d) + '\t' + str(u*xyz[1]+fac*d) + '\t' + str(u*xyz[2]+fac*d) + '\n'
                        elif to_displ_frac:
                            l = l_split[0] + '\t' + str(xyz[0]+fac*d/abc[0]) + '\t' + str(xyz[1]+fac*d/abc[1]) + '\t' + str(xyz[2]+fac*d/abc[2]) + '\n'
                    except ValueError:
                        raise CellError('Bad formatting in .cell file POSITION_* block')
            if "%block" in l_low:
                if "positions_abs" in l_low:
                    to_displ_abs = True
                    u = 1.0
                elif "positions_frac" in l_low:
                    to_displ_frac = True
        
        clines_displ.append(l)
    
    return clines_displ

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

def create_conv_folder(foldname, jobname, cut, kpn, prev_jobname=None):
    
    global ovwrite_files, stripped_cell, stripped_param, bool_par_vals, str_par_vals, pseudo_pots
    
    if not os.path.exists(foldname): 
        print "Creating folder " + foldname
        os.makedirs(foldname)
    elif not ovwrite_files and len(os.listdir(foldname)) > 0:
        to_del = raw_input("Warning: folder " + foldname + " already exists. \
        \nSome files might be overwritten. Continue (y/N/y-all)?")
        if to_del.lower() == 'y-all':
            ovwrite_files = True
        elif to_del.lower() != 'y':
            sys.exit("Aborting")
    
    print "Creating files for job " + jobname
    icell = open(os.path.join(foldname, jobname + '.cell'), 'w')
    iparam = open(os.path.join(foldname, jobname + '.param'), 'w')
    
    for l in stripped_cell:
        icell.write(l)
    
    # Write pseudopotentials block
    
    if pseudo_pots is not None and len(pseudo_pots) > 0:
        icell.write("%BLOCK SPECIES_POT\n")
        for p in pseudo_pots:
            icell.write(p[0] + '\t' + p[1] + '\n')
        icell.write("%ENDBLOCK SPECIES_POT\n")
    
    # Write kpoint grid
    
    icell.write("\nkpoint_mp_grid " + kgrid(kpn) + "\n")
    icell.close()
    
    iparam.write("task:\tSinglePoint\n")
    iparam.write("cut_off_energy:\t" + str(cut) + " eV\n")
    for l in stripped_param: 
        iparam.write(l)
    if prev_jobname is not None and bool_par_vals["sruse"]:
        iparam.write("reuse:\t" + prev_jobname + ".check\n")
    if bool_par_vals["cnvstr"]:
        iparam.write("calculate_stress:\ttrue\n")
    iparam.close()

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
                    to_del = raw_input("Warning: folder " + pp_folder + " already exists. \
                    \nSome files might be overwritten. Continue (y/N/y-all)?")
                    if to_del.lower() == 'y-all':
                        ovwrite_files = True
                    elif to_del.lower() != 'y':
                        sys.exit("Aborting")
                shutil.copy2(pseudo_pots[i][1], os.path.join(pp_folder, pseudo_pots[i][1]))
                pseudo_pots[i][1] = os.path.join('..', pp_folder, pseudo_pots[i][1])
            else:
                raise PotError("Pseudo potential file " + pseudo_pots[i][1] + " could not be found")

###### -- MAIN PROGRAM -- ######

if (sys.version_info[0] < 2 or sys.version_info[1] < 6):
    sys.exit("ERROR - Python version 2.6 or higher required to run the script")

seedname, cmdline_task = parse_cmd_args()

if seedname is None:
    sys.exit("ERROR - <seedname> is a required argument")

# PHASE 0 - Check for existence of all required files and read the necessary information

print "Reading " + seedname + ".conv"

try:
    job_convfile = open(seedname + ".conv", 'r')
except IOError:
    print "WARNING - Convergence parameter file for job " + seedname + " not found. Using default parameters"
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
        job_cellfile = open(seedname + ".cell", 'r')
        cellfile_lines = job_cellfile.readlines()
    except IOError:
        sys.exit("ERROR - .cell file for job " + seedname + " not found")
    
    print "Reading " + seedname + ".param"

    try:
        job_paramfile = open(seedname + ".param", 'r')
        paramfile_lines = job_paramfile.readlines()
    except IOError:
        print("WARNING - .param file for job " + seedname + " not found")
        paramfile_lines = None
    
else:
    
    cellfile_lines = None
    paramfile_lines = None

# Override task in .conv file with command line options

if cmdline_task is not None:
    str_par_vals['ctsk'] = {'c': 'clear', 'i': 'input', 'ir': 'inputrun', 'o': 'output', 'a': 'all'}[cmdline_task]

if (str_par_vals['ctsk'] not in ("clear", "input", "inputrun", "output", "all")):
    sys.exit("ERROR - Invalid convergence_task parameter")
    
# PHASE 0.5 - CLEAR
# Clear all files and folders from previous jobs

if (str_par_vals['ctsk'] in ("clear")):
    
    # Create a list of files and folders to delete
    
    to_del_fold = glob.glob(seedname + "_cut_*_kpn_*")
    
    if os.path.exists(seedname + "_conv"):
        to_del_fold += [seedname + "_conv"]
    if os.path.exists(seedname + "_pspot"):
        to_del_fold += [seedname + "_pspot"]
    
    # Create a list of files to delete
    
    to_del_files = []
    
    if os.path.isfile(seedname + ".conv_tab"):
        to_del_files += [seedname + ".conv_tab"]
    
    if os.path.isfile(seedname + "_cut_conv.dat"):
        to_del_files += [seedname + "_cut_conv.dat"]
    if os.path.isfile(seedname + "_kpn_conv.dat"):
        to_del_files += [seedname + "_kpn_conv.dat"]
        
    if os.path.isfile(seedname + "_cut_conv.gp"):
        to_del_files += [seedname + "_cut_conv.gp"]
    if os.path.isfile(seedname + "_cut_str_conv.gp"):
        to_del_files += [seedname + "_cut_str_conv.gp"]
    if os.path.isfile(seedname + "_kpn_conv.gp"):
        to_del_files += [seedname + "_kpn_conv.gp"]
    if os.path.isfile(seedname + "_kpn_str_conv.gp"):
        to_del_files += [seedname + "_kpn_str_conv.gp"]
    
    if os.path.isfile(seedname + "_cut_conv.agr"):
        to_del_files += [seedname + "_cut_conv.agr"]
    if os.path.isfile(seedname + "_cut_str_conv.agr"):
        to_del_files += [seedname + "_cut_str_conv.agr"]
    if os.path.isfile(seedname + "_kpn_conv.agr"):
        to_del_files += [seedname + "_kpn_conv.agr"]
    if os.path.isfile(seedname + "_kpn_str_conv.agr"):
        to_del_files += [seedname + "_kpn_str_conv.agr"]
    
    print "The following folders will be deleted:"
    
    for f in to_del_fold:
        sys.stdout.write(f + ' ')
    
    print "\nThe following files will be deleted:"
    
    for f in to_del_files:
        sys.stdout.write(f + ' ')
    
    to_del = raw_input("\nContinue (y/N)?")
    
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
        stripped_cell, abc_len, kpn_base, pseudo_pots = strip_cellfile(cellfile_lines)
    except CellError as CE:
        sys.exit("ERROR - " + str(CE))
    if paramfile_lines is not None:
        stripped_param = strip_paramfile(paramfile_lines)
    else:
        stripped_param = []
    
    # Check the positions and existence of pseudopotentials
    
    try:
        find_pseudopots(seedname, pseudo_pots)
    except PotError as PE:
        sys.exit("ERROR - " + str(PE))
    
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
        old_cutrange, old_kpnrange = parse_conv_tab_file(open(seedname + ".conv_tab", 'r'))
    else:
        old_cutrange = []
        old_kpnrange = []
    
    # Open a .conv_tab file to keep track of the created files and folders. Will be read if output is done as a separate operation
    
    if os.path.isfile(seedname + ".conv_tab") and not ovwrite_files:
        to_del = raw_input("Warning: " + seedname + ".conv_tab already exists. This file will be overwritten. Continue (y/N/y-all)?")
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
    kpnrange = [tuple([(int_par_vals["kpnmin"] + i * int_par_vals["kpnstep"]) * e for e in kpn_base]) for i in range(0, kpn_n)]
    
    if str_par_vals["rmode"] == "parallel":
        
        print "Creating folders for parallel convergence run"
        
        conv_tab_file.write("cutoff:\t")
        
        for i, cut in enumerate(cutrange):
            
            conv_tab_file.write(str(cut) + " eV\t")
            
            foldname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
            
            # If we're reusing data, skip recreating the folder only if we already have some results
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, foldname):
                        continue
                except JobError:
                    pass
                
            create_conv_folder(foldname, foldname, cut, kpnrange[0])
            
        conv_tab_file.write("\n")
            
        conv_tab_file.write("kpoint_n:\t" + kgrid(kpnrange[0]) + "\t|\t")
        
        for i, kpn in enumerate(kpnrange[1:]):
            
            conv_tab_file.write(kgrid(kpn) + "\t|\t")
            
            foldname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, foldname):
                        continue
                except JobError:
                    pass
            
            create_conv_folder(foldname, foldname, cutrange[0], kpn)
            
        conv_tab_file.close()
        
    elif str_par_vals["rmode"] == "serial":
                
        if str_par_vals["ctsk"] == "input":
            sys.exit("ERROR - Impossible to carry out a SERIAL convergence test for an INPUT task. SERIAL tests require the execution step to be automated")
        
        foldname = seedname + "_conv"
        
        print "Creating folder " + foldname + " for serial convergence run"
        
        if not os.path.exists(foldname): 
            os.makedirs(foldname)
        elif not ovwrite_files:
            to_del = raw_input("Warning: folder " + foldname + " already exists. \
            \nSome files might be overwritten. Continue (y/N/y-all)?")
            if to_del.lower() == 'y-all':
                ovwrite_files = True
            elif to_del.lower() != 'y':
                sys.exit("Aborting")
                
        conv_tab_file.write("cutoff:\t")
        
        prev_jobname = None
        
        for i, cut in enumerate(cutrange):
            
            conv_tab_file.write(str(cut) + " eV\t")
            
            jobname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, jobname):
                        continue
                except JobError:
                    pass
            
            create_conv_folder(foldname, jobname, cut, kpnrange[0], prev_jobname)
            
            prev_jobname = jobname
            
        
        conv_tab_file.write("\n")
            
        conv_tab_file.write("kpoint_n:\t" + kgrid(kpnrange[0]) + "\t|\t")
        
        for i, kpn in enumerate(kpnrange[1:]):
            
            conv_tab_file.write(kgrid(kpn) + "\t|\t")
            
            jobname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, jobname):
                        continue
                except JobError:
                    pass
            
            create_conv_folder(foldname, jobname, cutrange[0], kpn, prev_jobname)
            
            prev_jobname = jobname
            
        
        conv_tab_file.close()
    else:
        sys.exit("ERROR - Invalid value for the running_mode parameter")
    
# PHASE 2 - EXECUTION
# Run through the files (or create them as we go if the task is serial) and obtain the results

if (str_par_vals["ctsk"] in ("all", "inputrun")):
    
    if str_par_vals["rmode"] == "parallel":
        
        print "Running parallel convergence jobs"
        
        running_jobs = []
        
        for i, cut in enumerate(cutrange):
            
            foldname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
            
            # If we're reusing data, skip running only if we already have some results
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, foldname):
                        continue
                except JobError:
                    pass
            
            print "Running job " + foldname
            
            os.chdir(foldname)
            
            if os.path.isfile(foldname + ".castep") or os.path.isfile(foldname + ".check") or os.path.isfile(foldname + ".0001.err"):
                print "Removing output files from previous jobs for " + foldname
                if os.path.isfile(foldname + ".castep"):
                    os.remove(foldname + ".castep")
                if os.path.isfile(foldname + ".check"):
                    os.remove(foldname + ".check")
                if os.path.isfile(foldname + ".0001.err"):
                    os.remove(foldname + ".0001.err")
            
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(foldname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = foldname
            # Note: subprocess.Popen opens a subprocess without waiting for it to finish.
            # In fact, if we don't take care to check that all files are closed (see later) they might as well not be and we might end with a crash
            try:
                sp.Popen(cmd_line)
            except OSError:
                sys.exit("ERROR - Command:\n>\t" + cmd_line[0] + "\ndoes not exist on this system")
            os.chdir("..")
            running_jobs.append([foldname, foldname])
            
            if int_par_vals["maxjobs"] > 0:
                while len(running_jobs) >= int_par_vals["maxjobs"]:
                    for job_ind in range(len(running_jobs)-1, -1, -1):
                        try:
                            if jobfinish_check(running_jobs[job_ind][0], running_jobs[job_ind][1]):
                                del running_jobs[job_ind]
                        except JobError as JE:
                            del running_jobs[job_ind]
                            print "WARNING - " + str(JE)
        
        for i, kpn in enumerate(kpnrange[1:]):
            
            foldname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, foldname):
                        continue
                except JobError:
                    pass
            
            print "Running job " + foldname
            
            os.chdir(foldname)
            
            if os.path.isfile(foldname + ".castep") or os.path.isfile(foldname + ".check") or os.path.isfile(foldname + ".0001.err"):
                print "Removing output files from previous jobs for " + foldname
                if os.path.isfile(foldname + ".castep"):
                    os.remove(foldname + ".castep")
                if os.path.isfile(foldname + ".check"):
                    os.remove(foldname + ".check")
                if os.path.isfile(foldname + ".0001.err"):
                    os.remove(foldname + ".0001.err")
            
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(foldname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = foldname
            try:
                sp.Popen(cmd_line)
            except OSError:
                sys.exit("ERROR - Command:\n>\t" + cmd_line[0] + "\ndoes not exist on this system")
            os.chdir("..")
            running_jobs.append([foldname, foldname])
            
            if int_par_vals["maxjobs"] > 0:
                while len(running_jobs) >= int_par_vals["maxjobs"]:
                    for job_ind in range(len(running_jobs)-1, -1, -1):
                        try:
                            if jobfinish_check(running_jobs[job_ind][0], running_jobs[job_ind][1]):
                                del running_jobs[job_ind]
                        except JobError as JE:
                            del running_jobs[job_ind]
                            print "WARNING - " + str(JE)
            
        # If we're running an ALL job, wait for everyone to finish
        
        if str_par_vals["ctsk"] == "all":
            
            print "Waiting for all jobs to finish \
            WARNING: The program can be terminated with Ctrl+C, but that could terminate also the running jobs"
            
            for i, cut in enumerate(cutrange):
                foldname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
                try:
                    jobfinish_wait(foldname, foldname)
                except JobError as JE:
                    print "WARNING - " + str(JE)
                
            for i, kpn in enumerate(kpnrange[1:]):
                foldname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
                try:
                    jobfinish_wait(foldname, foldname)
                except JobError as JE:
                    print "WARNING - " + str(JE)
            
            print "\n -- All jobs finished. Proceeding to analyze output --"
        
    elif str_par_vals["rmode"] == "serial":
        
        print "Running serial convergence jobs"
        
        foldname = seedname + "_conv"
        
        for i, cut in enumerate(cutrange):
            
            jobname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, jobname):
                        continue
                except JobError:
                    pass
                
            print "Running job with " + str(cut) + " eV cutoff"
            
            os.chdir(foldname)
            
            if os.path.isfile(jobname + ".castep") or os.path.isfile(jobname + ".check") or os.path.isfile(jobname + ".0001.err"):
                print "Removing output files from previous jobs for " + jobname
                if os.path.isfile(jobname + ".castep"):
                    os.remove(jobname + ".castep")
                if os.path.isfile(jobname + ".check"):
                    os.remove(jobname + ".check")
                if os.path.isfile(jobname + ".0001.err"):
                    os.remove(jobname + ".0001.err")
            
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(jobname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = jobname
            try:
                sp.Popen(cmd_line)
            except OSError:
                sys.exit("ERROR - Command:\n>\t" + cmd_line[0] + "\ndoes not exist on this system")
            os.chdir("..")
            
            try:
                jobfinish_wait(foldname, jobname)
            except JobError as JE:
                print "WARNING - " + str(JE)
                    
        for i, kpn in enumerate(kpnrange[1:]):
            
            jobname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
            
            if bool_par_vals["rcalc"]:
                try:
                    if jobfinish_check(foldname, foldname):
                        continue
                except JobError:
                    pass
            
            print "Running job with kpoint grid " + kgrid(kpn)
            
            os.chdir(foldname)
            
            if os.path.isfile(jobname + ".castep") or os.path.isfile(jobname + ".check") or os.path.isfile(jobname + ".0001.err"):
                print "Removing output files from previous jobs for " + jobname
                if os.path.isfile(jobname + ".castep"):
                    os.remove(jobname + ".castep")
                if os.path.isfile(jobname + ".check"):
                    os.remove(jobname + ".check")
                if os.path.isfile(jobname + ".0001.err"):
                    os.remove(jobname + ".0001.err")
            
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(jobname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = jobname
            try:
                sp.Popen(cmd_line)
            except OSError:
                sys.exit("ERROR - Command:\n>\t" + cmd_line[0] + "\ndoes not exist on this system")
            os.chdir("..")
            
            try:
                jobfinish_wait(foldname, jobname)
            except JobError as JE:
                print "WARNING - " + str(JE)
        
        if str_par_vals["ctsk"] == "all":
            print "\n -- All jobs finished. Proceeding to analyze output --\n"
        
    else:
        sys.exit("ERROR - Invalid value for the running_mode parameter")

if (str_par_vals["ctsk"] in ("all", "output")):
    
    # If this is an output job, try opening a .conv_tab file and rebuild cutrange and kpnrange from that - otherwise just use the parameters from the .conv file
    
    if (cutrange is None or kpnrange is None):
        cut_n = int(math.ceil(float_par_vals["cutmax"]-float_par_vals["cutmin"])/float_par_vals["cutstep"])+1
        cutrange = [float_par_vals["cutmin"] + i * float_par_vals["cutstep"] for i in range(0, cut_n)]
        
        kpn_n = int(math.ceil(int_par_vals["kpnmax"]-int_par_vals["kpnmin"])/int_par_vals["kpnstep"])+1
        kpnrange = [tuple([(int_par_vals["kpnmin"] + i * int_par_vals["kpnstep"]) * e for e in kpn_base]) for i in range(0, kpn_n)]
        
        if os.path.isfile(seedname + ".conv_tab"):
            cutrange, kpnrange = parse_conv_tab_file(open(seedname + ".conv_tab", 'r'))
    
    calc_str = bool_par_vals["cnvstr"]
    
    cutnrg = []
    cutfor = []
    cutstr = []
    kpnnrg = []
    kpnfor = []
    kpnstr = []
    
    # Try opening the various .castep files and collect energy and forces
    
    for i, cut in enumerate(cutrange):
        
        jobname = seedname + "_cut_" + str(cut) + "_kpn_" + str(min(kpnrange[0]))
        
        if (str_par_vals["rmode"] == "serial"):
            foldname = seedname + "_conv"
        else:
            foldname = jobname
                
        filepath = os.path.join(foldname, jobname + '.castep')
        
        if not os.path.isfile(filepath):
            sys.exit("ERROR - File " + filepath + " not found. Please check the .conv file and that all calculations have actually finished")
        
        castepfile = open(filepath, 'r').readlines()
                
        atom_n = None
        i_nrg = None
        i_for = None
        i_str = None
        cut_check = None
        kpn_check = None
        
        start_l = rindex_cont(castepfile, __CASTEP_HEADER__)
        
        for j, l in enumerate(castepfile[start_l:]):
                                    
            try:
                if __CASTEP_ATOMN__ in l and atom_n is None:
                    atom_n = int(l.split()[7])
                elif __CASTEP_CUTOFF__ in l:
                    cut_check = float(l.split()[6])
                elif __CASTEP_KPOINTS__ in l:
                    kpn_check = tuple([int(x) for x in l.split()[7:]])
                elif __CASTEP_ENERGY__ in l and not has_fix_occ:
                    i_nrg = float(l.split()[4])
                elif __CASTEP_ENERGY_FIX__ in l and has_fix_occ:
                    i_nrg = float(l.split()[3])
                elif __CASTEP_FORCES__ in l or __CASTEP_FORCES_ALT__ in l:
                    i_for = parse_forces(castepfile[start_l+j:])
                elif calc_str and (__CASTEP_STRESSES__ in l or __CASTEP_STRESSES_ALT__ in l):
                    i_str = parse_stresses(castepfile[start_l+j:])
            except ValueError:
                sys.exit("ERROR - Corrupted " + filepath + " file detected")
            except CastepError as CE:
                sys.exit("ERROR - " + str(CE) + " in file " + filepath)
            
            
            if atom_n is not None and i_nrg is not None and i_for is not None and cut_check is not None and kpn_check is not None and ((calc_str and i_str is not None) or not calc_str):
                break
        
        if cut_check is None or kpn_check is None or cut_check != cut or kpn_check != kpnrange[0]:
            sys.exit("ERROR - Simulation parameters in " + filepath + " do not correspond to expected values")
        
        if i_nrg is None or i_for is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing energy or forces)")
        
        if calc_str and i_str is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing stresses)")
        
        if atom_n is None:
            sys.exit("ERROR - Corrupted " + filepath + " file")
        
        cutnrg.append(i_nrg)
        cutfor.append(i_for)
        if calc_str:
            cutstr.append(i_str)
    
    kpnnrg.append(cutnrg[0])
    kpnfor.append(cutfor[0])
    
    if calc_str:
        kpnstr.append(cutstr[0])
    
    for i, kpn in enumerate(kpnrange[1:]):
        
        jobname = seedname + "_cut_" + str(cutrange[0]) + "_kpn_" + str(min(kpn))
        
        if (str_par_vals["rmode"] == "serial"):
            foldname = seedname + "_conv"
        else:
            foldname = jobname
        
        filepath = os.path.join(foldname, jobname + '.castep')
        
        if not os.path.isfile(filepath):
            sys.exit("ERROR - File " + filepath + " not found. Please check the .conv file and that all calculations have actually finished")
        
        castepfile = open(filepath, 'r').readlines()
        
        i_nrg = None
        i_for = None
        i_str = None
        cut_check = None
        kpn_check = None
        
        start_l = rindex_cont(castepfile, __CASTEP_HEADER__)
        
        for j, l in enumerate(castepfile[start_l:]):
            
            try:
                if __CASTEP_CUTOFF__ in l:
                    cut_check = float(l.split()[6])
                elif __CASTEP_KPOINTS__ in l:
                    kpn_check = tuple([int(x) for x in l.split()[7:]])
                elif __CASTEP_ENERGY__ in l and not has_fix_occ:
                    i_nrg = float(l.split()[4])
                elif __CASTEP_ENERGY_FIX__ in l and has_fix_occ:
                    i_nrg = float(l.split()[3])
                elif __CASTEP_FORCES__ in l or __CASTEP_FORCES_ALT__ in l:
                    i_for = parse_forces(castepfile[start_l+j:])
                elif calc_str and (__CASTEP_STRESSES__ in l or __CASTEP_STRESSES_ALT__ in l):
                    i_str = parse_stresses(castepfile[start_l+j:])
            except ValueError:
                sys.exit("ERROR - Corrupted " + filepath + " file detected")
            except CastepError as CE:
                sys.exit("ERROR - " + str(CE) + " in file " + filepath)
            
            if i_nrg is not None and i_for is not None and cut_check is not None and kpn_check is not None and ((calc_str and i_str is not None) or not calc_str):
                break
        
        if cut_check is None or kpn_check is None or cut_check != cutrange[0] or kpn_check != kpn:
            sys.exit("ERROR - Simulation parameters in " + filepath + " do not correspond to expected values")
        
        if i_nrg is None or i_for is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing energy or forces)")
        
        if calc_str and i_str is None:
            sys.exit("ERROR - Incomplete simulation results in " + filepath + " (missing stresses)")
            
        kpnnrg.append(i_nrg)
        kpnfor.append(i_for)
        if calc_str:
            kpnstr.append(i_str)
    
    # Save the results in data files
    
    cutfile = open(seedname + "_cut_conv.dat", 'w')
    kpnfile = open(seedname + "_kpn_conv.dat", 'w')
    
    cutfile.write("Cutoff (eV)\tEnergy (eV)\tForces (eV/A)" + ("\tTotal stresses (GPa)" if calc_str else "") + "\n")
    kpnfile.write("kpoints (tot)\tEnergy (eV)\tForces (eV/A)" + ("\tTotal stresses (GPa)" if calc_str else "") + "\n")
    
    for i in range(0, len(cutrange)):
        cutfile.write(str(cutrange[i]) + '\t\t' + str(cutnrg[i]) + '\t\t' + str(cutfor[i]) + ('\t\t' + str(cutstr[i]) if calc_str else '') + '\n')
    for i in range(0, len(kpnrange)):
        kpnfile.write(str(kpnrange[i][0]*kpnrange[i][1]*kpnrange[i][2]) + '\t\t' + str(kpnnrg[i]) + '\t\t' + str(kpnfor[i]) + ('\t\t' + str(kpnstr[i]) if calc_str else '') + '\n')
    
    cutfile.close()
    kpnfile.close()
    
    # Make an estimate of the best values
        
    print "Convergence results:"
    
    if len(cutnrg) < 2:
        print "Impossible to give a convergence estimate with a single cutoff point"
    elif float_par_vals["nrgtol"] <= 0.0:
        print "Impossible to give a convergence estimate with a null or negative value for final_energy_delta"        
    else:
        
        delta_nrg = cutnrg[0]
        
        for i, nrg in enumerate(cutnrg[1:]):
            
            delta_nrg = abs(nrg - delta_nrg)
            if delta_nrg/atom_n < float_par_vals["nrgtol"]:
                print "Based on converging total energy to " + str(float_par_vals["nrgtol"]) + " eV per atom, minimum cutoff suggested is " + str(cutrange[i]) + " eV"
                break
            
            delta_nrg = nrg
            
            if i == len(cutnrg[1:])-1:
                print "Unable to converge cutoff with energy within given range.  Try increasing cutoff_max in " + seedname + ".conv"
        
    if len(cutfor) < 2:
        print "Impossible to give a convergence estimate with a single force point"
    elif float_par_vals["fortol"] <= 0.0:
        print "Impossible to give a convergence estimate with a null or negative value for forces_delta"        
    else:
        
        delta_for = cutfor[0]
        
        if delta_for < float_par_vals["fortol"]:
            print "WARNING - Maximum force is lower than " + str(float_par_vals["fortol"]) + " eV/Ang. A different atom displacement might be necessary to get meaningful results"
                    
        for i, force in enumerate(cutfor[1:]):
            
            delta_for = abs(force - delta_for)
            if delta_for < float_par_vals["fortol"]:
                print "Based on converging maximum force to " + str(float_par_vals["fortol"]) + " eV/Ang, minimum cutoff suggested is " + str(cutrange[i]) + " eV"
                break
            
            delta_for = force
            
            if i == len(cutfor[1:])-1:
                print "Unable to converge cutoff with forces within given range.  Try increasing cutoff_max in " + seedname + ".conv"
    
    if bool_par_vals["cnvstr"]:
        if len(cutstr) < 2:
            print "Impossible to give a convergence estimate with a single stresses point"
        elif float_par_vals["strtol"] <= 0.0:
            print "Impossible to give a convergence estimate with a null or negative value for stresses_delta"        
        else:
            
            delta_str = cutstr[0]
            
            if delta_str < float_par_vals["strtol"]:
                print "WARNING - Maximum stress is lower than " + str(float_par_vals["strtol"]) + " GPa. A different atom displacement might be necessary to get meaningful results"
                        
            for i, stress in enumerate(cutstr[1:]):
                
                delta_str = abs(stress - delta_str)
                if delta_str < float_par_vals["strtol"]:
                    print "Based on converging maximum stress to " + str(float_par_vals["fortol"]) + " GPa, minimum cutoff suggested is " + str(cutrange[i]) + " eV"
                    break
                
                delta_str = stress
                
                if i == len(cutstr[1:])-1:
                    print "Unable to converge cutoff with total stresses within given range.  Try increasing cutoff_max in " + seedname + ".conv"
    
    if len(kpnnrg) < 2:
        print "Impossible to give a convergence estimate with a single cutoff point"
    elif float_par_vals["nrgtol"] <= 0.0:
        print "Impossible to give a convergence estimate with a null or negative value for final_energy_delta"        
    else:
        
        delta_nrg = kpnnrg[0]
        
        for i, nrg in enumerate(kpnnrg[1:]):
            
            delta_nrg = abs(nrg - delta_nrg)
            if delta_nrg/atom_n < float_par_vals["nrgtol"]:
                print "Based on converging total energy to " + str(float_par_vals["nrgtol"]) + " eV per atom, minimum kpoint grid suggested is " + kgrid(kpnrange[i])
                break
            
            delta_nrg = nrg
            
            if i == len(kpnnrg[1:])-1:
                print "Unable to converge k-point grid with energy within given range.  Try increasing kpoint_n_max in " + seedname + ".conv"
        
    if len(kpnfor) < 2:
        print "Impossible to give a convergence estimate with a single force point"
    elif float_par_vals["fortol"] <= 0.0:
        print "Impossible to give a convergence estimate with a null or negative value for forces_delta"        
    else:
        
        delta_for = kpnfor[0]
        
        if delta_for < float_par_vals["fortol"]:
            print "WARNING - Maximum force is lower than " + str(float_par_vals["fortol"]) + " eV/Ang. A different atom displacement might be necessary to get meaningful results"
                    
        for i, force in enumerate(kpnfor[1:]):
            
            delta_for = abs(force - delta_for)
            if delta_for < float_par_vals["fortol"]:
                print "Based on converging maximum force to " + str(float_par_vals["fortol"]) + " eV/Ang, minimum kpoint grid suggested is " + kgrid(kpnrange[i])
                break
            
            delta_for = force
            
            if i == len(kpnfor[1:])-1:
                print "Unable to converge k-point grid with forces within given range.  Try increasing kpoint_n_max in " + seedname + ".conv"
    
    if bool_par_vals["cnvstr"]:
        if len(kpnstr) < 2:
            print "Impossible to give a convergence estimate with a single force point"
        elif float_par_vals["strtol"] <= 0.0:
            print "Impossible to give a convergence estimate with a null or negative value for stresses_delta"        
        else:
            
            delta_str = kpnstr[0]
            
            if delta_str < float_par_vals["strtol"]:
                print "WARNING - Maximum stress is lower than " + str(float_par_vals["strtol"]) + " GPa. A different atom displacement might be necessary to get meaningful results"
                        
            for i, stress in enumerate(kpnstr[1:]):
                
                delta_str = abs(stress - delta_str)
                if delta_str < float_par_vals["strtol"]:
                    print "Based on converging maximum stress to " + str(float_par_vals["strtol"]) + " GPa, minimum kpoint grid suggested is " + kgrid(kpnrange[i])
                    break
                
                delta_str = stress
                
                if i == len(kpnstr[1:])-1:
                    print "Unable to converge k-point grid with stresses within given range.  Try increasing kpoint_n_max in " + seedname + ".conv"
    
    if str_par_vals["outp"] == "gnuplot":
        
        out_file = open(seedname + "_cut_conv.gp", 'w')
        
        out_file.write("set xlabel \"Cutoff (eV)\"\n")
        out_file.write("set ylabel \"Final energy (eV)\"\n")
        out_file.write("set y2label \"Maximum force (eV/A)\"\n")
        out_file.write("set ytics nomirror\n")
        out_file.write("set y2tics\n")
        out_file.write("plot \"" + seedname + "_cut_conv.dat\" using 1:2 with linespoints pt 7 lc 1 ti \"Final energy\",")
        out_file.write("\"" + seedname + "_cut_conv.dat\" using 1:3 with linespoints pt 7 lc 2 axes x1y2 ti \"Force\"\n")
        out_file.write("pause -1 \"Hit return to continue\"\n")
        
        out_file.close()
        
        if bool_par_vals["cnvstr"]:
            
            out_file = open(seedname + "_cut_str_conv.gp", 'w')
            
            out_file.write("set xlabel \"Cutoff (eV)\"\n")
            out_file.write("set ylabel \"Final energy (eV)\"\n")
            out_file.write("set y2label \"Maximum stress (GPa)\"\n")
            out_file.write("set ytics nomirror\n")
            out_file.write("set y2tics\n")
            out_file.write("plot \"" + seedname + "_cut_conv.dat\" using 1:2 with linespoints pt 7 lc 1 ti \"Final energy\",")
            out_file.write("\"" + seedname + "_cut_conv.dat\" using 1:4 with linespoints pt 7 lc 3 axes x1y2 ti \"Stress\"\n")
            out_file.write("pause -1 \"Hit return to continue\"\n")
            
            out_file.close()
        
        out_file = open(seedname + "_kpn_conv.gp", 'w')
        
        out_file.write("set xlabel \"k-points\"\n")
        out_file.write("set ylabel \"Final energy (eV)\"\n")
        out_file.write("set y2label \"Total force (eV/A)\"\n")
        out_file.write("plot \"" + seedname + "_kpn_conv.dat\" using 1:2 with linespoints pt 7 ti \"Final energy\",")
        out_file.write("\"" + seedname + "_kpn_conv.dat\" using 1:3 with linespoints pt 7 axes x1y2 ti \"Forces\"\n")
        out_file.write("pause -1 \"Hit return to continue\"\n")
        
        out_file.close()
        
        if bool_par_vals["cnvstr"]:
            
            out_file = open(seedname + "_kpn_str_conv.gp", 'w')
            
            out_file.write("set xlabel \"k-points\"\n")
            out_file.write("set ylabel \"Final energy (eV)\"\n")
            out_file.write("set y2label \"Maximum stress (GPa)\"\n")
            out_file.write("set ytics nomirror\n")
            out_file.write("set y2tics\n")
            out_file.write("plot \"" + seedname + "_kpn_conv.dat\" using 1:2 with linespoints pt 7 lc 1 ti \"Final energy\",")
            out_file.write("\"" + seedname + "_kpn_conv.dat\" using 1:4 with linespoints pt 7 lc 3 axes x1y2 ti \"Stress\"\n")
            out_file.write("pause -1 \"Hit return to continue\"\n")
            
            out_file.close()
        
    elif str_par_vals["outp"] in ("xmgrace", "grace"):
        
        # Cutoff vs energy and forces
        
        out_file = open(seedname + "_cut_conv.agr", 'w')
        
        xrng  = (min(cutrange), max(cutrange))
        y1rng = (min(cutnrg)-0.1*(max(cutnrg)-min(cutnrg)), max(cutnrg)+0.1*(max(cutnrg)-min(cutnrg)))
        y2rng = (min(cutfor)-0.1*(max(cutfor)-min(cutfor)), max(cutfor)+0.1*(max(cutfor)-min(cutfor)))
        
        # Set up the graphics
        
        out_file.write("@version 50123\n")
        out_file.write("@title \"" + seedname + " - Energy and forces vs cutoff\"\n")
        
        out_file.write("@g0 on\n@g0 hidden false\n@with g0\n")
        out_file.write("@world " + str(xrng[0]) + ',' + str(y1rng[0]) + ',' + str(xrng[1]) + ',' + str(y1rng[1]) + '\n')
        out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
        out_file.write("@    xaxis label \"Cutoff (eV)\"\n")
        out_file.write("@    xaxis tick major " + str(float_par_vals["cutstep"]) + '\n')
        out_file.write("@    xaxis offset 0.0, 1.0\n")
        out_file.write("@    yaxis label \"Final energy (eV)\"\n")
        out_file.write("@    yaxis tick major " + str((y1rng[1]-y1rng[0])/8.0) + '\n')
        out_file.write("@    yaxis offset 0.0, 1.0\n")
        out_file.write("@    s0 hidden false\n@    s0 on\n")
        out_file.write("@    s0 legend \"Final energy\"\n")
        out_file.write("@    s0 line color 1\n")
        out_file.write("@    s0 symbol 1\n")
        out_file.write("@    s0 symbol size 0.7\n")
        out_file.write("@    s0 symbol color 1\n")
        out_file.write("@    s0 symbol fill color 1\n")
        out_file.write("@    s0 symbol fill pattern 1\n")
        
        out_file.write("@g1 on\n@g1 hidden false\n@with g1\n")
        out_file.write("@world " + str(xrng[0]) + ',' + str(y2rng[0]) + ',' + str(xrng[1]) + ',' + str(y2rng[1]) + '\n')
        out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
        out_file.write("@    xaxis label \"\"\n")
        out_file.write("@    xaxis tick off\n")
        out_file.write("@    xaxis tick major " + str(float_par_vals["cutstep"]) + '\n')
        out_file.write("@    yaxis label \"Force (eV/Ang)\"\n")
        out_file.write("@    yaxis tick major " + str((y2rng[1]-y2rng[0])/8.0) + '\n')
        out_file.write("@    yaxis ticklabel format exponential\n")
        out_file.write("@    yaxis ticklabel prec 1\n")
        out_file.write("@    yaxis offset 1.0, 0.0\n")
        out_file.write("@    yaxis label place opposite\n")
        out_file.write("@    yaxis tick place opposite\n")
        out_file.write("@    yaxis ticklabel place opposite\n")
        out_file.write("@    s0 hidden false\n@    s0 on\n")
        out_file.write("@    s0 legend \"Max force\"\n")
        out_file.write("@    s0 line color 2\n")
        out_file.write("@    s0 symbol 1\n")
        out_file.write("@    s0 symbol size 0.7\n")
        out_file.write("@    s0 symbol color 2\n")
        out_file.write("@    s0 symbol fill color 2\n")
        out_file.write("@    s0 symbol fill pattern 1\n")
        
        # Input the actual data
        
        out_file.write("@target G0.S0\n@type xy\n")
        for i, c in enumerate(cutrange):
            out_file.write(str(c) + '\t' + str(cutnrg[i]) + '\n')
        out_file.write('&\n')
        
        out_file.write("@target G1.S0\n@type xy\n")
        for i, c in enumerate(cutrange):
            out_file.write(str(c) + '\t' + str(cutfor[i]) + '\n')
        out_file.write('&\n')
        
        out_file.close()
        
        if bool_par_vals["cnvstr"]:
            
            # Cutoff vs energy and stresses
            
            out_file = open(seedname + "_cut_str_conv.agr", 'w')
            
            xrng  = (min(cutrange), max(cutrange))
            y1rng = (min(cutnrg)-0.1*(max(cutnrg)-min(cutnrg)), max(cutnrg)+0.1*(max(cutnrg)-min(cutnrg)))
            y2rng = (min(cutstr)-0.1*(max(cutstr)-min(cutstr)), max(cutstr)+0.1*(max(cutstr)-min(cutstr)))
            
            # Set up the graphics
            
            out_file.write("@version 50123\n")
            out_file.write("@title \"" + seedname + " - Energy and stresses vs cutoff\"\n")
            
            out_file.write("@g0 on\n@g0 hidden false\n@with g0\n")
            out_file.write("@world " + str(xrng[0]) + ',' + str(y1rng[0]) + ',' + str(xrng[1]) + ',' + str(y1rng[1]) + '\n')
            out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
            out_file.write("@    xaxis label \"Cutoff (eV)\"\n")
            out_file.write("@    xaxis tick major " + str(float_par_vals["cutstep"]) + '\n')
            out_file.write("@    xaxis offset 0.0, 1.0\n")
            out_file.write("@    yaxis label \"Final energy (eV)\"\n")
            out_file.write("@    yaxis tick major " + str((y1rng[1]-y1rng[0])/8.0) + '\n')
            out_file.write("@    yaxis offset 0.0, 1.0\n")
            out_file.write("@    s0 hidden false\n@    s0 on\n")
            out_file.write("@    s0 legend \"Final energy\"\n")
            out_file.write("@    s0 line color 1\n")
            out_file.write("@    s0 symbol 1\n")
            out_file.write("@    s0 symbol size 0.7\n")
            out_file.write("@    s0 symbol color 1\n")
            out_file.write("@    s0 symbol fill color 1\n")
            out_file.write("@    s0 symbol fill pattern 1\n")
            
            out_file.write("@g1 on\n@g1 hidden false\n@with g1\n")
            out_file.write("@world " + str(xrng[0]) + ',' + str(y2rng[0]) + ',' + str(xrng[1]) + ',' + str(y2rng[1]) + '\n')
            out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
            out_file.write("@    xaxis label \"\"\n")
            out_file.write("@    xaxis tick off\n")
            out_file.write("@    xaxis tick major " + str(float_par_vals["cutstep"]) + '\n')
            out_file.write("@    yaxis label \"Stress (GPa)\"\n")
            out_file.write("@    yaxis tick major " + str((y2rng[1]-y2rng[0])/8.0) + '\n')
            out_file.write("@    yaxis ticklabel format exponential\n")
            out_file.write("@    yaxis ticklabel prec 1\n")
            out_file.write("@    yaxis offset 1.0, 0.0\n")
            out_file.write("@    yaxis label place opposite\n")
            out_file.write("@    yaxis tick place opposite\n")
            out_file.write("@    yaxis ticklabel place opposite\n")
            out_file.write("@    s0 hidden false\n@    s0 on\n")
            out_file.write("@    s0 legend \"Max stress\"\n")
            out_file.write("@    s0 line color 2\n")
            out_file.write("@    s0 symbol 1\n")
            out_file.write("@    s0 symbol size 0.7\n")
            out_file.write("@    s0 symbol color 2\n")
            out_file.write("@    s0 symbol fill color 2\n")
            out_file.write("@    s0 symbol fill pattern 1\n")
            
            # Input the actual data
            
            out_file.write("@target G0.S0\n@type xy\n")
            for i, c in enumerate(cutrange):
                out_file.write(str(c) + '\t' + str(cutnrg[i]) + '\n')
            out_file.write('&\n')
            
            out_file.write("@target G1.S0\n@type xy\n")
            for i, c in enumerate(cutrange):
                out_file.write(str(c) + '\t' + str(cutstr[i]) + '\n')
            out_file.write('&\n')
            
            out_file.close()
        
        # K points vs energy and forces
        
        out_file = open(seedname + "_kpn_conv.agr", 'w')
        
        xrng  = (sum(kpnrange[0]), sum(kpnrange[-1]))
        y1rng = (min(kpnnrg)-0.1*(max(kpnnrg)-min(kpnnrg)), max(cutnrg)+0.1*(max(kpnnrg)-min(kpnnrg)))
        y2rng = (min(kpnfor)-0.1*(max(kpnfor)-min(kpnfor)), max(kpnfor)+0.1*(max(kpnfor)-min(kpnfor)))
        
        # Set up the graphics
        
        out_file.write("@version 50123\n")
        out_file.write("@title \"" + seedname + " - Energy and forces vs k points\"\n")
        
        out_file.write("@g0 on\n@g0 hidden false\n@with g0\n")
        out_file.write("@world " + str(xrng[0]) + ',' + str(y1rng[0]) + ',' + str(xrng[1]) + ',' + str(y1rng[1]) + '\n')
        out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
        out_file.write("@    xaxis label \"k points\"\n")
        out_file.write("@    xaxis tick major " + str(int_par_vals["kpnstep"]*3) + '\n')
        out_file.write("@    xaxis offset 0.0, 1.0\n")
        out_file.write("@    yaxis label \"Final energy (eV)\"\n")
        out_file.write("@    yaxis tick major " + str((y1rng[1]-y1rng[0])/8.0) + '\n')
        out_file.write("@    yaxis offset 0.0, 1.0\n")
        out_file.write("@    s0 hidden false\n@    s0 on\n")
        out_file.write("@    s0 legend \"Final energy\"\n")
        out_file.write("@    s0 line color 1\n")
        out_file.write("@    s0 symbol 1\n")
        out_file.write("@    s0 symbol size 0.7\n")
        out_file.write("@    s0 symbol color 1\n")
        out_file.write("@    s0 symbol fill color 1\n")
        out_file.write("@    s0 symbol fill pattern 1\n")
        
        out_file.write("@g1 on\n@g1 hidden false\n@with g1\n")
        out_file.write("@world " + str(xrng[0]) + ',' + str(y2rng[0]) + ',' + str(xrng[1]) + ',' + str(y2rng[1]) + '\n')
        out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
        out_file.write("@    xaxis label \"\"\n")
        out_file.write("@    xaxis tick off\n")
        out_file.write("@    xaxis tick major " + str(int_par_vals["kpnstep"]*3) + '\n')
        out_file.write("@    yaxis label \"Force (eV/Ang)\"\n")
        out_file.write("@    yaxis tick major " + str((y2rng[1]-y2rng[0])/8.0) + '\n')
        out_file.write("@    yaxis ticklabel format exponential\n")
        out_file.write("@    yaxis ticklabel prec 1\n")
        out_file.write("@    yaxis offset 1.0, 0.0\n")
        out_file.write("@    yaxis label place opposite\n")
        out_file.write("@    yaxis tick place opposite\n")
        out_file.write("@    yaxis ticklabel place opposite\n")
        out_file.write("@    s0 hidden false\n@    s0 on\n")
        out_file.write("@    s0 legend \"Max force\"\n")
        out_file.write("@    s0 line color 2\n")
        out_file.write("@    s0 symbol 1\n")
        out_file.write("@    s0 symbol size 0.7\n")
        out_file.write("@    s0 symbol color 2\n")
        out_file.write("@    s0 symbol fill color 2\n")
        out_file.write("@    s0 symbol fill pattern 1\n")
        
        # Input the actual data
        
        out_file.write("@target G0.S0\n@type xy\n")
        for i, k in enumerate(kpnrange):
            out_file.write(str(sum(k)) + '\t' + str(kpnnrg[i]) + '\n')
        out_file.write('&\n')
        
        out_file.write("@target G1.S0\n@type xy\n")
        for i, k in enumerate(kpnrange):
            out_file.write(str(sum(k)) + '\t' + str(kpnfor[i]) + '\n')
        out_file.write('&\n')
        
        out_file.close()
        
        if bool_par_vals["cnvstr"]:
            
            # K points vs energy and stresses
            
            out_file = open(seedname + "_kpn_str_conv.agr", 'w')
        
            xrng  = (sum(kpnrange[0]), sum(kpnrange[-1]))
            y1rng = (min(kpnnrg)-0.1*(max(kpnnrg)-min(kpnnrg)), max(cutnrg)+0.1*(max(kpnnrg)-min(kpnnrg)))
            y2rng = (min(kpnstr)-0.1*(max(kpnstr)-min(kpnstr)), max(kpnstr)+0.1*(max(kpnstr)-min(kpnstr)))
            
            # Set up the graphics
            
            out_file.write("@version 50123\n")
            out_file.write("@title \"" + seedname + " - Energy and forces vs k points\"\n")
            
            out_file.write("@g0 on\n@g0 hidden false\n@with g0\n")
            out_file.write("@world " + str(xrng[0]) + ',' + str(y1rng[0]) + ',' + str(xrng[1]) + ',' + str(y1rng[1]) + '\n')
            out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
            out_file.write("@    xaxis label \"k points\"\n")
            out_file.write("@    xaxis tick major " + str(int_par_vals["kpnstep"]*3) + '\n')
            out_file.write("@    xaxis offset 0.0, 1.0\n")
            out_file.write("@    yaxis label \"Final energy (eV)\"\n")
            out_file.write("@    yaxis tick major " + str((y1rng[1]-y1rng[0])/8.0) + '\n')
            out_file.write("@    yaxis offset 0.0, 1.0\n")
            out_file.write("@    s0 hidden false\n@    s0 on\n")
            out_file.write("@    s0 legend \"Final energy\"\n")
            out_file.write("@    s0 line color 1\n")
            out_file.write("@    s0 symbol 1\n")
            out_file.write("@    s0 symbol size 0.7\n")
            out_file.write("@    s0 symbol color 1\n")
            out_file.write("@    s0 symbol fill color 1\n")
            out_file.write("@    s0 symbol fill pattern 1\n")
            
            out_file.write("@g1 on\n@g1 hidden false\n@with g1\n")
            out_file.write("@world " + str(xrng[0]) + ',' + str(y2rng[0]) + ',' + str(xrng[1]) + ',' + str(y2rng[1]) + '\n')
            out_file.write("@    view 0.150000, 0.150000, 1.150000, 0.850000\n")
            out_file.write("@    xaxis label \"\"\n")
            out_file.write("@    xaxis tick off\n")
            out_file.write("@    xaxis tick major " + str(int_par_vals["kpnstep"]*3) + '\n')
            out_file.write("@    yaxis label \"Stress (GPa)\"\n")
            out_file.write("@    yaxis tick major " + str((y2rng[1]-y2rng[0])/8.0) + '\n')
            out_file.write("@    yaxis ticklabel format exponential\n")
            out_file.write("@    yaxis ticklabel prec 1\n")
            out_file.write("@    yaxis offset 1.0, 0.0\n")
            out_file.write("@    yaxis label place opposite\n")
            out_file.write("@    yaxis tick place opposite\n")
            out_file.write("@    yaxis ticklabel place opposite\n")
            out_file.write("@    s0 hidden false\n@    s0 on\n")
            out_file.write("@    s0 legend \"Max stress\"\n")
            out_file.write("@    s0 line color 2\n")
            out_file.write("@    s0 symbol 1\n")
            out_file.write("@    s0 symbol size 0.7\n")
            out_file.write("@    s0 symbol color 2\n")
            out_file.write("@    s0 symbol fill color 2\n")
            out_file.write("@    s0 symbol fill pattern 1\n")
            
            # Input the actual data
            
            out_file.write("@target G0.S0\n@type xy\n")
            for i, k in enumerate(kpnrange):
                out_file.write(str(sum(k)) + '\t' + str(kpnnrg[i]) + '\n')
            out_file.write('&\n')
            
            out_file.write("@target G1.S0\n@type xy\n")
            for i, k in enumerate(kpnrange):
                out_file.write(str(sum(k)) + '\t' + str(kpnstr[i]) + '\n')
            out_file.write('&\n')
            
            out_file.close()
        
    else:
        print "Only values currently supported for the output_type variable are gnuplot and [xm]grace"
