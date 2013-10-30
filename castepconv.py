#!/usr/bin/env python

# CASTEP convergence automation tool
# by Simone Sturniolo
#
# Copyright 2013 Science and Technology Facilities Council
# This software is distributed under the terms of the GNU General Public License (GNU GPL)

import sys, time, math, os
import subprocess as sp

class ConvError(Exception):
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

# Length units allowed by CASTEP. Internally used unit is always Angstroms 

length_units = {'ang':1.0, 'm':1.0e10, 'cm':1.0e8, 'nm':10.0, 'bohr':0.529, 'a0':0.529}

# Parameters from CONV files - names and values

str_par_names = {
"convergence_task": "ctsk",
"running_mode": "rmode",
"output_type":  "outp",
"running_command": "rcmd",
}

str_par_vals = {
"ctsk"   : "input",                         # Can be INPUT, INPUTRUN, OUTPUT or ALL
"rmode"  : "parallel",                      # Can be PARALLEL or SERIAL
"outp"   : "gnuplot",                       # Can be GNUPLOT or XMGRACE
"rcmd"   : "castep <seedname> -dryrun",
}

float_par_names = {
"cutoff_min": "cutmin",
"cutoff_max": "cutmax",
"cutoff_step": "cutstep",
"displace_atoms": "displ",
"final_energy_delta": "nrgtol",
"forces_delta": "fortol",
"stresses_delta": "strtol"
}

float_par_vals = { 
"cutmin" : 400.0,           # eV
"cutmax" : 800.0,           # eV
"cutstep": 100.0,           # eV
"displ"  : 0.0,             # Ang
"nrgtol" : 0.01,            # eV/atom
"fortol" : 0.001,           # eV/Ang
"strtol" : 0.01             # GPa
}

int_par_names = {
"kpoint_n_min": "kpnmin",
"kpoint_n_max": "kpnmax",
"kpoint_n_step": "kpnstep",
}

int_par_vals = {
"kpnmin"  : 1,
"kpnmax"  : 4,
"kpnstep" : 1,
}

bool_par_names = {
"converge_stress" : "cnvstr"
}

bool_par_vals = {
"cnvstr" : False
}

# CELL and PARAM files

cellfile_lines = None
paramfile_lines = None

cutrange = None
kpnrange = None

abc_len = None
kpn_base = (1, 1, 1)

__CASTEP_HEADER__       = "+-------------------------------------------------+"
__CASTEP_TIME__         = "Total time          ="
__CASTEP_ATOMN__        = "Total number of ions in cell = "
__CASTEP_CUTOFF__       = "plane wave basis set cut-off                   :"
__CASTEP_KPOINTS__      = "MP grid size for SCF calculation is"
__CASTEP_ENERGY__       = "Final energy, E             ="
__CASTEP_FORCES__       = "***************** Symmetrised Forces *****************"
__CASTEP_FORCES_END__   = "*                                                    *"
__CASTEP_STRESSES__     = "*********** Symmetrised Stress Tensor ***********"
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
            cline = cline.lower().split(':')
        if (len(cline) < 2):
            raise ConvError("Bad formatting in .conv file at line " + str(i))
        par_name = cline[0].strip()
        if (par_name in str_par_names):
            str_par_vals[str_par_names[par_name]] = cline[1].strip()
        elif (par_name in float_par_names):
            float_par_vals[float_par_names[par_name]] = float(cline[1])
        elif (par_name in int_par_names):
            int_par_vals[int_par_names[par_name]] = int(cline[1])
        elif (par_name in bool_par_names):
            bool_par_vals[bool_par_names[par_name]] = (cline[1].lower().strip() == "true")
        else:
            raise ConvError("Unrecognized option in .conv file at line " + str(i))
        
# Strip .cell file from unnecessary lines to get only the ones we care for (i.e. remove all reference to kpoints, we're going to put those in ourselves)
# Also read cell parameters and construct the proper kpn_base

def strip_cellfile(clines):
    
    stripped = []
    
    to_strip = False
    to_read_abc = False
    to_read_cart = False
    
    abc = None
    u = 1.0
    kbase = (1, 1, 1)
    
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
                        raise CellError('Bad formatting in .cell file LATTICE_ABC block')
                    to_read_cart = (min(abc) < 0.0)
            if "kpoints_mp_grid" in l_low:
                continue
            if "kpoints_mp_spacing" in l_low:
                continue
            if "%block" in l_low:
                if "kpoint_list" in l_low:
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
            if to_strip:
                if "%endblock" in l_low:
                    to_strip = False
                continue
        
        stripped.append(l)
        
    if abc is None:
        raise CellError('.cell file does not contain a LATTICE_* block')
    
    max_e = max(abc)
    kbase = tuple([int(max_e/x) for x in abc])
    
    return stripped, tuple(abc), kbase

# Strip .param file from unnecessary lines to get only the ones we care for (i.e. remove all reference to task and cutoff)

def strip_paramfile(plines):
    
    stripped = []
    
    for l in plines:
        
        l_split = l.strip().lower().split()
        
        # Skip empty lines
        
        if len(l_split) == 0:
            continue
        if "task" in l_split[0]:
            continue
        if "cut_off_energy" in l_split[0]:
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
        if not os.path.isfile(foldname + r'/' + jobname + '.castep'):
            time.sleep(1)
            continue
        # In case of errors, abort
        if os.path.isfile(foldname + r'/' + jobname + '.0001.err'):
            is_finished = True
            break
        cast_file = open(foldname + r'/' + jobname + '.castep', 'r')
        cast_lines = cast_file.readlines()
        cast_lines.reverse()
        cast_file.close()
        if cast_lines is None:
            time.sleep(1)
            continue
        for l in cast_lines:
            if __CASTEP_HEADER__ in l:
                break
            elif __CASTEP_TIME__ in l:
                is_finished = True
                break
        time.sleep(1)

# Read forces and sum their modules up in .castep file

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

def parse_stresses(cfile):
    
    global __CASTEP_STRESSES_END__
    
    tot_stress = 0.0
    
    for l in cfile[6:]:
        
        if __CASTEP_STRESSES_END__ in l:
            return tot_stress/3.0
        
        try:
            tot_stress += math.sqrt(sum([float(x)**2.0 for x in l.split()[2:5]]))
        except ValueError:
            pass
    
    raise CastepError("Corrupted stresses block")

###### -- MAIN PROGRAM -- ######

if (len(sys.argv) < 2):
    sys.exit("Not enough arguments provided - seedname is required")

if (sys.version_info[0] < 2 or sys.version_info[1] < 6):
    sys.exit("ERROR - Python version 2.6 or higher required to run the script")

seedname = sys.argv[1]

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
        sys.exit(".cell file for job " + seedname + " not found")
    
    print "Reading " + seedname + ".param"

    try:
        job_paramfile = open(seedname + ".param", 'r')
        paramfile_lines = job_paramfile.readlines()
    except IOError:
        print("WARNING: .param file for job " + seedname + " not found")
        paramfile_lines = None
    
else:
    
    cellfile_lines = None
    paramfile_lines = None
    
# PHASE 1 - INPUT
# Produce a batch of folders and files with the various configurations

if (str_par_vals['ctsk'] not in ("input", "inputrun", "output", "all")):
    sys.exit("ERROR - Invalid convergence_task parameter")

if (str_par_vals['ctsk'] in ("input", "inputrun", "all")):
    
    # First create a "stripped" version of cell and param files, to use as a template for the new files to generate
    
    try:
        stripped_cell, abc_len, kpn_base = strip_cellfile(cellfile_lines)
    except CellError as CE:
        sys.exit("ERROR - " + str(CE))
    if paramfile_lines is not None:
        stripped_param = strip_paramfile(paramfile_lines)
    else:
        stripped_param = []
    
    # Apply displacements to .cell atoms if needed to have non-zero forces
    
    if float_par_vals["displ"] > 0.0:
        print "Displacing atoms in .cell file of " + str(float_par_vals["displ"]) + " Angstroms"
        try:
            stripped_cell = displace_cell_atoms(stripped_cell, abc_len, float_par_vals["displ"])
        except CellError as CE:
            sys.exit("ERROR - " + str(CE))
            
    # Open a .conv_tab file to keep track of the created files and folders. Will be read if output is done as a separate operation
    
    conv_tab_file = open(seedname + ".conv_tab", 'w')
    
    if (float_par_vals["cutmin"] <= 0.0 or float_par_vals["cutstep"] <= 0.0 or float_par_vals["cutmax"] < float_par_vals["cutmin"]):
        sys.exit("ERROR - Invalid cutoff range defined in .conv file")
    if (int_par_vals["kpnmin"] <= 0 or int_par_vals["kpnstep"] <= 0 or int_par_vals["kpnmax"] < int_par_vals["kpnmin"]):
        sys.exit("ERROR - Invalid k-points range defined in .conv file")
    
    cut_n = int(math.ceil(float_par_vals["cutmax"]-float_par_vals["cutmin"])/float_par_vals["cutstep"])+1
    cutrange = [float_par_vals["cutmin"] + i * float_par_vals["cutstep"] for i in range(0, cut_n)]
        
    kpn_n = int(math.ceil(int_par_vals["kpnmax"]-int_par_vals["kpnmin"])/int_par_vals["kpnstep"])+1
    kpnrange = [tuple([(int_par_vals["kpnmin"] + i * int_par_vals["kpnstep"]) * e for e in kpn_base]) for i in range(0, kpn_n)]
    
    ovwrite_files = False
    
    if str_par_vals["rmode"] == "parallel":
        
        print "Creating folders for parallel convergence run"
        
        conv_tab_file.write("cutoff:\t")
        
        for i, cut in enumerate(cutrange):
            
            foldname = seedname + "_cut_" + str(i+1) + "_kpn_1"
            
            print "Creating folder " + foldname
            
            if not os.path.exists(foldname): 
                os.makedirs(foldname)
            elif not ovwrite_files:
                to_del = raw_input("Warning: folder " + foldname + " already exists. \
                \nSome files might be overwritten. Continue (y/N/y-all)?")
                if to_del == 'N':
                    sys.exit("Aborting")
                elif to_del == 'y-all':
                    ovwrite_files = True
            
            icell = open(foldname + r'/' + foldname + '.cell', 'w')
            iparam = open(foldname + r'/' + foldname + '.param', 'w')
            
            for l in stripped_cell: 
                icell.write(l)        
            icell.write("\nkpoint_mp_grid " + kgrid(kpnrange[0]) + "\n")
            icell.close()
            
            iparam.write("task:\tSinglePoint\n")
            iparam.write("cut_off_energy:\t" + str(cut) + " eV\n")
            for l in stripped_param: 
                iparam.write(l)
            iparam.close()
            
            conv_tab_file.write(str(cut) + " eV\t")
        
        conv_tab_file.write("\n")
            
        conv_tab_file.write("kpoint_n:\t" + kgrid(kpnrange[0]) + "\t|\t")
        
        for i, kpn in enumerate(kpnrange[1:]):
            
            foldname = seedname + "_cut_1_kpn_" + str(i+2)
            
            print "Creating folder " + foldname
            
            if not os.path.exists(foldname): 
                os.makedirs(foldname)
            elif not ovwrite_files:
                to_del = raw_input("Warning: folder " + foldname + " already exists. \
                \nSome files might be overwritten. Continue (y/N/y-all)?")
                if to_del == 'N':
                    sys.exit("Aborting")
                elif to_del == 'y-all':
                    ovwrite_files = True
                        
            icell = open(foldname + r'/' + foldname + '.cell', 'w')
            iparam = open(foldname + r'/' + foldname + '.param', 'w')
            
            for l in stripped_cell: 
                icell.write(l)        
            icell.write("\nkpoint_mp_grid " + kgrid(kpn) + "\n")
            icell.close()
            
            iparam.write("task:\tSinglePoint\n")
            iparam.write("cut_off_energy:\t" + str(cutrange[0]) + " eV\n")
            for l in stripped_param: 
                iparam.write(l)
            iparam.close()
            
            conv_tab_file.write(kgrid(kpn) + "\t|\t")
        
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
            if to_del == 'N':
                sys.exit("Aborting")
            elif to_del == 'y-all':
                ovwrite_files = True
                
        conv_tab_file.write("cutoff:\t")
        
        prev_jobname = None
        
        for i, cut in enumerate(cutrange):
            
            jobname = seedname + "_cut_" + str(i+1) + "_kpn_1"
            
            print "Creating files for job " + jobname
            
            icell = open(foldname + r'/' + jobname + '.cell', 'w')
            iparam = open(foldname + r'/' + jobname + '.param', 'w')
            
            for l in stripped_cell: 
                icell.write(l)        
            icell.write("\nkpoint_mp_grid " + kgrid(kpnrange[0]) + "\n")
            icell.close()
            
            iparam.write("task:\tSinglePoint\n")
            iparam.write("cut_off_energy:\t" + str(cut) + " eV\n")
            for l in stripped_param: 
                iparam.write(l)
            if prev_jobname is not None:
                iparam.write("reuse:\t" + prev_jobname + ".check")
            iparam.close()
            
            prev_jobname = jobname
            
            conv_tab_file.write(str(cut) + " eV\t")
        
        conv_tab_file.write("\n")
            
        conv_tab_file.write("kpoint_n:\t" + kgrid(kpnrange[0]) + "\t|\t")
        
        for i, kpn in enumerate(kpnrange[1:]):
            
            jobname = seedname + "_cut_1_kpn_" + str(i+2)
            
            print "Creating files for job " + jobname
            
                        
            icell = open(foldname + r'/' + jobname + '.cell', 'w')
            iparam = open(foldname + r'/' + jobname + '.param', 'w')
            
            for l in stripped_cell: 
                icell.write(l)
            icell.write("\nkpoint_mp_grid " + kgrid(kpn) + "\n")
            icell.close()
            
            iparam.write("task:\tSinglePoint\n")
            iparam.write("cut_off_energy:\t" + str(cutrange[0]) + " eV\n")
            for l in stripped_param: 
                iparam.write(l)
            if prev_jobname is not None:
                iparam.write("reuse:\t" + prev_jobname + ".check")
            iparam.close()
            
            prev_jobname = jobname
            
            conv_tab_file.write(kgrid(kpn) + "\t|\t")
        
        conv_tab_file.close()
    else:
        sys.exit("ERROR - Invalid value for the running_mode parameter")
    
# PHASE 2 - EXECUTION
# Run through the files (or create them as we go if the task is serial) and obtain the results

if (str_par_vals["ctsk"] in ("all", "inputrun")):
    
    if str_par_vals["rmode"] == "parallel":
        
        print "Running parallel convergence jobs"
        
        for i, cut in enumerate(cutrange):
            
            foldname = seedname + "_cut_" + str(i+1) + "_kpn_1"
            
            print "Running job " + foldname
            
            os.chdir(foldname)
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(foldname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = foldname
            # Note: subprocess.Popen opens a subprocess without waiting for it to finish.
            # In fact, if we don't take care to check that all files are closed (see later) they might as well not be and we might end with a crash
            sp.Popen(cmd_line)
            os.chdir("..")
                    
        for i, kpn in enumerate(kpnrange[1:]):
            
            foldname = seedname + "_cut_1_kpn_" + str(i+2)
            
            print "Running job " + foldname
            
            os.chdir(foldname)
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(foldname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = foldname
            sp.Popen(cmd_line)
            os.chdir("..")
        
        # If we're running an ALL job, wait for everyone to finish
        
        if str_par_vals["ctsk"] == "all":
            
            print "Waiting for all jobs to finish \
            WARNING: The program can be terminated with Ctrl+C, but that could terminate also the running jobs"
            
            for i, cut in enumerate(cutrange):
                foldname = seedname + "_cut_" + str(i+1) + "_kpn_1"
                jobfinish_wait(foldname, foldname)
                
            for i, kpn in enumerate(kpnrange[1:]):
                foldname = seedname + "_cut_1_kpn_" + str(i+2)
                jobfinish_wait(foldname, foldname)
            
            print "All jobs finished. Proceeding to analyze output"
        
    elif str_par_vals["rmode"] == "serial":
        
        print "Running serial convergence jobs"
        
        foldname = seedname + "_conv"
        
        for i, cut in enumerate(cutrange):
            
            jobname = seedname + "_cut_" + str(i+1) + "_kpn_1"
            
            print "Running job with " + str(cut) + " eV cutoff"
            
            os.chdir(foldname)
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(jobname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = jobname
            sp.Popen(cmd_line)
            os.chdir("..")
            
            jobfinish_wait(foldname, jobname)
                    
        for i, kpn in enumerate(kpnrange[1:]):
            
            jobname = seedname + "_cut_1_kpn_" + str(i+2)
            
            print "Running job with kpoint grid " + kgrid(kpn)
            
            os.chdir(foldname)
            cmd_line = str_par_vals["rcmd"].split()
            if not "<seedname>" in cmd_line:
                cmd_line.append(jobname)
            else:
                for j, l in enumerate(cmd_line):
                    if l == "<seedname>":
                        cmd_line[j] = jobname
            sp.Popen(cmd_line)
            os.chdir("..")
            
            jobfinish_wait(foldname, jobname)
        
        if str_par_vals["ctsk"] == "all":
            print "All jobs finished. Proceeding to analyze output"
        
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
            ctab_file = open(seedname + ".conv_tab", 'r').readlines()
            if len(ctab_file) == 2:
                try:
                    cutline = ctab_file[0].split(':')[1].split('eV')
                    kpnline = ctab_file[1].split(':')[1].split('|')
                    if (cutline is not None and kpnline is not None):
                        cutrange = [float(cut.strip()) for cut in cutline[:-1]]
                        kpnrange = [tuple([int(k) for k in kpn.strip().split()]) for kpn in kpnline[:-1]]
                except IndexError:
                    cutline = None
                    kpnline = None
    
    calc_str = bool_par_vals["cnvstr"]
    
    cutnrg = []
    cutfor = []
    cutstr = []
    kpnnrg = []
    kpnfor = []
    kpnstr = []
    
    # Try opening the various .castep files and collect energy and forces
    
    for i, cut in enumerate(cutrange):
        
        jobname = seedname + "_cut_" + str(i+1) + "_kpn_1"
        
        if (str_par_vals["rmode"] == "serial"):
            foldname = seedname + "_conv"
        else:
            foldname = jobname
                
        filepath = foldname + r'/' + jobname + '.castep'
        
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
                elif __CASTEP_ENERGY__ in l:
                    i_nrg = float(l.split()[4])
                elif __CASTEP_FORCES__ in l:
                    i_for = parse_forces(castepfile[start_l+j:])
                elif calc_str and __CASTEP_STRESSES__ in l:
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
        
        jobname = seedname + "_cut_1_kpn_" + str(i+2)
        
        if (str_par_vals["rmode"] == "serial"):
            foldname = seedname + "_conv"
        else:
            foldname = jobname
        
        filepath = foldname + r'/' + jobname + '.castep'
        
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
                elif __CASTEP_ENERGY__ in l:
                    i_nrg = float(l.split()[4])
                elif __CASTEP_FORCES__ in l:
                    i_for = parse_forces(castepfile[start_l+j:])
                elif calc_str and __CASTEP_STRESSES__ in l:
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
        
    print "\nConvergence results:"
    
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
                print "WARNING - Total stresses are lower than " + str(float_par_vals["strtol"]) + " GPa. A different atom displacement might be necessary to get meaningful results"
                        
            for i, stress in enumerate(cutstr[1:]):
                
                delta_str = abs(stress - delta_str)
                if delta_str < float_par_vals["strtol"]:
                    print "Based on converging total stresses to " + str(float_par_vals["fortol"]) + " GPa, minimum cutoff suggested is " + str(cutrange[i]) + " eV"
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
                print "WARNING - Total stresses are lower than " + str(float_par_vals["strtol"]) + " GPa. A different atom displacement might be necessary to get meaningful results"
                        
            for i, stress in enumerate(kpnstr[1:]):
                
                delta_str = abs(stress - delta_str)
                if delta_str < float_par_vals["strtol"]:
                    print "Based on converging total stresses to " + str(float_par_vals["strtol"]) + " GPa, minimum kpoint grid suggested is " + kgrid(kpnrange[i])
                    break
                
                delta_str = stress
                
                if i == len(kpnstr[1:])-1:
                    print "Unable to converge k-point grid with stresses within given range.  Try increasing kpoint_n_max in " + seedname + ".conv"
    
    if str_par_vals["outp"] == "gnuplot":
        
        out_file = open(seedname + "_conv.gp", 'w')
        
        out_file.write("set xlabel \"Cutoff (eV)\"\n")
        out_file.write("set ylabel \"Final energy (eV)\"\n")
        out_file.write("set y2label \"Total force (eV/A)\"\n")
        out_file.write("set ytics nomirror\n")
        out_file.write("set y2tics\n")
        out_file.write("plot \"" + seedname + "_cut_conv.dat\" using 1:2 with linespoints pt 7 ti \"Final energy\",")
        out_file.write("\"" + seedname + "_cut_conv.dat\" using 1:3 with linespoints pt 7 axes x1y2 ti \"Forces\"\n")
        out_file.write("pause -1 \"Hit return to continue\"\n")
        
        out_file.write("set xlabel \"k-points\"\n")
        out_file.write("set ylabel \"Final energy (eV)\"\n")
        out_file.write("set y2label \"Total force (eV/A)\"\n")
        out_file.write("plot \"" + seedname + "_kpn_conv.dat\" using 1:2 with linespoints pt 7 ti \"Final energy\",")
        out_file.write("\"" + seedname + "_kpn_conv.dat\" using 1:3 with linespoints pt 7 axes x1y2 ti \"Forces\"\n")
        out_file.write("pause -1 \"Hit return to continue\"\n")
        
        out_file.close()
    else:
        print "Only graphical output currently supported is gnuplot"
    
    
