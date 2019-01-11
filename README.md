# castepconv

A script to run converge tests with CASTEP. To install, clone or download and then use pip:

    pip install ./castepconv --user
    
The `--user` option is not necessary but it's recommended as it will not require admin access to the system to work.
To test, having CASTEP installed on your system, visit the directory `example` and run:

    castepconv.py Si

This will use the default command `castep.serial`. If your system uses a custom command for CASTEP, edit the file `Si.conv`
to replace it. A full guide on how to use it and how to run the example is found in the PDF manual.

## CASTEPconv 2.0

This new version has minimal changes in terms of interface, but a deep reworking of the code itself. The main changes are as follow:
* the "no dependencies" philosophy was dropped. CASTEPconv now depends on both Numpy and the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/). The latter provides the bulk of the interface to Castep and must be up to date (version 3.17 or higher).
* the code is now Python 3 compatible, in light of the imminent discontinuation of support for Python 2.
* CASTEPconv runs will now store input files and processed results in separated folders, allowing for easier retrieval or deletion.

Some new keywords also were introduced. Refer to the documentation for more details.