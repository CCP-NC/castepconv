# castepconv

A script to run converge tests with CASTEP. To install, clone or download and then use pip:

    pip install ./castepconv
    
To test, having CASTEP installed on your system, visit the directory `example` and run:

    castepconv.py Si

This will use the default command `castep.serial`. If your system uses a custom command for CASTEP, edit the file `Si.conv`
to replace it. A full guide on how to use it and how to run the example is found in the PDF manual.
