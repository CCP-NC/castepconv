***********************************************************
CASTEPconv v. 0.9.3

by Simone Sturniolo
Copyright 2014 Science and Technology Facilities Council

***********************************************************

EXAMPLE RUN - Si

This folder contains an example case you can run to get started with CASTEPconv. The folder contains the following files:

*Si.cell
*Si.conv

The Si.cell file only contains the lattice information on silicon and doesn't need to be modified in any way.
The Si.conv file contains only two instructions - task (set to 'all') and running_command. The current value assumes that on your system the command to run castep is castep.serial - if it isn't, modify it before running the test.
A number of other lines are provided and commented out (being preceeded by a '#'). Remove the # to make use of them and experiment with different settings.
Run the job simply by entering the folder and typing:

castepconv.py Si

The jobs should be easy enough to run fairly quickly on any average modern computer. If you want to delete everything and have a clean start just use the command:

castepconv.py -t c Si

which will run the 'clean' task.