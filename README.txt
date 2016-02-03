Author: Jonathan Mackey
Version: 1.0
Date: 2016.01.24

-------------------------------
--- IMPORTANT INFORMATION   ---
-------------------------------
Please read the licence file.  Using this code implies acceptance of
the terms of the licence for pion.  Sharing and/or publishing the
source code can also only be done according to the licence
conditions.

For more detailed information about using the code, see the Users Guide in
pion/docs/.

-------------------------------
---- COMPILING CODE:       ----
-------------------------------
Stages:
(1) compile external libraries (including manually downloading the
SUNDIALS library).
(2) run serial code compilation script.
(3) run parallel code compilation script.
(4) copy the executables to wherever you want to run the code.

(1)
First you need to compile extra libraries.
cd to "./extra_libraries"
run "bash install_all_libs.sh" to install.
You need to manually download the Sundials library (v.2.4.0 please)
from their website -- the install script will complain and tell you
the link to the library.
The script needs internet access with "wget" to work.

(2)
Then cd to "./bin_serial" and run "bash compile_code.sh"
if you are lucky it will detect that you are using a standard
workstation and will just compile.  If you are running Microsoft
Windows it almost certainly won't compile, and I have no expertise
to help you get it working, so "good luck".  On a UNIX/Linux
workstation or OS X (with x-code installed) it should compile ok.

If there are linking errors to silo (references to PDB and/or H5) or
fits (references to ff**** functions) or cvode/nvector then the 
libraries must have failed to compile.

If all went well, there should be no error messages and there should
be two executable files in ../ (i.e. the pion root directory):
 - "icgen_serial" for generating initial conditions, and
 - "pion_serial" for running the code.
Most code options are decided at run-time not compile-time, so you
shouldn't have to change any #define statements in the code.

If there are errors before linking, then there must be a problem with
the source code, maybe due to some c++ header files not being
installed.

Either way, probably email me at jmackey@astro.uni-bonn.de and I'll
see what I can do to help.  Quote the exact error message in the
email.

If you are not using Linux/UNIX then the Makefile and/or
compile_code.sh script may need editing.


(3)
cd to ../bin_parallel/ and run "bash compile_code.sh".
Again, if you are lucky it will just compile.  There will probably
be warning messages about pmpio.h from functions that are defined but
never used; this is not a problem.  If compilation is successful
then the files ../pion_parallel and ../icgen_parallel will
exist.

In the event that it doesn't compile follow the same procedure as for
the serial version of the code.


-------------------------------
----- RUNNING CODE:       -----
-------------------------------
----- TESTS:
-------------------------------

First you want to make sure the code tests run ok.  I'm still working
on this, but there is a subdirectory called
"uniform_grid_code/trunk/test_problems".  Here you can run
"./run_all_tests.sh" and it will run at least some tests, hopefully
without bugging out.  (some of the tests have directories hard-coded
to my desktop at AIfA).
First edit run_all_tests.sh so that the data_dir variable is set to 
some directory on your computer (and not my AIfA harddisk).

The Double-Mach-Reflection is a good test, and should run fine, and
if you have "eog" (eye-of-gnome) on your system, it will pop up jpeg
figures automatically which you can compare to
http://www.astro.uni-bonn.de/~jmackey/jmac/node10.html figures.  They
should look indistinguishable to the naked eye.  If not something is
definitely wrong.

Shock-tube tests will probably take a long time to run (a few hours),
so you may as well set that running overnight.

----------------------------------------
----- MORE COMPLICATED SIMULATIONS -----
----------------------------------------

cd to uniform_grid_code/trunk/bin
Here there are the two executables: "icgen_serial" and "pion_serial".

"icgen_serial" will read a parameter file and based on the values read
in, it will set up a grid with the appropriate geometry and size and
populate it with the requested data, also setting the correct boundary
conditions for the problem.  It then writes this to a initial
conditions file, which is identical to a restart/checkpointing file.

"pion_serial" will read an initial conditions file, or a
restart/checkpoint/output file (they are all the same) and
start/restart the simulation and run until it gets to "finishtime", at
which point it will output data and stop.  You can get a list of
command-line options by typing "./pion_serial" with no arguments.

A typical simulation run is as follows (copied from bin_serial/run.sh)
$ ./icgen_serial pf_test_winds.txt silo
$ ./pion_serial IC_wind_test1.silo 5 1 op_criterion=1 opfreq_time=1.58e10 cooling=0 redirect=test1 outfile=/export/aibn214_1/jmackey/testing/results/wind_test1 cfl=0.1

In icgen_serial, the first argument is the parameter-file, and the second
tells it to write a silo file (rather than fits).

for main serial, the arguments refer to:
REQUIRED ARGUMENTS:
(1) IC_wind_test1.silo	       initial conditions filename.
(2) 5			       That the IC file is silo format(fits=2)
(3) 1 			       That I'm using a uniform grid (always!)
OPTIONAL ARGUMENTS:
(4) op_criterion=1 	       output every n-time units (not n-steps)
(5) opfreq_time=1.58e10	       output frequency (in whatever units)
(6) cooling=0		       no cooling
(7) redirect=test1	       standard output to file ./test1info.txt
(8) outfile=/export/aibn214_1/jmackey/testing/results/wind_test1
			       output file name with path (will be
			       appended with step number)
(9) cfl=0.1		       Courant number (<1 in 1D, <0.5 2D, <0.35
    			       in 3D) and less with cooling etc.

---------------------
-- Parameter files --
---------------------
Have a look in source/ics/pfiles for some parameter
files.  Some of these are out of date, but if you look at the more
recent ones they should work ok.

You need to edit a parameter file for the problem you want to run.
Probably copy the template file from ics/pfiles/ to bin/ since
the templates are under version control and should only be updated if
they become obselete because of new code features.  a '#' at the start
of a line means it is ignored.

The most important parameter is "ics" which determines the correct
grid-setup routine to call.  Then down at the bottom there is a
section for "Parameters specific to a given problem"; these are the
physical properties of what you want to put into the simulation:
ambient medium density, dense clumps with some radius, shocks, etc.

The main way things fail for me is that I set the boundary conditions
string incorrectly.  "icgen_serial" doesn't check that it is correct
(it should, and will eventually).  There should be exactly two entries
for each dimension, and possibly "internal" extra boundaries.  If you
are running in 3D with only 4 boundaries the initial condition
generator will run, but pion_serial will bug out when it tries to
setup the grid.

The "BC" string is a sequence of 6-character boundary specifiers.  The
first two characters give the direction: XN=x-negative, XP=x-positive,
same for YN,YP,ZN,ZP, and IN=internal/special boundary.
The next three characters give the type of boundary: inflow, outflow,
fixed, one-way-outflow, periodic, reflecting.

"coordinates" = 
cylindrical (z,R) (2d only),
cartesian (x,y,z) (1d,2d,3d),
spherical (r)     (1d only)

"eqn" determines the equations to solve:
'hd' or 'euler' will choose the Euler equations.
'mhd' or 'glm-mhd' will choose the Ideal MHD equations.

"solver" is important: this is the solver to use for the flux
calculation.  For the Euler equations (no magnetic fields) you can use
any of 1-6.  I recommend 3, 4 or 6. (1 is bad, 2 is slow, 6 is pretty
good).
1=approximate linear Riemann solver
2=exact Riemann solver (slow!)
3=hybrid approximate/exact Riemann solver
4=Roe conserved variable Riemann solver
5=Roe primitive variable Riemann solver (may have problems!)
6=van Leer's flux vector splitting


-----------------------------
---- RUNNING IN PARALLEL ----
-----------------------------

The parallel compilation is in uniform_grid_code/trunk/bin_parallel
Again try "bash compile_code.sh", and the same caveats for OS X and
the Makefile apply for parallel as serial code.
Executables should be in trunk/bin/ called pion_parallel and
icgen_parallel. 

This should generate a lot of warnings about "PMPIO" but no errors.
If there are errors you should make sure that an MPI installation is
present on the system.  If you can compile and run an MPI 'hello
world' program, then something strange is going on and probably the
best thing is to contact me.

Parallel code should always be run with the "redirect=/path/to/file"
argument included, otherwise every process will output its info to the
console which is messy.  A typical example script is quoted below:

********************* run.sh ****************************
#!/bin/bash
# run code with 4 MPI processes.

sim_dir=/export/aibn129_1/jmackey/data_etc/stellar_winds/test_moving_src
mkdir ${sim_dir}
mkdir ${sim_dir}/run_log

mpirun -np 4 ./icgen_parallel pf_MSwind_Md1em6_v250_noI_adv000.txt silo redirect=${sim_dir}/run_log/ic_WIND_Md1em6_v250_noI_adv000
mpirun -np 4 ./icgen_parallel pf_MSwind_Md1em6_v250_noI_adv100.txt silo redirect=${sim_dir}/run_log/ic_WIND_Md1em6_v250_noI_adv100
rsync -vt IC_WIND_Md1em6_v250_noI_adv*.silo ${sim_dir}/
rm IC_WIND_Md1em6_v250_noI_adv*.silo

#
# First run the sims for a short time:
#
mpirun -np 4 ./pion_parallel ${sim_dir}/IC_WIND_Md1em6_v250_noI_adv000_0000.silo 5 1 \
outfile=${sim_dir}/WIND_Md1em6_v250_noI_adv000 \
redirect=${sim_dir}/run_log/msg_WIND_Md1em6_v250_noI_adv000 \
finishtime=3.16e11 opfreq_time=1.58e10 artvisc=0.15

mpirun -np 4 ./pion_parallel ${sim_dir}/IC_WIND_Md1em6_v250_noI_adv100_0000.silo 5 1 \
outfile=${sim_dir}/WIND_Md1em6_v250_noI_adv100 \
redirect=${sim_dir}/run_log/msg_WIND_Md1em6_v250_noI_adv100 \
finishtime=3.16e11 opfreq_time=1.58e10 artvisc=0.15

********************* run.sh ****************************

The first time you run parallel code there will probably be an error
about no "mpd" host running.  Run 'mpd&' and try again.

---------------------------------------------------------------------
Written by Jonathan Mackey (C) 2006-2013  jmackey@astro.uni-bonn.de
---------------------------------------------------------------------

