Author: Jonathan Mackey
Version: 1.0
Date: 2017.09.14

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
(1) compile external libraries.
(2) run serial code compilation script, and/or
(3) run parallel code compilation script.
(4) copy the executables to wherever you want to run the code.

(1)
First you need to compile extra libraries.
On debian/Ubuntu Linux distributions you should be able to install the
packages:
libsundials-serial-dev libsundials-serial 
libsiloh5-0 libsilo-bin libsilo-dev
libcfitsio-dev libcfitsio-bin

Because of a recent silo update, you must copy the silo_exports.h file directly into usr/include/, and give the file appropriate permissions (I use 777 so I don't have to think about it) using chmod. You can get silo_exports.h by running extra_libraries, as below, and it should be in the silo library.

On other systems you may have to compile Silo, Sundials, and Cfitsio
from source.  Here is how:
- cd to "./extra_libraries"
- From the command line, run "bash install_all_libs.sh" to install.
The script needs internet access to download the source code for the
libraries.  It also needs cmake, c, c++, fortran compilers to be
installed on the system.


(2) ** compile SERIAL version of the code **
Then cd to "./bin_serial" and run "bash compile_code.sh"
if you are lucky it will detect that you are using a standard
workstation and will just compile.  If you are running Microsoft
Windows it almost certainly won't compile, and I have no expertise
to help you get it working, so "good luck".  On a UNIX/Linux
workstation or OS X (with x-code installed) it should compile ok.

If there are linking errors to silo (references to PDB and/or H5) or
fits (references to ff**** functions) or cvode/nvector then the 
libraries didn't link and may have failed to compile.

If all went well, there should be no error messages and there should
be two executable files in ../ (i.e. the PION root directory):
 - "icgen_serial" for generating initial conditions, and
 - "pion_serial" for running the code.
Most code options are decided at run-time not compile-time, so you
shouldn't have to change any #define statements in the code.

If there are errors before linking, then there must be a problem with
the source code, maybe due to some c++ header files not being
installed.

Either way, probably email me at jmackey@cp.dias.ie and I'll
see what I can do to help.  Quote the exact error message in the
email.  Or open an "issue" on the bitbucket PION pages:
https://bitbucket.org/jmackey/pion/issues

If you are not using Linux/UNIX then the Makefile and/or
compile_code.sh script may need editing.


(3) ** compile PARALLEL version of the code **
cd to ../bin_parallel/ and run "bash compile_code.sh".
Again, if you are lucky it will just compile.
If compilation is successful then the files ../pion_parallel and 
../icgen_parallel will exist.

In the event that it doesn't compile follow the same procedure as for
the serial version of the code.


-------------------------------
----- RUNNING CODE:       -----
-------------------------------
----- TESTS:
-------------------------------

First you want to make sure the code tests run ok.  I'm still working
on this, but there is a subdirectory called "test_problems".  Here you
can run some test problems with known solutions.

The Double-Mach-Reflection is a good test, and should run fine with:
 "bash run_double_Mach_reflection_test.sh <RES> <DIR>"
where:
 - <RES> is 130, 260, or 520, for the spatial resolution of the
grid that you want to use,
 - <DIR> is the path to the directory to save the data.  If <DIR> is
 left blank then it saves in data_YYYYMMDD/ in the current directory.

The script runs the test with the given grid resolution for a number
of different Riemann solvers, saving a few snapshots in <DIR>.

The script will generate jpeg figures that you can compare with
https://homepages.dias.ie/jmackey/jmac/node10.html figures.
Note that you need to install VisIt for this script to work.  We are
working on a python script to do the plotting instead.
The figs should look indistinguishable to the naked eye.  If not
something is probably wrong.

[There are currently a lot of test problems in the 
test_problems/untested directory; these used to work, but need some
work to get them working again.]

----------------------------------------
----- MORE COMPLICATED SIMULATIONS -----
----------------------------------------

cd to the pion root directory.
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

A typical simulation run is as follows:
$ ./icgen_serial params_test_winds.txt silo
$ ./pion_serial wind_test1.00000000.silo op_criterion=1 \
 opfreq_time=1.58e10 cooling=0 redirect=test1 \
 outfile=/path/to/results/wind_test1 cfl=0.1

In icgen_serial, the first argument is the parameter-file, and the second
tells it to write a silo file (rather than fits).

for main serial, the arguments refer to:
REQUIRED ARGUMENTS:
(1) wind_test1.00000000.silo   initial conditions filename (set by user!).
OPTIONAL ARGUMENTS:
(2) op_criterion=1 	       output every n-time units (not n-steps)
(3) opfreq_time=1.58e10	       output frequency (in code seconds)
(4) cooling=0		       no cooling
(5) redirect=test1	       redirect stdout to file ./test1info.txt
(6) outfile=/path/to/results/wind_test1
			       output file name with path (will be
			       appended with step number)
(7) cfl=0.1		       Courant number (<1 in 1D, <0.5 2D, <0.35
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
of a line means it is ignored.  Alternatively you can get parameter
files from the test_problems directory.

The most important parameter is "ics" which determines the correct
grid-setup routine to call.  Then down at the bottom there is a
section for "Parameters specific to a given problem"; these are the
physical properties of what you want to put into the simulation:
ambient medium density, dense clumps with some radius, shocks, etc.

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
2=exact Riemann solver (slow!) (HD only)
3=hybrid approximate/exact Riemann solver (HD only)
4=Roe conserved variable Riemann solver
5=Roe primitive variable Riemann solver (may have problems!)
6=van Leer's flux vector splitting (HD only)
7=HLLD (MHD only)

-----------------------------
---- RUNNING IN PARALLEL ----
-----------------------------

The parallel compilation is in bin_parallel/
Again try "bash compile_code.sh", and the same caveats for OS X and
the Makefile apply for parallel as serial code.
Executables should be in the PION root directory called pion_parallel
and icgen_parallel. 

If there are errors you should make sure that an MPI installation is
present on the system.  If you can compile and run an MPI 'hello
world' program, then something strange is going on and probably the
best thing is to contact me.

A typical example script is quoted below:

********************* run.sh ****************************
#!/bin/bash
# run code with 4 MPI processes.

sim_dir=/mnt/jmackey/data_etc/stellar_winds/test_moving_src
mkdir -p ${sim_dir}
mkdir -p ${sim_dir}/run_log

mpirun -np 4 ./icgen_parallel pf_MSwind_Md1em6_v250_noI_adv000.txt silo \
 redirect=${sim_dir}/run_log/ic_WIND_Md1em6_v250_noI_adv000
mv WIND_Md1em6_v250_noI_adv*.silo ${sim_dir}/

#
# Run the sims for a short time:
#
mpirun -np 4 ./pion_parallel \
 ${sim_dir}/WIND_Md1em6_v250_noI_adv000_0000.00000000.silo \
 outfile=${sim_dir}/WIND_Md1em6_v250_noI_adv000 \
 redirect=${sim_dir}/log_WIND_Md1em6_v250_noI_adv000 \
 finishtime=3.16e11 opfreq_time=1.58e10

********************* run.sh ****************************

---------------------------------------------------------------------
Written by Jonathan Mackey (C) 2006-2017  jmackey@cp.dias.ie
---------------------------------------------------------------------

