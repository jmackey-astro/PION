
Installing PyPion on Kay
===========================

+ Pull latest PION commit to *openmp* branch in directory `<BASE>/pion/`.
+ `$ cd pion/extra_libraries`
+ `$ bash install_all_libs.sh` should create a new sub-dir called `python` containing shared libs, including `Silo.so`
+ Pull latest PyPion commits into directory `<BASE>/pion_python`


Running the test simulation to the end
========================================

This only needs to be run once to generate the sequence of snapshots.

+ `$ cd pion/test_problems/OpenMP`
+ Edit file `kay.Ostar2_B010_d2l1n0128.txt` to use correct paths, accounts, usernames, etc.
+ `$ sbatch kay.Ostar2_B010_d2l1n0128.txt` to run the 2D simulation on 32 cores for approx 15 mins.


Plotting results
=====================

+ `$ module load intel/2018u4`
+ `$ module load conda`
+ `$ source activate`
+ `$ pip install --user astropy`
+ `$ mkdir img`
+ `$ python PlotData.py /ichec/work/EuroCC-AF-2/jm/Test2D/ Ostar2_d2l1n0256_np032_0000.001 ./img Ostar2_d2l1n0256 2D`
+ `$ display img/Ostar2_d2l1n0256_00117*.png`


Restarting from 2nd last snapshot
=================================

For e.g. running with fewer MPI processes:

+ Restart with single MPI process from 2nd last file with `$ sbatch kay.restart_2d_sim.txt`
+ Check out the log with e.g. `$ tail -f /ichec/work/EuroCC-AF-2/jm/Test2D/log_restart_test_0_info.txt`
+ Plot the new results `$ python PlotData.py /ichec/work/EuroCC-AF-2/jm/Test2D/ Ostar2_d2l1n0256_np001_0000 ./img TEST_np001 2D`
+ View new image with `$ display img/TEST_np001_00117*.png`

This takes about 3 mins with a single MPI process and no multithreading.


Compare two snapshots with "silocompare"
==========================================

+ `$ cd pion/analysis/silocompare`
+ `$ module load gsl/intel/2.5`
+ `$ module load intel/2018u4`
+ `$ MAKE_UNAME="KAY" make -j12`
+ `$ ./silocompare /ichec/work/EuroCC-AF-2/jm/Test2D/ Ostar2_d2l1n0256_np032_0000.00117  /ichec/work/EuroCC-AF-2/jm/Test2D/ Ostar2_d2l1n0256_np001_0000.00117 0 cmp 2`

This displays the diffs on screen, and saves a file "cmp.00000.silo" which is an image of the diff on a cell-by-cell basis.
The last command-line argument can be 0=absolute value of diff, 1=relative diff, 2=calculate L1,L2 error.
Options 0 and 1 will save a difference image, option 2 will not.



