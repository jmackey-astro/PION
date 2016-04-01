#!/bin/bash

mpirun -np 1 ../../icgen_parallel params_DTE_D2_TTI_n00016.txt silo redirect=iclog
mpirun -np 1 ../../pion_parallel IC_DTE_D2_TTI_n00016_0000.silo 5 1 op_criterion=0 opfreq=10 outfile=DTE_D2_TTI_n00016_np01 redirect=log_np01

mpirun -np 2 ../../icgen_parallel params_DTE_D2_TTI_n00016.txt silo redirect=iclog
mpirun -np 2 ../../pion_parallel IC_DTE_D2_TTI_n00016_0000.silo 5 1 op_criterion=0 opfreq=10 outfile=DTE_D2_TTI_n00016_np02 redirect=log_np02

mpirun -np 4 ../../icgen_parallel params_DTE_D2_TTI_n00016.txt silo redirect=iclog
mpirun -np 4 ../../pion_parallel IC_DTE_D2_TTI_n00016_0000.silo 5 1 op_criterion=0 opfreq=10 outfile=DTE_D2_TTI_n00016_np04 redirect=log_np04

mpirun -np 8 ../../icgen_parallel params_DTE_D2_TTI_n00016.txt silo redirect=iclog
mpirun -np 8 ../../pion_parallel IC_DTE_D2_TTI_n00016_0000.silo 5 1 op_criterion=0 opfreq=10 outfile=DTE_D2_TTI_n00016_np08 redirect=log_np08

mpirun -np 16 ../../icgen_parallel params_DTE_D2_TTI_n00016.txt silo redirect=iclog
mpirun -np 16 ../../pion_parallel IC_DTE_D2_TTI_n00016_0000.silo 5 1 op_criterion=0 opfreq=10 outfile=DTE_D2_TTI_n00016_np16 redirect=log_np16

cd ../../analysis/silocompare
MAKE_UNAME=standard make clean
MAKE_UNAME=standard make -j8
cd -

../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00016_np01_0000.00 DOUBLE DTE_D2_TTI_n00016_np02_0000.00 DOUBLE cmp_DTE_D2_TTI_n00016_np01np02 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00016_np01_0000.00 DOUBLE DTE_D2_TTI_n00016_np04_0000.00 DOUBLE cmp_DTE_D2_TTI_n00016_np01np04 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00016_np01_0000.00 DOUBLE DTE_D2_TTI_n00016_np08_0000.00 DOUBLE cmp_DTE_D2_TTI_n00016_np01np08 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00016_np01_0000.00 DOUBLE DTE_D2_TTI_n00016_np16_0000.00 DOUBLE cmp_DTE_D2_TTI_n00016_np01np16 2


