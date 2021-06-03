#!/bin/bash

mpirun -np  1 ../../icgen-ug params_DTE_D2_TTI_n00256.txt silo redirect=iclog
mpirun -np  1 ../../pion-ug IC_DTE_D2_TTI_n00256_0000.silo 5 1 outfile=DTE_D2_TTI_n00256_np01 redirect=logDTE_D2_TTI_n00256_np01 &

sleep 5

mpirun -np  2 ../../icgen-ug params_DTE_D2_TTI_n00256.txt silo redirect=iclog
mpirun -np  2 ../../pion-ug IC_DTE_D2_TTI_n00256_0000.silo 5 1 outfile=DTE_D2_TTI_n00256_np02 redirect=logDTE_D2_TTI_n00256_np02 &

sleep 5

mpirun -np  4 ../../icgen-ug params_DTE_D2_TTI_n00256.txt silo redirect=iclog
mpirun -np  4 ../../pion-ug IC_DTE_D2_TTI_n00256_0000.silo 5 1 outfile=DTE_D2_TTI_n00256_np04 redirect=logDTE_D2_TTI_n00256_np04 &

wait

mpirun -np  8 ../../icgen-ug params_DTE_D2_TTI_n00256.txt silo redirect=iclog
mpirun -np  8 ../../pion-ug IC_DTE_D2_TTI_n00256_0000.silo 5 1 outfile=DTE_D2_TTI_n00256_np08 redirect=logDTE_D2_TTI_n00256_np08

mpirun -np 16 ../../icgen-ug params_DTE_D2_TTI_n00256.txt silo redirect=iclog
mpirun -np 16 ../../pion-ug IC_DTE_D2_TTI_n00256_0000.silo 5 1 outfile=DTE_D2_TTI_n00256_np16 redirect=logDTE_D2_TTI_n00256_np16

cd ../../analysis/silocompare
MAKE_UNAME=standard make clean
MAKE_UNAME=standard make -j8
cd -

../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00256_np01_0000.00 DOUBLE DTE_D2_TTI_n00256_np02_0000.00 DOUBLE cmp_DTE_D2_TTI_n00256_np02_np01np02 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00256_np01_0000.00 DOUBLE DTE_D2_TTI_n00256_np04_0000.00 DOUBLE cmp_DTE_D2_TTI_n00256_np02_np01np04 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00256_np01_0000.00 DOUBLE DTE_D2_TTI_n00256_np08_0000.00 DOUBLE cmp_DTE_D2_TTI_n00256_np02_np01np08 2
../../analysis/silocompare/silocompare . . DTE_D2_TTI_n00256_np01_0000.00 DOUBLE DTE_D2_TTI_n00256_np16_0000.00 DOUBLE cmp_DTE_D2_TTI_n00256_np02_np01np16 2


