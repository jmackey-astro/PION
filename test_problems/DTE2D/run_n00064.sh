#!/bin/bash

mpirun -np  1 ../../icgen_NG_parallel params_DTE_NG_D2_TTI_n00064.txt silo redirect=iclog

mpirun -np  1 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np01 redirect=logDTE_D2_TTI_n00064_np01 &
mpirun -np  2 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np02 redirect=logDTE_D2_TTI_n00064_np02 &
mpirun -np  4 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np04 redirect=logDTE_D2_TTI_n00064_np04 &
mpirun -np  8 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np08 redirect=logDTE_D2_TTI_n00064_np08 &
mpirun -np 16 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np16 redirect=logDTE_D2_TTI_n00064_np16 &
mpirun -np 32 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np32 redirect=logDTE_D2_TTI_n00064_np32 &
wait

mpirun -np 64 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=DTE_D2_TTI_n00064_np64 redirect=logDTE_D2_TTI_n00064_np64 &

cd ../../analysis/silocompare
MAKE_UNAME=ubuntu16 make clean
MAKE_UNAME=ubuntu16 make -j8
cd -

../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np02_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np02 2
../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np04_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np04 2
../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np08_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np08 2
../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np16_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np16 2
../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np32_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np32 2
../../analysis/silocompare/silocompare . DTE_D2_TTI_n00064_np01_level00_0000.00 \
 . DOUBLE DTE_D2_TTI_n00064_np64_level00_0000.00 0 cmp_DTE_D2_TTI_n00064_np01np64 2


