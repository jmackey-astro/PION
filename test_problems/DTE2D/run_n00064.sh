#!/bin/bash


mpirun -np  1 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np01 redirect=logOA2_np01 &
mpirun -np  2 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np02 redirect=logOA2_np02 &
mpirun -np  4 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np04 redirect=logOA2_np04 &
mpirun -np  8 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np08 redirect=logOA2_np08 &
mpirun -np 16 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np16 redirect=logOA2_np16 &
mpirun -np 32 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np32 redirect=logOA2_np32 &
wait

mpirun -np 64 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 outfile=OA2_np64 redirect=logOA2_np64

cd ../../analysis/silocompare
MAKE_UNAME=ubuntu16 make clean
MAKE_UNAME=ubuntu16 make -j8
cd -

../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np02_level00_0000.00 0 cmp_OA2_np01np02 2
../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np04_level00_0000.00 0 cmp_OA2_np01np04 2
../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np08_level00_0000.00 0 cmp_OA2_np01np08 2
../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np16_level00_0000.00 0 cmp_OA2_np01np16 2
../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np32_level00_0000.00 0 cmp_OA2_np01np32 2
../../analysis/silocompare/silocompare . OA2_np01_level00_0000.00 \
 . OA2_np64_level00_0000.00 0 cmp_OA2_np01np64 2

exit

mpirun -np  1 ../../icgen_NG_parallel params_DTE_NG_D2_TTI_n00064.txt silo redirect=iclog

mpirun -np  1 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np01 redirect=log_OA1_np01 &
mpirun -np  2 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np02 redirect=log_OA1_np02 &
mpirun -np  4 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np04 redirect=log_OA1_np04 &
mpirun -np  8 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np08 redirect=log_OA1_np08 &
mpirun -np 16 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np16 redirect=log_OA1_np16 &
mpirun -np 32 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np32 redirect=log_OA1_np32 &
wait

mpirun -np 64 ../../pion_NG_parallel DTE_NG_D2_TTI_n00064_level00_0000.00000000.silo \
 ooa=1 outfile=OA1_np64 redirect=log_OA1_np64

cd ../../analysis/silocompare
MAKE_UNAME=ubuntu16 make clean
MAKE_UNAME=ubuntu16 make -j8
cd -

../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np02_level00_0000.00 0 cmp_OA1_np01np02 2
../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np04_level00_0000.00 0 cmp_OA1_np01np04 2
../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np08_level00_0000.00 0 cmp_OA1_np01np08 2
../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np16_level00_0000.00 0 cmp_OA1_np01np16 2
../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np32_level00_0000.00 0 cmp_OA1_np01np32 2
../../analysis/silocompare/silocompare . OA1_np01_level00_0000.00 \
 . OA1_np64_level00_0000.00 0 cmp_OA1_np01np64 2

exit

