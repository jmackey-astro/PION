#!/bin/bash

../../icgen-ng params_DTE_NG_D2_TTI_n0128.txt silo
mpirun -np  1 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np01oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
mpirun -np  2 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np02oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
mpirun -np  4 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np04oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
mpirun -np  8 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np08oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
mpirun -np 16 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np16oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
mpirun -np 32 ../../pion-ng DTE_NG_D2_TTI_n0128_level00_0000.00000000.silo outfile=np32oa2 solver=8 cfl=0.4 AVtype=1 EtaVisc=0.15 &
wait
exit

for RES in "0032" "0064" "0128" "0256"; do
#for RES in "0256"; do
  ../../icgen-ng params_DTE_NG_D2_TTI_n0${RES}.txt silo
  mpirun -np 8 ../../pion-ng \
      DTE_NG_D2_TTI_n0${RES}_level00_0000.00000000.silo \
      outfile=tests/DTE_NG_D2_TTI_n0${RES} \
      redirect=log_DTE_NG_D2_TTI_n0${RES}
done

