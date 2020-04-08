#!/bin/bash

for RES in "0032" "0064" "0128" "0256"; do
#for RES in "0256"; do
  ../../icgen-ng params_DTE_NG_D2_TTI_n0${RES}.txt silo
  mpirun -np 8 ../../pion-ng \
      DTE_NG_D2_TTI_n0${RES}_level00_0000.00000000.silo \
      outfile=tests/DTE_NG_D2_TTI_n0${RES} \
      redirect=log_DTE_NG_D2_TTI_n0${RES}
done

