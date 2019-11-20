#!/bin/bash


for d in "B010_n0128l3" "B010_n0128l4" "B010_n0192l3" "B100_n0128l3" "B100_n0128l4" "B100_n0192l3"; do
  mpirun -np 8 ../../icgen_NG_parallel params_Ostar3D_${d}.txt silo
  mpirun -np 8 ../../pion_NG_parallel Ostar3D_${d}_level00_0000.00000000.silo outfile=Ostar3D_${d} redirect=log_Ostar3D_${d} finishtime=3.156e11
done
exit


