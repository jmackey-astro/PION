#!/bin/bash


#for d in "B010_n0128l3" "B010_n0128l4" "B010_n0192l3" "B100_n0128l3" "B100_n0128l4" "B100_n0192l3"; do
#for d in "B100_n0128l3" "B100_n0128l4" "B100_n0192l3" "B030_n0128l3" "B030_n0128l4" "B030_n0192l3"; do
for d in "B100_n0128l4" "B100_n0192l3" "B030_n0128l3" "B030_n0128l4" "B030_n0192l3"; do
  mpirun -np 64 ../../icgen-ng params_Ostar3D_${d}.txt silo
  mpirun -np 64 ../../pion-ng Ostar3D_${d}_level00_0000.00000000.silo outfile=Ostar3D_${d} redirect=log_Ostar3D_${d} finishtime=3.156e11
done
exit


