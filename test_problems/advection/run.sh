#!/bin/bash

for res in "l1n128" "l1n256" "l1n512"; do
  ../../icgen-ug params_advection_v020t30_${res}.txt silo
done
mpirun -np  4 ../../pion-ug advection_v020_t30_l1n128_0000.00000000.silo
mpirun -np 16 ../../pion-ug advection_v020_t30_l1n256_0000.00000000.silo
mpirun -np 64 ../../pion-ug advection_v020_t30_l1n512_0000.00000000.silo

for res in "l2n128" "l2n256", "l2n512"; do
  ../../icgen-ng params_advection_v020t30_${res}.txt silo
done
mpirun -np  4 ../../pion-ng advection_v020_t30_l2n128_level00_0000.00000000.silo
mpirun -np 16 ../../pion-ng advection_v020_t30_l2n256_level00_0000.00000000.silo
mpirun -np 64 ../../pion-ng advection_v020_t30_l2n512_level00_0000.00000000.silo

#python plot_final_state.py 
exit

