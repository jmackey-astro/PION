#!/bin/bash

for res in "l1n128" "l1n256" "l1n512"; do
  ../../icgen-ug params_advection_v020t30_${res}.txt silo
  mpirun -np 4 ../../pion-ug advection_v020_t30_${res}_0000.00000000.silo
done

for res in "l2n128" "l2n256", "l2n512"; do
  ../../icgen-ng params_advection_v020t30_${res}.txt silo
  mpirun -np 4 ../../pion-ng advection_v020_t30_${res}_level00_0000.00000000.silo
done

python plot_final_state.py 
exit

