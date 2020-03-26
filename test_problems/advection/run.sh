#!/bin/bash

for res in "l1n128" "l1n256" "l1n512"; do
  ../../icgen_parallel params_advection_v020t30_${res}.txt silo
  mpirun -np 4 ../../pion_parallel advection_v020_t30_${res}_0000.00000000.silo
done

for res in "l2n128" "l2n256"; do
  ../../icgen_NG_parallel params_advection_v020t30_${res}.txt silo
  mpirun -np 4 ../../pion_NG_parallel advection_v020_t30_${res}_level00_0000.00000000.silo
done

python plot_final_state.py 
exit

