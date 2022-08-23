#!/bin/bash

../../../build/icgen-ng param_hydro_d1l2n0128.txt 
../../../build/pion-ng Galrot_test_1.700_0.0_d1l2n0128_level00_0000.00000000.silo 
../../../build/icgen-ng param_hydro_d1l3n0128.txt 
../../../build/pion-ng Galrot_test_1.700_0.0_d1l3n0128_level00_0000.00000000.silo 
../../../build/icgen-ng param_hydro_d1l4n0128.txt 
../../../build/pion-ng Galrot_test_1.700_0.0_d1l4n0128_level00_0000.00000000.silo 
#../../../build/icgen-ng param_hydro_d1l5n0128.txt 
#../../../build/pion-ng Galrot_test_1.700_0.0_d1l5n0128_level00_0000.00000000.silo 

../../../build/icgen-ug param_hydro_d1l1n0256.txt 
../../../build/pion-ug Galrot_test_1.700_0.0_d1l1n0256_0000.00000000.silo 
../../../build/icgen-ug param_hydro_d1l1n0512.txt 
../../../build/pion-ug Galrot_test_1.700_0.0_d1l1n0512_0000.00000000.silo 




python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l2n0128 img
python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l3n0128 img
python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l4n0128 img
python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l5n0128 img
python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l1n0256 img
python3 plot_ng-1D_state.py silo/ Galrot_test_1.700_0.0_d1l1n0512 img




