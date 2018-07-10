#!/bin/bash

../../icgen_nest_serial params_Wind2D_n0128.txt silo
../../pion_nest_serial fastwind_n0128_level00.00000000.silo outfile=OA2_l4 nlevels=4 ooa=2 opfreq=512

