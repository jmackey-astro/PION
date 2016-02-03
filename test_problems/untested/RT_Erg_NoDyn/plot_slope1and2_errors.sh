#!/bin/bash

#################################################################
### 2011.04.18 JM: Comparing new and old radiative transfer.  ###
### 2011.06.22 JM: Updated for models with energetics.        ###
### 2011.06.23 JM: Added fractional error in radius plots.    ###
### 2011.07.09 JM: Adapted for 1D simulations (spherical)     ###
#################################################################

test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_Erg_NoDyn
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_Erg_NoDyn
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp

file1=rtt_Sph1D_ERG_n128_nh3_slope1_dE10dX01
file2=rtt_Sph1D_ERG_n128_nh3_slope1_dE10dX03
file3=rtt_Sph1D_ERG_n128_nh3_slope1_dX03
file4=rtt_Sph1D_ERG_n128_nh3_slope1_dX05
file5=rtt_Sph1D_ERG_n128_nh3_slope1_dX10
file6=rtt_Sph1D_ERG_n128_nh3_slope1_dX33
file7=rtt_Sph1D_ERG_n128_nh3_slope1_dX03HiAc

FNAME=${test_dir}/COMP1D_radius_error_nh3_slope1_n128

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set output "${FNAME}.eps"

set yrange [0.03:1]
set xrange [0.0003:10]
set log x
set log y
set key right bottom
set ylabel "R_{IF} (pc)" 2.5,0
set xlabel "Time (kyr)"
#set title 'Radial profile {/Symbol r} = 1/r'
plot '${data_dir}/${file1}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt -1 lw 2 title "N=128, High time acc.", \
     '${data_dir}/${file4}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 2  lw 2 title 'N=128, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/${file5}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 3  lw 2 title 'N=128, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/${file6}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 1  lw 2 title 'N=128, {/Symbol d}t= 1.00/xdot'
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg
#exit

file1=rtt_Sph1D_ERG_n128_nh3_slope2_dE10dX01
file2=rtt_Sph1D_ERG_n128_nh3_slope2_dE10dX03
file3=rtt_Sph1D_ERG_n128_nh3_slope2_dX03
file4=rtt_Sph1D_ERG_n128_nh3_slope2_dX05
file5=rtt_Sph1D_ERG_n128_nh3_slope2_dX10
file6=rtt_Sph1D_ERG_n128_nh3_slope2_dX33
file7=rtt_Sph1D_ERG_n128_nh3_slope2_dX03HiAc

FNAME=${test_dir}/COMP1D_radius_error_nh3_slope2_n128

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set output "${FNAME}.eps"

set yrange [0.03:1]
set xrange [0.0003:10]
set log x
set log y
set key right bottom
set ylabel "R_{IF} (pc)" 2.5,0
set xlabel "Time (kyr)"
#set title 'Radial profile {/Symbol r} = 1/r^{2}'
plot '${data_dir}/${file1}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt -1 lw 2 title "N=128, High time acc.", \
     '${data_dir}/${file4}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 2  lw 2 title 'N=128, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/${file5}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 3  lw 2 title 'N=128, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/${file6}.txt' u (\$1/3.17e10):(\$2/3.086e18) w l lt 1  lw 2 title 'N=128, {/Symbol d}t= 1.00/xdot'
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg
#exit

