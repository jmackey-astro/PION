#!/bin/bash

##################################
## Overstable Radiative Shocks  ##
## 2009-12-17                   ##
##################################
#
# JM 2009-12-17: 
# This file runs RadShockPosn to get the shock position as a function of time 
# to show oscillatory behaviour.  It then plots the position using gnuplot
# for various shock velocities.
#

test_dir=$1
data_dir=$2
plotfile=$3
#test_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_OverstableShock
#data_dir=/mnt/local/jm/temp_sims/code_test_dir
#plotfile=${test_dir}/OverstableShock

cd ${test_dir}
make -f Makefile.RadShockPosn

./RadShockPosn v100noch.txt ${data_dir}/RSH2D_n128_v100_nochem silo 0 50
./RadShockPosn v120noch.txt ${data_dir}/RSH2D_n128_v120_nochem silo 0 50
./RadShockPosn v130noch.txt ${data_dir}/RSH2D_n128_v130_nochem silo 0 50
./RadShockPosn v140noch.txt ${data_dir}/RSH2D_n128_v140_nochem silo 0 50
./RadShockPosn v150noch.txt ${data_dir}/RSH2D_n128_v150_nochem silo 0 50

./RadShockPosn v100chemV2.txt ${data_dir}/RSH2D_n256_v100_chemV2 silo 0 50
./RadShockPosn v120chemV2.txt ${data_dir}/RSH2D_n256_v120_chemV2 silo 0 50
./RadShockPosn v130chemV2.txt ${data_dir}/RSH2D_n256_v130_chemV2 silo 0 50
./RadShockPosn v140chemV2.txt ${data_dir}/RSH2D_n256_v140_chemV2 silo 0 50
./RadShockPosn v150chemV2.txt ${data_dir}/RSH2D_n256_v150_chemV2 silo 0 50


cat << EOF  > gnu.plt.noch
set terminal postscript enhanced color eps
set output "${plotfile}_noch.eps"
#set size 0.5, 0.5
set xrange [0:3.16e11]
set xlabel "Time (s)"
set yrange [0:1.0e17]
set ylabel "Shock Posn (cm)"
plot 'v100noch.txt' u 1:2 w l title "100 km/s", \
     'v120noch.txt' u 1:2 w l title "120 km/s", \
     'v130noch.txt' u 1:2 w l title "130 km/s", \
     'v140noch.txt' u 1:2 w l lw 2 title "140 km/s", \
     'v150noch.txt' u 1:2 w l lw 2 title "150 km/s"
#pause -1
EOF
gnuplot gnu.plt.noch

cat << EOF  > gnu.plt.chem
set terminal postscript enhanced color eps
set output "${plotfile}_chem.eps"
#set size 0.5, 0.5
set xrange [0:10]
set xlabel "Time (kyr)"
set yrange [0:40]
set ylabel "Shock Posn (milli-pc)"
plot 'v100chemV2.txt' u (\$1/3.16e10):(\$2/3.086e15) w l title "100 km/s", \
     'v120chemV2.txt' u (\$1/3.16e10):(\$2/3.086e15) w l title "120 km/s", \
     'v130chemV2.txt' u (\$1/3.16e10):(\$2/3.086e15) w l title "130 km/s", \
     'v140chemV2.txt' u (\$1/3.16e10):(\$2/3.086e15) w l title "140 km/s", \
     'v150chemV2.txt' u (\$1/3.16e10):(\$2/3.086e15) w l title "150 km/s"
#pause -1
EOF
gnuplot gnu.plt.chem

convert -density 300 -quality 100 ${plotfile}_noch.eps ${plotfile}_noch.jpeg
convert -density 300 -quality 100 ${plotfile}_chem.eps ${plotfile}_chem.jpeg

exit
