#!/bin/bash

#################################################################
### 2011.04.18 JM: Comparing new and old radiative transfer.  ###
#################################################################

test_dir=$1
exe_dir=$2
code_dir=$3
data_dir=$4
SUFFIX=$5

#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_RT_nodyn
#exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp


cat << EOF  > gnu.plt
#set terminal postscript enhanced eps
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.003:0.3]
set xrange [0:10]
set title ""
set output "${test_dir}/COMP_photoncons_norec2d_nh1.eps"
set label "2D, n_H=10 cm^{-3}" at 2,0.08
plot '${data_dir}/rtt2D_n32_nh1_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 ps 0.5 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh1_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 ps 0.5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n32_nh1_norec_${SUFFIX}.txt'  u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 4 title '32^2, ${SUFFIX}', \
     '${data_dir}/rtt2D_n100_nh1_norec_${SUFFIX}.txt' u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 5 title '100^2, ${SUFFIX}'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_norec2d_nh1.eps ${test_dir}/COMP_photoncons_norec2d_nh1.jpeg

cat << EOF  > gnu.plt
#set terminal postscript enhanced eps
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.003:0.3]
set xrange [0:10]
set title ""
set output "${test_dir}/COMP_photoncons_norec2d_nh2.eps"
set label "2D, n_H=100 cm^{-3}" at 2,0.08
plot '${data_dir}/rtt2D_n32_nh2_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 ps 0.5 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh2_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 ps 0.5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n32_nh2_norec_${SUFFIX}.txt'  u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 4 title '32^2, ${SUFFIX}', \
     '${data_dir}/rtt2D_n100_nh2_norec_${SUFFIX}.txt' u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 5 title '100^2, ${SUFFIX}'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_norec2d_nh2.eps ${test_dir}/COMP_photoncons_norec2d_nh2.jpeg

cat << EOF  > gnu.plt
#set terminal postscript enhanced eps
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.003:0.3]
set xrange [0:10]
set title ""
set output "${test_dir}/COMP_photoncons_norec2d_nh3.eps"
set label "2D, n_H=1000 cm^{-3}" at 2,0.08
plot '${data_dir}/rtt2D_n32_nh3_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 ps 0.5 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh3_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 ps 0.5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n32_nh3_norec_${SUFFIX}.txt'  u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 4 title '32^2, ${SUFFIX}', \
     '${data_dir}/rtt2D_n100_nh3_norec_${SUFFIX}.txt' u (\$1/3.16e10):(abs(1-\$6/\$7))    w lp lt 3 pt 5 title '100^2, ${SUFFIX}'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_norec2d_nh3.eps ${test_dir}/COMP_photoncons_norec2d_nh3.jpeg


######################################
######## RECOMBINATIONS PLOTS ########
######################################

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.1:1.1]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_photoncons_rec2d_nh1.eps"
set label "2D, n_H=10 cm^{-3}" at 1,0.85
set label "Cell {/Symbol d}{/Symbol t}=1 for 100^2 runs" at 1,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=5.520e17        # lower for new model which includes 10% He
tH=3.861e11
tI=3.510e11
plot '${data_dir}/rtt2D_n32_nh1_rec_dt100.txt'  u (\$1/tH):(\$2/RH)    w lp lt 1 pt 1 ps 0.5 title '32^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh1_rec_dt100.txt' u (\$1/tH):(\$2/RH)    w lp lt 1 pt 2 ps 0.5 title '100^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n32_nh1_rec_${SUFFIX}.txt'  u (\$1/tI):(\$2/RI)    w lp lt 2 pt 4 title '32^2, {/Symbol d}t= ${SUFFIX}', \
     '${data_dir}/rtt2D_n100_nh1_rec_${SUFFIX}.txt' u (\$1/tI):(\$2/RI)    w lp lt 2 pt 5 title '100^2, {/Symbol d}t= ${SUFFIX}', \
    g(x) w l lt -1 title 'Theoretical Radius'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh1.eps ${test_dir}/COMP_photoncons_rec2d_nh1.jpeg

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.1:1.1]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_photoncons_rec2d_nh2.eps"
set label "2D, n_H=100 cm^{-3}" at 1,0.85
set label "Cell {/Symbol d}{/Symbol t}=1 for 100^2 runs" at 1,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=5.520e17        # lower for new model which includes 10% He
tH=3.861e10
tI=3.510e10
plot '${data_dir}/rtt2D_n32_nh2_rec_dt100.txt'  u (\$1/tH):(\$2/RH)    w lp lt 1 pt 1 ps 0.5 title '32^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh2_rec_dt100.txt' u (\$1/tH):(\$2/RH)    w lp lt 1 pt 2 ps 0.5 title '100^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n32_nh2_rec_${SUFFIX}.txt'  u (\$1/tI):(\$2/RI)    w lp lt 2 pt 4 title '32^2, {/Symbol d}t= ${SUFFIX}', \
     '${data_dir}/rtt2D_n100_nh2_rec_${SUFFIX}.txt' u (\$1/tI):(\$2/RI)    w lp lt 2 pt 5 title '100^2, {/Symbol d}t= ${SUFFIX}', \
    g(x) w l lt -1 title 'Theoretical Radius'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh2.eps ${test_dir}/COMP_photoncons_rec2d_nh2.jpeg


