#!/bin/bash

#################################################################
### 2011.04.18 JM: Comparing new and old radiative transfer.  ###
### 2011.06.22 JM: Updated for models with energetics.        ###
### 2011.06.23 JM: Added fractional error in radius plots.    ###
#################################################################

#test_dir=$1
#exe_dir=$2
#code_dir=$3
#data_dir=$4

test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/RT_Erg_NoDyn
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp

######################################
######## RECOMBINATIONS PLOTS ########
### I-FRONT TO THEORETICAL RADIUS  ###
######################################

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R/R_{S}(t)"
set yrange [0.94:1.05]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_radius_errror_nh1_n100.eps"
set label "2D, n_H=10 cm^{-3}" at 0.5,0.995
set label "Cell {/Symbol d}{/Symbol t}=0.5 for 100^2 runs" at 0.5,0.99
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e11
tI=3.861e11
f(x)=1.0
h(x)=1.0+8.0e15/(RI*g(x))
j(x)=1.0-8.0e15/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX05.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX10.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX33.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_radius_errror_nh1_n100.eps ${test_dir}/COMP_radius_errror_nh1_n100.jpeg
#exit

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R/R_{S}(t)"
set yrange [0.93:1.03]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_radius_errror_nh2_n100.eps"
set label "2D, n_H=100 cm^{-3}" at 0.5,0.972
set label "Cell {/Symbol d}{/Symbol t}=5 for 100^2 runs" at 0.5,0.965
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e10
tI=3.861e10
f(x)=1.0
h(x)=1.0+8.0e15/(RI*g(x))
j(x)=1.0-8.0e15/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX05.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX10.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX33.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_radius_errror_nh2_n100.eps ${test_dir}/COMP_radius_errror_nh2_n100.jpeg
#exit

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R/R_{S}(t)"
set yrange [0.85:1.05]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_radius_errror_nh2_n032.eps"
set label "2D, n_H=100 cm^{-3}" at 0.5,0.98
set label "Cell {/Symbol d}{/Symbol t}=15.75 for 32^2 runs" at 0.5,0.97
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e10
tI=3.861e10
f(x)=1.0
h(x)=1.0+2.5e16/(RI*g(x))
j(x)=1.0-2.5e16/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     '${data_dir}/rtt2d_ERG_n032_nh2_dX05.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 2  title '32^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n032_nh2_dX10.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 3  title '32^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n032_nh2_dX33.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 4  title '32^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_radius_errror_nh2_n032.eps ${test_dir}/COMP_radius_errror_nh2_n032.jpeg
#exit

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R/R_{S}(t)"
set yrange [0.92:1.03]
set log x
set xrange [0.01:10]
set title ""
set output "${test_dir}/COMP_radius_errror_nh3_n100.eps"
set label "2D, n_H=10^{3} cm^{-3}" at 0.5,0.97
set label "Cell {/Symbol d}{/Symbol t}=50 for 100^2 runs" at 0.5,0.963
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e09
tI=3.861e09
f(x)=1.0
h(x)=1.0+8.0e15/(RI*g(x))
j(x)=1.0-8.0e15/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX05.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX10.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX33.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_radius_errror_nh3_n100.eps ${test_dir}/COMP_radius_errror_nh3_n100.jpeg
#exit


#####################################
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
set label "2D, n_H=10 cm^{-3}" at 0.02,0.85
set label "Cell {/Symbol d}{/Symbol t}=0.5 for 100^2 runs" at 0.02,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e11
tI=3.861e11
plot g(x) w l lw 2 lt -1 title 'Theoretical Radius', \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX05.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX10.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh1_dX33.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
     
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh1.eps ${test_dir}/COMP_photoncons_rec2d_nh1.jpeg
#exit


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
set label "2D, n_H=100 cm^{-3}" at 0.02,0.85
set label "Cell {/Symbol d}{/Symbol t}=5 for 100^2 runs" at 0.02,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e10
tI=3.861e10
plot g(x) w l lw 2 lt -1 title 'Theoretical Radius', \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX05.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX10.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh2_dX33.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh2.eps ${test_dir}/COMP_photoncons_rec2d_nh2.jpeg

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
set output "${test_dir}/COMP_photoncons_rec2d_nh2_LOWRES.eps"
set label "2D, n_H=100 cm^{-3}" at 0.02,0.85
set label "Cell {/Symbol d}{/Symbol t}=15.75 for 32^2 runs" at 0.02,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e10
tI=3.861e10
plot g(x) w l lw 2 lt -1 title 'Theoretical Radius', \
'${data_dir}/rtt2d_ERG_n032_nh2_dX05.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 2  title '32^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n032_nh2_dX10.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 3  title '32^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n032_nh2_dX33.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 4  title '32^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh2_LOWRES.eps ${test_dir}/COMP_photoncons_rec2d_nh2_LOWRES.jpeg

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
set output "${test_dir}/COMP_photoncons_rec2d_nh3.eps"
set label "2D, n_H=10^{3} cm^{-3}" at 0.02,0.85
set label "Cell {/Symbol d}{/Symbol t}=50 for 100^2 runs" at 0.02,0.8
g(x)=sqrt(1.0-exp(-x))
RH=6.072e17
RI=6.072e17        # lower for new model which includes 10% He
tH=3.861e09
tI=3.861e09
plot g(x) w l lw 2 lt -1 title 'Theoretical Radius', \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX05.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 2  title '100^2, {/Symbol d}t= 0.15/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX10.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 3  title '100^2, {/Symbol d}t= 0.30/xdot', \
     '${data_dir}/rtt2d_ERG_n100_nh3_dX33.txt' u (\$1/tI):(\$2/RI)    w l lw 2 lt 4  title '100^2, {/Symbol d}t= 1.00/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/COMP_photoncons_rec2d_nh3.eps ${test_dir}/COMP_photoncons_rec2d_nh3.jpeg


