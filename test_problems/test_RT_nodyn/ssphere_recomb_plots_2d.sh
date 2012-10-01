#!/bin/bash

###############################################################
# 23/11/2009 Plots of Photon Conservation with recombinations #
###############################################################

test_dir=$1
code_dir=$2
data_dir=$3
cd ${test_dir}

cat << EOF  > gnu.plt
#set terminal postscript enhanced color eps
set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.4:1.1]
set title ""
set output "${test_dir}/photoncons_rec2d_nh1.eps"
set label "2D, n_H=10 cm^{-3}" at 6,0.85
set label "Cell {/Symbol d}{/Symbol t}=1 for 100^2 runs" at 6,0.8
g(x)=sqrt(1.0-exp(-x))
plot '${data_dir}/rtt2D_n32_nh1_rec_dt010.txt'  u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 1 pt 1 title '32^2, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt2D_n32_nh1_rec_dt100.txt'  u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 1 pt 2 title '32^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh1_rec_dt010.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 2 pt 5 title '100^2, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh1_rec_dt100.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 3 pt 4 title '100^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n256_nh1_rec_dt010.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 4 pt 5 title '256^2, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt2D_n256_nh1_rec_dt100.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 4 pt 6 title '256^2, {/Symbol d}t= 0.1 t_{rec}', \
    g(x) w l lt -1 title 'Theoretical Radius'
#     '${data_dir}/rtt2D_n32_nh1_rec_dt500.txt'  u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 1 pt 3 title '32^2, 500 steps', \
#     '${data_dir}/rtt2D_n100_nh1_rec_dt500.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 3 pt 6 title '100^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh1_rec_dt500.txt' u (\$1/3.861e11):(\$2/6.072e17)    w lp lt 4 pt 6 title '256^2, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_rec2d_nh1.eps ${test_dir}/photoncons_rec2d_nh1.jpeg

cat << EOF  > gnu.plt
#set terminal postscript enhanced color eps
set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.4:1.1]
set title ""
set output "${test_dir}/photoncons_rec2d_nh2.eps"
set label "2D, n_H=100 cm^{-3}" at 6,0.85
set label "Cell {/Symbol d}{/Symbol t}=10 for 100^2 runs" at 6,0.8
g(x)=sqrt(1.0-exp(-x))
plot '${data_dir}/rtt2D_n32_nh2_rec_dt010.txt'  u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 1 pt 1 title '32^2, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt2D_n32_nh2_rec_dt100.txt'  u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 1 pt 2 title '32^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh2_rec_dt010.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 2 pt 5 title '100^2, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt2D_n100_nh2_rec_dt100.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 3 pt 4 title '100^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n256_nh2_rec_dt100.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 3 pt 5 title '256^2, {/Symbol d}t= 0.1 t_{rec}', \
     '${data_dir}/rtt2D_n256_nh2_rec_dt010.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 4 pt 5 title '256^2, {/Symbol d}t= 1.0 t_{rec}', \
     g(x) w l lt -1 title 'Theoretical Radius'
#     '${data_dir}/rtt2D_n32_nh2_rec_dt500.txt'  u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 1 pt 3 title '32^2, 500 steps', \
#     '${data_dir}/rtt2D_n100_nh2_rec_dt500.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 3 pt 6 title '100^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh2_rec_dt500.txt' u (\$1/3.861e10):(\$2/6.072e17)    w lp lt 4 pt 6 title '256^2, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_rec2d_nh2.eps ${test_dir}/photoncons_rec2d_nh2.jpeg

cp ${test_dir}/photoncons_rec2d* ${data_dir}
exit

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right top
set log y
set ylabel "Radius (R_{S})"
set yrange [0.01:1.1]
set title ""
set output "${test_dir}/photoncons_rec2d_nh3.eps"
set label "2D, n_H=1000 cm^{-3}" at 2,0.08
g(x)=sqrt(1.0-exp(-x))
plot '${data_dir}/rtt2D_n32_nh3_rec_dt010.txt'  u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 1 pt 1 title '32^2, 010 steps', \
     '${data_dir}/rtt2D_n32_nh3_rec_dt100.txt'  u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 1 pt 2 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n32_nh3_rec_dt500.txt'  u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 1 pt 3 title '32^2, 500 steps', \
     '${data_dir}/rtt2D_n100_nh3_rec_dt010.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 4 title '100^2, 010 steps', \
     '${data_dir}/rtt2D_n100_nh3_rec_dt100.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh3_rec_dt500.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 6 title '100^2, 500 steps', \
     '${data_dir}/rtt2D_n256_nh3_rec_dt010.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 4 title '256^2, 010 steps', \
     '${data_dir}/rtt2D_n256_nh3_rec_dt100.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 5 title '256^2, 100 steps', \
     '${data_dir}/rtt2D_n256_nh3_rec_dt500.txt' u (\$1/3.861e9):(\$2/6.072e17)    w lp lt 3 pt 6 title '256^2, 500 steps', \
     g(x) w l lt -1 title 'Theoretical Radius'
#pause -1
EOF
gnuplot gnu.plt

