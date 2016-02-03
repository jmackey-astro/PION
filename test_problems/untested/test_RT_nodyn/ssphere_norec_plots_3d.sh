#!/bin/bash

###########################################
# 21/11/2009 Plots of Photon Conservation.#
###########################################

test_dir=$1
code_dir=$2
data_dir=$3
cd ${test_dir}

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.004:0.1]
set xrange [0:10]
set title ""
set output "${test_dir}/photoncons_norec3d_nh1.eps"
set label "3D, n_H=10 cm^{-3}" at 2,0.08
plot '${data_dir}/rtt3D_n32_nh1_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^3, 010 steps', \
     '${data_dir}/rtt3D_n32_nh1_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^3, 100 steps', \
     '${data_dir}/rtt3D_n64_nh1_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 4 title '64^3, 010 steps', \
     '${data_dir}/rtt3D_n64_nh1_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 5 title '64^3, 100 steps'
#     '${data_dir}/rtt3D_n32_nh1_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^3, 500 steps', \
#     '${data_dir}/rtt3D_n64_nh1_norec_dt100.txt' u (\$1/3.16e10):(\$2/5.503473e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#     '${data_dir}/rtt3D_n64_nh1_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '64^3, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_norec3d_nh1.eps ${test_dir}/photoncons_norec3d_nh1.jpeg

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.003:0.1]
set xrange [0:10]
set title ""
set label "3D, n_H=100 cm^{-3}" at 2,0.08
set output "${test_dir}/photoncons_norec3d_nh2.eps"
plot '${data_dir}/rtt3D_n32_nh2_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^3, 010 steps', \
     '${data_dir}/rtt3D_n32_nh2_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^3, 100 steps', \
     '${data_dir}/rtt3D_n64_nh2_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 4 title '64^3, 010 steps', \
     '${data_dir}/rtt3D_n64_nh2_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 5 title '64^3, 100 steps'
#     '${data_dir}/rtt3D_n32_nh2_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^3, 500 steps', \
#     '${data_dir}/rtt3D_n64_nh2_norec_dt100.txt' u (\$1/3.16e10):(\$2/5.503417e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#     '${data_dir}/rtt3D_n64_nh2_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '64^3, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_norec3d_nh2.eps ${test_dir}/photoncons_norec3d_nh2.jpeg

cp ${test_dir}/photoncons_norec3d* ${data_dir}
exit


cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "|1-N_{i}/N_{/Symbol g}|"
set yrange [0.001:0.1]
set title ""
set label "3D, n_H=1000 cm^{-3}" at 2,0.08
set output "${test_dir}/photoncons_norec3d_nh3.eps"
plot '${data_dir}/rtt3D_n32_nh3_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^3, 010 steps', \
     '${data_dir}/rtt3D_n32_nh3_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^3, 100 steps', \
     '${data_dir}/rtt3D_n64_nh3_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 4 title '64^3, 010 steps', \
     '${data_dir}/rtt3D_n64_nh3_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 5 title '64^3, 100 steps'
#     '${data_dir}/rtt3D_n32_nh3_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^3, 500 steps', \
#     '${data_dir}/rtt3D_n64_nh3_norec_dt100.txt' u (\$1/3.16e10):(\$2/5.510277e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#     '${data_dir}/rtt3D_n64_nh3_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '64^3, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
exit

