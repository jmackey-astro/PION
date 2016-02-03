#!/bin/bash

###############################################################
# 23/11/2009 Plots of Photon Conservation with recombinations #
###############################################################

test_dir=$1
code_dir=$2
data_dir=$3
cd ${test_dir}

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.4:1.1]
set title ""
set output "${test_dir}/photoncons_rec3d_nh1.eps"
set label "3D, n_H=10 cm^{-3}" at 6,0.85
set label "Cell {/Symbol d}{/Symbol t}=1 for 64^3 runs" at 6,0.8
#g(x)=6.515e17*sqrt(1.0-exp(-x))
g(x)=1.0*(1.0-exp(-x))**(1.0/3.0)
plot '${data_dir}/rtt3D_n32_nh1_rec_dt010.txt'  u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 1 pt 1 title '32^3, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt3D_n64_nh1_rec_dt010.txt' u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 2 pt 5 title '64^3, {/Symbol d}t= 1.0 t_{rec}', \
     g(x) w l lt -1 title 'Theoretical Radius'

#     '${data_dir}/rtt3D_n32_nh1_rec_dt100.txt'  u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 1 pt 2 title '32^3, {/Symbol d}t= 0.1 t_{rec}', \
#     '${data_dir}/rtt3D_n64_nh1_rec_dt100.txt' u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 3 pt 4 title '64^3, {/Symbol d}t= 0.1 t_{rec}', \
#     '${data_dir}/rtt3D_n32_nh1_rec_dt500.txt'  u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 1 pt 3 title '32^3, {/Symbol d}t= 0.01 t_{rec}', \
#     '${data_dir}/rtt3D_n64_nh1_rec_dt500.txt' u (\$1/3.861e11):(\$2/6.515e17)    w lp lt 3 pt 6 title '64^3, {/Symbol d}t= 1.0 t_{rec}', \
# pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_rec3d_nh1.eps ${test_dir}/photoncons_rec3d_nh1.jpeg

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Radius (R_{S})"
set yrange [0.4:1.1]
#set xrange [0:7]
set title ""
set output "${test_dir}/photoncons_rec3d_nh2.eps"
set label "3D, n_H=100 cm^{-3}" at 6,0.85
set label "Cell {/Symbol d}{/Symbol t}=10 for 64^3 runs" at 6,0.8
g(x)=1.0*(1.0-exp(-x))**(1.0/3.0)
plot '${data_dir}/rtt3D_n32_nh2_rec_dt010.txt'  u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 1 pt 1 title '32^3, {/Symbol d}t= 1.0 t_{rec}', \
     '${data_dir}/rtt3D_n64_nh2_rec_dt010.txt' u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 2 pt 5 title '64^3, {/Symbol d}t= 1.0 t_{rec}', \
     g(x) w l lt -1 title 'Theoretical Position'

#     '${data_dir}/rtt3D_n32_nh2_rec_dt100.txt'  u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 1 pt 2 title '32^3, {/Symbol d}t= 0.1 t_{rec}', \
#     '${data_dir}/rtt3D_n64_nh2_rec_dt100.txt' u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 3 pt 4 title '64^3, {/Symbol d}t= 0.1 t_{rec}', \
#     '${data_dir}/rtt3D_n32_nh2_rec_dt500.txt'  u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 1 pt 3 title '32^3, 500 steps', \
#     '${data_dir}/rtt3D_n64_nh2_rec_dt500.txt' u (\$1/3.861e10):(\$2/6.515e17)    w lp lt 3 pt 6 title '64^3, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_rec3d_nh2.eps ${test_dir}/photoncons_rec3d_nh2.jpeg

cp ${test_dir}/photoncons_rec3d* ${data_dir}

exit
