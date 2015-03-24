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
set yrange [0.001:0.1]
set xrange [0:10]
set title ""
set output "${test_dir}/photoncons_norec2d_nh1.eps"
set label "2D, n_H=10 cm^{-3}" at 2,0.08
plot '${data_dir}/rtt2D_n32_nh1_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^2, 010 steps', \
     '${data_dir}/rtt2D_n32_nh1_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh1_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 4 title '100^2, 010 steps', \
     '${data_dir}/rtt2D_n100_nh1_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n256_nh1_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 6 title '256^2, 010 steps', \
     '${data_dir}/rtt2D_n256_nh1_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 7 title '256^2, 100 steps'

#     '${data_dir}/rtt2D_n32_nh1_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh1_norec_dt500.txt' u (\$1/3.16e10):(\$2/5.503473e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#     '${data_dir}/rtt2D_n100_nh1_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '100^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh1_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 9 title '256^2, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_norec2d_nh1.eps ${test_dir}/photoncons_norec2d_nh1.jpeg

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.001:0.1]
set xrange [0:10]
set title ""
set label "2D, n_H=100 cm^{-3}" at 2,0.07
set output "${test_dir}/photoncons_norec2d_nh2.eps"
plot '${data_dir}/rtt2D_n32_nh2_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^2, 010 steps', \
     '${data_dir}/rtt2D_n32_nh2_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh2_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 4 title '100^2, 010 steps', \
     '${data_dir}/rtt2D_n100_nh2_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n256_nh2_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 6 title '256^2, 010 steps', \
     '${data_dir}/rtt2D_n256_nh2_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 7 title '256^2, 100 steps'

#     '${data_dir}/rtt2D_n100_nh2_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '100^2, 500 steps', \
#     '${data_dir}/rtt2D_n32_nh2_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh2_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 9 title '256^2, 500 steps', \
#     '${data_dir}/rtt2D_n256_nh2_norec_dt500.txt' u (\$1/3.16e10):(\$2/5.503417e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_norec2d_nh2.eps ${test_dir}/photoncons_norec2d_nh2.jpeg

cat << EOF  > gnu.plt
set terminal postscript enhanced eps
#set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time"
set key right top
set log y
set ylabel "1 - N_{i} /N_{/Symbol g}"
set yrange [0.001:0.1]
set xrange [0:10]
set title ""
set label "2D n_H=1000 cm^{-3}" at 2,0.08
set output "${test_dir}/photoncons_norec2d_nh3.eps"
plot '${data_dir}/rtt2D_n32_nh3_norec_dt010.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 1 title '32^2, 010 steps', \
     '${data_dir}/rtt2D_n32_nh3_norec_dt100.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 2 title '32^2, 100 steps', \
     '${data_dir}/rtt2D_n100_nh3_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 2 pt 4 title '100^2, 010 steps', \
     '${data_dir}/rtt2D_n100_nh3_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 2 pt 5 title '100^2, 100 steps', \
     '${data_dir}/rtt2D_n256_nh3_norec_dt010.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 6 title '256^2, 010 steps', \
     '${data_dir}/rtt2D_n256_nh3_norec_dt100.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 7 title '256^2, 100 steps'

#     '${data_dir}/rtt2D_n256_nh3_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 4 pt 8 title '256^2, 500 steps'
#     '${data_dir}/rtt2D_n256_nh3_norec_dt500.txt' u (\$1/3.16e10):(\$2/5.510277e+19)  w l lt -1 title 'I-front radius r(t)/[100r(t_f)]'
#     '${data_dir}/rtt2D_n100_nh3_norec_dt500.txt' u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 3 pt 6 title '100^2, 500 steps', \
#     '${data_dir}/rtt2D_n32_nh3_norec_dt500.txt'  u (\$1/3.16e10):(1-\$6/\$7)    w lp lt 1 pt 3 title '32^2, 500 steps', \
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${test_dir}/photoncons_norec2d_nh3.eps ${test_dir}/photoncons_norec2d_nh3.jpeg

cp ${test_dir}/photoncons_norec2d* ${data_dir}
exit


########################################
######## OLD STUFF FROM 2008 ###########
########################################

# plot the angular errors...
set terminal postscript enhanced color eps
set output "azmthangle_n65d3const_nh100_weights.eps"
set title "Azimuthal Angle Dependence of I-Front Radius, n=65^3, constant density, 1/Tau Wt. (min:0.7), t_f=1e12"
set title ""
unset key
unset log y
set xrange [0:360]
dx=1.98462e17
set yrange [2.3e18:4.1e18]
set key right bottom
plot 'd3_const_n65/n65nh100ss3e46_op5_dt1e10_phi.txt' u 1:($3-1.4e18) lt 1     title "tau=125, op5", \
     'd3_const_n65/n65nh100ss3e46_op5_dt1e10_phi.txt' u 1:($7-1.4e18) w l lt 1 title "", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10_phi.txt' u 1:($3-1.1e18) lt 2     title "tau=125, op6", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10_phi.txt' u 1:($7-1.1e18) w l lt 2 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:($3-0.8e18) lt 3     title "tau=125, op7", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:($7-0.8e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op8_dt1e10_phi.txt' u 1:($3-0.5e18) lt 4     title "tau=125, op8", \
     'd3_const_n65/n65nh100ss3e46_op8_dt1e10_phi.txt' u 1:($7-0.5e18) w l lt 4 title "", \
     'd3_const_n65/n65nh100ss3e46_op9_dt1e10_phi.txt' u 1:($3-0.2e18) lt 0     title "tau=125, op9", \
     'd3_const_n65/n65nh100ss3e46_op9_dt1e10_phi.txt' u 1:($7-0.2e18) w l lt 0 title ""
pause -1
set title "Azimuthal Angle Dependence of I-Front Radius, n=65^3, constant density, 1/Tau Wt. (min:0.7), t_f=1e12"
set output "polarangle_n65d3const_nh100_weights.eps"
set title ""
unset key
unset log y
set xrange [0:180]
dx=1.98462e17
set yrange [2.3e18:4.1e18]
set key right bottom
plot 'd3_const_n65/n65nh100ss3e46_op5_dt1e10_theta.txt' u 1:($3-1.4e18) lt 1     title "tau=125, op5", \
     'd3_const_n65/n65nh100ss3e46_op5_dt1e10_theta.txt' u 1:($7-1.4e18) w l lt 1 title "", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10_theta.txt' u 1:($3-1.1e18) lt 2     title "tau=125, op6", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10_theta.txt' u 1:($7-1.1e18) w l lt 2 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:($3-0.8e18) lt 3     title "tau=125, op7", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:($7-0.8e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op8_dt1e10_theta.txt' u 1:($3-0.5e18) lt 4     title "tau=125, op8", \
     'd3_const_n65/n65nh100ss3e46_op8_dt1e10_theta.txt' u 1:($7-0.5e18) w l lt 4 title "", \
     'd3_const_n65/n65nh100ss3e46_op9_dt1e10_theta.txt' u 1:($3-0.2e18) lt 0     title "tau=125, op9", \
     'd3_const_n65/n65nh100ss3e46_op9_dt1e10_theta.txt' u 1:($7-0.2e18) w l lt 0 title ""
pause -1
quit

#*************** POLAR ANGLE PLOTS *************************
set term postscript enhanced color

set output "d3_const_n65/polar_meanrad_geowt.eps"
#set title "Polar Angle Dependence of I-Front Radius, n=65^3, constant density, Geometric Wt., t_f=1e12"
set title ""
set xlabel "Polar Angle (degrees)"
set ylabel "Mean I-Front Radius (points offset for clarity)"
set ylabel "Mean IF Radius (cm)"
set key right bottom
unset key
unset log y
dx=1.98462e17
set yrange [2.0e18:4.5e18]
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_theta.txt'   u 1:($3-0.0e18) lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_theta.txt'   u 1:($7-0.0e18) w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_theta.txt'   u 1:($3-0.6e18) lt 2     title "tau=4", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_theta.txt'   u 1:($7-0.6e18) w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt'  u 1:($3-1.0e18) lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt'  u 1:($7-1.0e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_theta.txt' u 1:($3-1.4e18) lt 4     title "tau=125", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_theta.txt' u 1:($7-1.4e18) w l lt 4 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt'  u 1:($7-1.0e18+dx) w l lt 0 title "+/- dx", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt'  u 1:($7-1.0e18-dx) w l lt 0 title ""
#pause -1
set output "d3_const_n65/polar_meanrad_invtauwt.eps"
#set title "Polar Angle Dependence of I-Front Radius, n=65^3, constant density, 1/Tau Wt. (min:0.7), t_f=1e12"
set title ""
set key right bottom
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_theta.txt'   u 1:($3-0.0e18) lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh1ss3e44_op7_dt1e10_theta.txt'   u 1:($7-0.0e18) w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_theta.txt'   u 1:($3-0.6e18) lt 2     title "tau=4", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_theta.txt'   u 1:($7-0.6e18) w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt'  u 1:($3-1.0e18) lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt'  u 1:($7-1.0e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:($3-1.4e18) lt 4     title "tau=125", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:($7-1.4e18) w l lt 4 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt'  u 1:($7-1.0e18+dx) w l lt 0 title "+/- dx", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt'  u 1:($7-1.0e18-dx) w l lt 0 title ""
#pause -1

set output "d3_const_n65/polar_numvarrad_geowt.eps"
set title "Numerical Variance of I-Front Radius: Polar Angle Dependence, n=65^3, rho=const., Geo.Wt."
set title ""
set ylabel "Max.-Min. Radius (units of cell size)"
#set label 1 "Note that in this case the I-Front is resolved, so the variance is real!"  at first 5, first 1.75e18
#set label 2 "Zero values are for angle bins which had 0 or 1 grid point(s), or all the same radius"  at first 5, first 1.6e18
set xrange [0:230]
set yrange [-1:13]
set key top right
dx=1.98462e+17
g=0.0
h=3.0
i=6.0
j=9.0
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_theta.txt' u 1:(($5-$6)/dx+j) lt 1     title "tau=1.25", \
      j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_theta.txt' u 1:(($5-$6)/dx+i) lt 2     title "tau=4.00", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt' u 1:(($5-$6)/dx+h) lt 3     title "tau=12.5", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_theta.txt' u 1:(($5-$6)/dx+g) lt 4     title "tau=125.", \
     g w l lt 4 title ""
#pause -1
set output "d3_const_n65/polar_numvarrad_invtauwt.eps"
set title "Numerical Variance of I-Front Radius: Polar Angle Dependence, n=65^3, rho=const., 1/Tau Wt."
set title ""
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_theta.txt' u 1:(($5-$6)/dx+j) lt 1     title "tau=1.25", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_theta.txt' u 1:(($5-$6)/dx+i) lt 2     title "tau=4.00", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt' u 1:(($5-$6)/dx+h) lt 3     title "tau=12.5", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:(($5-$6)/dx+g) lt 4     title "tau=125.", \
     g w l lt 4 title ""
#pause -1

set output "d3_const_n65/polar_MaxErrRad_geowt.eps"
set title "Max. Error in I-Front Radius: Polar Angle Dependence, n=65^3, rho=const., Geo.Wt."
set title ""
set ylabel "Analytic-Max/Min Radius (cell size)"
#set label 1 "In this case the I-Front is resolved, so variance is real."  at 10,1.95e18
#unset label 2
set xrange [0:240]
set yrange [-1:14]
g=0.0
h=4.0
i=7.0
j=11.0
set key top right
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_theta.txt' u 1:(($7-$5)/dx+j) lt 1 pt 1 title "tau=1.25, max", \
     'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_theta.txt' u 1:(($7-$6)/dx+j) lt 1 pt 2 title "tau=1.25, min", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_theta.txt' u 1:(($7-$5)/dx+i) lt 2 pt 1 title "tau=4.00, max", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_theta.txt' u 1:(($7-$6)/dx+i) lt 2 pt 2 title "tau=4.00, min", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt' u 1:(($7-$5)/dx+h) lt 3 pt 1 title "tau=12.5, max", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_theta.txt' u 1:(($7-$6)/dx+h) lt 3 pt 2 title "tau=12.5, min", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_theta.txt' u 1:(($7-$5)/dx+g) lt 4 pt 1 title "tau=125,  max", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_theta.txt' u 1:(($7-$6)/dx+g) lt 4 pt 2 title "tau=125,  min", \
     g w l lt 4 title ""
#pause -1

set output "d3_const_n65/polar_MaxErrRad_invtauwt.eps"
#set title "Max. Error in I-Front Radius: Polar Angle Dependence, n=65^3, rho=const., 1/Tau Wt.(max=0.7)"
set title ""
set ylabel "Analytic-Max/Min Radius (cell size)"
#set label 1 "In this case the I-Front is resolved, so variance is real."  at 10,1.95e18
#set label 2 "Note that Max. Error is always less than one cell for 1/Tau weighting" at 10,0.3e18
set xrange [0:240]
set yrange [-1:14]
g=0.0
h=4.0
i=7.0
j=11.0
set key top right
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_theta.txt' u 1:(($7-$5)/dx+j) lt 1 pt 1 title "tau=1.25, max", \
     'd3_const_n65/n65nh1ss3e44_op7_dt1e10_theta.txt' u 1:(($7-$6)/dx+j) lt 1 pt 2 title "tau=1.25, min", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_theta.txt' u 1:(($7-$5)/dx+i) lt 2 pt 1 title "tau=4.00, max", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_theta.txt' u 1:(($7-$6)/dx+i) lt 2 pt 2 title "tau=4.00, min", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt' u 1:(($7-$5)/dx+h) lt 3 pt 1 title "tau=12.5, max", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_theta.txt' u 1:(($7-$6)/dx+h) lt 3 pt 2 title "tau=12.5, min", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:(($7-$5)/dx+g) lt 4 pt 1 title "tau=125,  max", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_theta.txt' u 1:(($7-$6)/dx+g) lt 4 pt 2 title "tau=125,  min", \
     g w l lt 4 title ""
#pause -1
#unset label 1
#unset label 2
#*************** AZIMUTHAL ANGLE PLOTS *************************

set output "d3_const_n65/azimuth_meanrad_geowt.eps"
set title "Azimuthal Angle Dependence of I-Front Radius, n=65^3, constant density, Geometric Wt., t_f=1e12"
set title ""
set xlabel "Azimuthal Angle (degrees)"
set ylabel "Mean I-Front Radius (points offset for clarity)"
set ylabel "Mean IF Radius (cm)"
set key right bottom
unset key
unset log y
set xrange [0:360]
dx=1.98462e17
set yrange [2.0e18:4.5e18]
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_phi.txt'   u 1:($3-0.0e18) lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_phi.txt'   u 1:($7-0.0e18) w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_phi.txt'   u 1:($3-0.6e18) lt 2     title "tau=4", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_phi.txt'   u 1:($7-0.6e18) w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt'  u 1:($3-1.0e18) lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt'  u 1:($7-1.0e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_phi.txt' u 1:($3-1.4e18) lt 4     title "tau=125", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_phi.txt' u 1:($7-1.4e18) w l lt 4 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt'  u 1:($7-1.0e18+dx) w l lt 0 title "+/- dx", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt'  u 1:($7-1.0e18-dx) w l lt 0 title ""
#pause -1
set output "d3_const_n65/azimuth_meanrad_invtauwt.eps"
set title "Azimuthal Angle Dependence of I-Front Radius, n=65^3, constant density, 1/Tau Wt. (min:0.7), t_f=1e12"
set title ""
set key right bottom
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_phi.txt'   u 1:($3-0.0e18) lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh1ss3e44_op7_dt1e10_phi.txt'   u 1:($7-0.0e18) w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_phi.txt'   u 1:($3-0.6e18) lt 2     title "tau=4", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_phi.txt'   u 1:($7-0.6e18) w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt'  u 1:($3-1.0e18) lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt'  u 1:($7-1.0e18) w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:($3-1.4e18) lt 4     title "tau=125", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:($7-1.4e18) w l lt 4 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt'  u 1:($7-1.0e18+dx) w l lt 0 title "+/- dx", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt'  u 1:($7-1.0e18-dx) w l lt 0 title ""
#pause -1

set output "d3_const_n65/azimuth_numvarrad_geowt.eps"
set title "Numerical Variance of I-Front Radius: Azimuthal Angle Dependence, n=65^3, rho=const., Geo.Wt."
set title ""
set ylabel "Numerical Max.-Min. Radius (cell size)"
#set label 1 "Note that in this case the I-Front is resolved, so the variance is real!"  at first 10, first 1.95e18
#set label 2 "Zero values are for angle bins which had 0 or 1 grid point(s), or all the same radius"  at first 10, first 1.6e18
set xrange [0:450]
set yrange [-1:13]
set key top right
dx=1.98462e+17
g=0.0
h=3.0
i=6.0
j=9.0
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_phi.txt' u 1:(($5-$6)/dx+j) lt 1     title "tau=1.25", \
      j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_phi.txt' u 1:(($5-$6)/dx+i) lt 2     title "tau=4.00", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt' u 1:(($5-$6)/dx+h) lt 3     title "tau=12.5", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_phi.txt' u 1:(($5-$6)/dx+g) lt 4     title "tau=125.", \
     g w l lt 4 title ""
#pause -1
set output "d3_const_n65/azimuth_numvarrad_invtauwt.eps"
set title "Numerical Variance of I-Front Radius: Azimuthal Angle Dependence, n=65^3, rho=const., 1/Tau Wt."
set title ""
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_phi.txt' u 1:(($5-$6)/dx+j) lt 1     title "tau=1.25", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_phi.txt' u 1:(($5-$6)/dx+i) lt 2     title "tau=4.00", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt' u 1:(($5-$6)/dx+h) lt 3     title "tau=12.5", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:(($5-$6)/dx+g) lt 4     title "tau=125.", \
     g w l lt 4 title ""
#pause -1

set output "d3_const_n65/azimuth_MaxErrRad_geowt.eps"
set title "Max. Error in I-Front Radius: Azimuthal Angle Dependence, n=65^3, rho=const., Geo.Wt."
set title ""
set ylabel "Analytic-Max/Min Radius (cell size)"
#set label 1 "In this case the I-Front is resolved, so variance is real."  at 10,1.95e18
#unset label 2
set xrange [0:500]
set yrange [-2:14]
g=0.0
h=4.0
i=7.0
j=11.0
set key top right
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_phi.txt' u 1:(($7-$5)/dx+j) lt 1 pt 1 title "tau=1.25, max", \
     'd3_const_n65/n65nh1ss3e44_geowt_dt1e10_phi.txt' u 1:(($7-$6)/dx+j) lt 1 pt 2 title "tau=1.25, min", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_phi.txt' u 1:(($7-$5)/dx+i) lt 2 pt 1 title "tau=4.00, max", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10_phi.txt' u 1:(($7-$6)/dx+i) lt 2 pt 2 title "tau=4.00, min", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt' u 1:(($7-$5)/dx+h) lt 3 pt 1 title "tau=12.5, max", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10_phi.txt' u 1:(($7-$6)/dx+h) lt 3 pt 2 title "tau=12.5, min", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_phi.txt' u 1:(($7-$5)/dx+g) lt 4 pt 1 title "tau=125,  max", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10_phi.txt' u 1:(($7-$6)/dx+g) lt 4 pt 2 title "tau=125,  min", \
     g w l lt 4 title ""
#pause -1

set output "d3_const_n65/azimuth_MaxErrRad_invtauwt.eps"
set title "Max. Error in I-Front Radius: Azimuthal Angle Dependence, n=65^3, rho=const., 1/Tau Wt.(max=0.7)"
set title ""
set ylabel "Analytic-Max/Min Radius (cell size)"
#set label 1 "In this case the I-Front is resolved, so variance is real."  at 10,1.95e18
#set label 2 "Note that Max. Error is always less than one cell for 1/Tau weighting" at 10,0.3e18
set xrange [0:500]
set yrange [-2:14]
g=0.0
h=4.0
i=7.0
j=11.0
set key top right
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10_phi.txt' u 1:(($7-$5)/dx+j) lt 1 pt 1 title "tau=1.25, max", \
     'd3_const_n65/n65nh1ss3e44_op7_dt1e10_phi.txt' u 1:(($7-$6)/dx+j) lt 1 pt 2 title "tau=1.25, min", \
     j w l lt 1 title "", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_phi.txt' u 1:(($7-$5)/dx+i) lt 2 pt 1 title "tau=4.00, max", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10_phi.txt' u 1:(($7-$6)/dx+i) lt 2 pt 2 title "tau=4.00, min", \
     i w l lt 2 title "", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt' u 1:(($7-$5)/dx+h) lt 3 pt 1 title "tau=12.5, max", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10_phi.txt' u 1:(($7-$6)/dx+h) lt 3 pt 2 title "tau=12.5, min", \
     h w l lt 3 title "", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:(($7-$5)/dx+g) lt 4 pt 1 title "tau=125,  max", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10_phi.txt' u 1:(($7-$6)/dx+g) lt 4 pt 2 title "tau=125,  min", \
     g w l lt 4 title ""
#pause -1
quit


##############################################################################
# 3D 3D 3D!!! CONSTANT DENSITY PROFILE, N=129, RHO=1.67e-22, p=2.76e-12, ss=3e46 phot/s
# These lines make plots of photon conservation and asphericity for 3D pt.src
# raytracing.
##############################################################################
#set terminal postscript enhanced color
#set output "n129d3_const_photon_cons.eps"
set term x11

set title "Azimuthal Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.7, t_f=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op7_dt5e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op7_dt2e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e09_phi.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Azimuthal Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=1/sqrt(3), t_f=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_rt3_dt5e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt2e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e09_phi.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Azimuthal Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.5, t_f=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op5_dt5e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op5_dt2e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e09_phi.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Azimuthal Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.6, t_f=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op6_dt5e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op6_dt2e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e09_phi.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e09_phi.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

#
# now for polar angle dependence
#
set title "Polar Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.7, t_f=1e12"
set xlabel "Polar Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op7_dt5e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op7_dt2e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op7_dt1e09_theta.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Polar Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=1/sqrt(3), t_f=1e12"
set xlabel "Polar Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_rt3_dt5e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt2e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_rt3_dt1e09_theta.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Polar Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.5, t_f=1e12"
set xlabel "Polar Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op5_dt5e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op5_dt2e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op5_dt1e09_theta.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Polar Angle Dependence of I-Front Radius, n=129^3, constant density, Min.Tau=0.6, t_f=1e12"
set xlabel "Polar Angle (degrees)"
set ylabel "I-Front Radius"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_op6_dt5e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=5e10", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e10", \
     'd3_const_n129/n129nh100ss3e46_op6_dt2e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=2e09", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e09_theta.txt' u 1:3:5:6 w yerrorbars title "dt=1e09", \
     'd3_const_n129/n129nh100ss3e46_op6_dt1e09_theta.txt' u 1:7 w l title "Analytic Radius", \
     fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

# end of polar angle dependence





set title "Azimuthal Angle Dependence of I-Front Radius, n=129^3, constant density, dt=1e10, t=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "Radius(at Angle)"
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [4.0e18:4.3e18]
plot 'd3_const_n129/n129nh100ss3e46_cart3d_op7_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars pt 2 lt 2 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Polar Angle Dependence of I-Front Radius, n=129^3, constant density, dt=1e10, t=1e12" 
set xlabel "Polar Angle (degrees)"
set ylabel "Radius of I-Front"
set key right bottom
unset log y
f(x)=4.152831e18
fp(x)=4.252831e18
fm(x)=4.052831e18
set yrange [3.8e18:4.5e18]
plot 'd3_const_n129/n129nh100ss3e46_cart3d_op7_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars pt 2 lt 2 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set title "Azimuthal Angle Dependence of I-Front Radius, n=32^3, constant density, dt=1e10, t=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "Radius(at Angle)"
fp(x)=4.442831e18
fm(x)=3.862831e18
set yrange [3.8e18:4.5e18]
plot 'd3_const_n32/n32nh100ss3e46_cart3d_root3_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars pt 3 lt 3 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1


set title "Polar Angle Dependence of I-Front Radius, n=32^3, constant density, dt=1e10, t=1e12" 
set xlabel "Polar Angle (degrees)"
set ylabel "Radius of I-Front"
fp(x)=4.442831e18
fm(x)=3.862831e18
set grid
set yrange [3.8e18:4.5e18]
plot  'd3_const_n32/n32nh100ss3e46_cart3d_root3_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars pt 3 lt 3 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set yrange [*:*]
quit

set output "n129d2_rm2_dRoverR.eps"
set xlabel "Radius (cm)"
set key right bottom
set log y
f(x)=1.0e17/x
set ylabel "I-Front Width (dR/R)"
set title "Photon Conservation n=129^2, 2d Cartesian, Point Source, 1/r^2 profile"
plot 'RRadial2_dt3e10.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=3e10' lt 0 pt 1, \
     'RRadial2_dt1e10.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'RRadial2_dt3e9.txt'  u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=3e09' lt 2 pt 2, \
     'RRadial2_dt1e9.txt'  u 2:(($3-$4)/$2) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     'RRadial2_dt3e8.txt'  u 2:(($3-$4)/$2) w l title 'Tau-Weighted, dt=3e08' lt 4, \
     'RRadial2_dt1e8.txt'  u 2:(($3-$4)/$2) w l title 'Tau-Weighted, dt=1e08' lt 5, \
     f(x) title "Cell Size"
pause -1
set output "n129d2_rm2_dR.eps"
set xlabel "Radius (cm)"
set key right bottom
set log y
f(x)=1.0e17
set ylabel "I-Front Width (dR)"
set title "Photon Conservation n=129^2, 2d Cartesian, Point Source, 1/r^2 profile"
plot 'RRadial2_dt1e10.txt' u 1:(($3-$4)) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'RRadial2_dt3e9.txt'  u 1:(($3-$4)) w lp title 'Tau-Weighted, dt=3e09' lt 2 pt 2, \
     'RRadial2_dt1e9.txt'  u 1:(($3-$4)) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     'RRadial2_dt3e8.txt'  u 1:(($3-$4)) w l title 'Tau-Weighted, dt=3e08' lt 4, \
     'RRadial2_dt1e8.txt'  u 1:(($3-$4)) w l title 'Tau-Weighted, dt=1e08' lt 5, \
     f(x) title "Cell Size"
pause -1
quit

##############################################################################
# These file are moved to d2_const_n256/*.txt now
# These lines make plots of photon conservation and asphericity for 2D pt.src
# raytracing.
##############################################################################
set terminal postscript enhanced color
set output "n256d2_photon_cons.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=256^2, 2d Cartesian, Point Source"
plot 'test/n256dt5e10.txt' u 1:(abs(1-$6/$7)) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n256dt1e10.txt' u 1:(abs(1-$6/$7)) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n256dt2e9.txt'  u 1:(abs(1-$6/$7)) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n256dt5e10_new.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n256dt1e10_new.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n256dt2e9_new.txt'  u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=2e09' lt 4
pause -1

set output "n256d2_const_photon_cons2.eps"
set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation n=256^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n256/n256dt5e10_new.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd2_const_n256/n256dt1e10_new.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'd2_const_n256/n256dt2e9_new.txt' u 1:($6/$7) w l title 'Tau-Weighted, dt=2e09' lt 2, \
     'd2_const_n256/n256dt5e8_new.txt' u 1:($6/$7) w l title 'Tau-Weighted, dt=5e08' lt 3, \
     f(x) title "x=1"
pause -1

set terminal postscript enhanced color
set output "n256d2_fittoradius.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set yrange [1.e-4:0.1]
set ylabel "abs(1-R(fit)/R(Analytic))"
set title "n=256^2, 2d Cartesian; Fit Radius (mean radius of cells with 0.1<x(HII)<0.9) vs. Analytic Radius"
plot 'test/n256dt5e10.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n256dt1e10.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n256dt2e9.txt'  u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n256dt5e10_new.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n256dt1e10_new.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n256dt2e9_new.txt'  u 1:(abs(1-$2/sqrt($7/100/3.14159))) w l title 'Tau-Weighted, dt=2e09' lt 4
pause -1

set terminal postscript enhanced color
set output "n256d2_asphericity.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set yrange [1.e16:1.e18]
f(x)=5.02e16
set ylabel "I-Front Width"
set title "n=256^2, 2d Cartesian; Asphericity: Max. Radius Minus Min. Radius"
plot 'test/n256dt5e10.txt' u 1:($3-$4) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n256dt1e10.txt' u 1:($3-$4) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n256dt2e9.txt'  u 1:($3-$4) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n256dt5e10_new.txt' u 1:($3-$4) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n256dt1e10_new.txt' u 1:($3-$4) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n256dt2e9_new.txt'  u 1:($3-$4) w l title 'Tau-Weighted, dt=2e09' lt 4, \
     f(x) title "Cell Size" lt 0 lw 4
pause -1
quit

###############################################################################
# The following plots are for showing various conservation properties of code.#
###############################################################################
set terminal postscript enhanced color
set output "radius_error.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, n=129^2, No Recombinations, Photon Conservation as Function of Timestep"
set ylabel "sqrt(num.ions/num.photons)"
unset log y
set yrange [*:*]
f(x)=1
plot 'profile_dt1.txt'  u 1:(sqrt($5/$6)) w lp title "dt=1e11", \
     'profile_dt03.txt' u 1:(sqrt($5/$6)) w lp title "dt=3e10", \
     'profile_dt01.txt'  u 1:(sqrt($5/$6)) w l title "dt=1e10", \
     'profile_dt003.txt' u 1:(sqrt($5/$6)) w l title "dt=3e9", \
     'profile_dt001.txt'  u 1:(sqrt($5/$6)) w l title "dt=1e9", \
     'profile_dt0001.txt' u 1:(sqrt($5/$6)) w l title "dt=1e8", f(x) title ""
pause -1

set output "photons_dt.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, n=129^2, No Recombinations, Photon Conservation as Function of Timestep"
set ylabel "ion number/photon number"
f(x)=1
plot 'profile_dt1.txt'  u 1:(abs($5/$6)) w lp title "dt=1e11", \
     'profile_dt03.txt' u 1:(abs($5/$6)) w lp title "dt=3e10", \
     'profile_dt01.txt'  u 1:(abs($5/$6)) w l title "dt=1e10", \
     'profile_dt003.txt' u 1:(abs($5/$6)) w l title "dt=3e9", \
     'profile_dt001.txt'  u 1:(abs($5/$6)) w l title "dt=1e9", \
     'profile_dt0001.txt' u 1:(abs($5/$6)) w l title "dt=1e8", f(x) title ""
pause -1

set output "photons_dx.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, dt=1e10, No Recombinations, Photon Conservation as Function of Resolution"
set ylabel "ion number/photon number"
f(x)=1
plot 'profile_dt01_n32.txt'  u 1:(abs($5/$6)) w lp title "n=32^2", \
     'profile_dt01_n65.txt'  u 1:(abs($5/$6)) w lp title "n=65^2", \
     'profile_dt01_n129.txt' u 1:(abs($5/$6)) w lp title "n=129^2", \
     'profile_dt01_n256.txt' u 1:(abs($5/$6)) w lp title "n=256^2", f(x) title ""
pause -1

set output "photons_dxdt.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, No Recombinations, Photon Conservation as Function of Resolution and Timestep"
set ylabel "ion number/photon number"
f(x)=1
set yrange [0.95:1.02]
plot 'profile_dt02_n32.txt'  u 1:(abs($5/$6)) w lp title "n=32^2, dt=2e10", \
     'profile_dt01_n65.txt'  u 1:(abs($5/$6)) w lp title "n=65^2, dt=1e10", \
     'profile_dt005_n129.txt' u 1:(abs($5/$6)) w lp title "n=129^2, dt=5e9", \
     'profile_dt0025_n256.txt' u 1:(abs($5/$6)) w lp title "n=256^2, dt=2.5e9", f(x) title ""
pause -1

set output "ifront_posn.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, No Recombinations, I-Front Position"
set ylabel "Radius (cm)"
set grid
set key left
f(x)=4.e12*sqrt(x)
set yrange [4e17:6e18]
unset log y
unset log x
plot 'profile_dt02_n32.txt'  u 1:(1.0*$2) w lp title "n=32^2, dt=2e10, fit", \
     'profile_dt02_n32.txt'  u 1:4 w l title "n=32^2, dt=2e10, nions", \
     'profile_dt0025_n256.txt'  u 1:(sqrt($6/100/3.14159)) w l title "Analytic", \
     f(x) title "K*sqrt(x)", \
     'profile_dt01_n65.txt'  u 1:2 w lp title "n=65^2, dt=1e10, fit", \
     'profile_dt005_n129.txt' u 1:2 w lp title "n=129^2, dt=5e9, fit", \
     'profile_dt0025_n256.txt' u 1:2 w lp title "n=256^2, dt=2.5e9, fit"
pause -1



set term x11
#n=100.0
#s=1.e49
#alpha=2.59e-13
#trec=1/n/alpha
#rs=(3.*s/4/3.14159/alpha/n/n)**(1./3.)
#ri=(3.*s/4/3.14159/n)**(1./3.)
#f(x)=ri*x**(1./3.)
#g(x)=rs*(1.-exp(-x/trec))**(1./3.)
#set xrange [1.e7:1.e11]
#set xlabel "Time (seconds)"
#set ylabel "I-Front Radius (cm)"
#plot 'ss_axi_128_recomb3.txt' u 1:($2*2.5/3.086e18) title "Fit to Simulation", \
#     'ss_axi_128_recomb3.txt' u 1:($4*2.5/3.086e18) title "Alt. Fit", \
#     f(x)/3.086e18 title "No recombs", \
#     g(x)/3.086e18 title "With Recombs"
#pause -1

