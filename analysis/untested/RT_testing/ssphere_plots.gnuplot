##############################################################################
# ********** Monday 19/5/08 ************* photon conservation plots
# 2D with recombinations, plotting actual vs. theoretical radius.
##############################################################################
set terminal postscript enhanced color eps
set xlabel "Time (seconds)"
set key right bottom
unset log y
set title "Photon Conservation n=33^2, 2d Cartesian, Point Source, with Recombinations"
set title ""
set grid
#############################
#### first nh=10^3 plots ####
#############################
set title "dt=0.05 nh=100"
set title ""
set yrange [0.9:1.01]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh3ss5e31_dt020.eps"
plot 'd2_const_rr_n257/nh3ss5e31_dt20_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh3ss5e31_dt20_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh3ss5e31_dt20_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh3ss5e31_dt20_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh3ss5e31_dt20_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh3ss5e31_dt20_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.01 nh=100"
set title ""
set yrange [0.99:1.005]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh3ss5e31_dt100.eps"
plot 'd2_const_rr_n257/nh3ss5e31_dt100_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh3ss5e31_dt100_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh3ss5e31_dt100_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh3ss5e31_dt100_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh3ss5e31_dt100_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh3ss5e31_dt100_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.002 nh=100"
set title ""
set yrange [0.99:1.005]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh3ss5e31_dt500.eps"
plot 'd2_const_rr_n257/nh3ss5e31_dt500_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh3ss5e31_dt500_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh3ss5e31_dt500_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh3ss5e31_dt500_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh3ss5e31_dt500_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh3ss5e31_dt500_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1

#############################
#### first nh=10^2 plots ####
#############################
set title "dt=0.05 nh=100"
set title ""
set yrange [0.9:1.01]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh2ss5e29_dt020.eps"
plot 'd2_const_rr_n257/nh2ss5e29_dt20_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh2ss5e29_dt20_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh2ss5e29_dt20_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh2ss5e29_dt20_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh2ss5e29_dt20_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh2ss5e29_dt20_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.01 nh=100"
set title ""
set yrange [0.98:1.01]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh2ss5e29_dt100.eps"
plot 'd2_const_rr_n257/nh2ss5e29_dt100_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh2ss5e29_dt100_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh2ss5e29_dt100_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh2ss5e29_dt100_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh2ss5e29_dt100_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh2ss5e29_dt100_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.002 nh=100"
set title ""
set yrange [0.98:1.01]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh2ss5e29_dt500.eps"
plot 'd2_const_rr_n257/nh2ss5e29_dt500_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh2ss5e29_dt500_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh2ss5e29_dt500_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh2ss5e29_dt500_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh2ss5e29_dt500_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh2ss5e29_dt500_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1

################################
#### Now nh=10^1 plots  ########
################################
set title "dt=0.05 nh=10"
set title ""
set yrange [0.95:1.02]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh1ss5e27_dt020.eps"
plot 'd2_const_rr_n257/nh1ss5e27_dt20_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh1ss5e27_dt20_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh1ss5e27_dt20_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh1ss5e27_dt20_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh1ss5e27_dt20_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh1ss5e27_dt20_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.01 nh=10"
set title ""
set yrange [0.99:1.02]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh1ss5e27_dt100.eps"
plot 'd2_const_rr_n257/nh1ss5e27_dt100_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh1ss5e27_dt100_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh1ss5e27_dt100_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh1ss5e27_dt100_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh1ss5e27_dt100_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh1ss5e27_dt100_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.002 nh=10"
set title ""
set yrange [0.99:1.02]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh1ss5e27_dt500.eps"
plot 'd2_const_rr_n257/nh1ss5e27_dt500_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh1ss5e27_dt500_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh1ss5e27_dt500_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh1ss5e27_dt500_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh1ss5e27_dt500_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh1ss5e27_dt500_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
################################
#### Now nh=10^0 plots  ########
################################
set title "dt=0.05 nh=1"
set title ""
set yrange [0.98:1.06]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh0ss25e24_dt020.eps"
plot 'd2_const_rr_n257/nh0ss5e25_dt20_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh0ss5e25_dt20_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh0ss5e25_dt20_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh0ss5e25_dt20_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh0ss5e25_dt20_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh0ss5e25_dt20_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.01 nh=1"
set title ""
set yrange [0.99:1.06]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh0ss25e24_dt100.eps"
plot 'd2_const_rr_n257/nh0ss5e25_dt100_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh0ss5e25_dt100_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh0ss5e25_dt100_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh0ss5e25_dt100_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh0ss5e25_dt100_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh0ss5e25_dt100_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set title "dt=0.002 nh=1"
set title ""
set yrange [0.99:1.06]
set ylabel "Numerical/Analytic Radius"
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nh0ss25e24_dt500.eps"
plot 'd2_const_rr_n257/nh0ss5e25_dt500_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nh0ss5e25_dt500_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nh0ss5e25_dt500_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nh0ss5e25_dt500_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nh0ss5e25_dt500_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nh0ss5e25_dt500_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
################################
#### Now nh=10^-1 plots  ########
################################
set yrange [0.95:1.4]
set title "dt=0.05 nh=0.1"
set title ""
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nhm1ss25e22_dt020.eps"
plot 'd2_const_rr_n257/nhm1ss5e23_dt20_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nhm1ss5e23_dt20_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nhm1ss5e23_dt20_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nhm1ss5e23_dt20_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nhm1ss5e23_dt20_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nhm1ss5e23_dt20_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set yrange [0.95:1.4]
set title "dt=0.01 nh=0.1"
set title ""
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nhm1ss25e22_dt100.eps"
plot 'd2_const_rr_n257/nhm1ss5e23_dt100_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nhm1ss5e23_dt100_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nhm1ss5e23_dt100_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nhm1ss5e23_dt100_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nhm1ss5e23_dt100_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nhm1ss5e23_dt100_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
set yrange [0.95:1.4]
set title "dt=0.002 nh=0.1"
set title ""
set output "d2_const_rr_n257/photconsRR_wts_d2crt_n257const_nhm1ss25e22_dt500.eps"
plot 'd2_const_rr_n257/nhm1ss5e23_dt500_op1.txt' u 1:($2/$6) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_rr_n257/nhm1ss5e23_dt500_op5.txt' u 1:($2/$6) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_rr_n257/nhm1ss5e23_dt500_op6.txt' u 1:($2/$6) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_rr_n257/nhm1ss5e23_dt500_op7.txt' u 1:($2/$6) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_rr_n257/nhm1ss5e23_dt500_op8.txt' u 1:($2/$6) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_rr_n257/nhm1ss5e23_dt500_geo.txt' u 1:($2/$6) w lp lt 5 pt 5 ps 0.5 title 'Geometric'
pause -1
quit



##############################################################################
# ********** Friday 9/5/08 ************* photon conservation plots
# first a 2d plot, then 3d ones.
# ********** Thursday 15/5/08 ************* more 2d photon conservation plots
##############################################################################
set terminal postscript enhanced color eps
set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=257^2, 2d Cartesian, Point Source"
set title ""
set output "d2_const_n33/photoncons_geoop7_n33_variousdt.eps"
plot 'd2_const_n33/n33nh2ss1e28_geo_dt5e10.txt' u 1:(abs(1-$6/$7))    w lp lt 1 pt 1 title 'Geometric, 257^2, dt=5e10', \
     'd2_const_n33/n33nh2ss1e28_geo_dt1e10.txt' u 1:(abs(1-$6/$7))    w lp lt 1 pt 2 title 'Geometric, 257^2, dt=1e10', \
     'd2_const_n33/n33nh2ss1e28_geo_dt2e09.txt' u 1:(abs(1-$6/$7))    w lp lt 1 pt 3 title 'Geometric, 257^2, dt=2e09', \
     'd2_const_n33/n33nh2ss1e28_op7_dt5e10.txt' u 1:(abs(1-$6/$7))    w lp lt 3 pt 4 title 'Inv. Tau,  257^2, dt=5e10', \
     'd2_const_n33/n33nh2ss1e28_op7_dt1e10.txt' u 1:(abs(1-$6/$7))    w lp lt 3 pt 5 title 'Inv. Tau,  257^2, dt=1e10', \
     'd2_const_n33/n33nh2ss1e28_op7_dt2e09.txt' u 1:(abs(1-$6/$7))    w lp lt 3 pt 6 title 'Inv. Tau,  257^2, dt=2e09'
pause -1
set ylabel "N(ions)/N(photons)"
unset log y
set grid
set yrange [0.8:1.07]
#set yrange [*:*]
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh2ss1e28_dt1e10.eps"
plot 'd2_const_n33/n33nh2ss1e28_op5_dt1e10.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh2ss1e28_op6_dt1e10.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh2ss1e28_op7_dt1e10.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh2ss1e28_op8_dt1e10.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh2ss1e28_op9_dt1e10.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh2ss1e28_op1_dt1e10.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh2ss1e28_geo_dt1e10.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
set yrange [0.91:1.04]
#set yrange [*:*]
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh2ss1e28_dt5e10.eps"
plot 'd2_const_n33/n33nh2ss1e28_op5_dt5e10.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh2ss1e28_op6_dt5e10.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh2ss1e28_op7_dt5e10.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh2ss1e28_op8_dt5e10.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh2ss1e28_op9_dt5e10.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh2ss1e28_op1_dt5e10.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh2ss1e28_geo_dt5e10.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
set yrange [0.8:1.075]
#set yrange [*:*]
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh2ss1e28_dt2e09.eps"
plot 'd2_const_n33/n33nh2ss1e28_op5_dt2e09.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh2ss1e28_op6_dt2e09.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh2ss1e28_op7_dt2e09.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh2ss1e28_op8_dt2e09.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh2ss1e28_op9_dt2e09.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh2ss1e28_op1_dt2e09.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh2ss1e28_geo_dt2e09.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
set title "3.2 particles per cc"
set yrange [0.88:1.04]
#set yrange [*:*]
set title ""
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh0p5ss32e25_dt1e10.eps"
plot 'd2_const_n33/n33nh0p5ss32e25_op5_dt1e10.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh0p5ss32e25_op6_dt1e10.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh0p5ss32e25_op7_dt1e10.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh0p5ss32e25_op8_dt1e10.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh0p5ss32e25_op9_dt1e10.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh0p5ss32e25_op1_dt1e10.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh0p5ss32e25_geo_dt1e10.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
set yrange [0.91:1.04]
#set yrange [*:*]
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh0p5ss32e25_dt5e10.eps"
plot 'd2_const_n33/n33nh0p5ss32e25_op5_dt5e10.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh0p5ss32e25_op6_dt5e10.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh0p5ss32e25_op7_dt5e10.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh0p5ss32e25_op8_dt5e10.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh0p5ss32e25_op9_dt5e10.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh0p5ss32e25_op1_dt5e10.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh0p5ss32e25_geo_dt5e10.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
set yrange [0.88:1.04]
#set yrange [*:*]
set output "d2_const_n33/photcons_weights_d2crt_n33const_nh0p5ss32e25_dt2e09.eps"
plot 'd2_const_n33/n33nh0p5ss32e25_op5_dt2e09.txt' u 1:($6/$7) w lp lt 1 pt 1 ps 0.5 title 't-min=0.5', \
     'd2_const_n33/n33nh0p5ss32e25_op6_dt2e09.txt' u 1:($6/$7) w lp lt 2 pt 2 ps 0.5 title 't-min=0.6', \
     'd2_const_n33/n33nh0p5ss32e25_op7_dt2e09.txt' u 1:($6/$7) w lp lt 3 pt 3 ps 0.5 title 't-min=0.7', \
     'd2_const_n33/n33nh0p5ss32e25_op8_dt2e09.txt' u 1:($6/$7) w lp lt 4 pt 4 ps 0.5 title 't-min=0.8', \
     'd2_const_n33/n33nh0p5ss32e25_op9_dt2e09.txt' u 1:($6/$7) w lp lt 5 pt 5 ps 0.5 title 't-min=0.9', \
     'd2_const_n33/n33nh0p5ss32e25_op1_dt2e09.txt' u 1:($6/$7) w lp lt 7 pt 7 ps 0.5 title 't-min=0.1', \
     'd2_const_n33/n33nh0p5ss32e25_geo_dt2e09.txt' u 1:($6/$7) w lp lt 0 pt 6 ps 1 title 'Geometric'
pause -1
quit

##############################################################################
# ********** Tuesday 6/5/08 ************
# 3D 3D 3D!!! CONSTANT DENSITY PROFILE, N=65^3, RHO=var, p=var, ss=var phot/s
# These lines make plots of photon conservation and asphericity for 3D pt.src
# raytracing.
# This is testing asphericity/conservation for geometric weighting as a funtion
# of the optical depth of the cells.
# RHO= 1.67e-24*{1, 3.2, 10, 100} (so tau~{1.25,4,12.5,125} per cell)
# P_G= 2.76e-14*{1, 3.2, 10, 100} (so T0=200K)
# ss = 3.00e+44*{1, 3.33, 10, 100} (so equal volume ionised in equal times)
#
# ********** Friday 9/5/08 ************* photon conservation plots
##############################################################################
set xlabel "Time (seconds)"
set ylabel "abs(1-Nions/Nphotons)"
set key right bottom
set title "3D Photon Conservation, 65^3, Point Source, Constant Density, 1/Tau Weighting"
set title ""
set grid
unset log y
set ylabel "Nions/Nphotons"
set yrange [0.7:1.05]
set output "photcons_n65d3const_nh100_weights.eps"
plot 'd3_const_n65/n65nh100ss3e46_op5_dt1e10.txt' u 1:(abs($6/$7)) w l lt 1     title "tau=125, dt=1e10, tm=0.5", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10.txt' u 1:(abs($6/$7)) w l lt 2     title "tau=125, dt=1e10, tm=0.6", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10.txt' u 1:(abs($6/$7)) w l lt 3     title "tau=125, dt=1e10, tm=0.7", \
     'd3_const_n65/n65nh100ss3e46_op8_dt1e10.txt' u 1:(abs($6/$7)) w l lt 4     title "tau=125, dt=1e10, tm=0.8", \
     'd3_const_n65/n65nh100ss3e46_op9_dt1e10.txt' u 1:(abs($6/$7)) w l lt 5     title "tau=125, dt=1e10, tm=0.9", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10.txt' u 1:(abs($6/$7)) w l lt 0   title "tau=125, dt=1e10, geo.wt."
pause -1
unset log y
set key right bottom
set ylabel "Nions/Nphotons"
set yrange [0.9:1.02]
set output "photcons_n65d3const_nh3p2_weights.eps"
plot 'd3_const_n65/n65nh3ss1e45_op5_dt1e10.txt' u 1:(abs($6/$7)) w l lt 1     title "tau=4, dt=1e10, tm=0.5", \
     'd3_const_n65/n65nh3ss1e45_op6_dt1e10.txt' u 1:(abs($6/$7)) w l lt 2     title "tau=4, dt=1e10, tm=0.6", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10.txt' u 1:(abs($6/$7)) w l lt 3     title "tau=4, dt=1e10, tm=0.7", \
     'd3_const_n65/n65nh3ss1e45_op8_dt1e10.txt' u 1:(abs($6/$7)) w l lt 4     title "tau=4, dt=1e10, tm=0.8", \
     'd3_const_n65/n65nh3ss1e45_op9_dt1e10.txt' u 1:(abs($6/$7)) w l lt 5     title "tau=4, dt=1e10, tm=0.9", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10.txt' u 1:(abs($6/$7)) w l lt 0   title "tau=4, dt=1e10, geo.wt."
pause -1

set output "null.eps"
set terminal x11
set yrange [*:*]
unset grid
unset log y
f(x)=1.0
#
set title "3D Photon Conservation, 65^3, Point Source, Constant Density, Geometric Weighting."
set ylabel "N(ions)/N(photons)"
f(x)=1.0
unset log y
plot 'd3_const_n65/n65nh1ss3e44_geowt_dt1e10.txt'   u 1:($6/$7) w lp lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh3ss1e45_geowt_dt1e10.txt'   u 1:($6/$7) w lp lt 2     title "tau=4", \
     'd3_const_n65/n65nh10ss3e45_geowt_dt1e10.txt'  u 1:($6/$7) w lp lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh100ss3e46_geowt_dt1e10.txt' u 1:($6/$7) w lp lt 4     title "tau=125", \
     f(x) w l lt 0 title ""
pause -1
set title "3D Photon Conservation, 65^3, Point Source, Constant Density, 1/Tau Weighting, dt=1e10s."
set ylabel "N(ions)/N(photons)"
f(x)=1.0
unset log y
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10.txt'   u 1:($6/$7) w lp lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10.txt'   u 1:($6/$7) w lp lt 2     title "tau=4", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10.txt'  u 1:($6/$7) w lp lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10.txt' u 1:($6/$7) w lp lt 4     title "tau=125", \
     f(x) w l lt 0 title ""
pause -1
set log y
set ylabel "abs(1-Nions/Nphotons)"
plot 'd3_const_n65/n65nh1ss3e44_op7_dt1e10.txt'   u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=1.25", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10.txt'   u 1:(abs(1-$6/$7)) w lp lt 2     title "tau=4", \
     'd3_const_n65/n65nh10ss3e45_op7_dt1e10.txt'  u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=12.5", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp lt 4     title "tau=125"
pause -1
#
set title "3D Photon Conservation, 65^3, Point Source, Constant Density, 1/Tau Weighting"
set log y
set yrange [0.003:0.03]
set ylabel "abs(1-Nions/Nphotons)"
plot 'd3_const_n65/n65nh100ss3e46_op7_dt5e10.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=5e10", \
     'd3_const_n65/n65nh100ss3e46_op7_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=1e10", \
     'd3_const_n65/n65nh100ss3e46_op7_dt2e09.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=2e09", \
     'd3_const_n65/n65nh3ss1e45_op7_dt5e10.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=5e10", \
     'd3_const_n65/n65nh3ss1e45_op7_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=1e10", \
     'd3_const_n65/n65nh3ss1e45_op7_dt2e09.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=2e09"
pause -1
set yrange [0.0001:0.03]
set key left bottom
plot 'd3_const_n65/n65nh100ss3e46_op6_dt5e10.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=5e10", \
     'd3_const_n65/n65nh100ss3e46_op6_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=1e10", \
     'd3_const_n65/n65nh100ss3e46_op6_dt2e09.txt' u 1:(abs(1-$6/$7)) w lp lt 1     title "tau=125, dt=2e09", \
     'd3_const_n65/n65nh3ss1e45_op6_dt5e10.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=5e10", \
     'd3_const_n65/n65nh3ss1e45_op6_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=1e10", \
     'd3_const_n65/n65nh3ss1e45_op6_dt2e09.txt' u 1:(abs(1-$6/$7)) w lp lt 3     title "tau=4, dt=2e09"
pause -1

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

set title "Azimuthal Angle Dependence of I-Front Radius, n=33^3, constant density, dt=1e10, t=1e12"
set xlabel "Azimuthal Angle (degrees)"
set ylabel "Radius(at Angle)"
fp(x)=4.442831e18
fm(x)=3.862831e18
set yrange [3.8e18:4.5e18]
plot 'd3_const_n33/n33nh100ss3e46_cart3d_root3_dt1e10_phi.txt' u 1:3:5:6 w yerrorbars pt 3 lt 3 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1


set title "Polar Angle Dependence of I-Front Radius, n=33^3, constant density, dt=1e10, t=1e12" 
set xlabel "Polar Angle (degrees)"
set ylabel "Radius of I-Front"
fp(x)=4.442831e18
fm(x)=3.862831e18
set grid
set yrange [3.8e18:4.5e18]
plot  'd3_const_n33/n33nh100ss3e46_cart3d_root3_dt1e10_theta.txt' u 1:3:5:6 w yerrorbars pt 3 lt 3 title "IF Radius", \
     f(x) title "Analytic Radius", fp(x) title "R+dR", fm(x) title "R-dR"
pause -1

set yrange [*:*]
quit

set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation 3d Cartesian, Point Source, Constant Density"
plot 'd3_const_n33/n33nh100ss3e46_3d_root3.txt'       u 1:($6/$7) w l title 'Tau-1/root3, dt=1e10, n=33^3' lt 2, \
     'd3_const_n33/n33nh100ss3e46_cart3d_root3_dt1e09.txt' u 1:($6/$7) w l title 'Tau-1/root3, dt=1e09, n=33^3' lt 6, \
     'd3_const_n129/n129nh100ss3e46_cart3d_root3_dt1e10.txt' u 1:($6/$7) w l title 'Tau-1/root3, dt=1e10, n=129^3' lt 3, \
     'd3_const_n129/n129nh100ss3e46_cart3d_root3_dt2e09.txt' u 1:($6/$7) w l title 'Tau-1/root3, dt=2e09, n=129^3' lt 4, \
     'd3_const_n129/n129nh100ss3e46_cart3d_op7_dt1e10.txt' u 1:($6/$7) w l title 'Tau-0.7, dt=1e10, n=129^3' lt 5, \
     f(x) lt 0 title ""
pause -1

set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=33^3, 3d Cartesian, Point Source, Constant Density"
plot 'd3_const_n33/n33nh100ss3e46_3d_geo.txt' u 1:(abs(1-$6/$7)) w l title 'Geo-Weighted, dt=1e10' lt 0, \
     'd3_const_n33/n33nh100ss3e46_3d_root2.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-1/root2, dt=1e10' lt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root3.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-1/root3, dt=1e10' lt 2, \
     'd3_const_n33/n33nh100ss3e46_3d_p6.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-0.6, dt=1e10' lt 3, \
     'd3_const_n33/n33nh100ss3e46_3d_p5.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-0.5, dt=1e10' lt 4, \
     'd3_const_n33/n33nh100ss3e46_3d_p5_dt1e9.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-0.5, dt=1e9' lt 5
pause -1

#set output "n33d3_const_photon_cons2.eps"
set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation n=33^3, 3d Cartesian, Point Source, Constant Density"
plot 'd3_const_n33/n33nh100ss3e46_3d_geo.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root2.txt' u 1:($6/$7) w lp title 'Tau-1/root2, dt=1e10' lt 1 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root3.txt' u 1:($6/$7) w lp title 'Tau-1/root3, dt=1e10' lt 2 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root3_dt1e9.txt' u 1:($6/$7) w l title 'Tau-1/root3, dt=1e09' lt 6, \
     'd3_const_n33/n33nh100ss3e46_3d_p6.txt' u 1:($6/$7) w l title 'Tau-0.6, dt=1e10' lt 3, \
     'd3_const_n33/n33nh100ss3e46_3d_p5.txt' u 1:($6/$7) w l title 'Tau-0.5, dt=1e10' lt 4, \
     'd3_const_n33/n33nh100ss3e46_3d_p5_dt1e9.txt' u 1:($6/$7) w l title 'Tau-0.5, dt=1e09' lt 5, f(x) title ""
pause -1

#set output "n33d3_const_dRoverR.eps"
set xlabel "Radius (cm)"
set key right top
unset log y
f(x)=3.91e17/x
set ylabel "I-Front Width (dR/R)"
set title "Photon Conservation n=33^3, 3d Cartesian, Point Source, Constant Density"
plot 'd3_const_n33/n33nh100ss3e46_3d_geo.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root2.txt' u 2:(($3-$4)/$2) w lp title 'Tau-1/root2, dt=1e10' lt 1 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root3.txt' u 2:(($3-$4)/$2) w lp title 'Tau-1/root3, dt=1e10' lt 2 pt 1, \
     f(x) title "Cell Size"
pause -1
#set output "n33d3_const_dR.eps"
set xlabel "Radius (cm)"
set key left top
unset log y
f(x)=3.91e17
g(x)=1.0e17
set ylabel "I-Front Width (dR)"
set title "Photon Conservation n=33^3, 3d Cartesian, Point Source, Constant Density"
plot 'd3_const_n33/n33nh100ss3e46_3d_geo.txt' u 2:(($3-$4)) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root2.txt' u 2:(($3-$4)) w lp title 'Tau-1/root2, dt=1e10' lt 1 pt 1, \
     'd3_const_n33/n33nh100ss3e46_3d_root3.txt' u 2:(($3-$4)) w lp title 'Tau-1/root3, dt=1e10' lt 2 pt 1, \
     'd3_const_n129/n129nh100ss3e46_cart3d_root3_dt1e10.txt' u 2:($3-$4) w l title 'Tau-1/root3, dt=1e10, n=129^3' lt 3, \
     'd3_const_n129/n129nh100ss3e46_cart3d_root3_dt2e09.txt' u 2:($3-$4) w l title 'Tau-1/root3, dt=2e09, n=129^3' lt 4, \
     f(x) title "Cell Size 33^3", g(x) title "Cell Size 129^3"
pause -1
quit

##############################################################################
# CONSTANT DENSITY PROFILE, N=33, RHO=1.67e-22, p=2.76e-12, ss=1e28 phot/s
# These lines make plots of photon conservation and asphericity for 2D pt.src
# raytracing.
##############################################################################
set terminal postscript enhanced color
set output "n33d2_const_photon_cons.eps"
#set term x11
set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=33^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n33/n33n100ss28_dt5e10.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd2_const_n33/n33n100ss28_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'd2_const_n33/n33n100ss28_dt2e09.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=2e09' lt 2 pt 2, \
     'd2_const_n33/n33n100ss28_dt1e09.txt' u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=1e09' lt 3
pause -1

set output "n33d2_const_photon_cons2.eps"
set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation n=33^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n33/n33n100ss28_dt5e10.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd2_const_n33/n33n100ss28_dt1e10.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'd2_const_n33/n33n100ss28_dt2e09.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=2e09' lt 2 pt 2, \
     'd2_const_n33/n33n100ss28_dt1e09.txt' u 1:($6/$7) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     f(x) title "x=1"
pause -1

set output "n33d2_const_dRoverR.eps"
set xlabel "Radius (cm)"
set key right bottom
set log y
f(x)=3.91e17/x
set ylabel "I-Front Width (dR/R)"
set title "Photon Conservation n=33^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n33/n33n100ss28_dt5e10.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd2_const_n33/n33n100ss28_dt1e10.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'd2_const_n33/n33n100ss28_dt2e09.txt' u 2:(($3-$4)/$2) w lp title 'Tau-Weighted, dt=2e09' lt 2 pt 2, \
     'd2_const_n33/n33n100ss28_dt1e09.txt' u 2:(($3-$4)/$2) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     f(x) title "Cell Size"
pause -1
set output "n33d2_const_dR.eps"
set xlabel "Radius (cm)"
set key right bottom
set log y
f(x)=3.91e17
set ylabel "I-Front Width (dR)"
set title "Photon Conservation n=33^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n33/n33n100ss28_dt5e10.txt' u 2:(($3-$4)) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'd2_const_n33/n33n100ss28_dt1e10.txt' u 2:(($3-$4)) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'd2_const_n33/n33n100ss28_dt2e09.txt' u 2:(($3-$4)) w l title 'Tau-Weighted, dt=2e09' lt 3, \
     'd2_const_n33/n33n100ss28_dt1e09.txt' u 2:(($3-$4)) w l title 'Tau-Weighted, dt=1e09' lt 4, \
     f(x) title "Cell Size"
pause -1
quit

##############################################################################
# 1/R^2 DENSITY PROFILE
# These lines make plots of photon conservation and asphericity for 2D pt.src
# raytracing.
##############################################################################
set terminal postscript enhanced color
set output "n129d2_rm2_photon_cons.eps"
#set term x11
set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=129^2, 2d Cartesian, Point Source, 1/r^2 profile"
plot 'RRadial2_dt3e10.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=3e10' lt 0 pt 1, \
     'RRadial2_dt1e10.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'RRadial2_dt3e9.txt'  u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=3e09' lt 2 pt 2, \
     'RRadial2_dt1e9.txt'  u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     'RRadial2_dt3e8.txt'  u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=3e08' lt 4, \
     'RRadial2_dt1e8.txt'  u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=1e08' lt 5
pause -1

set output "n129d2_rm2_photon_cons2.eps"
set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation n=129^2, 2d Cartesian, Point Source, 1/r^2 profile"
plot 'RRadial2_dt3e10.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=3e10' lt 0 pt 1, \
     'RRadial2_dt1e10.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'RRadial2_dt3e9.txt'  u 1:($6/$7) w lp title 'Tau-Weighted, dt=3e09' lt 2 pt 2, \
     'RRadial2_dt1e9.txt'  u 1:($6/$7) w l title 'Tau-Weighted, dt=1e09' lt 3, \
     'RRadial2_dt3e8.txt'  u 1:($6/$7) w l title 'Tau-Weighted, dt=3e08' lt 4, \
     'RRadial2_dt1e8.txt'  u 1:($6/$7) w l title 'Tau-Weighted, dt=1e08' lt 5, \
     f(x) title "x=1"
pause -1

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
# These file are moved to d2_const_n257/*.txt now
# These lines make plots of photon conservation and asphericity for 2D pt.src
# raytracing.
##############################################################################
set terminal postscript enhanced color
set output "n257d2_photon_cons.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set ylabel "abs(1-Nions/Nphotons)"
set title "Photon Conservation n=257^2, 2d Cartesian, Point Source"
plot 'test/n257dt5e10.txt' u 1:(abs(1-$6/$7)) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n257dt1e10.txt' u 1:(abs(1-$6/$7)) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n257dt2e9.txt'  u 1:(abs(1-$6/$7)) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n257dt5e10_new.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n257dt1e10_new.txt' u 1:(abs(1-$6/$7)) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n257dt2e9_new.txt'  u 1:(abs(1-$6/$7)) w l title 'Tau-Weighted, dt=2e09' lt 4
pause -1

set output "n257d2_const_photon_cons2.eps"
set xlabel "Time (seconds)"
set key right bottom
unset log y
f(x)=1.0
set ylabel "N(ions)/N(photons)"
set title "Photon Conservation n=257^2, 2d Cartesian, Point Source, Constant Density"
plot 'd2_const_n257/n257dt5e10_new.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=5e10' lt 0 pt 1, \
     'd2_const_n257/n257dt1e10_new.txt' u 1:($6/$7) w lp title 'Tau-Weighted, dt=1e10' lt 1 pt 1, \
     'd2_const_n257/n257dt2e9_new.txt' u 1:($6/$7) w l title 'Tau-Weighted, dt=2e09' lt 2, \
     'd2_const_n257/n257dt5e8_new.txt' u 1:($6/$7) w l title 'Tau-Weighted, dt=5e08' lt 3, \
     f(x) title "x=1"
pause -1

set terminal postscript enhanced color
set output "n257d2_fittoradius.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set yrange [1.e-4:0.1]
set ylabel "abs(1-R(fit)/R(Analytic))"
set title "n=257^2, 2d Cartesian; Fit Radius (mean radius of cells with 0.1<x(HII)<0.9) vs. Analytic Radius"
plot 'test/n257dt5e10.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n257dt1e10.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n257dt2e9.txt'  u 1:(abs(1-$2/sqrt($7/100/3.14159))) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n257dt5e10_new.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n257dt1e10_new.txt' u 1:(abs(1-$2/sqrt($7/100/3.14159))) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n257dt2e9_new.txt'  u 1:(abs(1-$2/sqrt($7/100/3.14159))) w l title 'Tau-Weighted, dt=2e09' lt 4
pause -1

set terminal postscript enhanced color
set output "n257d2_asphericity.eps"
set xlabel "Time (seconds)"
set key right bottom
set log y
set yrange [1.e16:1.e18]
f(x)=5.02e16
set ylabel "I-Front Width"
set title "n=257^2, 2d Cartesian; Asphericity: Max. Radius Minus Min. Radius"
plot 'test/n257dt5e10.txt' u 1:($3-$4) title 'Geometric, dt=5e10' w lp lt 3 pt 1, \
     'test/n257dt1e10.txt' u 1:($3-$4) title 'Geometric, dt=1e10' w lp lt 3 pt 2, \
     'test/n257dt2e9.txt'  u 1:($3-$4) title 'Geometric, dt=2e09' w lp lt 3 pt 4, \
     'test/n257dt5e10_new.txt' u 1:($3-$4) w lp title 'Tau-Weighted, dt=5e10' lt 1 pt 1, \
     'test/n257dt1e10_new.txt' u 1:($3-$4) w lp title 'Tau-Weighted, dt=1e10' lt 2 pt 2, \
     'test/n257dt2e9_new.txt'  u 1:($3-$4) w l title 'Tau-Weighted, dt=2e09' lt 4, \
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
plot 'profile_dt01_n33.txt'  u 1:(abs($5/$6)) w lp title "n=33^2", \
     'profile_dt01_n65.txt'  u 1:(abs($5/$6)) w lp title "n=65^2", \
     'profile_dt01_n129.txt' u 1:(abs($5/$6)) w lp title "n=129^2", \
     'profile_dt01_n257.txt' u 1:(abs($5/$6)) w lp title "n=257^2", f(x) title ""
pause -1

set output "photons_dxdt.eps"
set xlabel "Time (seconds)"
set title "2D Point Source, No Recombinations, Photon Conservation as Function of Resolution and Timestep"
set ylabel "ion number/photon number"
f(x)=1
set yrange [0.95:1.02]
plot 'profile_dt02_n33.txt'  u 1:(abs($5/$6)) w lp title "n=33^2, dt=2e10", \
     'profile_dt01_n65.txt'  u 1:(abs($5/$6)) w lp title "n=65^2, dt=1e10", \
     'profile_dt005_n129.txt' u 1:(abs($5/$6)) w lp title "n=129^2, dt=5e9", \
     'profile_dt0025_n257.txt' u 1:(abs($5/$6)) w lp title "n=257^2, dt=2.5e9", f(x) title ""
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
plot 'profile_dt02_n33.txt'  u 1:(1.0*$2) w lp title "n=33^2, dt=2e10, fit", \
     'profile_dt02_n33.txt'  u 1:4 w l title "n=33^2, dt=2e10, nions", \
     'profile_dt0025_n257.txt'  u 1:(sqrt($6/100/3.14159)) w l title "Analytic", \
     f(x) title "K*sqrt(x)", \
     'profile_dt01_n65.txt'  u 1:2 w lp title "n=65^2, dt=1e10, fit", \
     'profile_dt005_n129.txt' u 1:2 w lp title "n=129^2, dt=5e9, fit", \
     'profile_dt0025_n257.txt' u 1:2 w lp title "n=257^2, dt=2.5e9, fit"
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

