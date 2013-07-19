#!/bin/bash

DDIR=/vol/klaipeda3/scratch/jmackey/HII1Dv2

n00128=`ls ${DDIR}/HII1Dv2_n00128.*.txt | head -n2 | tail -n1`
n00256=`ls ${DDIR}/HII1Dv2_n00256.*.txt | head -n2 | tail -n1`
n00512=`ls ${DDIR}/HII1Dv2_n00512.*.txt | head -n2 | tail -n1`
n01024=`ls ${DDIR}/HII1Dv2_n01024.*.txt | head -n2 | tail -n1`
n02048=`ls ${DDIR}/HII1Dv2_n02048.*.txt | head -n2 | tail -n1`
n04096=`ls ${DDIR}/HII1Dv2_n04096.*.txt | head -n2 | tail -n1`
n08192=`ls ${DDIR}/HII1Dv2_n08192.*.txt | head -n2 | tail -n1`
n16384=`ls ${DDIR}/HII1Dv2_n16384.*.txt | head -n2 | tail -n1`
n32768=`ls ${DDIR}/HII1Dv2_n32768.*.txt | head -n2 | tail -n1`

echo $n08192
grep time $n00128
grep time $n08192

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times-New-Roman,22" \
 fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set lmargin 5.5
set tmargin 1.5
set rmargin 1.25
set bmargin 2.5
set key top left 
set yrange [*:*]
unset log y
unset log x
set xlabel "r (pc)" offset 0,0.6

set output "HII1D_Density_t0p1Myr_res.eps"
set xrange [7.25:7.8]
set key top left 
set yrange [0:80]
set ylabel "n_{H} (cm^{-3})" offset 2,0
set title "Density field at t = 0.1 Myr" offset 0,-0.6
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 9 lt 1 lw 2 title "N_x=32768"



set output "HII1D_Velocity_t0p1Myr_res.eps"
set key top left 
set ylabel "Velocity (km/s)" offset 2.5,0
set title "Velocity field at t = 0.1 Myr" offset 0,-0.6
set yrange [-1:10]
set lmargin 4.5
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$4*1e-5) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$4*1e-5) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$4*1e-5) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$4*1e-5) w l lc 9 lt 1 lw 2 title "N_x=32768"


set output "HII1D_IonFraction_t0p1Myr_res.eps"
set key top left 
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.1 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
set lmargin 7.5
var=7
plot "${n00128}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):var w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):var w l lc 9 lt 1 lw 2 title "N_x=32768"

set output "HII1D_Pressure_t0p1Myr_res.eps"
set key bottom left 
set ylabel "P (10^{-12} dyne/cm^2)" offset 2,0
set title "Gas Pressure at t = 0.1 Myr" offset 0,-0.6
unset log y
set yrange [*:*]
set lmargin 5.5
plot "${n00128}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 9 lt 1 lw 2 title "N_x=32768"

EOF



n00128=`ls ${DDIR}/HII1Dv2_n00128.*.txt | head -n11 | tail -n1`
n00256=`ls ${DDIR}/HII1Dv2_n00256.*.txt | head -n11 | tail -n1`
n00512=`ls ${DDIR}/HII1Dv2_n00512.*.txt | head -n11 | tail -n1`
n01024=`ls ${DDIR}/HII1Dv2_n01024.*.txt | head -n11 | tail -n1`
n02048=`ls ${DDIR}/HII1Dv2_n02048.*.txt | head -n11 | tail -n1`
n04096=`ls ${DDIR}/HII1Dv2_n04096.*.txt | head -n11 | tail -n1`
n08192=`ls ${DDIR}/HII1Dv2_n08192.*.txt | head -n11 | tail -n1`
n16384=`ls ${DDIR}/HII1Dv2_n16384.*.txt | head -n11 | tail -n1`
n32768=`ls ${DDIR}/HII1Dv2_n32768.*.txt | head -n11 | tail -n1`

grep time $n08192

cat << EOF >> temp.gp
set output "HII1D_Density_t1Myr_res.eps"
set lmargin 5.5
set xrange [13:15.3]
set key top left 
set yrange [0:50]
set ylabel "n_{H} (cm^{-3})" offset 2,0
set title "Density field at t = 1.0 Myr" offset 0,-0.6
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 9 lt 1 lw 2 title "N_x=32768"



set output "HII1D_Velocity_t1Myr_res.eps"
set key bottom left 
set ylabel "Velocity (km/s)" offset 1.5,0
set title "Velocity field at t = 1.0 Myr" offset 0,-0.6
set yrange [-1:7]
set lmargin 4.5
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$4*1e-5) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$4*1e-5) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$4*1e-5) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$4*1e-5) w l lc 9 lt 1 lw 2 title "N_x=32768"


set output "HII1D_IonFraction_t1Myr_res.eps"
set key top left 
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 1.0 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
set lmargin 7.5
var=7
plot "${n00128}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):var w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):var w l lc 9 lt 1 lw 2 title "N_x=32768"

set output "HII1D_Pressure_t1Myr_res.eps"
set key bottom left 
set ylabel "P (10^{-12} dyne/cm^2)" offset 1.5,0
set title "Gas Pressure at t = 1.0 Myr" offset 0,-0.6
unset log y
set yrange [1:8]
set lmargin 5.5
plot "${n00128}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 9 lt 1 lw 2 title "N_x=32768"

EOF

n00128=`ls ${DDIR}/HII1Dv2_n00128.*.txt | head -n26 | tail -n1`
n00256=`ls ${DDIR}/HII1Dv2_n00256.*.txt | head -n26 | tail -n1`
n00512=`ls ${DDIR}/HII1Dv2_n00512.*.txt | head -n26 | tail -n1`
n01024=`ls ${DDIR}/HII1Dv2_n01024.*.txt | head -n26 | tail -n1`
n02048=`ls ${DDIR}/HII1Dv2_n02048.*.txt | head -n26 | tail -n1`
n04096=`ls ${DDIR}/HII1Dv2_n04096.*.txt | head -n26 | tail -n1`
n08192=`ls ${DDIR}/HII1Dv2_n08192.*.txt | head -n26 | tail -n1`
n16384=`ls ${DDIR}/HII1Dv2_n16384.*.txt | head -n26 | tail -n1`
n32768=`ls ${DDIR}/HII1Dv2_n32768.*.txt | head -n26 | tail -n1`

grep time $n08192

cat << EOF >> temp.gp
set output "HII1D_Density_t2p5Myr_res.eps"
set xrange [18:25]
set key top left 
set lmargin 5.5
set yrange [0:30]
set ylabel "n_{H} (cm^{-3})" offset 2,0
set title "Density field at t = 2.5 Myr" offset 0,-0.6
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 9 lt 1 lw 2 title "N_x=32768"



set output "HII1D_Velocity_t2p5Myr_res.eps"
set key top left 
set ylabel "Velocity (km/s)" offset 2.5,0
set title "Velocity field at t = 2.5 Myr" offset 0,-0.6
set yrange [-0.1:4]
set lmargin 5.5
unset log y
plot "${n00128}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$4*1e-5) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$4*1e-5) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$4*1e-5) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$4*1e-5) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$4*1e-5) w l lc 9 lt 1 lw 2 title "N_x=32768"


set output "HII1D_IonFraction_t2p5Myr_res.eps"
set key top left 
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 2.5 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
set lmargin 7.5
var=7
plot "${n00128}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):var w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):var w l lc 9 lt 1 lw 2 title "N_x=32768"

set output "HII1D_Pressure_t2p5Myr_res.eps"
set key top left 
set ylabel "P (10^{-12} dyne/cm^2)" offset 2.5,0
set title "Gas Pressure at t = 2.5 Myr" offset 0,-0.6
unset log y
set yrange [1:5]
set lmargin 5.5
plot "${n00128}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 4 lw 2 title "N_x=00128", \
     "${n00256}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 4 lw 2 title "N_x=00256", \
     "${n00512}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 3 lt 4 lw 2 title "N_x=00512", \
     "${n01024}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 0 lt 1 lw 2 title "N_x=01024", \
     "${n02048}" u (\$1/3.086e18):(\$3*1.0e12) w lp lc 2 lt 1 lw 2 title "N_x=02048", \
     "${n04096}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 3 lt 1 lw 2 title "N_x=04096", \
     "${n08192}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 5 lt 1 lw 2 title "N_x=08192", \
     "${n16384}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 8 lt 1 lw 2 title "N_x=16384", \
     "${n32768}" u (\$1/3.086e18):(\$3*1.0e12) w l lc 9 lt 1 lw 2 title "N_x=32768"

quit
EOF

gnuplot temp.gp
~/active/bin/eps2jpeg.sh HII1D

exit


