#!/bin/bash

DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_ascii

n0080=`ls ${DDIR}/S1D_n0080.*.txt | head -n101 | tail -n1`
n0160=`ls ${DDIR}/S1D_n0160.*.txt | head -n101 | tail -n1`
n0320=`ls ${DDIR}/S1D_n0320.*.txt | head -n101 | tail -n1`
n0640=`ls ${DDIR}/S1D_n0640.*.txt | head -n101 | tail -n1`
n1280=`ls ${DDIR}/S1D_n1280.*.txt | head -n101 | tail -n1`
n2560=`ls ${DDIR}/S1D_n2560.*.txt | head -n101 | tail -n1`
n5120=`ls ${DDIR}/S1D_n5120.*.txt | head -n101 | tail -n1`

echo $n5120
grep time $n0080
grep time $n5120

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times-New-Roman,22" 
#\
# fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set lmargin 7.5
set tmargin 1.5
set rmargin 1.25
set bmargin 2.5
set key top left 
set yrange [*:*]
unset log y
unset log x
set xlabel "r (pc)" offset 0,0.6

set output "Density_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set yrange [*:*]
set ylabel "Density (g/cm^{3})" offset 2.5,0
set title "Density field at t = 0.1 Myr" offset 0,-0.6
set log y
var=2
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

set output "Velocity_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.1 Myr" offset 0,-0.6
unset log y
var=4
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

set output "IonFraction_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.1 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

EOF

n0080=`ls ${DDIR}/S1D_n0080.*.txt | head -n51 | tail -n1`
n0160=`ls ${DDIR}/S1D_n0160.*.txt | head -n51 | tail -n1`
n0320=`ls ${DDIR}/S1D_n0320.*.txt | head -n51 | tail -n1`
n0640=`ls ${DDIR}/S1D_n0640.*.txt | head -n51 | tail -n1`
n1280=`ls ${DDIR}/S1D_n1280.*.txt | head -n51 | tail -n1`
n2560=`ls ${DDIR}/S1D_n2560.*.txt | head -n51 | tail -n1`
n5120=`ls ${DDIR}/S1D_n5120.*.txt | head -n51 | tail -n1`

grep time $n5120

cat << EOF >> temp.gp
set output "Density_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set yrange [*:*]
set ylabel "Density (g/cm^{3})" offset 2.5,0
set title "Density field at t = 0.05 Myr" offset 0,-0.6
set log y
var=2
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"


set output "Velocity_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.05 Myr" offset 0,-0.6
unset log y
var=4
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

set output "IonFraction_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.05 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"


EOF

n0080=`ls ${DDIR}/S1D_n0080.*.txt | head -n26 | tail -n1`
n0160=`ls ${DDIR}/S1D_n0160.*.txt | head -n26 | tail -n1`
n0320=`ls ${DDIR}/S1D_n0320.*.txt | head -n26 | tail -n1`
n0640=`ls ${DDIR}/S1D_n0640.*.txt | head -n26 | tail -n1`
n1280=`ls ${DDIR}/S1D_n1280.*.txt | head -n26 | tail -n1`
n2560=`ls ${DDIR}/S1D_n2560.*.txt | head -n26 | tail -n1`
n5120=`ls ${DDIR}/S1D_n5120.*.txt | head -n26 | tail -n1`

grep time $n5120

cat << EOF >> temp.gp
set output "Density_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set yrange [*:*]
set ylabel "Density (g/cm^{3})" offset 2.5,0
set title "Density field at t = 0.025 Myr" offset 0,-0.6
set log y
var=2
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"


set output "Velocity_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.025 Myr" offset 0,-0.6
unset log y
var=4
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

set output "IonFraction_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.025 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n0080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=0080", \
     "${n0160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=0160", \
     "${n0320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=0320", \
     "${n0640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=0640", \
     "${n1280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=1290", \
     "${n2560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=2560", \
     "${n5120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=5120"

quit
EOF

gnuplot temp.gp

exit


