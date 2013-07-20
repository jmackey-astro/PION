#!/bin/bash

DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_asciiV2

n00080=`ls ${DDIR}/S1D_n00080.*.txt | head -n101 | tail -n1`
n00160=`ls ${DDIR}/S1D_n00160.*.txt | head -n101 | tail -n1`
n00320=`ls ${DDIR}/S1D_n00320.*.txt | head -n101 | tail -n1`
n00640=`ls ${DDIR}/S1D_n00640.*.txt | head -n101 | tail -n1`
n01280=`ls ${DDIR}/S1D_n01280.*.txt | head -n101 | tail -n1`
n02560=`ls ${DDIR}/S1D_n02560.*.txt | head -n101 | tail -n1`
n05120=`ls ${DDIR}/S1D_n05120.*.txt | head -n101 | tail -n1`
n10240=`ls ${DDIR}/S1D_n10240.*.txt | head -n101 | tail -n1`
n20480=`ls ${DDIR}/S1D_n20480.*.txt | head -n101 | tail -n1`
n40960=`ls ${DDIR}/S1D_n40960.*.txt | head -n101 | tail -n1`

echo $n05120
grep time $n00080
grep time $n05120

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
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120", \
     "${n10240}" u (\$1/3.086e18):var w l lc 0 lt 4 lw 2 title "N_x=10240", \
     "${n20480}" u (\$1/3.086e18):var w l lc 2 lt 4 lw 2 title "N_x=20480", \
     "${n40960}" u (\$1/3.086e18):var w l lc 3 lt 4 lw 2 title "N_x=40960"

set output "Velocity_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.1 Myr" offset 0,-0.6
unset log y
var=4
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"

set output "IonFraction_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.1 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"

EOF

n00080=`ls ${DDIR}/S1D_n00080.*.txt | head -n51 | tail -n1`
n00160=`ls ${DDIR}/S1D_n00160.*.txt | head -n51 | tail -n1`
n00320=`ls ${DDIR}/S1D_n00320.*.txt | head -n51 | tail -n1`
n00640=`ls ${DDIR}/S1D_n00640.*.txt | head -n51 | tail -n1`
n01280=`ls ${DDIR}/S1D_n01280.*.txt | head -n51 | tail -n1`
n02560=`ls ${DDIR}/S1D_n02560.*.txt | head -n51 | tail -n1`
n05120=`ls ${DDIR}/S1D_n05120.*.txt | head -n51 | tail -n1`
n10240=`ls ${DDIR}/S1D_n10240.*.txt | head -n51 | tail -n1`
n20480=`ls ${DDIR}/S1D_n20480.*.txt | head -n51 | tail -n1`
n40960=`ls ${DDIR}/S1D_n40960.*.txt | head -n51 | tail -n1`

grep time $n05120

cat << EOF >> temp.gp
set output "Density_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set yrange [*:*]
set ylabel "Density (g/cm^{3})" offset 2.5,0
set title "Density field at t = 0.05 Myr" offset 0,-0.6
set log y
var=2
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"


set output "Velocity_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.05 Myr" offset 0,-0.6
unset log y
var=4
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"

set output "IonFraction_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.05 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"


EOF

n00080=`ls ${DDIR}/S1D_n00080.*.txt | head -n26 | tail -n1`
n00160=`ls ${DDIR}/S1D_n00160.*.txt | head -n26 | tail -n1`
n00320=`ls ${DDIR}/S1D_n00320.*.txt | head -n26 | tail -n1`
n00640=`ls ${DDIR}/S1D_n00640.*.txt | head -n26 | tail -n1`
n01280=`ls ${DDIR}/S1D_n01280.*.txt | head -n26 | tail -n1`
n02560=`ls ${DDIR}/S1D_n02560.*.txt | head -n26 | tail -n1`
n05120=`ls ${DDIR}/S1D_n05120.*.txt | head -n26 | tail -n1`
n10240=`ls ${DDIR}/S1D_n10240.*.txt | head -n26 | tail -n1`
n20480=`ls ${DDIR}/S1D_n20480.*.txt | head -n26 | tail -n1`
n40960=`ls ${DDIR}/S1D_n40960.*.txt | head -n26 | tail -n1`

grep time $n05120

cat << EOF >> temp.gp
set output "Density_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set yrange [*:*]
set ylabel "Density (g/cm^{3})" offset 2.5,0
set title "Density field at t = 0.025 Myr" offset 0,-0.6
set log y
var=2
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"


set output "Velocity_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "Velocity (cm/s)" offset 2.5,0
set title "Velocity field at t = 0.025 Myr" offset 0,-0.6
unset log y
var=4
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"

set output "IonFraction_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.025 Myr" offset 0,-0.6
set log y
set yrange [5e-6:1.1]
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_x=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_x=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_x=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_x=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_x=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_x=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_x=05120"

quit
EOF

gnuplot temp.gp

exit


