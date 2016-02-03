#!/bin/bash

DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_asciiV2
#DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_siloV2

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

echo "**** T=0.1 Myr****"
grep time $n00080
grep time $n00160
grep time $n00320
grep time $n00640
grep time $n01280
grep time $n02560
grep time $n05120
grep time $n10240
grep time $n20480
grep time $n40960

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times-New-Roman,22" \
 fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set tmargin 1.5
set rmargin 1.25
set bmargin 2.5
set key top left 
set yrange [*:*]
unset log y
unset log x
set xlabel "r (pc)" offset 0,0.6

set style line 10 lt 1 lc rgb "gray20"    lw 3 pt 0 ps 1
set style line 11 lt 1 lc rgb "gray20"    lw 2 pt 1 ps 1
set style line 14 lt 1 lc rgb "gray20"    lw 2 pt 4 ps 1
set style line 20 lt 1 lc rgb "blue"      lw 3 pt 0 ps 1
set style line 21 lt 1 lc rgb "blue"      lw 2 pt 1 ps 1
set style line 24 lt 1 lc rgb "blue"      lw 2 pt 4 ps 1
set style line 30 lt 1 lc rgb "web-green" lw 3 pt 0 ps 1
set style line 31 lt 1 lc rgb "web-green" lw 2 pt 1 ps 1
set style line 34 lt 1 lc rgb "web-green" lw 2 pt 4 ps 1

set output "figS1D_Density_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set yrange [3e2:5e5]
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set title "Density field at t = 0.1 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set key top left 
plot "${n00080}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Density_t0p1Myr_res.eps"
unset title
set xrange [1.125:1.145]
plot "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 10 title "N_r=40960"

set output "figS1D_Velocity_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "Velocity (km/s)" offset 1.5,0
set title "Velocity field at t = 0.1 Myr" offset 0,-0.6
set yrange [*:*]
unset log y
set lmargin 4.5
plot "${n00080}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Velocity_t0p1Myr_res.eps"
unset title
set xrange [1.125:1.145]
plot "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 10 title "N_r=40960"

set output "figS1D_IonFraction_t0p1Myr_res.eps"
set xrange [1.03:1.15]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.1 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set yrange [5e-6:1.1]
set key bottom left
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_IonFraction_t0p1Myr_res.eps"
unset title
set xrange [1.128:1.1365]
set key bottom left
plot "${n05120}" u (\$1/3.086e18):var w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l ls 10 title "N_r=40960"

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

echo "**** T=0.05 Myr****"
grep time $n00080
grep time $n00160
grep time $n00320
grep time $n00640
grep time $n01280
grep time $n02560
grep time $n05120
grep time $n10240
grep time $n20480
grep time $n40960

cat << EOF >> temp.gp
set output "figS1D_Density_t0p05Myr_res.eps"
set xrange [0.74:0.85]
set yrange [3e2:5e5]
set title "Density field at t = 0.05 Myr" offset 0,-0.6
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set log y
set lmargin 7.5
set key top left 
plot "${n00080}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Density_t0p05Myr_res.eps"
unset title
set xrange [0.796:0.8047]
plot "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 10 title "N_r=40960"

set output "figS1D_Velocity_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "Velocity (km/s)" offset 2.5,0
set title "Velocity field at t = 0.05 Myr" offset 0,-0.6
set yrange [*:*]
unset log y
set lmargin 4.5
plot "${n00080}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Velocity_t0p05Myr_res.eps"
unset title
set key bottom left
set xrange [0.796:0.8047]
plot "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 10 title "N_r=40960"

set output "figS1D_IonFraction_t0p05Myr_res.eps"
set xrange [0.7:0.85]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.05 Myr" offset 0,-0.6
set key top left
set log y
set lmargin 7.5
set yrange [5e-6:1.1]
set key bottom left
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_IonFraction_t0p05Myr_res.eps"
unset title
set xrange [0.796:0.8025]
set key bottom left
plot "${n05120}" u (\$1/3.086e18):var w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l ls 10 title "N_r=40960"


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

echo "**** T=0.025 Myr****"
grep time $n00080
grep time $n00160
grep time $n00320
grep time $n00640
grep time $n01280
grep time $n02560
grep time $n05120
grep time $n10240
grep time $n20480
grep time $n40960

cat << EOF >> temp.gp
set output "figS1D_Density_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set yrange [3e2:5e5]
set title "Density field at t = 0.025 Myr" offset 0,-0.6
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set log y
set lmargin 7.5
set key top left 
plot "${n00080}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$2/1.67e-24) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Density_t0p025Myr_res.eps"
unset title
set xrange [0.58:0.59]
plot "${n05120}" u (\$1/3.086e18):(\$2/1.67e-24) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$2/1.67e-24) w l ls 10 title "N_r=40960"

set output "figS1D_Velocity_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "Velocity (km/s)" offset 2.5,0
set title "Velocity field at t = 0.025 Myr" offset 0,-0.6
set yrange [*:*]
unset log y
set lmargin 4.5
plot "${n00080}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):(\$4*1.0e-5) w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_Velocity_t0p025Myr_res.eps"
unset title
set key bottom left
set xrange [0.581:0.5887]
plot "${n05120}" u (\$1/3.086e18):(\$4*1.0e-5) w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):(\$4*1.0e-5) w l ls 10 title "N_r=40960"

set output "figS1D_IonFraction_t0p025Myr_res.eps"
set xrange [0.52:0.62]
set ylabel "H^{+} Fraction" offset 3.5,0
set title "H^{+} Fraction at t = 0.025 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set yrange [5e-6:1.1]
set key bottom left
var=7
plot "${n00080}" u (\$1/3.086e18):var w lp lc 0 lt 4 lw 2 title "N_r=00080", \
     "${n00160}" u (\$1/3.086e18):var w lp lc 2 lt 4 lw 2 title "N_r=00160", \
     "${n00320}" u (\$1/3.086e18):var w lp lc 3 lt 4 lw 2 title "N_r=00320", \
     "${n00640}" u (\$1/3.086e18):var w lp lc 0 lt 1 lw 2 title "N_r=00640", \
     "${n01280}" u (\$1/3.086e18):var w lp lc 2 lt 1 lw 2 title "N_r=01280", \
     "${n02560}" u (\$1/3.086e18):var w l lc 3 lt 1 lw 2 title "N_r=02560", \
     "${n05120}" u (\$1/3.086e18):var w l lc 5 lt 1 lw 2 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l lc 0 lt 4 lw 2 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l lc 2 lt 4 lw 2 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l lc 3 lt 4 lw 2 title "N_r=40960"

set output "figS1D_hires_IonFraction_t0p025Myr_res.eps"
unset title
set xrange [0.582:0.5882]
set key bottom left
plot "${n05120}" u (\$1/3.086e18):var w lp ls 31 title "N_r=05120", \
     "${n10240}" u (\$1/3.086e18):var w l ls 30 title "N_r=10240", \
     "${n20480}" u (\$1/3.086e18):var w l ls 20 title "N_r=20480", \
     "${n40960}" u (\$1/3.086e18):var w l ls 10 title "N_r=40960"

quit
EOF

gnuplot temp.gp
~/active/bin/eps2jpeg.sh figS1D
exit


