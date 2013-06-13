#!/bin/bash

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times-New-Roman,22" 
#\
# fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set lmargin 5.0
set tmargin 0.5
set rmargin 1.25
set bmargin 3.5
set key bottom right 
set xrange [0:4.3]
set yrange [*:*]
unset log y
unset log x
set xlabel "time (Myr)" offset 0,0
set ylabel "R_{IF} (pc)" offset 2,0

set output "S1D_Ifront_resolution.eps"
a=2
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "dx=0.2000 pc", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "dx=0.1000 pc", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "dx=0.0500 pc", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "dx=0.0250 pc", \
     "../analysis/S1Dn1280.txt" u 1:a w lp title "dx=0.0125 pc"

set output "S1D_Shock_resolution.eps"
set ylabel "R_{SF} (pc)" offset 2,0
a=3
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "dx=0.2000 pc", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "dx=0.1000 pc", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "dx=0.0500 pc", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "dx=0.0250 pc", \
     "../analysis/S1Dn1280.txt" u 1:a w lp title "dx=0.0125 pc"

exit
EOF
gnuplot temp.gp

~/Documents/active/bin/eps2jpeg.sh S1D

