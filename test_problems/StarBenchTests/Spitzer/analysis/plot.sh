#!/bin/bash

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times-New-Roman,22" \
 fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set lmargin 6.0
set tmargin 0.5
set rmargin 1.25
set bmargin 3.5
set key bottom right 
set xrange [0:0.12]
set yrange [*:*]
unset log y
unset log x
set xlabel "time (Myr)" offset 0,0
set ylabel "R_{IF} (pc)" offset 2,0

set output "S1D_Ifront_resolution.eps"
a=2
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "N_x=0080", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "N_x=0160", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "N_x=0320", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "N_x=0640", \
     "../analysis/S1Dn1280.txt" u 1:a w l lt 1 lc 5 lw 2 title "N_x=1280", \
     "../analysis/S1Dn2560.txt" u 1:a w l lt 1 lc 0 lw 2 title "N_x=2560", \
     "../analysis/S1Dn5120.txt" u 1:a w l lt 1 lc 3 lw 2 title "N_x=5120", \
     "early_phase.txt" u 2:3 w l title "Raga", "" u 2:4 w l title "Spitzer", "" u 2:6 w l lw 4 title "HI"


set output "S1D_Shock_resolution.eps"
set ylabel "R_{SF} (pc)" offset 2,0
a=3
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "N_x=0080", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "N_x=0160", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "N_x=0320", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "N_x=0640", \
     "../analysis/S1Dn1280.txt" u 1:a w l lt 1 lc 5 lw 2 title "N_x=1280", \
     "../analysis/S1Dn2560.txt" u 1:a w l lt 1 lc 0 lw 2 title "N_x=2560", \
     "../analysis/S1Dn5120.txt" u 1:a w l lt 1 lc 3 lw 2 title "N_x=5120", \
     "early_phase.txt" u 2:5 w l title "Raga-SF", "" u 2:6 w l title "HI"
#plot "../analysis/S1Dn0080.txt" u 1:a w lp title "dx=0.2000 pc", \
#     "../analysis/S1Dn0160.txt" u 1:a w lp title "dx=0.1000 pc", \
#     "../analysis/S1Dn0320.txt" u 1:a w lp title "dx=0.0500 pc", \
#     "../analysis/S1Dn0640.txt" u 1:a w lp title "dx=0.0250 pc", \
#     "../analysis/S1Dn1280.txt" u 1:a w l lt 1 lc 5 lw 2 title "dx=0.0125 pc", \
#     "../analysis/S1Dn2560.txt" u 1:a w l lt -1 title "dx=0.01125 pc"

set output "S1D_iMass_resolution.eps"
set ylabel "Ionised mass (Msun)" offset 2,0
a=4
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "N_x=0080", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "N_x=0160", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "N_x=0320", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "N_x=0640", \
     "../analysis/S1Dn1280.txt" u 1:a w l lt 1 lc 5 lw 2 title "N_x=1280", \
     "../analysis/S1Dn2560.txt" u 1:a w l lt 1 lc 0 lw 2 title "N_x=2560", \
     "../analysis/S1Dn5120.txt" u 1:a w l lt 1 lc 3 lw 2 title "N_x=5120"

set output "S1D_ShellMass_resolution.eps"
set ylabel "Shell mass (Msun)" offset 2,0
a=5
plot "../analysis/S1Dn0080.txt" u 1:a w lp title "N_x=0080", \
     "../analysis/S1Dn0160.txt" u 1:a w lp title "N_x=0160", \
     "../analysis/S1Dn0320.txt" u 1:a w lp title "N_x=0320", \
     "../analysis/S1Dn0640.txt" u 1:a w lp title "N_x=0640", \
     "../analysis/S1Dn1280.txt" u 1:a w l lt 1 lc 5 lw 2 title "N_x=1280", \
     "../analysis/S1Dn2560.txt" u 1:a w l lt 1 lc 0 lw 2 title "N_x=2560", \
     "../analysis/S1Dn5120.txt" u 1:a w l lt 1 lc 3 lw 2 title "N_x=5120"
exit
EOF
gnuplot temp.gp

#~/Documents/active/bin/eps2jpeg.sh S1D
~/active/bin/eps2jpeg.sh S1D

