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

#
# Hosokawa & Inutsuka (2006) shell radius.
#
HI(x)=0.313838*(1.0+84.7879*x)**(4.0/7.0)

#
# Thomas's ionised mass estimate (in Msun)
#
MI(x) = 56.694*(HI(x))**1.5

#
# Neutral mass in shell estimate (in Msun)
#
MN(x) = 322.463*(HI(x))**3.0 -MI(x)

set style line 10 lt 1 lc rgb "gray20"    lw 2 pt 0 ps 1
set style line 11 lt 1 lc rgb "gray20"    lw 2 pt 1 ps 1
set style line 14 lt 1 lc rgb "gray20"    lw 2 pt 4 ps 1
set style line 20 lt 1 lc rgb "blue"      lw 2 pt 0 ps 1
set style line 21 lt 1 lc rgb "blue"      lw 2 pt 1 ps 1
set style line 24 lt 1 lc rgb "blue"      lw 2 pt 4 ps 1
set style line 30 lt 1 lc rgb "web-green" lw 2 pt 0 ps 1
set style line 31 lt 1 lc rgb "web-green" lw 2 pt 1 ps 1
set style line 34 lt 1 lc rgb "web-green" lw 2 pt 4 ps 1
set style line 40 lt 1 lc rgb "cyan" lw 2 pt 0 ps 1
set style line 41 lt 1 lc rgb "cyan" lw 2 pt 1 ps 1
set style line 44 lt 1 lc rgb "cyan" lw 2 pt 4 ps 1

set yrange [0.3:1.3]
set output "S1D_Ifront_resolution_low.eps"
a=2
plot "../analysis/S1Dn00080.txt" u 1:a w l ls 10 title "N_x=00080", \
     "../analysis/S1Dn00160.txt" u 1:a w l ls 20 title "N_x=00160", \
     "../analysis/S1Dn00320.txt" u 1:a w l ls 30 title "N_x=00320", \
     "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "early_phase.txt" u 2:3 w l lc 4 title "Raga", "" u 2:4 w l lc 7 title "Spitzer", HI(x) w l lw 4 lc 0 title "HI"

set output "S1D_Ifront_resolution_med.eps"
plot "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "../analysis/S1Dn01280.txt" u 1:a w l ls 10 title "N_x=01280", \
     "../analysis/S1Dn02560.txt" u 1:a w l ls 20 title "N_x=02560", \
     "../analysis/S1Dn05120.txt" u 1:a w l ls 30 title "N_x=05120", \
     "early_phase.txt" u 2:3 w l lc 4 title "Raga", "" u 2:4 w l lc 7 title "Spitzer", HI(x) w l lw 4 lc 0 title "HI"

set output "S1D_Ifront_resolution_high.eps"
plot "../analysis/S1Dn05120.txt" u 1:a w lp ls 31 title "N_x=05120", \
     "../analysis/S1Dn10240.txt" u 1:a w l ls 40 title "N_x=10240", \
     "../analysis/S1Dn20480.txt" u 1:a w l ls 10 title "N_x=20480", \
     "../analysis/S1Dn40960.txt" u 1:a w l ls 20 title "N_x=40960", \
     "early_phase.txt" u 2:3 w l lc 4 title "Raga", "" u 2:4 w l lc 7 title "Spitzer", HI(x) w l lw 4 lc 0 title "HI"


set ylabel "R_{SF} (pc)" offset 2,0
a=3
set output "S1D_Shock_resolution_low.eps"
plot "../analysis/S1Dn00080.txt" u 1:a w l ls 10 title "N_x=00080", \
     "../analysis/S1Dn00160.txt" u 1:a w l ls 20 title "N_x=00160", \
     "../analysis/S1Dn00320.txt" u 1:a w l ls 30 title "N_x=00320", \
     "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "early_phase.txt" u 2:5 w l lc 4 title "Raga-SF", HI(x) w l lw 4 lc 0 title "HI"

set output "S1D_Shock_resolution_med.eps"
plot "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "../analysis/S1Dn01280.txt" u 1:a w l ls 10 title "N_x=01280", \
     "../analysis/S1Dn02560.txt" u 1:a w l ls 20 title "N_x=02560", \
     "../analysis/S1Dn05120.txt" u 1:a w l ls 30 title "N_x=05120", \
     "early_phase.txt" u 2:5 w l lc 4 title "Raga-SF", HI(x) w l lw 4 lc 0 title "HI"

set output "S1D_Shock_resolution_high.eps"
plot "../analysis/S1Dn05120.txt" u 1:a w lp ls 31 title "N_x=05120", \
     "../analysis/S1Dn10240.txt" u 1:a w l ls 40 title "N_x=10240", \
     "../analysis/S1Dn20480.txt" u 1:a w l ls 10 title "N_x=20480", \
     "../analysis/S1Dn40960.txt" u 1:a w l ls 20 title "N_x=40960", \
     "early_phase.txt" u 2:5 w l lc 4 title "Raga-SF", HI(x) w l lw 4 lc 0 title "HI"
     

set yrange [*:*]
set ylabel "Ionised mass (Msun)" offset 2,0
a=4
set output "S1D_iMass_resolution_low.eps"
plot "../analysis/S1Dn00080.txt" u 1:a w l ls 10 title "N_x=00080", \
     "../analysis/S1Dn00160.txt" u 1:a w l ls 20 title "N_x=00160", \
     "../analysis/S1Dn00320.txt" u 1:a w l ls 30 title "N_x=00320", \
     "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     MI(x) w l lw 4 lc 0 title "HI"

set output "S1D_iMass_resolution_med.eps"
plot "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "../analysis/S1Dn01280.txt" u 1:a w l ls 10 title "N_x=01280", \
     "../analysis/S1Dn02560.txt" u 1:a w l ls 20 title "N_x=02560", \
     "../analysis/S1Dn05120.txt" u 1:a w l ls 30 title "N_x=05120", \
     MI(x) w l lw 4 lc 0 title "HI"


set output "S1D_iMass_resolution_high.eps"
plot "../analysis/S1Dn05120.txt" u 1:a w lp ls 31 title "N_x=05120", \
     "../analysis/S1Dn10240.txt" u 1:a w l ls 40 title "N_x=10240", \
     "../analysis/S1Dn20480.txt" u 1:a w l ls 10 title "N_x=20480", \
     "../analysis/S1Dn40960.txt" u 1:a w l ls 20 title "N_x=40960", \
     MI(x) w l lw 4 lc 0 title "HI"

     

set ylabel "Shell mass (Msun)" offset 2,0
a=5
set output "S1D_ShellMass_resolution_low.eps"
plot "../analysis/S1Dn00080.txt" u 1:a w l ls 10 title "N_x=00080", \
     "../analysis/S1Dn00160.txt" u 1:a w l ls 20 title "N_x=00160", \
     "../analysis/S1Dn00320.txt" u 1:a w l ls 30 title "N_x=00320", \
     "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     MN(x) w l lw 4 lc 0 title "HI"


set output "S1D_ShellMass_resolution_med.eps"
plot "../analysis/S1Dn00640.txt" u 1:a w l ls 40 title "N_x=00640", \
     "../analysis/S1Dn01280.txt" u 1:a w l ls 10 title "N_x=01280", \
     "../analysis/S1Dn02560.txt" u 1:a w l ls 20 title "N_x=02560", \
     "../analysis/S1Dn05120.txt" u 1:a w l ls 30 title "N_x=05120", \
     MN(x) w l lw 4 lc 0 title "HI"


set output "S1D_ShellMass_resolution_high.eps"
plot "../analysis/S1Dn05120.txt" u 1:a w lp ls 31 title "N_x=05120", \
     "../analysis/S1Dn10240.txt" u 1:a w l ls 40 title "N_x=10240", \
     "../analysis/S1Dn20480.txt" u 1:a w l ls 10 title "N_x=20480", \
     "../analysis/S1Dn40960.txt" u 1:a w l ls 20 title "N_x=40960", \
     MN(x) w l lw 4 lc 0 title "HI"


exit
EOF
gnuplot temp.gp

#~/Documents/active/bin/eps2jpeg.sh S1D
~/active/bin/eps2jpeg.sh S1D

