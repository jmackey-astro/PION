#!/bin/bash
# This file plots a bunch of hard-coded infiles with columns of time and
# shock-position for the radiative shock stability test problem.
# see e.g. Stone & Norman (1993), among others.

if [ "$1" = "" ]
then
echo Usage: $0 outfile
outfile=jono

else
outfile=$1
fi
#g++ -Wall -DSERIAL RadShockPosn.cc ../testing/global.cc ../testing/uniformGrid.cc ../testing/dataio.cc -lreadline -lcfitsio
#./a.out rsh_v100n0064.txt ../results/RSh1D_n0064_v100 0 100 > tmp.txt
#./a.out rsh_v100n0128.txt ../results/RSh1D_n0128_v100 0 100 > tmp.txt
#./a.out rsh_v100n0256.txt ../results/RSh1D_n0256_v100 0 200 > tmp.txt
#./a.out rsh_v100n0512.txt ../results/RSh1D_n0512_v100 0 300 > tmp.txt
#./a.out rsh_v100n1024.txt ../results/RSh1D_n1024_v100 0 400 > tmp.txt
#./a.out rsh_v100n2048.txt ../results/RSh1D_n2048_v100 0 800 > tmp.txt
#
#./a.out rsh_v120n0064.txt ../results/RSh1D_n0064_v120 0 100 > tmp.txt
#./a.out rsh_v120n0128.txt ../results/RSh1D_n0128_v120 0 100 > tmp.txt
#./a.out rsh_v120n0256.txt ../results/RSh1D_n0256_v120 0 200 > tmp.txt
#./a.out rsh_v120n0512.txt ../results/RSh1D_n0512_v120 0 300 > tmp.txt
#./a.out rsh_v120n1024.txt ../results/RSh1D_n1024_v120 0 400 > tmp.txt
#
#./a.out rsh_v130n0064.txt ../results/RSh1D_n0064_v130 0 100 > tmp.txt
#./a.out rsh_v130n0128.txt ../results/RSh1D_n0128_v130 0 100 > tmp.txt
#./a.out rsh_v130n0256.txt ../results/RSh1D_n0256_v130 0 200 > tmp.txt
#./a.out rsh_v130n0512.txt ../results/RSh1D_n0512_v130 0 300 > tmp.txt
#./a.out rsh_v130n1024.txt ../results/RSh1D_n1024_v130 0 400 > tmp.txt
#./a.out rsh_v130n4096.txt ../results/RSh1D_n4096_v130 0 800 > tmp.txt
#
#./a.out rsh_v140n0064.txt ../results/RSh1D_n0064_v140 0 100 > tmp.txt
#./a.out rsh_v140n0128.txt ../results/RSh1D_n0128_v140 0 100 > tmp.txt
#./a.out rsh_v140n0256.txt ../results/RSh1D_n0256_v140 0 200 > tmp.txt
#./a.out rsh_v140n0512.txt ../results/RSh1D_n0512_v140 0 300 > tmp.txt
#
#./a.out rsh_v150n0064.txt ../results/RSh1D_n0064_v150 0 100 > tmp.txt
#./a.out rsh_v150n0128.txt ../results/RSh1D_n0128_v150 0 100 > tmp.txt
#./a.out rsh_v150n0256.txt ../results/RSh1D_n0256_v150 0 200 > tmp.txt
#./a.out rsh_v150n0512.txt ../results/RSh1D_n0512_v150 0 300 > tmp.txt
#
#./a.out rsh_v170n0064.txt ../results/RSh1D_n0064_v170 0 100 > tmp.txt
#./a.out rsh_v170n0128.txt ../results/RSh1D_n0128_v170 0 100 > tmp.txt
#./a.out rsh_v170n0256.txt ../results/RSh1D_n0256_v170 0 200 > tmp.txt
#./a.out rsh_v170n0512.txt ../results/RSh1D_n0512_v170 0 300 > tmp.txt
#
echo plotting...
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set size square
#set size 0.5,0.5
#set cntrparam levels 35
#set contour base
#show contour
#set grid front
#set cntrparam levels 100
#set samples 100
#set contour base
set xrange [:]
set yrange [:]
#set title 'Shock Position as Function of Time (128 cells)'
#set title 'Shock Position for 140km/s Shock as Function of Resolution (RK4)'
set ylabel "Shock Position (cm)"
set xlabel "Time (s)"
#unset key 
set title 'Shock Position for 100km/s Shock as Function of Resolution'
set output '$outfile.v100.eps' 
plot 'rsh_v100n0064.txt' u 1:2 w l lw 2, \
     'rsh_v100n0128.txt' u 1:2 w l lw 2, \
     'rsh_v100n0256.txt' u 1:2 w l lw 2, \
     'rsh_v100n0512.txt' u 1:2 w l lw 2, \
     'rsh_v100n1024.txt' u 1:2 w l lw 2, \
     'rsh_v100n2048.txt' u 1:2 w l lw 2
pause -1
#
set title 'Shock Position for 120km/s Shock as Function of Resolution'
set output '$outfile.v120.eps' 
plot 'rsh_v120n0064.txt' u 1:2 w l lw 2, \
     'rsh_v120n0128.txt' u 1:2 w l lw 2, \
     'rsh_v120n0256.txt' u 1:2 w l lw 2, \
     'rsh_v120n0512.txt' u 1:2 w l lw 2, \
     'rsh_v120n1024.txt' u 1:2 w l lw 2
pause -1
#
set title 'Shock Position for 130km/s Shock as Function of Resolution'
set output '$outfile.v130.eps'
set xrange [0:1.2e11]
plot 'rsh_v130n0064.txt' u 1:2 w l lw 1, \
     'rsh_v130n0128.txt' u 1:2 w l lw 1, \
     'rsh_v130n0256.txt' u 1:2 w l lw 1, \
     'rsh_v130n0512.txt' u 1:2 w l lw 1, \
     'rsh_v130n1024.txt' u 1:2 w l lw 1, \
     'rsh_v130n4096.txt' u 1:2 w l lw 1 lt 0

pause -1
#
set title 'Shock Position for 140km/s Shock as Function of Resolution'
set output '$outfile.v140.eps' 
set xrange [*:*]
plot 'rsh_v140n0064.txt' u 1:2 w l lw 2, \
     'rsh_v140n0128.txt' u 1:2 w l lw 2, \
     'rsh_v140n0256.txt' u 1:2 w l lw 2, \
     'rsh_v140n0512.txt' u 1:2 w l lw 2
pause -1
#
set title 'Shock Position for 150km/s Shock as Function of Resolution'
set output '$outfile.v150.eps' 
plot 'rsh_v150n0064.txt' u 1:2 w l lw 2, \
     'rsh_v150n0128.txt' u 1:2 w l lw 2, \
     'rsh_v150n0256.txt' u 1:2 w l lw 2, \
     'rsh_v150n0512.txt' u 1:2 w l lw 2
pause -1
#
set title 'Shock Position for 170km/s Shock as Function of Resolution'
set output '$outfile.v170.eps' 
plot 'rsh_v170n0064.txt' u 1:2 w l lw 2, \
     'rsh_v170n0128.txt' u 1:2 w l lw 2, \
     'rsh_v170n0256.txt' u 1:2 w l lw 2, \
     'rsh_v170n0512.txt' u 1:2 w l lw 2
pause -1
set title 'Shock Position as Function of Inflow Velocity (n=512)'
set output '$outfile.allspeeds.eps' 
set xrange [0:3.0e11]
plot 'rsh_v100n0512.txt' u 1:2 w l lw 1 title "100km/s", \
     'rsh_v120n0512.txt' u 1:2 w l lw 1 title "120km/s", \
     'rsh_v130n0512.txt' u 1:2 w l lw 1 title "130km/s", \
     'rsh_v140n0512.txt' u 1:2 w l lw 1 title "140km/s", \
     'rsh_v150n0512.txt' u 1:2 w l lw 1 title "150km/s", \
     'rsh_v170n0512.txt' u 1:2 w l lw 1 lt 0 title "170km/s"
pause -1
EOF
gnuplot gnu.plt
#
#plot 'v120.txt' u 1:2 w l, \
#     'v130.txt' u 1:2 w l, \
#     'v140.txt' u 1:2 w l, \
#     'v150.txt' u 1:2 w l
#
#plot 'n0064.txt'      u 1:2 w l lt 1, \
#     'n0064_noBC.txt' u 1:2 w p lt 1, \
#     'n0128.txt'      u 1:2 w l lt 2, \
#     'n0128_noBC.txt' u 1:2 w p lt 2, \
#     'n0256.txt'      u 1:2 w l lt 3, \
#     'n0256_noBC.txt' u 1:2 w p lt 3, \
#     'n0512.txt'      u 1:2 w l lt 4, \
#     'n0512_noBC.txt' u 1:2 w p lt 4, \
#     'n1024.txt'      u 1:2 w l lt 5, \
#     'n1024_noBC.txt' u 1:2 w p lt 5
#     'n0064_Adaptive.txt' u 1:2 w p 1, \
#     'n0128_Adaptive.txt' u 1:2 w p 2, \
#     'n0256_Adaptive.txt' u 1:2 w p 3, \
#     'n0512_Adaptive.txt' u 1:2 w p lt -1 pt 0
#
#plot 'v140n0064_noBC.txt' u 1:2 w l lt 1, \
#     'v140n0128_noBC.txt' u 1:2 w l lt 2, \
#     'v140n0256_noBC.txt' u 1:2 w l lt 3, \
#     'v140n0512_noBC.txt' u 1:2 w p lt 4, \
#     'v140n1024_noBC.txt' u 1:2 w p lt 5
#
#plot 'v140n0128_wBC.txt' u 1:2 w l lt 2, \
#     'v140n0256_wBC.txt' u 1:2 w l lt 3, \
#     'v140n0512_wBC.txt' u 1:2 w l lt 4, \
#     'v140n1024_wBC.txt' u 1:2 w l lt 5, \
#     'v140n2048_wBC.txt' u 1:2 w l lt 6, \
#     'v140n4096_wBC.txt' u 1:2 w l lt 7
#plot 'n0128v150_OF.txt' u 1:2 w l lt 2, \
#     'n0256v150_OF.txt' u 1:2 w l lt 3, \
#     'n0512v150_OF.txt' u 1:2 w l lt 4, \
#     'n1024v150_OF.txt' u 1:2 w l lt 5 #, \
#     'v150n2048_wBC.txt' u 1:2 w l lt 6, \
#     'v150n4096_wBC.txt' u 1:2 w l lt 7
#pause -1
#EOF
#gnuplot gnu.plt
