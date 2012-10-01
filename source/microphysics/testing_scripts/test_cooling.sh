#!/bin/bash


g++ ../mp_only_cooling.cc ../cooling_SD93_cie.cc ../../global.cc ../../cell_interface.cc -o a.out -lreadline
./a.out > cooling.txt

cat << EOF > gnu.plt
set xrange [4.0:8.5]
set yrange [9:16]
set grid
set title "Cooling time as fn of temperature (log-log)"
plot './cooling.txt' u (log10(\$2)):(log10(\$3)) w lp
pause -1
set xrange [*:*]
set yrange [*:*]
set title "Temperature as a function of time (lin-log)"
plot './cooling.txt' u 1:(log10(\$2)) w lp
pause -1
set log xy
set yrange [1.e-30:1.e-18]
set title "Temperature as a function of time (log-log)"
plot 'coolingcurve.txt' u 1:2 w lp lt -1
pause -1
quit
EOF

gnuplot gnu.plt

