#!/bin/bash

rm rtt1D_NoRec_* temp/rtt1D_NoRec_*
mkdir temp;

#  100 km/s model
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0100_n1e1_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0100_n1e2_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0100_n1e3_alg1.txt silo
OOA=C2RAY

#n1e1
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} redirect=temp/v0100_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} redirect=temp/v0100_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} redirect=temp/v0100_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e1_${OOA}_${ACC} 0 1 5 silo

#n1e2
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} redirect=temp/v0100_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} redirect=temp/v0100_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} redirect=temp/v0100_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0100_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e2_${OOA}_${ACC} 0 1 5 silo

#n1e3
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} redirect=temp/v0100_n1e3_${OOA}_${ACC}_ opfreq=1
./plot_radius temp/results_v0100_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} redirect=temp/v0100_n1e3_${OOA}_${ACC}_ opfreq=2
./plot_radius temp/results_v0100_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} 0 2 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} redirect=temp/v0100_n1e3_${OOA}_${ACC}_ opfreq=9
./plot_radius temp/results_v0100_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0100_n1e3_${OOA}_${ACC} 0 9 5 silo

#  300 km/s model
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0300_n1e1_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0300_n1e2_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v0300_n1e3_alg1.txt silo

#n1e1
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} redirect=./temp/v0300_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} redirect=./temp/v0300_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} redirect=./temp/v0300_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e1_${OOA}_${ACC} 0 1 5 silo

#n1e2
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} redirect=./temp/v0300_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} redirect=./temp/v0300_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} redirect=./temp/v0300_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v0300_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e2_${OOA}_${ACC} 0 1 5 silo

#n1e3
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} redirect=./temp/v0300_n1e3_${OOA}_${ACC}_ opfreq=1
./plot_radius temp/results_v0300_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} redirect=./temp/v0300_n1e3_${OOA}_${ACC}_ opfreq=2
./plot_radius temp/results_v0300_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} 0 2 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0300_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} redirect=./temp/v0300_n1e3_${OOA}_${ACC}_ opfreq=9
./plot_radius temp/results_v0300_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v0300_n1e3_${OOA}_${ACC} 0 9 5 silo

#  1000 km/s model
./icgen_serial pf_rtt_NoRec_cart1D_n128_v1000_n1e1_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v1000_n1e2_alg1.txt silo
./icgen_serial pf_rtt_NoRec_cart1D_n128_v1000_n1e3_alg1.txt silo

#n1e1
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} redirect=./temp/v1000_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} redirect=./temp/v1000_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e1_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} redirect=./temp/v1000_n1e1_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e1_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e1_${OOA}_${ACC} 0 1 5 silo

#n1e2
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} redirect=./temp/v1000_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} redirect=./temp/v1000_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} 0 1 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e2_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} redirect=./temp/v1000_n1e2_${OOA}_${ACC}_
./plot_radius temp/results_v1000_n1e2_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e2_${OOA}_${ACC} 0 1 5 silo

# n1e3
ACC=030tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} redirect=./temp/v1000_n1e3_${OOA}_${ACC}_ opfreq=1
./plot_radius temp/results_v1000_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} 0 1 5 silo
ACC=010tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} redirect=./temp/v1000_n1e3_${OOA}_${ACC}_ opfreq=2
./plot_radius temp/results_v1000_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} 0 2 5 silo
ACC=003tr
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v1000_n1e3_alg1.silo 5 1 optype=6 outfile=temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} redirect=./temp/v1000_n1e3_${OOA}_${ACC}_ opfreq=9
./plot_radius temp/results_v1000_n1e3_${OOA}_${ACC} temp/rtt1D_NoRec_v1000_n1e3_${OOA}_${ACC} 0 9 5 silo



### O1 PLOTS ###
BASE1=v0100_n1e1
BASE2=v0100_n1e2
BASE3=v0100_n1e3
FNAME=Accuracy_C2RAY_v0100
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (kyr)"
set ylabel "Fractional Error" 2.0,0.0
set log y
set yrange [1e-8:1.0e-3]
set xrange [0.0:20.0]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=100 cm^{-3}, Alg1" at 0.5,0.75
TS=3.16e10
plot 'temp/results_${BASE1}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.30 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.10 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.03 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.30 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.10 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.03 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.30 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.10 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.03 t_{rec}"
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg

### O1 PLOTS ###
BASE1=v0300_n1e1
BASE2=v0300_n1e2
BASE3=v0300_n1e3
FNAME=Accuracy_C2RAY_v0300
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (kyr)"
set ylabel "Fractional Error" 2.0,0.0
set log y
set yrange [1e-8:1.0e-3]
#set yrange [0.0003:0.03]
set xrange [0.0:6.0]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=100 cm^{-3}, Alg1" at 0.5,0.75
TS=3.16e10
plot 'temp/results_${BASE1}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.30 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.10 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.03 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.30 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.10 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.03 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.30 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.10 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.03 t_{rec}"
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg

### O1 PLOTS ###
BASE1=v1000_n1e1
BASE2=v1000_n1e2
BASE3=v1000_n1e3
FNAME=Accuracy_C2RAY_v1000
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (kyr)"
set ylabel "Fractional Error" 2.0,0.0
set log y
set yrange [1e-8:1.0e-3]
#set yrange [0.0005:0.05]
set xrange [0.0:2.0]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=100 cm^{-3}, Alg1" at 0.5,0.75
TS=3.16e10
plot 'temp/results_${BASE1}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.30 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.10 t_{rec}", \
     'temp/results_${BASE1}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^1, t=0.03 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.30 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.10 t_{rec}", \
     'temp/results_${BASE2}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^2, t=0.03 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_030tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.30 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_010tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.10 t_{rec}", \
     'temp/results_${BASE3}_C2RAY_003tr.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w lp lw 2 title "n=10^3, t=0.03 t_{rec}"
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg

##### ALL DONE ###

rm temp/rtt1D_NoRec_*C2RAY*
rmdir temp;

