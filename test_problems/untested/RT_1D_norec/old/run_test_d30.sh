#!/bin/bash

rm rtt1D_NoRec_*

./icgen_serial pf_rtt_NoRec_cart1D_n128_v0100_d30.txt silo


OOA=O2
ACC=DX10000
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX05000
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX02500
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX01250
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX00625
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo

OOA=O1
ACC=DX10000
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX05000
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX02500
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX01250
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo
ACC=DX00625
./main_serial_${OOA}_${ACC} IC_rtt_NoRec_cart1D_n128_v0100_d30.silo 5 1 optype=6 outfile=rtt1D_NoRec_v0100_d30_${OOA}_${ACC} redirect=./temp_
./plot_radius ./results_v0100_d30_${OOA}_${ACC} ./rtt1D_NoRec_v0100_d30_${OOA}_${ACC} 0 5 5 silo


### O1 PLOTS ###
SPEED=v0100_d30
FNAME=AccuracyO1_${SPEED}
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (kyr)"
set ylabel "Fractional Error" 2.0,0.0
set log y
set yrange [0.01:1.0]
set xrange [0.0:9.0]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=30 cm^{-3}, O1" at 0.5,0.75
TS=3.16e10
plot './results_${SPEED}_O1_DX10000.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=1.0000", \
     './results_${SPEED}_O1_DX05000.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.5000", \
     './results_${SPEED}_O1_DX02500.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.2500", \
     './results_${SPEED}_O1_DX01250.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.1250", \
     './results_${SPEED}_O1_DX00625.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.0625"
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg


### O2 ###
SPEED=v0100_d30
FNAME=AccuracyO2_${SPEED}
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set size 0.7071,0.7071
set xlabel "Time (kyr)"
set ylabel "Fractional Error" 3.0,0.0
set log y
set yrange [0.0001:0.8]
set xrange [0.0:9.0]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=30 cm^{-3}, O2" at 0.5,0.4
TS=3.16e10
plot './results_${SPEED}_O2_DX10000.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=1.0000", \
     './results_${SPEED}_O2_DX05000.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.5000", \
     './results_${SPEED}_O2_DX02500.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.2500", \
     './results_${SPEED}_O2_DX01250.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.1250", \
     './results_${SPEED}_O2_DX00625.txt' u (\$1/TS):(abs(1.0-\$6/\$7)) w l lw 2 title "K=0.0625"
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg


### ALL DONE ###

rm rtt1D_NoRec_*

