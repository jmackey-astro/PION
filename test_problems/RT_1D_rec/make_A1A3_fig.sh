#!/bin/bash

BASE_A1=$1
BASE_A2=$2
BOXSIZE=$3
TREC=$4
NX=$5
LABEL=$6

BASE02=${BASE}_dt02
BASE05=${BASE}_dt05
FNAME=${BASE}

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})" 0.0,0.5
set key right bottom
unset log y
set ylabel "Ratio R_{a}/R_{if}(t)" 2.0,0.0
set yrange [0.95:1.02]
set log x
set xrange [0.001:11]
set title ""
set label "$LABEL" at 0.0015,1.015
g(x)=(1.0-exp(-x))**(1.0/3.0)
RH=${BOXSIZE}
RI=0.75287*${BOXSIZE}
#RI=0.75588*${BOXSIZE}
tI=${TREC}
DX=${BOXSIZE}/${NX}
f(x)=1.0
h(x)=1.0+DX/(RI*g(x))
j(x)=1.0-DX/(RI*g(x))
set output "${FNAME}_A1A3.eps"
plot f(x) w l lw 2 lt -1 title 'Analytic Solution', \
     h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     '${BASE_A3}.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 3  title 'A3-dt02', \
     '${BASE_A1}.txt' u (\$1/tI):((\$2/RI)/g(\$1/tI)) w l lw 2 lt 1  title 'A1-dt05'

quit
EOF

gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}_A1A3.eps ${FNAME}_A1A3.jpeg
cp ${FNAME}_*.eps ${FNAME}_*.jpeg ./


