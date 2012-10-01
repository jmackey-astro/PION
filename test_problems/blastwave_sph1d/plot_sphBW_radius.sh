#!/bin/bash

#
# Plot the radius of a blastwave as a function of time, from a number of different files.
#
outfile=$1
F1=$2
F2=$3
F3=$4
F1title="$5"
F2title="$6"
F3title="$7"

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
set output '${outfile}.eps' 
set size 0.7071,0.7071
#set term x11 noraise size 1500,500
set log xy
set xrange [1.0:7.0e4]
set yrange [0.2:12.0]
set ylabel "Shock Radius (pc)" offset 2,0
set xlabel "Time (yr)" offset 0,0.5
set key bottom right
#
# f(x) is appropriate for the Sedov-Taylor solution.
#
f(x)=0.125*x**(0.4)
#
plot '$F1' u (\$1/3.16e7):(\$2/3.086e18) w p title "$F1title", \
     '$F2' u (\$1/3.16e7):(\$2/3.086e18) w p title "$F2title", \
     '$F3' u (\$1/3.16e7):(\$2/3.086e18) w p title "$F3title", \
     f(x) title "k*t^{2/5}"
#pause -1
quit
EOF

gnuplot gnu.plt


convert -density 300 -quality 100 ${outfile}.eps ${outfile}.jpeg
