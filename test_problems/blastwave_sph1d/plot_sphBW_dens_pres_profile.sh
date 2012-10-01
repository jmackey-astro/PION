#!/bin/bash

#
# Plot the density and pressure for blast waves at the final time (50kyr)
#
outfile=$1
F1=$2
F2=$3
F3=$4
F1title="$5"
F2title="$6"
F3title="$7"

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps size 4,2 "Times-Roman" 14
set output '${outfile}.eps' 
set size 1.0,1.0
set origin 0.0,0.0
set key bottom right

#set size 0.7071,0.7071
#set term x11 noraise size 1500,500

set multiplot
#set lmargin 15
#set rmargin 1
#set size 1.0, 0.35
#set origin 0.0,0.0
set size 1.0, 0.57
set origin 0.0,0.0
set xlabel "R (pc)" offset 0,1.0
set ylabel "Density" offset 2.0,0
set yrange [1.0e-25:2.0e-21]
set xrange [0.0:10.2]
set log y
plot '$F1' u (\$1/3.086e18):2 w l lw 2 title "$F1title", \
     '$F2' u (\$1/3.086e18):2 w l lw 2 title "$F2title", \
     '$F3' u (\$1/3.086e18):2 w l lw 2 title "$F3title"
#
# SECOND PLOT
#
#set origin 0.0,0.33
#set size 1.0,0.35
set key off
set origin 0.0,0.5
set size 1.0,0.5
set ylabel "Pressure" offset 2.0,0
set yrange [1.0e-11:1.0e-7]
set format x ""
set xlabel ""
set label 1 "t = 50 kyr" at screen 0.2, screen 0.6
plot '$F1' u (\$1/3.086e18):3 w l lw 2 title "$F1title, P_g", \
     '$F2' u (\$1/3.086e18):3 w l lw 2 title "$F2title, P_g", \
     '$F3' u (\$1/3.086e18):3 w l lw 2 title "$F3title, P_g"

#
# THIRD PLOT
#
#set origin 0.0,0.66
#set size 1.0,0.35
#set key on
#set yrange [1.0e4:1.0e10]
#set ylabel "E_{int}" offset 0.7,0
#plot '$F1' u (\$1/3.086e18):4 w l title "$F1title, Eint", \
#     '$F2' u (\$1/3.086e18):4 w l title "$F2title, Eint", \
#     '$F3' u (\$1/3.086e18):4 w l title "$F3title, Eint"
#
unset multiplot
#
#
#pause -1
EOF



gnuplot gnu.plt

convert -density 300 -quality 100 ${outfile}.eps ${outfile}.jpeg
