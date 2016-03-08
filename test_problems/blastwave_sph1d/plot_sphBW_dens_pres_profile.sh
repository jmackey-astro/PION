#!/bin/bash

#
# Plot the density and pressure for blast waves at the final time (50kyr)
# - 2016.03.08 JM: updated for new pion.
#
outfile=$1
F1=$2
F2=$3
F3=$4
F1title="$5"
F2title="$6"
F3title="$7"

cat << EOF  > gnu.plt
set terminal postscript enhanced color eps size 4,2 "Helvetica" 14
set output '${outfile}.eps' 
set size 1.0,1.0
set origin 0.0,0.0
set key bottom right

#set size 0.7071,0.7071
#set term x11 noraise size 1500,500

set multiplot
set ytics -100,1,0
#set lmargin 15
#set rmargin 1
#set size 1.0, 0.35
#set origin 0.0,0.0
set size 1.0, 0.57
set origin 0.0,0.0
set xlabel "R (pc)" offset 0,1.0
set ylabel "log_{10} rho/(g/cm^3)" offset 1.0,0
set yrange [-26:-21]
set xrange [0.0:10.2]
unset log y
plot '$F1' u (\$1/3.086e18):(log10(\$2)) w l lw 2 title "$F1title", \
     '$F2' u (\$1/3.086e18):(log10(\$2)) w l lw 2 title "$F2title", \
     '$F3' u (\$1/3.086e18):(log10(\$2)) w l lw 2 title "$F3title"
#
# SECOND PLOT
#
#set origin 0.0,0.33
#set size 1.0,0.35
set key off
set origin 0.0,0.5
set size 1.0,0.5
set ylabel "log_{10} P/(dyne/cm^2)" offset 1.0,0
set yrange [-11:-7]
set format x ""
set xlabel ""
set label 1 "t = 50 kyr" at screen 0.2, screen 0.6
set ytics -100,1,0
plot '$F1' u (\$1/3.086e18):(log10(\$3)) w l lw 2 title "$F1title, P_g", \
     '$F2' u (\$1/3.086e18):(log10(\$3)) w l lw 2 title "$F2title, P_g", \
     '$F3' u (\$1/3.086e18):(log10(\$3)) w l lw 2 title "$F3title, P_g"

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

convert -density 400 ${outfile}.eps -background white -flatten ${outfile}.png
