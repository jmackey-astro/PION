#!/bin/bash

#
# This file should generate a gnu.plt file which gnuplot can use to
# read data and generate a figure.  I have hardwired this script for
# the Alfven Wave test problem, so i plot rho, vx, vy, vz
#
infile1=$1
infile2=$2
infile3=$3
outfile=$4

cat << EOF  > gnu.plt
#set terminal postscript enhanced eps
set terminal postscript enhanced color eps
set output "${outfile}.eps"
#set size 3.0, 1.0
set size 1.5, 0.5
set origin 0.0,0.0
set multiplot
#
# first plot
#
#set size 1.0, 1.0
set size 0.5, 0.5
set size square
set origin -0.1,0.0
#set yrange [0:1.05]
set xlabel "Position"
set ylabel "Density" offset 1.3,0.0
set xtics 0,.2
set key off
#set key noautotitles
#set key 0.45,0.2
set yrange [0.9:1.1]
plot '${infile3}'  using 1:2 w lp lt 5 pt 5 ps 0.5, \
     '${infile1}'  using 1:2 w lp lt 1 pt 1 title 'density', \
     '${infile2}'  using 1:2 w l lt -1
#
# second plot
#
set origin 0.3,0.0
set yrange [-0.1:0.1]
set ylabel "v_x" offset 1.3,0.0
plot '${infile3}'  using 1:4 w lp lt 5 pt 5 ps 0.5, \
     '${infile1}'  using 1:4 w lp lt 2 pt 2 , \
     '${infile2}'  using 1:4 w l lt -1

# THIRD PLOT
set origin 0.675,0.0
set yrange [-1.1:1.1]
set ylabel "v_y" offset 1.3,0.0
plot '${infile3}'  using 1:5 w lp lt 5 pt 5 ps 0.5, \
     '${infile1}'  using 1:5 with lp  lt 3 pt 3, \
     '${infile2}'  using 1:5 w l lt -1

# Fourth PLOT
set origin 1.05,0.0
set yrange [-0.1:2.1]
set ylabel "v_z" offset 1.3,0.0
plot '${infile3}'  using 1:6 w lp lt 5 pt 5 ps 0.5, \
     '${infile1}'  using 1:6 with lp  lt 4 pt 4, \
     '${infile2}'  using 1:6 w l lt -1


unset multiplot
#pause -1
EOF

gnuplot gnu.plt
convert -density 300 -quality 100 ${outfile}.eps ${outfile}.jpeg
convert -adaptive-resize 900x210 ${outfile}.jpeg ${outfile}.jpeg

