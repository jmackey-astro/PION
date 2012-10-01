#!/bin/bash
#
# 2009-12-25 JM: changed point size.
#
#
# This file should generate a gnu.plt file which gnuplot can use to
# read data and generate a figure.  I have hardwired this script for
# the Alfven Wave test problem, so i plot B_t as a function of resolution
#
# infile1 = n016
# infile2 = n032
# infile3 = n064
# infile4 = n128
#
infile1=$1
infile2=$2
infile3=$3
infile4=$4
outfile=$5
ix=$6

cat << EOF  > gnu.plt
#set terminal postscript enhanced eps
set terminal postscript enhanced color eps
set output "${outfile}.eps"
set size 0.7071, 0.7071
set yrange [-0.1:0.1]
set xlabel "Position"
set ylabel "B_t" offset 1.3,0.0
#set xtics 0,.2
set key top right
plot '${infile1}' using 1:(\$6) w p lt 1 pt 1 ps 0.5 title "N_x=016", \
     '${infile2}' using 1:(\$6) w p lt 2 pt 2 ps 0.5 title "N_x=032", \
     '${infile3}' using 1:(\$6) w p lt 3 pt 3 ps 0.5 title "N_x=064", \
     '${infile4}' using 1:(\$6) w p lt 4 pt 4 ps 0.5 title "N_x=128"
#pause -1
EOF

gnuplot gnu.plt
convert -density 300 -quality 100 ${outfile}.eps ${outfile}.jpeg

