#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename first-element-name
file=out_0001
elh=1	
else
file=$1
elP=3
elV=4
elT=11
echo Using $file with elements
fi

cat << EOF  > gnu.plt
#set terminal postscript eps color  # no enhanced, as it turns underscore into subscript!
#set output 'rsh_Cooling_$file.log.eps' 
set size 1.0,1.0
set origin 0.0,0.0
set label 1 "Shock-Tube Cooling Tests (1D, 81.4km/s) - $file" at screen 0.2, screen 0.98
show label 1
unset key 
#
set multiplot
#
set xrange [0e15:5e15]
set nokey
#set xtics 10
#set mxtics 10
set lmargin 15
set rmargin 1
# first plot
set size 0.9, 0.27
#set size square
set origin 0.0,0.0
set xlabel "x (cm)" 0,0.2
set ylabel "Density" 0.7
set key noautotitles
set log y
#set yrange [0.99999e-22:1.00001e-22]
plot '$file'  using 1:2 w l #with points pointsize 0.5 pt 2  #, '$3' using 1:2 w l
#plot '$file'  using 1:2 with lines 
#
# SECOND PLOT
set origin 0.0,0.23
set size 0.9, 0.27
set yrange [*:*]
#set yrange [-0.1:1.1]
set log y
#set yrange [1.e-12:*]
#set yrange [-0.1:1.1]
#set nolog y
set format x ""
set xlabel ""
set ylabel "Pressure" 0.7
plot '$file'  using 1:3 with l #points ps 0.5  #, '$3' using 1:$elh w l
#plot '$file'  using 1:$elh w l
#
# THIRD PLOT
set origin 0.0,0.46
set nolog y
set yrange [*:*]
#set yrange [-7.000003e6:-6.999997e6]
set ylabel "Velocity" 0.7
plot '$file'  using 1:4 w l #p ps 0.5  #, '$3' using 1:$elhe w l
#plot '$file'  using 1:$elhe w l
#
# FOURTH PLOT
set origin 0.0,0.69
#set nolog y
set log y
set yrange [*:*]
set ylabel "T (K)" 0.7
#set yrange [5e10:2.0e11]
plot '$file'  using 1:$elT w l #with points ps 0.5  #, '$3' using 1:$elhepp w l
#plot '$file'  using 1:$elhepp w l
#
unset multiplot
#
#
pause -1
EOF
gnuplot gnu.plt

