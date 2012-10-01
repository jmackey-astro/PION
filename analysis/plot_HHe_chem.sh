#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename first-element-name
file=out_0001
elh=1	
else
file=$1
elh=$2
elhe=$elh+2
elhepp=$elhe+1
echo Using $file with elements $elh $elhe $elhepp
fi

cat << EOF  > gnu.plt
set terminal postscript eps color  # no enhanced, as it turns underscore into subscript!
set output 'ST_chemistry_$file.lin.eps' 
set size 1.0,1.0
set origin 0.0,0.0
set label 1 "Shock-Tube Microphysics Tests (1D) - $file" at screen 0.2, screen 0.98
show label 1
unset key 
#
set multiplot
#
set xrange [6.0e16:9.0e16]
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
set yrange [0.75e-22:4.25e-22]
plot '$file'  using 1:2 with points pointsize 0.5 pt 2  #, '$3' using 1:2 w l
#plot '$file'  using 1:2 with lines 
#
# SECOND PLOT
set origin 0.0,0.23
set size 0.9, 0.27
set yrange [*:*]
#set yrange [-0.1:1.1]
#set log y
#set yrange [1.e-12:*]
#set nolog y
set format x ""
set xlabel ""
set ylabel "x(HII)" 0.3
plot '$file'  using 1:$elh with points ps 0.5  #, '$3' using 1:$elh w l
#plot '$file'  using 1:$elh w l
#
# THIRD PLOT
set origin 0.0,0.46
set yrange [*:*]
#set yrange [1.e-14:*]
set ylabel "x(HeII)" 0.5
plot '$file'  using 1:$elhe w p ps 0.5  #, '$3' using 1:$elhe w l
#plot '$file'  using 1:$elhe w l
#
# FOURTH PLOT
set origin 0.0,0.69
set yrange [*:*]
set ylabel "x(HeIII)" 0.0
#set yrange [0:2.0e-12]
plot '$file'  using 1:$elhepp with points ps 0.5  #, '$3' using 1:$elhepp w l
#plot '$file'  using 1:$elhepp w l
#
unset multiplot
#
#
pause -1
EOF
gnuplot gnu.plt
