#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename1 filename2
file=out_0001
file2=twodoned
	
else
file=$1
file2=$2
fi

cat << EOF  > gnu.plt
set terminal postscript enhanced color 
set output 'briowu_DensByUint.eps' 
set style line 1 lt -1 lw 1 pt 1
set style line 2 lt 1 lw 1 pt 1 ps 0.5


set size 1.0,0.5
set origin 0.0,0.0
set label 1 "Brio and Wu Test Problem for 1D - MHD" at screen 0.35, screen 0.45
show label 1
set multiplot
# first plot
set size 0.44, 0.44
set size square
set origin 0.0,0.0
set yrange [0:1.05]
set xlabel "x"
set ylabel "Density" 1.3
#set label 1 "Density" at -0.25,0.4 rotate
set xtics 0,.2
set key noautotitles
#set key 0.45,0.2
plot '$file'  using 1:2 with points  ls 2
# second plot
set origin 0.33,0.0
#set yrange [0:1.05]
set yrange [-1.05:1.05]
set ylabel "Tangential Field, B_y" 1.3
#set label 1 "Pressure" at -0.25,0.4 rotate
#plot '$file'  using 1:3 with lines  ls 2 
plot '$file'  using 1:9 with points  ls 2 

# THIRD PLOT
set origin 0.66,0.0
set ylabel "Internal Energy" 1.3
set yrange [0.5:2.5]
#plot '$file'  using 1:7 with lines  ls 1
plot '$file'  using 1:7 with points  ls 2

unset multiplot
#pause -1
EOF
gnuplot gnu.plt
