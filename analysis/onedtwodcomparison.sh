#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]||  [ "$3" = "" ]||  [ "$4" = "" ]
then
echo Usage: $0 filename1 filename2 element-name ylabel
file=out_0001
element=1
	
else
file=$1
file2=$2
#export arg1=`printf "%04d" $1`
#file=test.$1
element=$3
yl=$4
fi

cat << EOF  > gnu.plt
set terminal postscript enhanced color 
set output '1d2d_$yl.eps' 
set size 0.5,0.5
#set size square
#set cntrparam levels 35
#set samples 100
#set contour base
#show contour
set grid front
set xrange [:]
set yrange [0:1.1]
set xlabel "x"
set ylabel "$yl"
#set key left bottom
set key 0.45,0.2
# define a new line style with terminal-independent color cyan,
# linewidth 3, and associated point type 6 (a circle with a dot in it).
#    set style line 5 lt rgb "cyan" lw 3 pt 6
set style line 1 lt -1 lw 1 pt 1
set style line 2 lt 2 lw 1 pt 3 ps 0.5
set title "Toro (1999) Test 1, Ch. 6.4"
plot '$file'  using 1:$element with lines  ls 1 title   "1D 4000x1  grid", \
     '$file2' using 1:$element with points ls 2 title "2D 100x100 grid"
#pause -1
EOF
gnuplot gnu.plt
