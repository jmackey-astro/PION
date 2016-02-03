#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename element-name
file=out_0001
element=1
	
else
file=$1
#export arg1=`printf "%04d" $1`
#file=test.$1
element=$2
fi

cat << EOF  > gnu.plt
#set terminal postscript enhanced color 
#set output 'file$arg1.ps' 
#set size square
set grid front
#set cntrparam levels 100
#set samples 100
set cntrparam levels 8
set contour base
show contour
set xlabel "X"
set ylabel "Y"
set xrange [:]
set yrange [:]
set zrange [:]
#set zrange [-4.e-14:4.e-14]
splot '$file' using 1:2:$element with lines
#plot 'Output.88' using 1 with lines
#plot 'Output.218' using 3
pause -1
EOF
gnuplot gnu.plt
