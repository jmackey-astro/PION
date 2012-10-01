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
#set terminal postscript eps enhanced color 
#set output 'jono$arg1.eps' 
#set size square
#set size 0.5,0.5
#set cntrparam levels 35
#set contour base
#show contour
#set grid front
#set cntrparam levels 100
#set samples 100
#set contour base
set xrange [:]
set yrange [:]
set title "Shock Tube Problem $file"
set ylabel "Variable"
set xlabel "X"
unset key 
#set yrange [0.95:1.05]
plot '$file' using 1:$element with lines 
#plot 'Output.88' using 1 with lines
#plot 'Output.218' using 3
pause -1
EOF
gnuplot gnu.plt
