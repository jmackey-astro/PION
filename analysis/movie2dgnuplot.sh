#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
 echo Usage: $0 filename step element-no
 file=out_0001
 element=1
else
 file=$1
 step=$2
 element=$3
fi
counter=0

for i in $( ls *$file* ); 
do
src=$file.$counter
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
#set xrange [0:1]
#set yrange [0:1]
set xrange [:]
set yrange [:]
#set zrange [0.1249:0.1251]
set zrange [:]
splot '$src' using 1:2:$element with lines
#plot 'Output.88' using 1 with lines
#plot 'Output.218' using 3
pause 1
EOF
gnuplot gnu.plt
rm gnu.plt
echo The counter is $counter
let counter=counter+$step
done
