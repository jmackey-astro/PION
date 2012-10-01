#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ] || [ "$3" = "" ]
then
echo Usage: $0 filename1 filename2 nelements element
file=out_0001
element=1
	
else
file1=$1
#export arg1=`printf "%04d" $1`
#file=test.$1
file2=$2
nels=$3
element=$4
#el2=0
let el2=element+nels
echo el2 is $el2
fi

paste $file1 $file2 > temp.txt
for i in `seq 2 $nels`;
do
  echo $i
  let element=i
  let el2=element+nels
  echo el2 is $el2 and element is $element
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
set title "Difference between two shocktubes"
set ylabel "Variable"
set xlabel "X"
unset key 
#set yrange [0.95:1.05]
plot 'temp.txt' using 1:(\$$element-\$$el2) with points 
#plot 'Output.88' using 1 with lines
#plot 'Output.218' using 3
pause -1
EOF
gnuplot gnu.plt

done
rm temp.txt
