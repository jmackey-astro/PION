#!/bin/bash
#
# 2010.10.04, author: Jonathan Mackey
# Script to plot a sequence of files on the same scale with gnuplot
#
# Assumes x-position is column 1 of the file.
# Makes most sense when each file is a snapshot in time of a 1D line.
#
if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ] ||  [ "$4" = "" ] ||  [ "$5" = "" ] ||  [ "$6" = "" ] ||  [ "$7" = "" ] 
then
 echo Usage: $0 File1Base File2Base Var.Column var.min var.max var.name var.log?
 exit
fi

F1BASE=$1
F2BASE=$2
VAR=$3
VMIN=$4
VMAX=$5
VNAME=$6
VLOG=$7

cat << EOF  > gnu.plt
#set terminal postscript enhanced color 
#set output 'file$arg1.ps' 
set size 1,1
set term x11 noraise size 1500,500
#set xrange [0:1.0e19]
set yrange [${VMIN}:${VMAX}]
set ylabel "${VNAME}"
#
# f(x) is an empirically normalised fit to the pressure profile in a test
# problem.
f(x)=7.e38*x**(-10./3.)
#
EOF

if [ "$VLOG" = "yes" ]
then 
    cat << EOF  >> gnu.plt
set log y
EOF
fi

LIST1=(`ls ${F1BASE}*.txt`)
LIST2=(`ls ${F2BASE}*.txt`)
#
# These are two plot commands: one also plots f(x), the other doesn't.  Pick
# one and put it in the plot line below; also the pause time between plots is
# in seconds.
#
# plot '${LIST1[$ii]}' u 1:${VAR} w l lt 1 lw 3, '${LIST2[$ii]}' u 1:${VAR} w l lt 3 lw 1, f(x) lt 4
# plot '${LIST1[$ii]}' u 1:${VAR} w l lt 2 lw 3, '${LIST2[$ii]}' u 1:${VAR} w l lt 3 lw 2
#
for ii in $(seq 0 $((${#LIST1[@]} - 1)))
do
    cat << EOF  >> gnu.plt
plot '${LIST1[$ii]}' u 1:${VAR} w l lt 2 lw 3, '${LIST2[$ii]}' u 1:${VAR} w l lt 3 lw 2
pause -1
EOF
done

gnuplot gnu.plt
exit

#plot '${LIST1[$ii]}' u 1:(abs(1.0-(\$3/\$13))) w l lt 2 lw 3, '${LIST2[$ii]}' u 1:(\$3*1.0e6) w l lt 3 lw 2
#plot '${LIST1[$ii]}' u 1:(abs(1.0-(\$7/\$17))) w l lt 2 lw 3, '${LIST2[$ii]}' u 1:7 w l lt 3 lw 2
#set yrange [1.e-10:1.3]
#set grid x
#set grid y


