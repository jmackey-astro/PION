#!/bin/bash
#
# 2010.10.04, author: Jonathan Mackey
# Script to plot a sequence of files on the same scale with gnuplot
#
# Assumes x-position is column 1 of the file.
# Makes most sense when each file is a snapshot in time of a 1D line.
#
if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ] ||  [ "$4" = "" ] ||  [ "$5" = "" ] ||  [ "$6" = "" ] ||  [ "$7" = "" ] ||  [ "$8" = "" ] 
then
 echo Usage: $0 File1Base File2Base File3Base Var.Column var.min var.max var.name var.log?
 exit
fi

F1BASE=$1
F2BASE=$2
F3BASE=$3
VAR=$4
VMIN=$5
VMAX=$6
VNAME=$7
VLOG=$8

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

LIST1=(`ls ${F1BASE}*`)
LIST2=(`ls ${F2BASE}*`)
LIST3=(`ls ${F3BASE}*`)
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
plot '${LIST1[$ii]}' u 1:${VAR} w lp lt 2 lw 3, \
     '${LIST2[$ii]}' u 1:${VAR} w lp lt 3 lw 2, \
     '${LIST3[$ii]}' u 1:${VAR} w lp lt 4 lw 2
pause -1
EOF
done

gnuplot gnu.plt
exit

