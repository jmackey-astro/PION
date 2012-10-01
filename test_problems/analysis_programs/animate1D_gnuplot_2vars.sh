#!/bin/bash
#
# 2010.10.04, author: Jonathan Mackey
# Script to plot a sequence of files on the same scale with gnuplot
# This version plots the difference between two variables from the same file.
#
# Assumes x-position is column 1 of the file.
# Makes most sense when each file is a snapshot in time of a 1D line.
#
if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ] ||  [ "$4" = "" ] ||  [ "$5" = "" ] ||  [ "$6" = "" ] 
then
 echo Usage: $0 FileBase Var.Column var.min var.max var.name var.log?
 exit
fi

FBASE=$1
VAR1=$2
VAR2=$3
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
EOF

if [ "$VLOG" = "yes" ]
then 
    cat << EOF  >> gnu.plt
set log y
EOF
fi

LIST=(`ls ${FBASE}*.txt`)
for ii in $(seq 0 $((${#LIST[@]} - 1)))
do
    cat << EOF  >> gnu.plt
plot '${LIST[$ii]}' u 1:(abs(1.0-(\$${VAR1}/\$${VAR2}))) w lp lt 2
pause 0.2
EOF
done

gnuplot gnu.plt
exit
#plot '${LIST[$ii]}' u 1:(abs(1.0-(\$7/\$17))) w lp lt 2
#set yrange [1.e-8:1.0]
#set grid x
#set grid y

