#!/bin/bash
#
# 2010.10.04, author: Jonathan Mackey
# Script to plot a sequence of files on the same scale with gnuplot
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
VAR=$2
VMIN=$3
VMAX=$4
VNAME=$5
VLOG=$6

cat << EOF  > gnu.plt
#set terminal postscript enhanced color 
#set output 'file$arg1.ps' 
set size 1,1
set term x11 noraise size 1500,500
#set xrange [:]
set yrange [${VMIN}:${VMAX}]
set ylabel "${VNAME}"
EOF

if [ "$VLOG" = "yes" ]
then 
    cat << EOF  >> gnu.plt
set log y
EOF
fi

LIST=(`ls ${FBASE}*`)
for ii in $(seq 0 $((${#LIST[@]} - 1)))
do
    cat << EOF  >> gnu.plt
plot '${LIST[$ii]}' u 1:${VAR} w lp lt 2
pause 0.15
#pause -1
EOF
done

gnuplot gnu.plt
exit

