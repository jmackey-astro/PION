#!/bin/bash
#
# Script to plot evolution of density, pressure and velocity on the 
# same graph for a 1D spherically symmetric stellar wind simulation.
# 
# Density is in g/cm3
# Pressure in dyne/cm2
# Velocity in cm/s

if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ] ||  [ "$4" = "" ] ||  [ "$5" = "" ] ||  [ "$6" = "" ] ||  [ "$7" = "" ]
then
 echo Usage: $0 FileBase x.min x.max rho_max pg_max vel_max T_max
 exit
fi

FBASE=$1
XMIN=$2
XMAX=$3
DMAX=$4
PMAX=$5
VMAX=$6
TMAX=$7

cat << EOF  > gnu.plt
#set terminal postscript enhanced color 
#set output 'file$arg1.ps' 
set size 1,1
#set term x11 noraise size 1500,1000
set xrange [${XMIN}:${XMAX}]
set yrange [1.0e-8:3.0]
set ylabel "rho,P,v"
#set log x
set log y
set grid
EOF

LIST=`ls ${FBASE}*.txt`
for ii in $LIST
do
    TITLE=`grep "time =" $ii`
    cat << EOF  >> gnu.plt
set title "$ii ${TITLE}"
plot '$ii' u (\$1*3.24e-19):(\$2/$DMAX) w lp lt -1 title "Density ($DMAX g/cm3)", \
                 '' u (\$1*3.24e-19):(\$3/$PMAX) w l lt 2 lw 2     title "Pressure ($PMAX dyne/cm2", \
                 '' u (\$1*3.24e-19):(abs(\$4/$VMAX)+1.0e-15) w l lt 3 lw 2 title "Velocity ($VMAX cm/s)", \
                 '' u (\$1*3.24e-19):(\$3*1.45e-8/(\$2*$TMAX)) w l lt 4 lw 2 title "~Temperature ($TMAX K)", \
                 '' u (\$1*3.24e-19):7 w l lt 1 lw 2 title "Ion fraction, x"
#pause 0.1
pause -1
EOF
done

gnuplot gnu.plt

exit
