#!/bin/bash
#
# Script to make a movie from 1D simulations of HII region and wind around a massive star.
#

FDIR=$1
FBASE=$2
NTR=$3
MAG=$4

IMGDIR=${FDIR}/img2
if [ ! -d "${IMGDIR}" ]; then
  mkdir $IMGDIR
fi

#
# 1D data, so first cols are <pos> <rho> <p> <vx> <vy> <vz>
# If MHD then we have another 3 columns, and then tracers,
# and finally temperature.

IonFracVar=7
ColourVar=8
TempVar=7

TempVar=$[TempVar+$NTR]

if [ "$MAG" = "yes" ]; then
  IonFracVar=$[IonFracVar+3]
  ColourVar=$[ColourVar+3]
  TempVar=$[TempVar+3]
fi

cat << EOF  > gnu.plt
set terminal png truecolor notransparent linewidth 2 font "times,22" enhanced size 1024,768
set xrange [*:*]
#set xrange [0:1.2]
set yrange [1.0e-3:20000.0]
set ylabel "rho,T,v" offset 3,0
set xlabel "r (pc)" offset 0,0.5
set rmargin 1.0
set lmargin 7.0
set tmargin 2.0
set bmargin 3.0
#set log x
set log y
set grid
Myr=3.156e13
PC=3.086e18
EOF
#
# Get list of files, and write code to make image for each one.
#
LIST=(`ls ${FDIR}/${FBASE}*.txt`)
for ii in $(seq 0 $((${#LIST[@]} - 1)))
do
    TITLE=`grep "time =" ${LIST[$ii]}`
    set -- $TITLE
    Time=$4
    printf -v Time "%f" "$Time"
    Myr=3.156e13
    printf -v Myr "%f" "$Myr"
    Time=$(echo "scale=6; $Time/$Myr" | bc);
    printf -v Time "%f" "$Time";
    TITLE="Time = $Time Myr"
    Time=${Time:2:6}
    printf -v num "%03d" "$ii"
    OUTFILE="${FBASE}_${num}.png"  #`echo ${LIST[$ii]} | sed -e "s/txt/png/"`
    #
    cat << EOF  >> gnu.plt
!echo "${IMGDIR}/${OUTFILE}"
set output "${IMGDIR}/${OUTFILE}"
set title "${TITLE}" offset 0,-0.5
plot '${LIST[$ii]}' u (\$1/PC):(\$2*4.277e21)      w l lt 1 lc -1 lw 2 title "n(H) (100 cm-3)", \
                 '' u (\$1/PC):(\$$TempVar*1.0e-4) w l lt 1 lc 1 lw 2 title "Temperature (10^{4}K)", \
                 '' u (\$1/PC):(\$4/1.0e5)    w l lt 3 lc 3 lw 2 title "Velocity (km/s)", \
                 '' u (\$1/PC):((-1)*\$4/1.0e5)    w l lt 3 lc 3 lw 1 notitle, \
                 '' u (\$1/PC):(\$$IonFracVar)     w l lt 1 lc 4 lw 2 title "Ion fraction", \
                 '' u (\$1/PC):(\$$ColourVar)      w l lt 5 lc 0 lw 1 title "Contact"

EOF
done

gnuplot gnu.plt

ffmpeg -y -r 4.0 -f image2 -i ${IMGDIR}/${FBASE}_%03d.png -q:v 0 -s 1024x768 -pix_fmt yuv420p -vcodec h264 ${IMGDIR}/${FBASE}.mp4




