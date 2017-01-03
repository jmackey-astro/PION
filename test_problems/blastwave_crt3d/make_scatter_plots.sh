#!/bin/bash


DDIR=$1
RES=$2
echo $DDIR
echo $RES

for solver in "HYB_FKJav01" "RPV_FKJav01" "FVS_FKJav01" "RCV_FKJav01" "RCV_Hcorr"; do
  bash ./make_scatterplot.sh \
    $DDIR BWcrt3Dpll_Octant_NR${RES}_${solver}_0000 \
    $DDIR/IMG Density_BWcrt3Dpll_Octant_NR${RES}_${solver} \
    "BW3D ${solver} N=${RES}" 1
done

exit

