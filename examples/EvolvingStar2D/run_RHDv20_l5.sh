#!/bin/bash

# Nova
#PION=../nested_pion/pion-ng
#DDIR=/mnt/data/jm/WRNeb
#mkdir -p $DDIR

# local computer for testing:
PION=../../pion-ng
ICGEN=../../icgen-ng
DDIR=./WRNeb
FBASE=35Msun_RHDv20_l5n256
NC=8



mkdir -p $DDIR

mpirun -np $NC $ICGEN params_${FBASE}.txt silo

mpirun -np $NC $PION ${FBASE}_level00_0000.00000000.silo \
 outfile=${DDIR}/${FBASE} \
 redirect=${DDIR}/log_${FBASE}_v1 \
 opfreq=1024 cfl=0.25 finishtime=1.50031e14

RESTART=`ls ${DDIR}/${FBASE}_level00_0000.*.silo | tail -n1`
mpirun -np $NC $PION $RESTART \
 outfile=${DDIR}/${FBASE} \
 redirect=${DDIR}/log_${FBASE}_v2 \
 opfreq=1024 cfl=0.05 finishtime=1.50033e14

RESTART=`ls ${DDIR}/${FBASE}_level00_0000.*.silo | tail -n1`
mpirun -np $NC $PION $RESTART \
 outfile=${DDIR}/${FBASE} \
 redirect=${DDIR}/log_${FBASE}_v3 \
 opfreq=1024 cfl=0.25 finishtime=1.520e14



########## restart b/c of crash ################
#RESTART=`ls ${DDIR}/${FBASE}_level00_0000.*.silo | tail -n1`
#mpirun -np $NC $PION $RESTART \
# redirect=${DDIR}/log_${FBASE}_v4 \
# finishtime=1.600e14
#exit
########## restart b/c of crash ################




