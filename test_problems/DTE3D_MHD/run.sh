#!/bin/bash

mpirun -np 8 ../../icgen_parallel params_DTE_D3_KSG07_n00256.txt silo

DDIR=/mnt/data/jm/scratch/KSG07

mv DTE_D3_KSG07_n0256_0000.00000000.silo $DDIR

mpirun -np 8 ../../pion_parallel \
 ${DDIR}/DTE_D3_KSG07_n0256_0000.00000000.silo \
 outfile=${DDIR}/DTE_D3_KSG07_n0256_S4AV1 \
 redirect=${DDIR}/log_DTE_D3_KSG07_n0256_S4AV1 \
 solver=4


