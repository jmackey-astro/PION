#!/bin/bash
#
# 2018.01.25 JM: Test calculation
#
# Need to have file xray-table.txt in the path where you are running
# the projection code.

DATA_DIR=~/DIAS/comp_astrophysics/3d_test
OP_DIR=~/DIAS/comp_astrophysics/3d_test
mkdir -p $OP_DIR

SNAPSHOT=test_simulation_0000.0
VTKFILE=proj3d_ZY_test_simulation
mpirun -np 1 ./projection3D $DATA_DIR $SNAPSHOT ${OP_DIR}/${VTKFILE}_t15 3 1 ZN YP 25 20 0

exit


