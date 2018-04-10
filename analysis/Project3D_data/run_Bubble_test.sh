#!/bin/bash
#
# 2018.01.25 JM: Test for Sam's Simulations
#

DATA_DIR=~/DIAS/comp_astrophysics/3d_test
OP_DIR=~/DIAS/comp_astrophysics/3d_test

DATA_DIR=/Users/jm/data
OP_DIR=/Users/jm/data
mkdir -p $OP_DIR

#MAKE_UNAME=OSX make clean
MAKE_UNAME=OSX make -j8

#
# projection <data_dir> <input-file> <output-path+file> <ftype=1(FITS)> 
#            <multi-files=1> <LOS> <VERT> <ANGLE-to-LOS> <What2int>
#            [<Nbins> <v_min> <v_max> <smooth?> <how-much>]
#
#

SNAPSHOT=bubble3D_1b_medres_0000.0
VTKFILE=proj3d_ZY_bubble3D_1b_medres
mpirun -np 1 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t15 3 1 ZN YP 25 7 0
#mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t30 3 1 ZN YP 30 7 0
#mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t45 3 1 ZN YP 45 7 0
#mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t00 3 1 ZN YP 00 7 0

exit


