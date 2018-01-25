#!/bin/bash
#
# 2018.01.25 JM: Test for Sam's Simulations
#

DATA_DIR=/mnt/data/jm/BubbleNeb/
OP_DIR=/mnt/data/jm/BubbleNeb/proj3d
mkdir $OP_DIR

MAKE_UNAME=standard make clean
MAKE_UNAME=standard make -j8

#
# projection <data_dir> <input-file> <output-path+file> <ftype=1(FITS)> 
#            <multi-files=1> <LOS> <VERT> <ANGLE-to-LOS> <What2int>
#            [<Nbins> <v_min> <v_max> <smooth?> <how-much>]
#
#

SNAPSHOT=bubble3D_1b_medres_0000.0
VTKFILE=proj3d_ZY_bubble3D_1b_medres
mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t15 3 1 ZN YP 15 7
mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t30 3 1 ZN YP 30 7
mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t45 3 1 ZN YP 45 7
mpirun -np 8 ./projection ${DATA_DIR} $SNAPSHOT ${OP_DIR}/${VTKFILE}_t00 3 1 ZN YP 00 7

exit


