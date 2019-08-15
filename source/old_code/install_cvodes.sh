#!/bin/bash
#
# Install a basic CVODES library
# Download must be done manually, by free registration.  Please put the file
# sundials-2.4.0.tar.gz in this directory.
# Available under a BSD license.
#
echo "********************************"
echo "*** INSTALLING CVODES LIBRARY ****"
echo "********************************"

#################################
# Change these for new versions:
FILE=sundials-2.4.0.tar.gz
SRC_DIR=sundials-2.4.0
REMOTE_URL=https://computation.llnl.gov/casc/sundials/download/download.html
#################################

if [ -e $FILE ]; then
	echo "*** File downloaded, moving on... ***"
else 
	echo "***** File does not exist ******"
	echo "********************************"
	echo "*** Automatic download is not permitted so please: "
        echo "*** 1) Download $FILE from $REMOTE_URL"
        echo "*** 2) Save to current directory."
        echo "*** 3) Re-run this script."
	echo "********************************"
	exit
fi 

export CC=gcc
export CXX=g++
export F77=gfortran
export CFLAGS='-O3'

#################################
### EDIT FOR JUROPA #############
#################################
#export CC=icc
#export CXX=icpc
#export FC=ifort
#################################
### EDIT FOR JUROPA #############
#################################


#################################
### HACK FOR PHALANX ICC/ICPC ###
#################################
if [ "${HOST}" = 'phalanx.star.ucl.ac.uk' ]; then
    source /opt/intel/Compiler/11.1/073/bin/ifortvars.sh intel64
    source /opt/intel/Compiler/11.1/073/bin/iccvars.sh intel64
    export SGIMPT=/opt/sgi/mpt/mpt-1.26
    export PATH=$SGIMPT/bin:/opt/sgi/perfcatcher/bin:$PATH
    export LD_LIBRARY_PATH=$SGIMPT/lib:$LD_LIBRARY_PATH
    export CC=icc
    export CXX=icpc
    export FC=ifort
    echo "***** PHALANX $FC $CC $CXX "
fi
#################################
### HACK FOR PHALANX ICC/ICPC ###
#################################
 

echo "***********************************"
echo "*** EXTRACTING SUNDIALS LIBRARY ***"
echo "***********************************"
tar zxf $FILE
echo "***********************************"
echo "*** RUNNING CONFIGURE ***"
echo "***********************************"
BASE_PATH=`pwd`
echo "***Path = $BASE_PATH ***"
cd $SRC_DIR
make distclean
./configure --prefix=${BASE_PATH} --disable-shared --disable-ida \
 --disable-idas --disable-kinsol --disable-cpodes --disable-mpi --disable-fcmix \
 --disable-blas --disable-lapack
echo "********************************"
echo "*** RUNNING MAKE ***"
echo "********************************"
make
echo "*********************************"
echo "*** INSTALLING CVODES LIBRARY ***"
echo "*********************************"
make install
cd -
echo "********************************"
echo "*** FINISHED! ***"
echo "********************************"
