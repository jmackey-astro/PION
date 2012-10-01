#!/bin/bash
#
# Install a basic silo library, without HDF5 or silex (needs Qt)
# silo from https://wci.llnl.gov/codes/silo/index.html
# Available under a BSD license.
#
# 2011.04.14 JM: changed install path
#
echo "********************************"
echo "*** INSTALLING SILO LIBRARY ****"
echo "********************************"

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
 

#################################
# Change these for new versions:
FILE=silo-4.8-bsd.tar.gz
SRC_DIR=silo-4.8-bsd
REMOTE_URL=https://wci.llnl.gov/codes/silo/silo-4.8/silo-4.8-bsd.tar.gz
#################################

if [ -e $FILE ]; then
	echo "*** File exists, no need to download ***"
else 
	echo "***** File does not exist ******"
	echo "********************************"
	echo "*** DOWNLOADING SILO LIBRARY ***"
	echo "********************************"
	wget --no-check-certificate $REMOTE_URL
fi 
echo "********************************"
echo "*** EXTRACTING SILO LIBRARY ***"
echo "********************************"
tar zxf $FILE
echo "********************************"
echo "*** RUNNING CONFIGURE ***"
echo "********************************"
BASE_PATH=`pwd`
echo "***Path = $BASE_PATH ***"
cd $SRC_DIR
make distclean
./configure --prefix=${BASE_PATH} --enable-hdf5=no --without-readline --enable-fortran=no --enable-silex=no
echo "********************************"
echo "*** RUNNING MAKE ***"
echo "********************************"
make
echo "********************************"
echo "*** RUNNING TESTS ***"
echo "********************************"
#cd tests/
#make check
#cd ..
echo "********************************"
echo "*** INSTALLING SILO LIBRARY ***"
echo "********************************"
make install
cd -
echo "********************************"
echo "*** FINISHED! ***"
echo "********************************"
