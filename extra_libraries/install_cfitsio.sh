#!/bin/bash
#
# 2010.10.11 JM: Script to download, configure, compile and install cfitsio
# 2011.04.14 JM: Changed install path.
#
echo "*******************************"
echo "*** INSTALLING FITS LIBRARY ***"
echo "*******************************"
#################################
# Change these for new versions:
FILE=cfitsio3250.tar.gz
SRC_DIR=cfitsio
REMOTE_URL=ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3250.tar.gz
#################################

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
### HACK FOR PHALANX ICC/ICPC ###
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
### HACK FOR PHALANX ICC/ICPC ###
#################################


if [ -e $FILE ]; then
	echo "*** File exists, no need to download ***"
else 
	echo "***** File does not exist ******"
	echo "*******************************"
	echo "*** DOWNLOADING FITS LIBRARY ***"
	echo "*******************************"
	wget $REMOTE_URL
fi
echo "*******************************"
echo "*** EXTRACTING FITS LIBRARY ***"
echo "*******************************"
tar zxf $FILE
echo "*******************************"
echo "*** RUNNING CONFIGURE ***"
echo "*******************************"
BASE_PATH=`pwd`
echo "***Path = $BASE_PATH ***"
cd $SRC_DIR
make clean
./configure --prefix=${BASE_PATH}
echo "*******************************"
echo "*** RUNNING MAKE ***"
echo "*******************************"
make
echo "*******************************"
echo "*** INSTALLING CFITSIO ***"
echo "*******************************"
make install
echo "*******************************"
echo "*** FINISHED! ***"
echo "*******************************"


