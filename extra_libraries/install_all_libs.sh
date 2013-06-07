#!/bin/bash
#
# - 2012.02.22 JM/HD: Added options for compiling on UCL-dougal
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.02.07 JM: Updated for new library versions (all).


mkdir include
mkdir bin
mkdir lib

MAKE_UNAME=standard
NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran


#################################
### TEST FOR Dougal ICC/ICPC ###
#################################
if [ "${HOST}" = 'dougal.hpc.phys.ucl.ac.uk' ]; then
    source /opt/intel/Compiler/11.1/046/bin/ifortvars.sh intel64
    source /opt/intel/Compiler/11.1/046/bin/iccvars.sh intel64
    export CC=icc
    export CXX=icpc
    export FC=ifort
    echo "***** COMPILING ON ${HOST}: COMPILERS ARE $CC $CXX "  
    MAKE_UNAME=dougal
    NCORES=8
fi
#################################

#################################
### TEST FOR PHALANX ICC/ICPC ###
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
    echo "***** COMPILING WITH PHALANX: COMPILERS ARE $CC $CXX "  
    MAKE_UNAME=phalanx
    NCORES=8
fi
#################################

#######################
### TEST FOR JUROPA ###
#######################
case $HOST in
  jj[0-9][0-9]l[0-9][0-9])
    echo "Compiling on JUROPA"
    module purge
    #module load mkl/10.2.2.025 intel/11.1.059 sundials/2.4.0 parastation/mpi2-intel-5.0.25-2
    module load mkl/10.2.5.035 intel/11.1.072 parastation/mpi2-intel-5.0.26-1
    export CC=icc
    export CXX=icpc
    export FC=ifort
    MAKE_UNAME=JUROPA
    NCORES=8
    ;;
esac
#######################

#######################
## TEST FOR SuperMUC ##
#######################
case $HOST in
  login[0-9][0-9])
    echo "Compiling on SuperMUC"
    export CC=icc
    export CXX=icpc
    export FC=ifort
    MAKE_UNAME=SUPERMUC
    NCORES=8
    ;;
esac
#######################

###################################
### TEST FOR DIRAC-2-COMPLEXITY ###
###################################
case $HOSTNAME in
  dirac[0-9][0-9])
    echo "Compiling on DIRAC-Complexity"
    module list
    module load intel/compilers/13.0.0 intel/impi/4.1.0 intel/mkl/11.0.0 cfitsio
    MAKE_UNAME=DIRAC
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export CC=icc
    export CXX=icpc
    export FC=ifort
  ;;
esac
#######################



export MAKE_UNAME
export NCORES
echo "COMPILING WITH MACHINE: ${MAKE_UNAME}. Compilers: CC=$CC FC=$FC CXX=$CXX"
CURDIR=`pwd`

##################################
##########     SILO     ##########
##################################
echo "********************************"
echo "*** INSTALLING SILO LIBRARY ****"
echo "********************************"
#################################
# Change these for new versions:
FILE=silo-4.9-bsd.tar.gz
SRC_DIR=silo-4.9-bsd
REMOTE_URL=https://wci.llnl.gov/codes/silo/silo-4.9/silo-4.9-bsd.tar.gz
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
./configure --prefix=${BASE_PATH} \
--enable-browser \
--disable-fortran \
--disable-silex --with-readline

#--with-readline \
#--enable-hdf5=no
echo "********************************"
echo "*** RUNNING MAKE ***"
echo "********************************"
make -j$NCORES
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
cd $CURDIR
echo "********************************"
echo "*** FINISHED! ***"
echo "********************************"

##################################
##########   SUNDIALS   ##########
##################################

#################################
# Change these for new versions:
FILE=sundials-2.5.0.tar.gz
SRC_DIR=sundials-2.5.0
REMOTE_URL=https://computation.llnl.gov/casc/sundials/download/download.html
echo "********************************"
echo "*** INSTALLING CVODES LIBRARY FILE=${FILE}****"
echo "********************************"
#################################
if [ -e $FILE ]; then
	echo "*** File downloaded, moving on... ***"
else 
	echo "***** File does not exist ******"
	echo "********************************"
	echo "*** Automatic download is not permitted so please: "
        echo "***  -Download ${FILE} from ${REMOTE_URL}"
        echo "***  -Save to current directory."
        echo "***  -Re-run this script."
	echo "********************************"
	exit
fi 
export CFLAGS='-O3'
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
make -j$NCORES
echo "*********************************"
echo "*** INSTALLING CVODES LIBRARY ***"
echo "*********************************"
make install
cd $CURDIR
echo "********************************"
echo "*** FINISHED! ***"
echo "********************************"

##################################
##########   CFITSIO    ##########
##################################
echo "*******************************"
echo "*** INSTALLING FITS LIBRARY ***"
echo "*******************************"
#################################
# Change these for new versions:
FILE=cfitsio3310.tar.gz
SRC_DIR=cfitsio
REMOTE_URL=ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3310.tar.gz
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
make -j$NCORES
echo "*******************************"
echo "*** INSTALLING CFITSIO ***"
echo "*******************************"
make install
echo "*******************************"
echo "*** FINISHED! ***"
echo "*******************************"
cd $CURDIR

