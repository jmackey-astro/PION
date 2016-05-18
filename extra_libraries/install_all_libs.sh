#!/bin/bash
#
# - 2012.02.22 JM/HD: Added options for compiling on UCL-dougal
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.02.07 JM: Updated for new library versions (all).
# - 2014.04.14 JM: Added section for JUDGE at JSC.
# - 2015.01.14 JM: Section for Juropatest system at JSC.
# - 2016.04.29 JM: updated for sundials 2.6.2 and cfitsio 3390.
# - 2016.05.04 JM: Added FIONN to list of machines.

mkdir include
mkdir bin
mkdir lib

MAKE_UNAME=standard
NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran


#################################
### TEST FOR Dougal ICC/ICPC ###
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
### TEST FOR PHALANX ICC/ICPC ###
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
### TEST FOR JUROPA ###
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
### TEST FOR JUROPATEST ###
#######################
MACHINE=$(cat /etc/FZJ/systemname)
if test "${MACHINE}" = "juropatest"; then
    echo "Compiling on JUROPATEST"
    module purge
    module load intel-para
    export CC=icc
    export CXX=icpc
    export FC=ifort
    MAKE_UNAME=JUROPA
    NCORES=8
fi
#######################

#######################
## TEST FOR SuperMUC ##
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

#######################
### TEST FOR JUDGE ###
#######################
case $HOST in
  judgel[0-9])
    echo "Compiling on JUDGE"
    module purge
    module load intel/11.1.072 mkl/10.2.5.035 parastation/intel
    MAKE_UNAME=JUDGE
    NCORES=8
    export CC=icc
    export CXX=icpc
    export FC=ifort
  ;;
esac
#######################

#######################
### TEST FOR FIONN  ###
#######################
case $HOST in
  fionn[0-9])
    echo "Compiling on FIONN/ICHEC"
    source /usr/share/modules/init/bash
    module purge
    module load dev intel
    #module load dev cmake/intel/latest
    module load dev cmake/intel/3.0.2
    MAKE_UNAME=FIONN
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
FILE=silo-4.10.2-bsd.tar.gz
SRC_DIR=silo-4.10.2-bsd
REMOTE_URL=https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz
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
--disable-silex \
--with-readline \
--enable-pythonmodule

# Silex is broken because I can't get Qt working...
#--enable-silex \
#--with-Qt-include-dir=/usr/include/qt4 --with-Qt-bin-dir=/usr/lib/x86_64-linux-gnu/qt4/bin --with-Qt-lib-dir=/usr/lib/x86_64-linux-gnu --with-Qt-dir=/usr/lib/x86_64-linux-gnu/qt4

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
FILE=sundials-2.6.2.tar.gz
SRC_DIR=sundials-2.6.2
BLD_DIR=sundials_build
REMOTE_URL=http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers/download/sundials-2.6.2.tar.gz
echo "********************************"
echo "*** INSTALLING CVODES LIBRARY FILE=${FILE}****"
echo "********************************"
#################################
if [ -e $FILE ]; then
	echo "*** File downloaded, moving on... ***"
else 
	echo "***** File does not exist ******"
	echo "********************************"
	echo "*** Automatic download of ${FILE}"
        echo "from ${REMOTE_URL}"
        wget --no-check-certificate $REMOTE_URL
        echo "***  Downloaded."
	echo "********************************"
fi 
export CFLAGS='-O3'
echo "***********************************"
echo "*** EXTRACTING SUNDIALS LIBRARY ***"
echo "***********************************"
tar zxf $FILE
echo "***********************************"
echo "*** RUNNING CMAKE CONFIG ***"
echo "***********************************"
BASE_PATH=`pwd`
echo "***Path = $BASE_PATH ***"
mkdir -p $BLD_DIR
cd $BLD_DIR
cmake -DCMAKE_INSTALL_PREFIX=${BASE_PATH} \
 -DEXAMPLES_INSTALL_PATH=${BASE_PATH} -DEXAMPLES_INSTALL=OFF \
 ${BASE_PATH}/${SRC_DIR}
echo "********************************"
echo "*** RUNNING MAKE ***"
echo "********************************"
make -j$NCORES install
echo "*********************************"
echo "*** INSTALLED CVODES LIBRARY ***"
echo "*********************************"
cd $CURDIR
rm -rf $BLD_DIR
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
FILE=cfitsio3390.tar.gz
SRC_DIR=cfitsio
REMOTE_URL=http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3390.tar.gz
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

