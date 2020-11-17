#!/bin/bash
#
# Author: Jonathan Mackey
# Date:   2010.XX.XX
# Description:
#  Compilation script for the parallel version of the code.
# 
# Modifications:
# - 2010.11.09 JM: Added options for Phalanx at UCL
# - 2011.12.20 JM: Changed completely to work with new Makefile.
#     This script now makes the decisions about compilers and optimisation.
# - 2012.02.22 JM/HD: Added options for Dougal at UCL
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.01.14 JM: Added section for DIRAC/Complexity (it works now).
# - 2013.01.17 JM: Got rid of readline/ncurses from link line in
#    production version of pion.
# - 2013.02.27 JM: Added extensions for contributed code.
#    Added NEW_METALLICITY flag for testing the new microphysics
#    classes.
# - 2013.09.20 JM: New Intel compilers for Juropa.
# - 2014.04.14 JM: Section for Judge at JSC.
# - 2015.01.14 JM: Section for Juropatest system at JSC.
# - 2016.05.04 JM: Added FIONN to list of machines

#
# We first need to set MAKE_UNAME which is an identifier for the computer
# we are compiling on.  If it is not a known computer it is just set to
# "standard" which is a standard linux system, assumed to have 8 cores (for
# compiling code only, not running).
#
MAKE_UNAME=standard
#MAKE_UNAME=locallibs
NCORES=8
#
# Production code options:
#
export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS"
export PION_OPTIMISE=HIGH
export CXX=mpicxx
#
# Debugging code options, using text files to communicate data:
#
#export PION_OPTIONS="-DPARALLEL -DUSE_FILE_COMMS -DSILO -DFITS"
#export PION_OPTIMISE=LOW
#NCORES=1
#export CXX=g++

#################################
## First try debian-based linux #
#################################
LINUX=""
id=""
ver=""
code=""
nc=""
if ! [ -x "$(command -v lsb_release)" ]; then
  LINUX="NO"
else
  LINUX="YES"
  id=`lsb_release -s -i`
  ver=`lsb_release -s -r`
  code=`lsb_release -s -c`
  nc=`nproc --all`
fi

if [ "$LINUX" == "YES" ]; then

  if [ "$id" == "Ubuntu" ] && [ "$ver" == "20.04" ]; then
    echo "Detected Ubuntu 20.04: Note system Silo library has a bug, you must install yourself using the script in PION/extra_libraries"
    MAKE_UNAME=debian10
    export CXX=mpicxx
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE3"
    export PION_OPTIMISE=HIGH
    NCORES=$nc
    #NCORES=1
  elif [ "$id" == "Ubuntu" ] && [ "$ver" == "18.04" ]; then
    echo "Detected Ubuntu 18.04: Note system Silo library has a bug, you must install yourself using the script in PION/extra_libraries"
    MAKE_UNAME=ubuntu18
    export CXX=mpicxx
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    NCORES=$nc
    #NCORES=1
  elif [  "$id" == "Ubuntu" ] && [ "$ver" == "16.04" ]; then
    echo "Detected Ubuntu 16.04 (xenial): compiling extra libraries"
    MAKE_UNAME=ubuntu16
    export CXX=mpicxx
    export CC=mpicc
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    NCORES=$nc
  elif [ "$id" == "Debian" ] && [ "$code" == "stretch" ]; then
    echo "Detected Debian 9 (stretch), using system libs for SILO, FITS, GSL, SUNDIALS"
    MAKE_UNAME=debian9
    export CXX=mpicxx
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE2"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    NCORES=$nc
  elif [ "$id" == "Debian" ] && [ "$code" == "buster" ]; then
    echo "Detected Debian 10 (buster), using system libs for SILO, FITS, GSL, SUNDIALS"
    MAKE_UNAME=debian10
    export CXX=mpicxx
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE3"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    NCORES=$nc
  elif [ "$id" == "ManjaroLinux" ]; then
    echo "Detected ManjaroLinux, using local libs for SILO FITS SUNDIALS"
    echo "Assuming using system libs for MPI and GSL"
    MAKE_UNAME=ManjaroLinux
    export CXX=mpicxx
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    NCORES=$nc
  else
    echo "Failed to find a known version of Linux: checking for other OS types."
    MAKE_UNAME=standard
  fi

fi
#################################


################### --- KAY at ICHEC.IE ---######################
# Options for kay.ichec.ie
######################################################################
case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC: USING SELF-COMPILED LIBRARIES"
    source /usr/share/Modules/init/bash
    module load intel/2018u4
    module load gsl/intel/2.5
    module list
    MAKE_UNAME=KAY
    NCORES=8
    export CC=mpiicc
    export CXX=mpiicpc
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DINTEL -DCVODE5"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    PION_PATH=`pwd`
    PION_PATH=${PION_PATH}/../extra_libraries/lib
    export LD_LIBRARY_PATH=${PION_PATH}${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH:-}
    echo "***** COMPILING WITH KAY: COMPILERS ARE $CC $CXX "  
    ;;
esac
################### --- KAY at ICHEC.IE ---######################


#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  if ! [ -x "$(command -v /opt/local/bin/port)" ]; then
    SL="Homebrew"
  else
    SL="MacPorts"
  fi
  echo "detected support software ${SL} and assuming that is what you are using"
  if [ "$SL" == "MacPorts" ]; then
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE3"
    export PION_OPTIMISE=HIGH
    export CXX=mpicxx
    export CC=mpicc
    echo "*** COMPILING WITH OS-X: host ${HOST}: libs=${SL}: COMPILERS ARE $CC $CXX"
    echo "Make sure you installed silo, sundials, cfitsio with MacPorts"
    MAKE_UNAME=OSX-MP
  elif [ "$SL" == "Homebrew" ]; then
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    export CXX=mpicxx
    export CC=mpicc
    echo "*** COMPILING WITH OS-X: host ${HOST}: libs=${SL}: COMPILERS ARE $CC $CXX"
    echo "Make sure you installed open-mpi, sundials, cfitsio with Homebrew"
    MAKE_UNAME=OSX-HB
  else
    echo "ERROR: Need Macports or Homebrew installed on OSX"
    exit
  fi
  NCORES=`sysctl -n hw.ncpu`
fi
#################################

#########################
### TEST FOR SuperMUC ###
#########################
case $HOST in
  login[0-9][0-9])
    echo "Compiling on SuperMUC: WARNING THIS CONFIG IS >4 YEARS OLD"
    echo "MKL LIBS: $MKL_SHLIB"
    export CC=mpicc
    export CXX=mpiCC
    export FC=mpif90
    MAKE_UNAME=SUPERMUC
    NCORES=8
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL -DCVODE2"
    export PION_OPTIMISE=HIGH
  ;;
esac
#######################

###################################
### TEST FOR DIRAC-2-COMPLEXITY ###
###################################
case $HOSTNAME in
  dirac[0-9][0-9])
    echo "Compiling on DIRAC-Complexity: WARNING THIS CONFIG IS >4 YEARS OLD"
    #module list
    module load intel/compilers/13.0.0 intel/mkl/11.0.0 intel/impi/4.1.0
    module list
    MAKE_UNAME=DIRAC
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    #export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL -DTESTING"
    #export PION_OPTIMISE=LOW
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL -DCVODE2"
    export PION_OPTIMISE=HIGH
    export CXX=mpiicpc
  ;;
esac
#######################


#####################################################################
# For testing/debugging, we need to add -DTESTING to the compile   ##
# flags, and also to add readline and maybe ncurses to the linker. ##
#####################################################################
if [ $PION_OPTIMISE == LOW ]
  then
  echo "LDFLAGS= $DLFLAGS"
  export LDFLAGS=" -lreadline -lncurses "
  echo "LDFLAGS= $LDFLAGS"
  export PION_OPTIONS="$PION_OPTIONS -DTESTING"
else
  export LDFLAGS=""
fi
#####################################################################


#####################################################################
############ EXTRA CODE EXTENSIONS FROM CONTRIBUTED CODE ############
#####################################################################
# Harpreet Dhanoa's chemistry/microphysics module
#PION_OPTIONS+=" -DHARPREETS_CODE_EXT"

# Read in turbulence simulations provided by Blakesley Burkhart
#PION_OPTIONS+=" -DBBTURBULENCE_CODE_EXT"

#PION_OPTIONS="$PION_OPTIONS -DLEGACY_CODE"
#export PION_OPTIONS

PION_OPTIONS+=" -DCODE_EXT_SBII"
export PION_OPTIONS
echo PION_OPTIONS: $PION_OPTIONS

#####################################################################


#####################################################################
#### now compile the code:
#####################################################################
export MAKE_UNAME
echo "COMPILING WITH MACHINE: $MAKE_UNAME"
make -j${NCORES} -f Makefile uni
make -j${NCORES} -f Makefile NG

#####################################################################

exit

