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
MAKE_UNAME=locallibs
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


################### --- KAY at ICHEC.IE ---######################
# Options for fionn.ichec.ie
######################################################################
case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    source /usr/share/Modules/init/bash
    #module purge
    module load intel
    module list
    MAKE_UNAME=KAY
    NCORES=8
    export CC=mpiicc
    export CXX=mpiicpc
    # -DINTEL means the code uses the intel math headers instead of gnu.
    #export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL"
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DINTEL"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    PION_PATH=`pwd`
    PION_PATH=${PION_PATH}/../extra_libraries/lib
    export LD_LIBRARY_PATH=${PION_PATH}${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH:-}
    echo $LD_LIBRARY_PATH
    echo "***** COMPILING WITH KAY: COMPILERS ARE $CC $CXX "  
    ;;
esac
################### --- KAY at ICHEC.IE ---######################


#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
#  export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS"
#  export PION_OPTIMISE=HIGH
  export CXX=mpicxx
  export CC=mpicc
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=OSX
  #NCORES=1
  path=`pwd`
fi
#################################

#########################
### TEST FOR SuperMUC ###
#########################
case $HOST in
  login[0-9][0-9])
    echo "Compiling on SuperMUC"
    echo "MKL LIBS: $MKL_SHLIB"
    export CC=mpicc
    export CXX=mpiCC
    export FC=mpif90
    MAKE_UNAME=SUPERMUC
    NCORES=8
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL"
    export PION_OPTIMISE=HIGH
  ;;
esac
#######################

###################################
### TEST FOR DIRAC-2-COMPLEXITY ###
###################################
case $HOSTNAME in
  dirac[0-9][0-9])
    echo "Compiling on DIRAC-Complexity"
    #module list
    module load intel/compilers/13.0.0 intel/mkl/11.0.0 intel/impi/4.1.0
    module list
    MAKE_UNAME=DIRAC
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    #export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL -DTESTING"
    #export PION_OPTIMISE=LOW
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL"
    export PION_OPTIMISE=HIGH
    export CXX=mpiicpc
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
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export PION_OPTIONS="-DPARALLEL -DUSE_MPI -DSILO -DFITS -DINTEL"
    export PION_OPTIMISE=HIGH
    export CXX=mpicxx
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

