# Makefile for silo2text code, Jonathan Mackey
#
# 2011.12.14 JM: Adapted from silocompare Makefile.
#

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
srcdir = .

# --------------------
# Compile-time options
# --------------------
# Tell it that we're using MPI
OPT  = -DPARALLEL
#IOMODE=USE_FILE_COMMS
IOMODE=USE_MPI
OPT += -D$(IOMODE)
# if using Silo
OPT += -DSILO

#MAKE_UNAME=JUROPA
MAKE_UNAME=standard
#MAKE_UNAME=aibn129
#MAKE_UNAME=epona
#MAKE_UNAME=stokes

# All machines and options need readline and maths
LDFLAGS = -lreadline
# display all warnings
CXXFLAGS = -O3 
#CXXFLAGS += -Wall

######################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)

ifeq ($(IOMODE), USE_FILE_COMMS)
CXX = g++
endif
ifeq ($(IOMODE), USE_MPI)
CXX = mpicxx
endif

#CXXFLAGS = -O3
#CXXFLAGS = -O3 -g
CXXFLAGS = -O0 -Wall -g
CXXFLAGS+= -Wno-unused-function

MYLIB = ../../extra_libraries

# Includes for FITS/SILO
INC = -I${MYLIB}/include
# for Silo, without HDF5 here
INC += -I${MYLIB}/include

LDFLAGS =  -lm  -lreadline

## If we want static linking this is the flag for it!
#LDFLAGS += -static
# Libraries for fits, silo
LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
endif

######################## -- JUROPA -- ################################
# Options for the Juelich Supercomputer JUROPA
# Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), JUROPA)
CXX = mpicxx
CXXFLAGS = -O3  # optimised code
#CXXFLAGS = -O0 -Wall -g -wd981 -wd383   
# -wd981 disables a warning about order of evaluation of operands:
#  see http://software.intel.com/en-us/forums/showthread.php?t=62309
#  and also http://software.intel.com/en-us/articles/cdiag981/
# -wd383 also diasables a warning about passing a temporary object to a function.
#CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -fp-stack-check -fstack-security-check

# CURDIR may only work with GNU Make!  It gives absolute paths, so that
# the executable can always link to shared libs.
MYLIB = $(CURDIR)/../../extra_libraries

# Includes for FITS/SILO/CVODES
INC = -I${MYLIB}/include
#LDFLAGS =  -lm  -lreadline -lncurses
LDFLAGS = -L${MYLIB}/lib
# Libraries for fits, silo
LDFLAGS += -lcfitsio
LDFLAGS += -lsilo 
# libs for CVODES
#LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/cvodes_lib/lib -lsundials_cvodes -lsundials_nvecserial
#LDFLAGS += -llapack -lblas
LDFLAGS += ${MYLIB}/lib/libsundials_cvodes.a ${MYLIB}/lib/libsundials_nvecserial.a

OPT += -DINTEL

# Libraries for MPI/MKL
LDFLAGS += -lmkl_sequential -lmkl_core

endif
##################### -- JUROPA -- ###################################

#############################################
# Specific Options for STOKES - ICHEC cluster
# Intel Compilers! (OR gcc -- need to choose which and comment out section below)
#############################################
# module load mvapich2-intel intel-cc intel-mkl valgrind ddt
#############################################
ifeq ($(MAKE_UNAME), stokes)
ifeq ($(IOMODE), USE_FILE_COMMS)
CXX = icpc
endif
ifeq ($(IOMODE), USE_MPI)
CXX = mpicxx
endif
CXXFLAGS = -O3  # optimised code
#CXXFLAGS = -O0
#CXXFLAGS += -Wall
#CXXFLAGS += -g   # for debugging with ddt. 
#CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -fp-stack-check -fstack-security-check
CXXFLAGS +=  -DINTEL -axSSE4.2,SSSE3
MYLIB = /ichec/work/dsast002b/jm_code/extra_libs
# gcc
#INC  = -I$(MYLIB)/cfitsio_gcc/include
#LDFLAGS += -L$(MYLIB)/cfitsio_gcc/lib/ -lcfitsio
#INC += -I$(MYLIB)/silo_gcc_nohdf/include
#LDFLAGS += -L$(MYLIB)/silo_gcc_nohdf/lib/ -lsilo
#LDFLAGS += -lm
#intel
INC  = -I$(MYLIB)/cfitsio/include
LDFLAGS += -L$(MYLIB)/cfitsio/lib/ -lcfitsio
INC += -I$(MYLIB)/silo_intel_nohdf/include
LDFLAGS += -L$(MYLIB)/silo_intel_nohdf/lib/ -lsilo
LDFLAGS += -lmkl_sequential -lmkl_core
endif


##################################################
##################################################

CXXFLAGS += $(OPT)
VPATH = ../../source

OBJECTS = silo2text.o uniformGrid.o dataio.o readparams.o global.o \
	dataio_silo.o uniformGridMPI.o comm_mpi.o comm_files.o \
	VectorOps.o VectorOps_spherical.o \
	stellar_wind_BC.o cell_interface.o

silo2text: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o silo2text $(LDFLAGS)

# make has rules for making objects from .cc files, so I don't need them
# explicitly, unless I am doing something unusual.
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<

silo2text.o:	silo2text.cc global.h uniformGrid.h dataio.h
global.o: 	  global.cc global.h
uniformGrid.o: 	  uniformGrid.cc uniformGrid.h global.h stellar_wind_BC.h \
		VectorOps.h VectorOps_spherical.h
uniformGridMPI.o: uniformGridMPI.cc uniformGrid.cc uniformGrid.h global.h comms.h
dataio.o:       dataio.cc dataio.h global.h solver_eqn_base.h readparams.h
#dataio.o:	dataio.cc dataio.h global.h solver.h
readparams.o:	readparams.cc readparams.h global.h
dataio_silo.o:  dataio_silo.cc dataio_silo.h dataio.h global.h solver_eqn_base.h
#dataio_silo.o:	dataio_silo.cc dataio_silo.h dataio.h global.h
comm_mpi.o: comm_mpi.cc comm_mpi.h comms.h
comm_files.o: comm_files.cc comm_files.h comms.h global.h

stellar_wind_BC.o: stellar_wind_BC.cc stellar_wind_BC.h global.h
VectorOps.o:           VectorOps.cc VectorOps.h global.h uniformGrid.h
VectorOps_spherical.o: VectorOps_spherical.cc VectorOps_spherical.h VectorOps.h \
		       global.h uniformGrid.h

clean:
	rm -f *.o silo2text

