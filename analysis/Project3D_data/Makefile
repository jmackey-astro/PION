# Makefile for sim_projection code, Jonathan Mackey
#
# 2010.12.13 JM: Changed aibn129 to standard workstation.
#   Updated Vpath to point to the tagged release used for MHD code.
#

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cxx .o
srcdir = .

GLOBAL = global.h functionality_flags.h testing_flags.h


# --------------------
# Compile-time options
# --------------------
# Tell it that we're using MPI
OPTION  = -DPARALLEL
#IOMODE=USE_FILE_COMMS
IOMODE=USE_MPI
OPTION += -D$(IOMODE)
# if using Silo
OPTION += -DSILO
OPTION += -DFITS

#
# 2010.12.13 JM: New calculation of projected field! (Thanks to Will Henney!)
#
OPTION += -DNEW_STOKES_CALC

#
# This is to multi-thread the pixel calculations.
#
#MT=USE_THREADS
#OPTION += -DTHREADS -D__WXGTK__
MT=NO_THREADS

#MAKE_UNAME=JUROPA
MAKE_UNAME=standard
#MAKE_UNAME=imac
#MAKE_UNAME=epona
#MAKE_UNAME=stokes
#MAKE_UNAME=aurora
#MAKE_UNAME=furfur

# All machines and options need readline and maths
LDFLAGS = -lreadline
ifeq ($(MT), USE_THREADS)
LDFLAGS += -pthread
endif

###########################################
# Options for standard workstation (e.g. Ubuntu 10.04 64 bit)
# make sure cfitsio and silo are compiled in ../extra_libraries/
###########################################
ifeq ($(MAKE_UNAME), standard)

  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif
  #CXXFLAGS = -O3
  CXXFLAGS = -O0 -Wall -g
  CXXFLAGS+= -Wno-unused-function


  MYLIB = $(CURDIR)/../../extra_libraries

  INC = -I${MYLIB}/include
  LDFLAGS =  -lm  -lreadline -lncurses
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES

  VPATH=../../source:../../source/defines
  #VPATH=../../../uniform_grid_code/tags/ET_MHD2010_r215/source
endif

######################## -- JUROPA -- ################################
# Options for the Juelich Supercomputer JUROPA.  Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), JUROPA)
  # intel math-header
  OPTION += -DINTEL

  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = icc
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpic++
  endif
  CXXFLAGS = -O3
  #CXXFLAGS = -O0 -Wall -g -wd981 -wd383
  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  #CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check

  MYLIB = ../../extra_libraries
  INC  = -I${MYLIB}/include
  LDFLAGS  = -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -lmkl_sequential -lmkl_core
  LDFLAGS += -lreadline

  VPATH=../../source:../../source/defines
endif
##################### -- JUROPA -- ###################################



###########################################
# Specific Options for my imac
###########################################
ifeq ($(MAKE_UNAME), imac)
ifeq ($(IOMODE), USE_FILE_COMMS)
CXX = g++
endif
ifeq ($(IOMODE), USE_MPI)
CXX = mpicxx
endif

#CXXFLAGS = -O3
CXXFLAGS = -O0
#CXXFLAGS += -Wall
CXXFLAGS += -g    # enable debugging with gdb
CXXFLAGS += -arch x86_64

MYLIB = /Users/jmackey/Library/my_local_libs
# for Silo, which can use HDF5.
INC += -I$(MYLIB)/silo_gcc64/include
INC += -I$(MYLIB)/cfitsio64/include

# Libraries for fits, silo (with hdf5)
LDFLAGS += -L$(MYLIB)/cfitsio64/lib/ -lcfitsio
LDFLAGS += -L$(MYLIB)/silo_gcc64/lib/ -lsilo

endif
#############################################

###########################################
# Options for klaipedaX.astro.uni-bonn.de (Debian 64 bit)
###########################################
ifeq ($(MAKE_UNAME), klaipeda)
ifeq ($(IOMODE), USE_FILE_COMMS)
CXX = g++
endif
ifeq ($(IOMODE), USE_MPI)
CXX = mpicxx
endif
#CXXFLAGS = -O3
CXXFLAGS = -O3 -Wall -g

MYLIB = /scratch/jmackey/extra_libraries

# Includes for FITS/SILO
INC = -I${MYLIB}/JM_cfitsio_gcc/include
INC += -I${MYLIB}/JM_silo_gcc/include

LDFLAGS =  -lm  -lreadline

## If we want static linking this is the flag for it!
#LDFLAGS += -static
# Libraries for fits, silo (with system hdf5)
LDFLAGS += -L${MYLIB}/JM_cfitsio_gcc/lib -lcfitsio
LDFLAGS += -L${MYLIB}/JM_silo_gcc/lib -lsilo
#LDFLAGS += -L${MYLIB}/hdf5/1.8.4/linux-i686_gcc-4.4/lib -lhdf5
#LDFLAGS += -L${MYLIB}/szip/2.1/linux-i686_gcc-4.4/lib -lsz
LDFLAGS += -lz -lpthread
endif


#############################################
# Specific Options for STOKES - ICHEC cluster
# Intel Compilers! (OR gcc -- need to choose which and comment out section below)
#############################################
# module load mvapich2-intel intel-cc intel-mkl valgrind ddt
#############################################
ifeq ($(MAKE_UNAME), stokes)
CXX = mpicxx
CXXFLAGS = -O3  # optimised code
#CXXFLAGS = -O0 -Wall -g -wd981 -wd383
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -fp-stack-check -fstack-security-check
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

CXXFLAGS += $(OPTION)

OBJECTS = uniformGrid.o dataio.o readparams.o global.o                          \
	dataio_utility.o dataio_silo.o dataio_fits.o uniformGridMPI.o           \
	comm_mpi.o comm_files.o utility_fits_class.o                            \
	main_projection.o sim_projection.o image_io.o dataio_silo_MPI.o         \
        VectorOps.o VectorOps_spherical.o cell_interface.o                      \
        stellar_wind_BC.o 
        

ifeq ($(MT), USE_THREADS)
OBJECTS += logmessages.o threadpool.o threads.o
VPATH +=:./andys_threads:./andys_threads/threadpool
THREAD_INCLUDE=threads.h msvc_constants.h constants.h
endif

projection: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o projection $(LDFLAGS)

# make has rules for making objects from .cc files, so I don't need them
# explicitly, unless I am doing something unusual.
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<

################################
# Andy's multi-threading files #
################################
logmessages.o: logmessages.cxx logmessages.h constants.h msvc_constants.h threads.h
threadpool.o: threadpool.cxx threadpool.h constants.h msvc_constants.h logmessages.h threads.h
threads.o: threads.cxx msvc_constants.h constants.h threads.h

#################################
# Stuff From my Sim Code Source #
#################################
global.o: 	  global.cc ${GLOBAL}

cell_interface.o: cell_interface.cc cell_interface.h ${GLOBAL}

uniformGrid.o: 	  uniformGrid.cc uniformGrid.h ${GLOBAL} stellar_wind_BC.h

uniformGridMPI.o: uniformGridMPI.cc uniformGrid.h ${GLOBAL} comms.h

VectorOps.o: 	VectorOps.cc VectorOps.h ${GLOBAL} uniformGrid.h

VectorOps_spherical.o: VectorOps_spherical.cc VectorOps_spherical.h VectorOps.h \
		       ${GLOBAL} uniformGrid.h


dataio.o:       dataio.cc dataio.h ${GLOBAL} solver_eqn_base.h readparams.h \
		stellar_wind_BC.h
dataio_silo.o:  dataio_silo.cc dataio_silo.h dataio.h ${GLOBAL} solver_eqn_base.h file_status.h
dataio_fits.o:  dataio_fits.cc dataio_fits.h dataio.h ${GLOBAL} solver_eqn_base.h \
		solver_eqn_hydro_adi.h solver_eqn_hydro_iso.h file_status.h
readparams.o:	readparams.cc readparams.h ${GLOBAL}
dataio_silo_MPI.o:    dataio_silo_MPI.cc dataio_silo.h dataio.h ${GLOBAL} comms.h file_status.h
utility_fits_class.o: utility_fits_class.cc dataio_fits.h dataio.h ${GLOBAL}

comm_mpi.o: comm_mpi.cc comm_mpi.h comms.h ${GLOBAL}
comm_files.o: comm_files.cc comm_files.h comms.h ${GLOBAL}

#dataio_utility.o:	dataio_utility.cc dataio_utility.h dataio_silo.h dataio.h ${GLOBAL}
stellar_wind_BC.o: stellar_wind_BC.cc stellar_wind_BC.h ${GLOBAL}

###############################
# Stuff Specific to this code #
###############################
sim_projection.o:	sim_projection.cc sim_projection.h ${GLOBAL} uniformGrid.h dataio.h
main_projection.o:	main_projection.cc sim_projection.h ${GLOBAL} uniformGrid.h dataio.h ${THREAD_INCLUDE} 
image_io.o:		image_io.cc image_io.h dataio.h dataio_fits.h ${GLOBAL}

clean:
	rm -f *.o projection
