# Makefile for silo2text code, Jonathan Mackey
#
# 2011.12.14 JM: Adapted from silocompare Makefile.
# 2013.09.05 JM: updated for pion.
# 2016.03.03 JM: updated for pion_dev: uses MCMD, SimSetup, constants.h


############################
### BASIC INFO AND PATHS ###
############################

SHELL = /bin/bash
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cxx .o
.SUFFIXES: .cpp .o

#######################
## COMPILATION RULES ##
#######################
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################


SRCDIR = ${CURDIR}/../../source

bindir=.

CODEDIR=${CURDIR}/../..

VPATH=${SRCDIR}:\
${SRCDIR}/grid:\
${SRCDIR}/defines:\
${SRCDIR}/dataIO:\
${SRCDIR}/comms:\
${SRCDIR}/coord_sys:\
${SRCDIR}/equations:\
${SRCDIR}/flux_calc:\
${SRCDIR}/spatial_solvers:\
${SRCDIR}/time_integrators:\
${SRCDIR}/microphysics:\
${SRCDIR}/raytracing:\
${SRCDIR}/contrib\
${SRCDIR}/tools



GLOBAL = functionality_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h


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
#MAKE_UNAME=standard
#MAKE_UNAME=stokes

# All machines and options need readline and maths
LDFLAGS = -lreadline
# display all warnings
CXXFLAGS = -O3 
#CXXFLAGS += -Wall

######################################################################
# Options for a "standard" linux workstation.
######################################################################
ifeq ($(MAKE_UNAME), standard)

  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  CXXFLAGS = -O3
  #CXXFLAGS = -O3 -g
  #CXXFLAGS = -O0 -Wall -g
  CXXFLAGS+= -Wno-unused-function

  MYLIB = ${CODEDIR}/extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}

  LDFLAGS =  -lm  -lreadline
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################

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



###########################################
# Specific Options for my imac
###########################################
ifeq ($(MAKE_UNAME), OSX)
  ifeq ($(IOMODE), USE_FILE_COMMS)
  CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
  CXX = mpicxx
  endif

  CXXFLAGS = -O3
#  CXXFLAGS = -O0
#CXXFLAGS += -Wall
#  CXXFLAGS += -g    # enable debugging with gdb
  #CXXFLAGS += -arch x86_64

  MYLIB = $(CURDIR)/../../extra_libraries
  #MYLIB = /Users/jmackey/Library/my_local_libs
  INC += -I$(MYLIB)/include
  INC += -I${SRCDIR}

  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  # libs for CVODES
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
#############################################


##################################################
##################################################

CXXFLAGS += $(OPT)

#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################

OBJECTS = \
  silo2text.o                 \
  constants.o \
  timer.o \
  mem_manage.o \
  reporting.o \
  command_line_interface.o \
  interpolate.o \
  sim_params.o \
  MCMD_control.o \
  setup_fixed_grid.o \
  setup_fixed_grid_MPI.o \
  uniform_grid.o              \
  dataio.o                    \
  readparams.o                \
  dataio_silo.o               \
  dataio_silo_utility.o       \
  dataio_silo_MPI.o           \
  uniform_grid_pllel.o        \
  comm_mpi.o                  \
  comm_files.o                \
  VectorOps.o                 \
  VectorOps_spherical.o       \
  stellar_wind_BC.o           \
  cell_interface.o            \
  file_status.o               \
  cooling_SD93_cie.o          \
  hydrogen_recomb_Hummer94.o  \
  hydrogen_mp.o               \
  hydrogen_photoion.o         \
  mp_explicit_H.o             \
  cvode_integrator.o          \
  mp_implicit_H.o             \
  mp_v2_aifa.o                \
  integrator.o                \
  cooling.o                   \
  microphysics_base.o         \
  microphysics.o              \
  microphysics_v1.o           \
  microphysics_lowZ.o         \
  mp_only_cooling.o           \
  mpv5_molecular.o            \
  mpv6_PureH.o                \
  mpv7_TwoTempIso.o           \
  mpv8_StarBench_heatcool.o   \
  raytracer_base.o            \
  raytracer_SC_pllel.o        \
  raytracer_SC.o

##########################
### -- EXECTUTABLES -- ###
##########################

silo2text: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o silo2text $(LDFLAGS)


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################
silo2text.o:	silo2text.cpp uniform_grid.h dataio.h ${GLOBAL}
global.o: global.cpp ${GLOBAL}

uniform_grid.o: uniform_grid.cpp uniform_grid.h ${GLOBAL} stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} comms.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h
dataio.o: dataio.cpp dataio.h ${GLOBAL} solver_eqn_base.h \
    readparams.h stellar_wind_BC.h file_status.h
dataio_silo.o:  dataio_silo.cpp dataio_silo.h dataio.h ${GLOBAL} \
    solver_eqn_base.h file_status.h
readparams.o:	readparams.cpp readparams.h ${GLOBAL}

comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h ${GLOBAL}
comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}

stellar_wind_BC.o: stellar_wind_BC.cpp stellar_wind_BC.h ${GLOBAL}

VectorOps.o: VectorOps.cpp VectorOps.h ${GLOBAL} uniform_grid.h
VectorOps_spherical.o: VectorOps_spherical.cpp VectorOps_spherical.h VectorOps.h \
 ${GLOBAL} uniform_grid.h


######################
# MICROPHYSICS STUFF #
######################

integrator.o: integrator.cpp integrator.h ${GLOBAL}
cooling.o:    cooling.cpp cooling.h ${GLOBAL}

cooling_SD93_cie.o: cooling_SD93_cie.cpp cooling_SD93_cie.h ${GLOBAL}

hydrogen_recomb_Hummer94.o: hydrogen_recomb_Hummer94.h hydrogen_recomb_Hummer94.cpp \
    ${GLOBAL}

mp_only_cooling.o: mp_only_cooling.cpp mp_only_cooling.h ${GLOBAL} \
   cooling_SD93_cie.h microphysics_base.h hydrogen_recomb_Hummer94.h

microphysics.o:	microphysics.cpp microphysics.h microphysics_base.h \
   ${GLOBAL} cooling.h integrator.h
microphysics_v1.o: microphysics_v1.cpp microphysics.h \
   microphysics_base.h ${GLOBAL} cooling.h integrator.h
microphysics_lowZ.o: microphysics_lowZ.cpp microphysics_lowZ.h \
   microphysics_base.h ${GLOBAL} integrator.h

hydrogen_mp.o: hydrogen_mp.cpp hydrogen_mp.h ${GLOBAL} hydrogen_recomb_Hummer94.h \
  hydrogen_photoion.h 

hydrogen_photoion.o: hydrogen_photoion.cpp hydrogen_photoion.h ${GLOBAL}

mp_v2_aifa.o: mp_v2_aifa.cpp mp_v2_aifa.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

cvode_integrator.o: cvode_integrator.cpp cvode_integrator.h ${GLOBAL}

mp_explicit_H.o: mp_explicit_H.cpp mp_explicit_H.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

mp_implicit_H.o: mp_implicit_H.cpp mp_implicit_H.h mp_explicit_H.h \
  ${GLOBAL} hydrogen_mp.h hydrogen_recomb_Hummer94.h hydrogen_photoion.h \
  cooling_SD93_cie.h microphysics_base.h

mpv5_molecular.o: mpv5_molecular.cpp mpv5_molecular.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

mpv6_PureH.o: mpv6_PureH.cpp mpv6_PureH.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

mpv7_TwoTempIso.o: mpv7_TwoTempIso.cpp mpv7_TwoTempIso.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

mpv8_StarBench_heatcool.o: mpv8_StarBench_heatcool.cpp \
  mpv8_StarBench_heatcool.h mp_explicit_H.h \
  microphysics_base.h ${GLOBAL}


clean:
	rm -f *.o silo2text

