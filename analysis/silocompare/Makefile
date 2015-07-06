# Makefile for silocompare code, Jonathan Mackey
#
# 2010.10.13 JM: Added new stellar_wind_BC.o object and new standard compile.
# 2010.12.07 JM: Added geometric grid support, and also cell_interface object.
# 2011.04.23 JM: now works on Juropa.
# 2013.09.05 JM: Updated for pion.
# 2015.03.26 JM: updated for pion v0.2.

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
srcdir = .
SRCDIR=${CURDIR}/../../source
GLOBAL = global.h functionality_flags.h testing_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h

#######################
## COMPILATION RULES ##
#######################
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################

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
#MAKE_UNAME=aibn129
#MAKE_UNAME=epona
#MAKE_UNAME=stokes

# All machines and options need readline and maths
LDFLAGS = -lreadline
# display all warnings
CXXFLAGS = -O3 
CXXFLAGS += -g -Wall

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

  CXXFLAGS = -O3 -Wno-unused-function
  #CXXFLAGS = -O3 -g
  #CXXFLAGS = -O0 -Wall -g
  #CXXFLAGS+= -Wno-unused-function

  MYLIB = $(CURDIR)/../../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}

  LDFLAGS =  -lm  -lreadline
  # Libraries for fits, silo
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################

###########################################
# Specific Options for OS-X
###########################################
ifeq ($(MAKE_UNAME), osx)
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  CXXFLAGS = -O3 -Wno-unused-function
  #CXXFLAGS = -Wall -g -Wno-unused-function

  # CURDIR may only work with GNU Make!  It gives absolute paths, so that
  # the executable can always link to shared libs.
  #MYLIB = /Users/jmackey/Library/JMs_stuff
  MYLIB = $(CURDIR)/../../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  # libs for CVODES
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
#################### --- OS X Darwin --- #############################

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
VPATH = \
../../source:\
../../source/ics:\
../../source/Riemann_solvers:\
../../source/defines:\
../../source/microphysics:\
../../source/raytracing:\
../../source/microphysics/Harpreets_code:\
../../source/grid:\
../../source/dataIO:\
../../source/comms:\
../../source/coord_sys:\
../../source/equations:\
../../source/flux_calc:\
../../source/spatial_solvers:\
../../source/time_integrators:\
../../source/tools


#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################

OBJECTS = file_status.o  \
        constants.o \
        timer.o \
        mem_manage.o \
        reporting.o \
        interpolate.o \
        MCMD_control.o \
        setup_fixed_grid.o \
        setup_fixed_grid_MPI.o \
        silocompare.o \
        readparams.o \
        global.o \
        dataio.o \
	dataio_silo.o \
        dataio_silo_MPI.o \
        dataio_silo_utility.o \
        cell_interface.o \
	stellar_wind_BC.o \
        uniform_grid.o \
        uniform_grid_pllel.o \
        comm_mpi.o \
        comm_files.o \
        VectorOps.o \
        VectorOps_spherical.o \
        microphysics_base.o \
	integrator.o \
        cooling.o \
        microphysics.o \
        microphysics_v1.o \
	mp_only_cooling.o \
	cooling_SD93_cie.o \
	hydrogen_recomb_Hummer94.o  \
        hydrogen_mp.o \
	hydrogen_photoion.o \
        mp_explicit_H.o \
	cvode_integrator.o \
	mp_implicit_H.o \
	mp_v2_aifa.o \
        raytracer_base.o  \
	raytracer_SC_pllel.o \
        raytracer_SC.o


OBJECTS += mpv5_molecular.o mpv6_PureH.o mpv7_TwoTempIso.o mpv8_StarBench_heatcool.o


##########################
### -- EXECTUTABLES -- ###
##########################

silocompare: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o silocompare $(LDFLAGS)



####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################

silocompare.o:	silocompare.cc global.h uniform_grid.h dataio.h ${GLOBAL}

reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
constants.o: constants.cpp ${GLOBAL}
timer.o: tools/timer.cpp tools/timer.h ${GLOBAL}
global.o: 	  global.cc  ${GLOBAL}

comm_mpi.o: comm_mpi.cc comm_mpi.h comms.h
comm_files.o: comm_files.cc comm_files.h comms.h ${GLOBAL}

uniform_grid.o: uniform_grid.cc uniform_grid.h ${GLOBAL} stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cc uniform_grid.h ${GLOBAL} comms.h
stellar_wind_BC.o: stellar_wind_BC.cc stellar_wind_BC.h ${GLOBAL}

VectorOps.o:           VectorOps.cc VectorOps.h uniform_grid.h ${GLOBAL}
VectorOps_spherical.o: VectorOps_spherical.cc VectorOps_spherical.h VectorOps.h \
		       uniform_grid.h ${GLOBAL}

MCMD_control.o: MCMD_control.cpp MCMD_control.h ${GLOBAL}

setup_fixed_grid.o: setup_fixed_grid.cpp setup_fixed_grid.h ${GLOBAL} \
    grid_base_class.h uniform_grid.h

setup_fixed_grid_MPI.o: setup_fixed_grid_MPI.cpp setup_fixed_grid_MPI.h \
    setup_fixed_grid.h ${GLOBAL} \
    grid_base_class.h uniform_grid.h

######################
# MICROPHYSICS STUFF #
######################
microphysics_base.o: microphysics_base.cpp microphysics_base.h ${GLOBAL}

integrator.o: integrator.cc integrator.h ${GLOBAL}
cooling.o:    cooling.cc cooling.h ${GLOBAL}

cooling_SD93_cie.o: cooling_SD93_cie.cc cooling_SD93_cie.h ${GLOBAL}

hydrogen_recomb_Hummer94.o: hydrogen_recomb_Hummer94.h hydrogen_recomb_Hummer94.cc \
    ${GLOBAL}

mp_only_cooling.o: mp_only_cooling.cc mp_only_cooling.h ${GLOBAL} \
   cooling_SD93_cie.h microphysics_base.h hydrogen_recomb_Hummer94.h

microphysics.o:	microphysics.cc microphysics.h microphysics_base.h \
   ${GLOBAL} cooling.h integrator.h
microphysics_v1.o: microphysics_v1.cc microphysics.h \
   microphysics_base.h ${GLOBAL} cooling.h integrator.h
microphysics_lowZ.o: microphysics_lowZ.cc microphysics_lowZ.h \
   microphysics_base.h ${GLOBAL} integrator.h

hydrogen_mp.o: hydrogen_mp.cc hydrogen_mp.h ${GLOBAL} hydrogen_recomb_Hummer94.h \
  hydrogen_photoion.h 

hydrogen_photoion.o: hydrogen_photoion.cc hydrogen_photoion.h ${GLOBAL}

mp_v2_aifa.o: mp_v2_aifa.cc mp_v2_aifa.h ${GLOBAL} hydrogen_mp.h \
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
  mpv8_StarBench_heatcool.h mp_explicit_H.h microphysics_base.h ${GLOBAL}


########################
# Raytracing
########################
raytracer_base.o: raytracer_base.cpp raytracer_base.h ${GLOBAL}

raytracer_SC.o: raytracer_SC.cc raytracer_SC.h ${GLOBAL} \
    uniform_grid.h raytracer_base.h

raytracer_SC_pllel.o : raytracer_SC_pllel.cc raytracer_SC.h \
    ${GLOBAL} uniform_grid.h comms.h raytracer_base.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h

dataio.o: dataio.cc dataio.h ${GLOBAL} solver_eqn_base.h \
    readparams.h stellar_wind_BC.h file_status.h

dataio_silo.o:  dataio_silo.cc dataio_silo.h dataio.h ${GLOBAL} \
    solver_eqn_base.h file_status.h

dataio_silo_MPI.o:    dataio_silo_MPI.cc dataio_silo.h dataio.h \
    ${GLOBAL} comms.h file_status.h

readparams.o:	readparams.cc readparams.h ${GLOBAL}

dataio_silo_MPI.o:    dataio_silo_MPI.cc dataio_silo.h dataio.h \
    ${GLOBAL} comms.h file_status.h

dataio_silo_utility.o: dataio_silo_utility.cpp dataio_silo_utility.h\
    ${GLOBAL} file_status.h dataio.h dataio_silo.h comms.h


clean:
	rm -f *.o silocompare

