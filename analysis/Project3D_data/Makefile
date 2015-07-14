# Makefile for sim_projection code, Jonathan Mackey
#
# 2010.12.13 JM: Changed aibn129 to standard workstation.
#   Updated Vpath to point to the tagged release used for MHD code.
# 2013.02.05 JM: updated to work with pion.
# 2013.10.15 JM: Updated to use microphysics to get temperature/n(H).
# 2015.07.03 JM: updated for pion_dev: uses MCMD, SimSetup, constants.h
#

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


srcdir = ${CURDIR}/../../source

bindir=.

CODEDIR=${CURDIR}/../..

VPATH=${srcdir}:\
${srcdir}/grid:\
${srcdir}/defines:\
${srcdir}/dataIO:\
${srcdir}/comms:\
${srcdir}/coord_sys:\
${srcdir}/equations:\
${srcdir}/flux_calc:\
${srcdir}/spatial_solvers:\
${srcdir}/time_integrators:\
${srcdir}/microphysics:\
${srcdir}/raytracing:\
${srcdir}/contrib\
${srcdir}/tools



GLOBAL = functionality_flags.h testing_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h


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
# This is to multi-thread the pixel calculations.
#
MT=USE_THREADS
OPTION += -DTHREADS -D__WXGTK__
#MT=NO_THREADS

#MAKE_UNAME=JUROPA
#MAKE_UNAME=standard
#MAKE_UNAME=OSX
#MAKE_UNAME=SUPERMUC

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
  CXXFLAGS = -O1 -Wall -g
  CXXFLAGS+= -Wno-unused-function


  MYLIB = ${CODEDIR}/extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${srcdir}
  LDFLAGS += -lm
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
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

  #VPATH=../../source:../../source/defines
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

#CXXFLAGS = -O3
  CXXFLAGS = -O0
#CXXFLAGS += -Wall
  CXXFLAGS += -g    # enable debugging with gdb
  #CXXFLAGS += -arch x86_64

  MYLIB = $(CURDIR)/../../extra_libraries
  #MYLIB = /Users/jmackey/Library/my_local_libs
  INC += -I$(MYLIB)/include
  INC += -I${srcdir}

  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  # libs for CVODES
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
#############################################

######################## -- SuperMUC -- ##############################
# Options for the LRZ Supercomputer SuperMUC.  Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), SUPERMUC)
  ifeq ($(IOMODE), USE_FILE_COMMS)
  CXX=icc
  endif
  ifeq ($(IOMODE), USE_MPI)
  CXX = mpiCC
  endif
  CXXFLAGS = -O1 -Wall -g #-wd981 -wd383
  #CXXFLAGE = -O1 -Wall
  CXXFLAGS += -DINTEL
  #CXXFLAGS += -DTESTING

  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  #CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check
  MYLIB = ${CODEDIR}/extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${srcdir}
  INC += -I/usr/local/sundials/2.4.0/include
  INC += ${MKL_INC} ${MPI_INC}
  LDFLAGS = -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += ${MKL_SHLIB} ${MPI_CXX_LIB}
  LDFLAGS += -lreadline
  ifeq ($(MT), USE_THREADS)
    LDFLAGS += -pthread
  endif


endif
##################### -- SuperMUC -- #################################


##################################################
##################################################

CXXFLAGS += $(OPTION)

OBJECTS = \
        constants.o \
        timer.o \
        mem_manage.o \
        reporting.o \
        command_line_interface.o \
        interpolate.o \
        MCMD_control.o \
        setup_fixed_grid.o \
        setup_fixed_grid_MPI.o \
sim_params.o \
cell_interface.o \
comm_mpi.o \
comm_files.o \
stellar_wind_BC.o \
uniform_grid.o \
uniform_grid_pllel.o \
file_status.o \
readparams.o \
dataio.o \
dataio_fits.o \
dataio_silo.o \
dataio_silo_utility.o \
dataio_silo_MPI.o \
utility_fits_class.o \
main_projection.o \
sim_projection.o \
image_io.o \
VectorOps.o \
VectorOps_spherical.o \
  \
  raytracer_base.o            \
  raytracer_SC_pllel.o        \
  raytracer_SC.o              \
  \
  microphysics_base.o         \
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
  microphysics.o              \
  microphysics_v1.o           \
  microphysics_lowZ.o         \
  mp_only_cooling.o           \
  mpv5_molecular.o            \
  mpv6_PureH.o                \
  mpv7_TwoTempIso.o           \
  mpv8_StarBench_heatcool.o

        

ifeq ($(MT), USE_THREADS)
OBJECTS += logmessages.o threadpool.o threads.o
VPATH +=:./andys_threads:./andys_threads/threadpool
THREAD_INCLUDE=./andys_threads/threads.h ./andys_threads/msvc_constants.h ./andys_threads/reefa_constants.h
endif

projection: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o projection $(LDFLAGS)


################################
# Andy's multi-threading files #
################################
logmessages.o: logmessages.cxx logmessages.h reefa_constants.h msvc_constants.h threads.h
threadpool.o: threadpool.cxx threadpool.h reefa_constants.h msvc_constants.h logmessages.h threads.h
threads.o: threads.cxx msvc_constants.h reefa_constants.h threads.h

#################################
# Stuff From my Sim Code Source #
#################################

reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
timer.o: tools/timer.cpp tools/timer.h ${GLOBAL}
command_line_interface.o: tools/command_line_interface.cpp ${GLOBAL}
constants.o: constants.cpp ${GLOBAL}
sim_params.o: 	  sim_params.cpp ${GLOBAL}

cell_interface.o: cell_interface.cc  cell_interface.h ${GLOBAL}
uniform_grid.o:   uniform_grid.cc  uniform_grid.h  ${GLOBAL}  stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cc uniform_grid.h ${GLOBAL} comms.h

VectorOps.o:  VectorOps.cc  coord_sys/VectorOps.h ${GLOBAL} uniform_grid.h
VectorOps_spherical.o:  VectorOps_spherical.cc  coord_sys/VectorOps_spherical.h  coord_sys/VectorOps.h \
		       ${GLOBAL} uniform_grid.h


file_status.o: dataIO/file_status.cpp dataIO/file_status.h
dataio.o:       dataIO/dataio.cc  dataIO/dataio.h  ${GLOBAL}  solver_eqn_base.h  dataIO/readparams.h \
		stellar_wind_BC.h dataIO/file_status.h
dataio_silo.o:  dataIO/dataio_silo.cc  dataIO/dataio_silo.h  dataIO/dataio.h ${GLOBAL}  solver_eqn_base.h  dataIO/file_status.h
dataio_fits.o:  dataIO/dataio_fits.cc  dataIO/dataio_fits.h  dataIO/dataio.h ${GLOBAL}  solver_eqn_base.h \
		solver_eqn_hydro_adi.h  solver_eqn_hydro_iso.h  dataIO/file_status.h
readparams.o:	dataIO/readparams.cc  dataIO/readparams.h  ${GLOBAL}
dataio_silo_MPI.o:    dataIO/dataio_silo_MPI.cc  dataIO/dataio_silo.h  dataIO/dataio.h ${GLOBAL}  comms.h  dataIO/file_status.h
utility_fits_class.o: dataIO/utility_fits_class.cc  dataIO/dataio_fits.h  dataIO/dataio.h ${GLOBAL}

comm_mpi.o:    comm_mpi.cc  comm_mpi.h  comms.h  ${GLOBAL}
comm_files.o:  comm_files.cc  comm_files.h  comms.h  ${GLOBAL}

dataio_silo_utility.o:  dataIO/dataio_silo_utility.cpp  dataIO/dataio_silo_utility.h  dataIO/dataio_silo.h  dataIO/dataio.h  ${GLOBAL}
stellar_wind_BC.o: stellar_wind_BC.cc  stellar_wind_BC.h ${GLOBAL}

######################
# MICROPHYSICS STUFF #
######################

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

mpv5_molecular.o: mpv5_molecular.cpp mpv5_molecular.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

mpv6_PureH.o: mpv6_PureH.cpp mpv6_PureH.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

mpv7_TwoTempIso.o: mpv7_TwoTempIso.cpp mpv7_TwoTempIso.h \
  mp_explicit_H.h microphysics_base.h ${GLOBAL}

cvode_integrator.o: cvode_integrator.cpp cvode_integrator.h ${GLOBAL}

mp_explicit_H.o: mp_explicit_H.cpp mp_explicit_H.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

mp_implicit_H.o: mp_implicit_H.cpp mp_implicit_H.h mp_explicit_H.h \
  ${GLOBAL} hydrogen_mp.h hydrogen_recomb_Hummer94.h hydrogen_photoion.h \
  cooling_SD93_cie.h microphysics_base.h

###############################
# Stuff Specific to this code #
###############################
sim_projection.o:	sim_projection.cc sim_projection.h ${GLOBAL} uniform_grid.h dataIO/dataio.h
main_projection.o:	main_projection.cc sim_projection.h ${GLOBAL} uniform_grid.h dataIO/dataio.h ${THREAD_INCLUDE} 
image_io.o:		image_io.cc image_io.h  dataIO/dataio.h  dataIO/dataio_fits.h ${GLOBAL}

clean:
	rm -f *.o projection
