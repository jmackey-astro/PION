# Makefile for projecting 2D axisymmetric data onto the plane of sky code, Jonathan Mackey
#
# 2012.01.09 JM Modified from pllel_analyse_data
# 2012.01.10 JM Debugged code and modified makefile a little.
# 2012.11.23 JM: Updated to use new code branch.

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cxx .o
.SUFFIXES: .cpp .o
srcdir = .

#######################
## COMPILATION RULES ##
#######################
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################


####### PATHS TO SOURCES #######
PION=${CURDIR}/../..
srcdir=${PION}/source
bindir=.

VPATH = \
${srcdir}:\
${srcdir}/ics:\
${srcdir}/Riemann_solvers:\
${srcdir}/defines:\
${srcdir}/microphysics:\
${srcdir}/raytracing:\
${srcdir}/grid:\
${srcdir}/dataIO:\
${srcdir}/comms:\
${srcdir}/coord_sys:\
${srcdir}/equations:\
${srcdir}/flux_calc:\
${srcdir}/spatial_solvers:\
${srcdir}/time_integrators:\
${srcdir}/contrib\
${srcdir}/tools



GLOBAL = functionality_flags.h testing_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h

####### PATHS TO SOURCES #######


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
ifeq ($(MT), USE_THREADS)
  OPTION += -DTHREADS -D__WXGTK__
endif

# This variable includes a #define for ABSORPTION, if requested
OPTION += ${ABS}

# All machines and options need readline and maths
LDFLAGS  = -lreadline
LDFLAGS += -lm
ifeq ($(MT), USE_THREADS)
  LDFLAGS += -pthread
endif

######################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  CXX=mpicxx
  #ifeq ($(PION_OPTIMISE), HIGH)
  CXXFLAGS = -O3 -Wno-unused-function
  #else
  #CXXFLAGS = -Wall -g -Wno-unused-function
  #endif

  MYLIB = ${PION}/extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${srcdir}
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
################### --- STANDARD WORKSTATION ---######################

######################################################################
# Options for an OS-X workstation.
######################################################################
ifeq ($(MAKE_UNAME), OSX)
  CXX = mpicxx
  #CXXFLAGS = -O3 -Wno-unused-function
  #CXXFLAGS = -O3 -g -Wno-unused-function
  CXXFLAGS = -O0 -Wall -g -Wno-unused-function

# Hardcode for my mbp2015
  MYLIB = ${PION}/extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${srcdir}

  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif

####################   OS X Options ############################


#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################

CXXFLAGS += $(OPTION)

OBJECTS = \
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
  uniform_grid.o             \
  readparams.o               \
  file_status.o              \
  dataio.o                   \
  dataio_silo_utility.o      \
  dataio_silo.o              \
  dataio_silo_MPI.o          \
  dataio_fits.o              \
  uniform_grid_pllel.o       \
  utility_fits_class.o       \
  VectorOps.o                \
  VectorOps_spherical.o      \
  stellar_wind_BC.o          \
  image_io.o                 \
  cell_interface.o           \
  raytracer_base.o            \
  raytracer_SC_pllel.o        \
  raytracer_SC.o              \
  microphysics_base.o         \
  cooling_SD93_cie.o         \
  hydrogen_recomb_Hummer94.o \
  hydrogen_mp.o              \
  hydrogen_photoion.o        \
  mp_explicit_H.o            \
  cvode_integrator.o         \
  mp_implicit_H.o            \
  mp_v2_aifa.o               \
  integrator.o               \
  cooling.o                  \
  microphysics.o             \
  microphysics_v1.o          \
  microphysics_lowZ.o        \
  mp_only_cooling.o          \
  mpv5_molecular.o           \
  mpv6_PureH.o               \
  mpv7_TwoTempIso.o          \
  mpv8_StarBench_heatcool.o  \
  xray_emission.o            \
  emission_absorption.o

ifeq ($(IOMODE), USE_FILE_COMMS)
  OBJECTS += comm_files.o 
else
ifeq ($(IOMODE), USE_MPI)
  OBJECTS += comm_mpi.o
endif
endif


ifeq ($(MT), USE_THREADS)
  OBJECTS += logmessages.o threadpool.o threads.o
  OBJECTS += project_axi2Ddata_thread.o
  VPATH +=:${srcdir}/tools/threads_AJL:${srcdir}/tools/threads_AJL/threadpool
  THREAD_INCLUDE=${srcdir}/tools/threads_AJL/threads.h ${srcdir}/tools/threads_AJL/msvc_constants.h ${srcdir}/tools/threads_AJL/reefa_constants.h
else
  OBJECTS += project_axi2Ddata.o
endif


##########################
### -- EXECTUTABLES -- ###
##########################
project_data2D: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o project_data2D $(LDFLAGS)

####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################

#################################
# Stuff From PION Source #
#################################
reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
timer.o: tools/timer.cpp tools/timer.h ${GLOBAL}
command_line_interface.o: tools/command_line_interface.cpp ${GLOBAL}
constants.o: constants.cpp ${GLOBAL}
sim_params.o: 	  sim_params.cpp ${GLOBAL}
interpolate.o: tools/interpolate.h tools/interpolate.cpp ${GLOBAL}

cell_interface.o: cell_interface.cpp  cell_interface.h ${GLOBAL}
uniform_grid.o:   uniform_grid.cpp  uniform_grid.h  ${GLOBAL}  stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} comms.h

VectorOps.o:  VectorOps.cpp  coord_sys/VectorOps.h ${GLOBAL} uniform_grid.h
VectorOps_spherical.o:  VectorOps_spherical.cpp  coord_sys/VectorOps_spherical.h \
    coord_sys/VectorOps.h ${GLOBAL} uniform_grid.h


file_status.o: dataIO/file_status.cpp dataIO/file_status.h
dataio.o:       dataIO/dataio.cpp  dataIO/dataio.h  ${GLOBAL}  solver_eqn_base.h  dataIO/readparams.h \
		stellar_wind_BC.h dataIO/file_status.h
dataio_silo.o:  dataIO/dataio_silo.cpp  dataIO/dataio_silo.h  dataIO/dataio.h ${GLOBAL}  solver_eqn_base.h  dataIO/file_status.h
dataio_fits.o:  dataIO/dataio_fits.cpp  dataIO/dataio_fits.h  dataIO/dataio.h ${GLOBAL}  solver_eqn_base.h \
		solver_eqn_hydro_adi.h  solver_eqn_hydro_iso.h  dataIO/file_status.h
readparams.o:	dataIO/readparams.cpp  dataIO/readparams.h  ${GLOBAL}
dataio_silo_MPI.o:    dataIO/dataio_silo_MPI.cpp  dataIO/dataio_silo.h  dataIO/dataio.h ${GLOBAL}  comms.h  dataIO/file_status.h
utility_fits_class.o: dataIO/utility_fits_class.cpp  dataIO/dataio_fits.h  dataIO/dataio.h ${GLOBAL}

comm_mpi.o:    comm_mpi.cpp  comm_mpi.h  comms.h  ${GLOBAL}
comm_files.o:  comm_files.cpp  comm_files.h  comms.h  ${GLOBAL}

dataio_silo_utility.o:  dataIO/dataio_silo_utility.cpp  dataIO/dataio_silo_utility.h  dataIO/dataio_silo.h  dataIO/dataio.h  ${GLOBAL}
stellar_wind_BC.o: stellar_wind_BC.cpp  stellar_wind_BC.h ${GLOBAL}

################################
# Andy's multi-threading files #
################################
logmessages.o: logmessages.cxx logmessages.h reefa_constants.h \
  msvc_constants.h threads.h
threadpool.o: threadpool.cxx threadpool.h reefa_constants.h \
  msvc_constants.h logmessages.h threads.h
threads.o: threads.cxx msvc_constants.h reefa_constants.h threads.h

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
project_axi2Ddata.o: project_axi2Ddata.cpp image_io.h ${GLOBAL}
image_io.o: image_io.cpp image_io.h ${GLOBAL}

clean:
	rm -f *.o project_data2D
