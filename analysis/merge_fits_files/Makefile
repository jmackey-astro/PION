#
# Makefile for code to stitch together a bunch of fits outputs from
# many subdomains into a single large file.
#
# Author: Jonathan Mackey
# Date: Sometime in 2007
# Modifications:
# - 2013.04.16 JM: Updated for new pion source code.


SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cpp .o
.SUFFIXES: .cpp .o
bindir = ..
SRCDIR=${CURDIR}/../../source

VPATH = \
../../source:\
../../source/ics:\
../../source/Riemann_solvers:\
../../source/defines:\
../../source/microphysics:\
../../source/raytracing:\
../../source/grid:\
../../source/dataIO:\
../../source/comms:\
../../source/coord_sys:\
../../source/equations:\
../../source/flux_calc:\
../../source/spatial_solvers:\
../../source/time_integrators:\
../../source/contrib\
../../source/tools


GLOBAL = functionality_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h


#######################
## COMPILATION RULES ##
#######################
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################

OPT = -DPARALLEL
OPT += -DUSE_MPI
OPT += -DFITS
CXXFLAGS = $(OPT)

#MAKE_UNAME=standard

#####################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  CXXFLAGS += -g -O3 -Wno-unused-function
  CXX=mpicxx

  MYLIB = $(CURDIR)/../../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -lm
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
################### --- STANDARD WORKSTATION ---######################


###########################################
# Options for OS X 
###########################################
ifeq ($(MAKE_UNAME), OSX)
  CXX = mpicxx

#CXXFLAGS = -O3
  CXXFLAGS += -O0
#CXXFLAGS += -Wall
  CXXFLAGS += -g    # enable debugging with gdb
  #CXXFLAGS += -arch x86_64

  MYLIB = $(CURDIR)/../../extra_libraries
  #MYLIB = /Users/jmackey/Library/my_local_libs
  INC += -I$(MYLIB)/include
  INC += -I${SRCDIR}

  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  #LDFLAGS += -lsilo 
  # libs for CVODES
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
#############################################



OBJECTS = merge_fits_files.o \
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
  uniform_grid_pllel.o        \
  readparams.o                \
  comm_mpi.o                  \
  comm_files.o                \
  VectorOps.o                 \
  VectorOps_spherical.o       \
  stellar_wind_BC.o           \
  cell_interface.o            \
  file_status.o               \
  dataio.o                    \
  dataio_fits.o utility_fits_class.o   \
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


########################
merge_fits_files: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o merge_fits_files $(LDFLAGS)
########################

########################
merge_fits_files.o: merge_fits_files.cpp ${GLOBAL} dataio.h
global.o: global.cpp ${GLOBAL}
cell_interface.o: cell_interface.cpp cell_interface.h ${GLOBAL}
########################

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h
readparams.o: readparams.cpp readparams.h ${GLOBAL}
dataio.o: dataio.cpp dataio.h ${GLOBAL} solver_eqn_base.h \
    readparams.h stellar_wind_BC.h file_status.h
dataio_fits.o: dataio_fits.cpp dataio_fits.h dataio.h ${GLOBAL}
utility_fits_class.o: utility_fits_class.cpp dataio_fits.h dataio.h \
    ${GLOBAL} file_status.h


uniform_grid.o: uniform_grid.cpp uniform_grid.h ${GLOBAL} stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} comms.h
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

#######################
## MPI/files comms classes
########################
comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}
comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h ${GLOBAL}

clean:
	rm -f *.o merge_fits_files
