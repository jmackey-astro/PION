#
# Makefile to convert silo outputs from pion into fits files.
#
# Author: Jonathan Mackey
# Date: 2013.04.18
# Modifications:
# - 2016.06.16 JM: updated to work with new PION code.


SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
.SUFFIXES: .cxx .o
bindir = .
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
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################

OPT = -DPARALLEL
OPT += -DUSE_MPI
OPT += -DFITS
OPT += -DSILO
CXXFLAGS = $(OPT)

MAKE_UNAME=standard
MAKE_UNAME=ubuntu18

#####################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  #CXXFLAGS += -O3 -Wno-unused-function -Wall
  CXXFLAGS += -g -Wno-unused-function
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

######################################################################
# Options for a debian 9 (stretch) workstation with system libs
######################################################################
ifeq ($(MAKE_UNAME), debian9)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  INC  = -I/usr/include
  INC += -I${srcdir}
  LDFLAGS += -lm
  LDFLAGS += -L/usr/lib
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  LDFLAGS += -lgsl
  LDFLAGS +=-Wl,--rpath -Wl,/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
################### ---  DEBIAN 9 (STRETCH)  ---######################

######################################################################
# Options for an Ubuntu 18 workstation with locally installed libraries.
######################################################################
ifeq ($(MAKE_UNAME), ubuntu18)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    #CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  MYLIB = $(CODEDIR)/extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${srcdir}
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -Wl,-rpath,${MYLIB}/lib
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  LDFLAGS += -lm
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  LDFLAGS += -lgsl -lgslcblas
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
################### --- UBUNTU 18 ---######################

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



OBJECTS = silo2fits.o   \
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
        comm_files.o \
        comm_mpi.o         \
        dataio.o dataio_fits.o utility_fits_class.o                 \
        file_status.o readparams.o                                  \
        dataio_silo.o dataio_silo_utility.o dataio_silo_MPI.o       \
        cell_interface.o      \
        uniform_grid_pllel.o uniform_grid.o                         \
        VectorOps.o VectorOps_spherical.o stellar_wind_BC.o         \
	stellar_wind_angle.o \
  microphysics_base.o         \
microphysics.o \
microphysics_v1.o \
integrator.o \
cooling_SD93_cie.o \
cooling.o \
hydrogen_recomb_Hummer94.o \
hydrogen_mp.o \
hydrogen_photoion.o \
mp_explicit_H.o \
cvode_integrator.o \
mp_implicit_H.o \
mp_v2_aifa.o \
microphysics_lowZ.o \
mp_only_cooling.o \
mpv5_molecular.o \
mpv6_PureH.o \
mpv7_TwoTempIso.o   \
  mpv8_StarBench_heatcool.o   \
  raytracer_base.o            \
  raytracer_SC_pllel.o        \
  raytracer_SC.o


silo2fits: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o SILO2FITS $(LDFLAGS)


silo2fits.o: silo2fits.cpp ${GLOBAL} dataio.h


cell_interface.o: cell_interface.cpp  cell_interface.h ${GLOBAL}
uniform_grid.o:   uniform_grid.cpp  uniform_grid.h  ${GLOBAL} \
    stellar_wind_BC.h stellar_wind_angle.h  \
    VectorOps.h VectorOps_spherical.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} \
     comms.h stellar_wind_angle.h stellar_wind_BC.h \
    VectorOps.h VectorOps_spherical.h
stellar_wind_BC.o: stellar_wind_BC.cpp  stellar_wind_BC.h ${GLOBAL}
stellar_wind_angle.o: stellar_wind_angle.cpp stellar_wind_angle.h stellar_wind_BC.h \
  tools/interpolate.h ${GLOBAL}

VectorOps.o:  VectorOps.cpp  coord_sys/VectorOps.h ${GLOBAL} uniform_grid.h
VectorOps_spherical.o:  VectorOps_spherical.cpp  coord_sys/VectorOps_spherical.h  coord_sys/VectorOps.h \
		       ${GLOBAL} uniform_grid.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h
readparams.o: readparams.cpp readparams.h ${GLOBAL}
dataio.o: dataio.cpp dataio.h ${GLOBAL}
dataio_fits.o: dataio_fits.cpp dataio_fits.h dataio.h ${GLOBAL}
utility_fits_class.o: utility_fits_class.cpp dataio_fits.h dataio.h \
    ${GLOBAL} file_status.h
dataio_silo_MPI.o:    dataio_silo_MPI.cpp dataio_silo.h dataio.h \
    ${GLOBAL} comms.h file_status.h

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
	rm -f *.o SILO2FITS
