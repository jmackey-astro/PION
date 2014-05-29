#
# Makefile to convert silo outputs from pion into fits files.
#
# Author: Jonathan Mackey
# Date: 2013.04.18
# Modifications:


SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
.SUFFIXES: .cxx .o
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
../../source/contrib

GLOBAL = global.h functionality_flags.h testing_flags.h

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

#####################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  #CXXFLAGS += -O3 -Wno-unused-function
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




OBJECTS = silo2fits.o global.o                                      \
        comm_files.o comm_mpi.o                                     \
        dataio.o dataio_fits.o utility_fits_class.o                 \
        file_status.o readparams.o                                  \
        dataio_silo.o dataio_silo_utility.o dataio_silo_MPI.o       \
        cell_interface.o                                            \
        uniform_grid_pllel.o uniform_grid.o                         \
        VectorOps.o VectorOps_spherical.o stellar_wind_BC.o         \
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
mpv7_TwoTempIso.o 


silo2fits: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o SILO2FITS $(LDFLAGS)


silo2fits.o: silo2fits.cpp ${GLOBAL} dataio.h
global.o: global.cc ${GLOBAL}
cell_interface.o: cell_interface.cc cell_interface.h ${GLOBAL}

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h
readparams.o: readparams.cc readparams.h ${GLOBAL}
dataio.o: dataio.cc dataio.h ${GLOBAL}
dataio_fits.o: dataio_fits.cc dataio_fits.h dataio.h ${GLOBAL}
utility_fits_class.o: utility_fits_class.cc dataio_fits.h dataio.h \
    ${GLOBAL} file_status.h
dataio_silo_MPI.o:    dataio_silo_MPI.cc dataio_silo.h dataio.h \
    ${GLOBAL} comms.h file_status.h

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

#######################
## MPI/files comms classes
########################
comm_files.o: comm_files.cc comm_files.h comms.h ${GLOBAL}
comm_mpi.o: comm_mpi.cc comm_mpi.h comms.h ${GLOBAL}

clean:
	rm -f *.o SILO2FITS
