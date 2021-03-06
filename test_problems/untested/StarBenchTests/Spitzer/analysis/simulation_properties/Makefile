# Makefile for HIIregion_KineticE.cpp

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cxx .o
.SUFFIXES: .cpp .o
srcdir=.
bindir=.

# make has rules for making objects from .cc files, so I don't need them
# explicitly, unless I am doing something unusual.
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<


CODEDIR=${CURDIR}/../../../../..

VPATH=${CODEDIR}/source:\
${CODEDIR}/source/defines:\
${CODEDIR}/source/ics:\
${CODEDIR}/source/microphysics:\
${CODEDIR}/source/raytracing:\
${CODEDIR}/source/grid:\
${CODEDIR}/source/dataIO:\
${CODEDIR}/source/coord_sys:\
${CODEDIR}/source/comms:\
${CODEDIR}/source/defines

GLOBAL = global.h defines/functionality_flags.h defines/testing_flags.h

# --------------------
# Compile-time options
# --------------------
# parallel code
OPT = -DPARALLEL
OPT += -DUSE_MPI
OPT += -DSILO
OPT += -DFITS
OPT += -DNEW_METALLICITY

MAKE_UNAME=standard
#MAKE_UNAME=OSX
#MAKE_UNAME=JUROPA
#


######################################################################
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  CXX = mpic++
  CXXFLAGS = -O3 -Wno-unused-function
  #CXXFLAGS = -O3 -g -Wno-unused-function
  #CXXFLAGS = -O0 -Wall -g -Wno-unused-function
  MYLIB = ${CODEDIR}/extra_libraries
  INC = -I${MYLIB}/include -I${CODEDIR}/source
  LDFLAGS = -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm  

endif
######################################################################
# Options for a "standard" workstation.
######################################################################

###########################################
# Specific Options for OS-X
###########################################
ifeq ($(MAKE_UNAME), OSX)
  CXX = mpic++
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif
  #CXXFLAGS += -Wall
  #CXXFLAGS += -g    # enable debugging with gdb
  #CXXFLAGS += -arch x86_64

  # CURDIR may only work with GNU Make!  It gives absolute paths, so that
  # the executable can always link to shared libs.
  MYLIB = ${CODEDIR}/extra_libraries
  INC = -I${MYLIB}/include -I${CODEDIR}/source
  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  # libs for CVODES
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
#################### --- OS X Darwin --- #############################


##################################################
##################################################

#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################
CXXFLAGS += $(OPT)

OBJECTS = cell_interface.o global.o get_sim_info.o \
	dataio.o readparams.o dataio_silo.o  dataio_silo_MPI.o \
        dataio_fits.o utility_fits_class.o dataio_silo_utility.o \
        uniform_grid_pllel.o uniform_grid.o \
        VectorOps.o VectorOps_spherical.o stellar_wind_BC.o  \
        cooling_SD93_cie.o hydrogen_recomb_Hummer94.o \
        hydrogen_mp.o hydrogen_photoion.o \
        mp_explicit_H.o cvode_integrator.o mp_implicit_H.o mp_v2_aifa.o        \
	comm_mpi.o comm_files.o                                                \
        integrator.o cooling.o microphysics.o microphysics_v1.o                \
        microphysics_lowZ.o mp_only_cooling.o calc_radius.o \
        file_status.o
     
ifneq (,$(findstring NEW_METALLICITY,$(OPT)))
  OBJECTS += mpv5_molecular.o mpv6_PureH.o mpv7_TwoTempIso.o mpv8_StarBench_heatcool.o
endif
#####################################################
##### END OF LIST OF OBJECT FILES               #####
#####################################################

##########################
### -- EXECTUTABLES -- ###
##########################
${bindir}/Radii : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o ${bindir}/Radii $(LDFLAGS)


##################################################
##################################################

####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################


global.o: global.cc ${GLOBAL}
cell_interface.o:  grid/cell_interface.cc grid/cell_interface.h \
    ${GLOBAL}
get_sim_info.o: ics/get_sim_info.h ics/get_sim_info.cc ${GLOBAL} \
   dataIO/readparams.h grid/stellar_wind_BC.h


file_status.o: file_status.cpp file_status.h
dataio.o: dataIO/dataio.cc dataIO/dataio.h ${GLOBAL} \
    spatial_solvers/solver_eqn_base.h dataIO/readparams.h \
    grid/stellar_wind_BC.h
dataio_silo.o:  dataIO/dataio_silo.cc dataIO/dataio_silo.h \
    dataIO/dataio.h ${GLOBAL} spatial_solvers/solver_eqn_base.h
dataio_fits.o:  dataIO/dataio_fits.cc dataIO/dataio_fits.h \
    dataIO/dataio.h ${GLOBAL} spatial_solvers/solver_eqn_base.h \
    spatial_solvers/solver_eqn_hydro_adi.h \
    spatial_solvers/solver_eqn_hydro_iso.h comms/comms.h
readparams.o:   dataIO/readparams.cc dataIO/readparams.h ${GLOBAL}
utility_fits_class.o: dataIO/utility_fits_class.cc \
    dataIO/dataio_fits.h dataIO/dataio.h ${GLOBAL}
dataio_silo_utility.o: dataIO/dataio_silo_utility.cpp dataIO/dataio_silo_utility.h \
    dataIO/dataio_silo.h dataIO/dataio.h ${GLOBAL}
dataio_silo_MPI.o: dataIO/dataio_silo_MPI.cc dataIO/dataio_silo.h \
    dataIO/dataio.h ${GLOBAL} comms/comms.h dataIO/file_status.h
VectorOps.o:  coord_sys/VectorOps.cc coord_sys/VectorOps.h \
    ${GLOBAL} grid/uniform_grid.h
VectorOps_spherical.o: coord_sys/VectorOps_spherical.cc \
    coord_sys/VectorOps_spherical.h coord_sys/VectorOps.h \
    ${GLOBAL} grid/uniform_grid.h
uniform_grid.o: 	grid/uniform_grid.cc grid/uniform_grid.h \
    ${GLOBAL} grid/stellar_wind_BC.h
uniform_grid_pllel.o: grid/uniform_grid_pllel.cc grid/uniform_grid.h \
    ${GLOBAL} comms/comms.h
stellar_wind_BC.o: grid/stellar_wind_BC.cc grid/stellar_wind_BC.h \
    ${GLOBAL}

integrator.o: microphysics/integrator.cc \
    microphysics/integrator.h ${GLOBAL}
cooling.o:    microphysics/cooling.cc \
    microphysics/cooling.h ${GLOBAL}
cooling_SD93_cie.o: microphysics/cooling_SD93_cie.cc \
    microphysics/cooling_SD93_cie.h ${GLOBAL}
hydrogen_recomb_Hummer94.o: microphysics/hydrogen_recomb_Hummer94.h \
    microphysics/hydrogen_recomb_Hummer94.cc \
    ${GLOBAL}

mp_only_cooling.o: microphysics/mp_only_cooling.cc \
    microphysics/mp_only_cooling.h ${GLOBAL} \
    microphysics/cooling_SD93_cie.h \
    microphysics/microphysics_base.h \
    microphysics/hydrogen_recomb_Hummer94.h
microphysics.o:	microphysics/microphysics.cc \
    microphysics/microphysics.h \
    microphysics/microphysics_base.h \
    ${GLOBAL} microphysics/cooling.h \
    microphysics/integrator.h
microphysics_v1.o: microphysics/microphysics_v1.cc \
    microphysics/microphysics.h \
    microphysics/microphysics_base.h ${GLOBAL} \
    microphysics/cooling.h microphysics/integrator.h
microphysics_lowZ.o: microphysics/microphysics_lowZ.cc \
    microphysics/microphysics_lowZ.h \
    microphysics/microphysics_base.h ${GLOBAL} \
    microphysics/integrator.h
hydrogen_mp.o: microphysics/hydrogen_mp.cc \
    microphysics/hydrogen_mp.h ${GLOBAL} \
    microphysics/hydrogen_recomb_Hummer94.h \
    microphysics/hydrogen_photoion.h 
hydrogen_photoion.o: microphysics/hydrogen_photoion.cc \
    microphysics/hydrogen_photoion.h ${GLOBAL}
mp_v2_aifa.o: microphysics/mp_v2_aifa.cc \
    microphysics/mp_v2_aifa.h ${GLOBAL} \
    microphysics/hydrogen_mp.h \
    microphysics/hydrogen_recomb_Hummer94.h \
    microphysics/hydrogen_photoion.h \
    microphysics/cooling_SD93_cie.h \
    microphysics/microphysics_base.h
cvode_integrator.o: microphysics/cvode_integrator.cpp \
    microphysics/cvode_integrator.h ${GLOBAL}
mp_explicit_H.o: microphysics/mp_explicit_H.cpp \
    microphysics/mp_explicit_H.h ${GLOBAL} \
    microphysics/hydrogen_mp.h \
    microphysics/hydrogen_recomb_Hummer94.h \
    microphysics/hydrogen_photoion.h \
    microphysics/cooling_SD93_cie.h \
    microphysics/microphysics_base.h
mp_implicit_H.o: microphysics/mp_implicit_H.cpp \
    microphysics/mp_implicit_H.h \
    microphysics/mp_explicit_H.h \
    ${GLOBAL} microphysics/hydrogen_mp.h \
    microphysics/hydrogen_recomb_Hummer94.h \
    microphysics/hydrogen_photoion.h \
    microphysics/cooling_SD93_cie.h \
    microphysics/microphysics_base.h

mpv5_molecular.o: microphysics/mpv5_molecular.cpp \
  microphysics/mpv5_molecular.h microphysics/mp_explicit_H.h \
  microphysics/microphysics_base.h ${GLOBAL}

mpv6_PureH.o: microphysics/mpv6_PureH.cpp microphysics/mpv6_PureH.h \
  microphysics/mp_explicit_H.h microphysics/microphysics_base.h \
  ${GLOBAL}

mpv7_TwoTempIso.o: microphysics/mpv7_TwoTempIso.cpp \
  microphysics/mpv7_TwoTempIso.h microphysics/mp_explicit_H.h \
  microphysics/microphysics_base.h ${GLOBAL}

mpv8_StarBench_heatcool.o: microphysics/mpv8_StarBench_heatcool.cpp \
  microphysics/mpv8_StarBench_heatcool.h microphysics/mp_explicit_H.h \
  microphysics/microphysics_base.h ${GLOBAL}




#######################
## MPI/files comms classes
########################
comm_mpi.o: comms/comm_mpi.cc comms/comm_mpi.h comms/comms.h \
  ${GLOBAL}
comm_files.o: comms/comm_files.cc comms/comm_files.h comms/comms.h \
  ${GLOBAL}



calc_radius.o: calc_radius.cpp ${GLOBAL}

clean:
	rm -f *.o ${bindir}/Radii

