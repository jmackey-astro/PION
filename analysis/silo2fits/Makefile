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
CODEDIR=${CURDIR}/../..

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
../../source/tools:\
../../source/decomposition:\
../../source/Riemann_solvers:\
../../source/spatial_solvers:\
../../source/boundaries


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
MAKE_UNAME=debian9

IOMODE = USE_MPI

######################################################################
# Options for a debian 9 (stretch) workstation with system libs
######################################################################
ifeq ($(MAKE_UNAME), debian9)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS += -O3 -Wno-unused-function
  else
    CXXFLAGS += -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  CXXFLAGS += -DCVODE2
  INC  = -I/usr/include
  INC += -I${SRCDIR}
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
    CXXFLAGS += -O3 -Wno-unused-function
    #CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS += -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  CXXFLAGS += -DCVODE5
  MYLIB = $(CODEDIR)/extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
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

OBJ_CHEM=integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
  hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
  microphysics_base.o mp_only_cooling.o \
  MPv3.o MPv5.o MPv6.o MPv7.o MPv10.o photo_xsections.o

OBJ_IO =  parameter_defs.o file_status.o readparams.o \
  dataio_silo.o dataio_fits.o utility_fits_class.o \
  dataio_base.o dataio_text.o dataio_fits_MPI.o \
  dataio_silo_MPI.o dataio_silo_utility.o

OBJ_GRID =  cell_interface.o uniform_grid.o setup_fixed_grid.o \
  stellar_wind_BC.o stellar_wind_angle.o MCMD_control.o \
  raytracer_SC.o setup_fixed_grid_MPI.o \
  uniform_grid_pllel.o raytracer_SC_pllel.o \
  setup_grid_NG_MPI.o setup_NG_grid.o

OBJ_BCS = assign_update_bcs.o periodic_boundaries.o \
  outflow_boundaries.o inflow_boundaries.o \
  reflecting_boundaries.o oneway_out_boundaries.o \
  fixed_boundaries.o jet_boundaries.o \
  jetreflect_boundaries.o double_Mach_ref_boundaries.o \
  stellar_wind_boundaries.o MCMD_boundaries.o \
  periodic_boundaries_MPI.o assign_update_bcs_MPI.o \
  RT_MPI_boundaries.o

OBJ_NG_BCS = assign_update_bcs_NG_MPI.o assign_update_bcs_NG.o \
  NG_coarse_to_fine_boundaries.o \
  NG_fine_to_coarse_boundaries.o \
  NG_MPI_coarse_to_fine_boundaries.o \
  NG_MPI_fine_to_coarse_boundaries.o NG_BC89flux.o \
  NG_MPI_BC89flux.o

OBJ_EQNS =  eqns_base.o eqns_hydro_adiabatic.o eqns_mhd_adiabatic.o 

OBJ_SOLVE = riemannMHD.o HLLD_MHD.o Roe_MHD_ConservedVar_solver.o \
  riemann.o findroot.o Riemann_FVS_hydro.o HLL_hydro.o \
  Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
  solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o

OBJECTS = silo2fits.o   \
  constants.o \
  timer.o \
  mem_manage.o \
  reporting.o \
  command_line_interface.o \
  interpolate.o \
  sim_params.o \
        comm_files.o \
        comm_mpi.o         \
        VectorOps.o VectorOps_spherical.o

OBJECTS += $(OBJ_CHEM) $(OBJ_IO) $(OBJ_GRID) $(OBJ_BCS) $(OBJ_NG_BCS) \
   $(OBJ_EQNS) $(OBJ_SOLVE)

silo2fits: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o SILO2FITS $(LDFLAGS)


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################

silo2fits.o: silo2fits.cpp ${GLOBAL} dataio_base.h


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
VectorOps_spherical.o:  VectorOps_spherical.cpp  coord_sys/VectorOps_spherical.h \
      coord_sys/VectorOps.h ${GLOBAL} uniform_grid.h

MCMD_control.o: MCMD_control.cpp MCMD_control.h ${GLOBAL}


###################
# RIEMANN SOLVERS #
###################

riemann.o: riemann.cpp riemann.h findroot.h ${GLOBAL} \
 eqns_hydro_adiabatic.h  eqns_base.h

findroot.o: findroot.cpp findroot.h ${GLOBAL}

HLL_hydro.o: HLL_hydro.cpp HLL_hydro.h \
 ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

Riemann_FVS_hydro.o: Riemann_FVS_hydro.cpp Riemann_FVS_hydro.h \
 ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

Roe_Hydro_ConservedVar_solver.o: Roe_Hydro_ConservedVar_solver.cpp \
 Roe_Hydro_ConservedVar_solver.h \
 ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

Roe_Hydro_PrimitiveVar_solver.o: Roe_Hydro_PrimitiveVar_solver.cpp \
   Roe_Hydro_PrimitiveVar_solver.h \
   ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

riemannMHD.o: riemannMHD.cpp riemannMHD.h ${GLOBAL} \
 eqns_mhd_adiabatic.h eqns_base.h 

Roe_MHD_ConservedVar_solver.o : Roe_MHD_ConservedVar_solver.cpp \
   Roe_MHD_ConservedVar_solver.h riemannMHD.h \
   ${GLOBAL} eqns_mhd_adiabatic.h eqns_base.h

HLLD_MHD.o : HLLD_MHD.cpp HLLD_MHD.h  \
   ${GLOBAL} eqns_mhd_adiabatic.h eqns_base.h

#####################
# BOUNDARY CLASSES  #
#####################

periodic_boundaries.o: periodic_boundaries.cpp periodic_boundaries.h\
  boundaries.h ${GLOBAL} grid_base_class.h

outflow_boundaries.o: outflow_boundaries.cpp outflow_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

oneway_out_boundaries.o: oneway_out_boundaries.cpp \
  oneway_out_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h outflow_boundaries.h

inflow_boundaries.o: inflow_boundaries.cpp inflow_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

reflecting_boundaries.o: reflecting_boundaries.cpp \
  reflecting_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

fixed_boundaries.o: fixed_boundaries.cpp fixed_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

jet_boundaries.o: jet_boundaries.cpp jet_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

jetreflect_boundaries.o: jetreflect_boundaries.cpp \
  jetreflect_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

double_Mach_ref_boundaries.o: double_Mach_ref_boundaries.cpp \
  double_Mach_ref_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h

stellar_wind_boundaries.o: stellar_wind_boundaries.cpp \
  stellar_wind_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h stellar_wind_BC.h

stellar_wind_angle.o: stellar_wind_angle.cpp stellar_wind_angle.h \
  ${GLOBAL} cell_interface.h grid_base_class.h 

stellar_wind_BC.o: stellar_wind_BC.cpp stellar_wind_BC.h \
  ${GLOBAL} cell_interface.h grid_base_class.h

assign_update_bcs.o: assign_update_bcs.cpp assign_update_bcs.h \
  periodic_boundaries.h outflow_boundaries.h \
  oneway_out_boundaries.h inflow_boundaries.h \
  reflecting_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

assign_update_bcs_NG.o: assign_update_bcs_NG.cpp \
  assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h solver_eqn_base.h

NG_coarse_to_fine_boundaries.o: NG_coarse_to_fine_boundaries.h \
  NG_coarse_to_fine_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h

NG_fine_to_coarse_boundaries.o: NG_fine_to_coarse_boundaries.h \
  NG_fine_to_coarse_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h

NG_MPI_coarse_to_fine_boundaries.o: \
  NG_MPI_coarse_to_fine_boundaries.h NG_coarse_to_fine_boundaries.h \
  NG_MPI_coarse_to_fine_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h MCMD_control.h

NG_MPI_fine_to_coarse_boundaries.o: \
  NG_MPI_fine_to_coarse_boundaries.h NG_fine_to_coarse_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h MCMD_control.h

assign_update_bcs_MPI.o: assign_update_bcs_MPI.cpp \
  assign_update_bcs_MPI.h periodic_boundaries_MPI.h \
  MCMD_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h solver_eqn_base.h MCMD_control.h 

assign_update_bcs_NG_MPI.o: assign_update_bcs_NG_MPI.cpp \
  assign_update_bcs_NG_MPI.h assign_update_bcs_MPI.h \
  assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.h \
  NG_MPI_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h periodic_boundaries_MPI.h \
  MCMD_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h solver_eqn_base.h MCMD_control.h 

periodic_boundaries_MPI.o: periodic_boundaries_MPI.cpp \
  periodic_boundaries_MPI.h boundaries.h ${GLOBAL} \
  grid_base_class.h MCMD_control.h MCMD_boundaries.h

MCMD_boundaries.o: MCMD_boundaries.cpp \
  MCMD_boundaries.h boundaries.h ${GLOBAL} \
  grid_base_class.h MCMD_control.h

RT_MPI_boundaries.o: RT_MPI_boundaries.cpp RT_MPI_boundaries.h \
  boundaries.h ${GLOBAL} \
  grid_base_class.h MCMD_control.h

NG_BC89flux.o: NG_BC89flux.cpp NG_BC89flux.h ${GLOBAL} \
  grid_base_class.h uniform_grid.h solver_eqn_base.h

NG_MPI_BC89flux.o: NG_MPI_BC89flux.cpp NG_MPI_BC89flux.h \
  NG_BC89flux.h ${GLOBAL} solver_eqn_base.h \
  grid_base_class.h uniform_grid.h MCMD_control.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h ${GLOBAL}

parameter_defs.o:  parameter_defs.cpp parameter_defs.h ${GLOBAL} \
    solver_eqn_base.h readparams.h \
    stellar_wind_angle.h stellar_wind_BC.h cell_interface.h

dataio_base.o:  dataio_base.cpp dataio_base.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h

readparams.o: readparams.cpp readparams.h ${GLOBAL}

dataio_text.o:  dataio_text.cpp dataio_text.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h grid_base_class.h dataio_base.h


dataio_silo.o:  dataio_silo.cpp dataio_silo.h ${GLOBAL} \
    parameter_defs.h  grid_base_class.h dataio_base.h \
    solver_eqn_base.h file_status.h

dataio_fits.o:  dataio_fits.cpp dataio_fits.h ${GLOBAL} \
    solver_eqn_base.h solver_eqn_hydro_adi.h  \
    parameter_defs.h  grid_base_class.h \
    file_status.h dataio_base.h

dataio_fits_MPI.o:  dataio_fits_MPI.cpp dataio_fits_MPI.h \
    dataio_fits.h ${GLOBAL} dataio_base.h \
    solver_eqn_base.h solver_eqn_hydro_adi.h  \
    parameter_defs.h  grid_base_class.h \
    file_status.h uniform_grid_pllel.h

readparams.o:	readparams.cpp readparams.h ${GLOBAL}

dataio_silo_MPI.o:    dataio_silo_MPI.cpp dataio_silo.h dataio_base.h \
    parameter_defs.h  grid_base_class.h solver_eqn_base.h \
    ${GLOBAL} comms.h file_status.h uniform_grid_pllel.h

dataio_silo_utility.o:    dataio_silo_utility.cpp dataio_silo_utility.h \
    dataio_silo.h dataio_base.h solver_eqn_base.h \
    parameter_defs.h  grid_base_class.h \
    ${GLOBAL} comms.h file_status.h uniform_grid_pllel.h

utility_fits_class.o: utility_fits_class.cpp dataio_fits.h dataio_base.h \
    parameter_defs.h  grid_base_class.h \
    ${GLOBAL} file_status.h uniform_grid_pllel.h


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

raytracer_SC.o: raytracer_SC.cpp raytracer_SC.h ${GLOBAL} \
    uniform_grid.h raytracer_base.h

raytracer_SC_pllel.o : raytracer_SC_pllel.cpp raytracer_SC.h \
    ${GLOBAL} uniform_grid.h comms.h raytracer_base.h

#######################
## MPI/files comms classes
########################
comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}
comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h ${GLOBAL}

clean:
	rm -f *.o SILO2FITS
