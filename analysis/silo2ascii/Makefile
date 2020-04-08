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
${SRCDIR}/Riemann_solvers:\
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
${SRCDIR}/tools:\
${SRCDIR}/decomposition:\
${SRCDIR}/boundaries



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
MAKE_UNAME=debian10

# All machines and options need readline and maths
LDFLAGS = -lreadline
# display all warnings
CXXFLAGS = -O3 
#CXXFLAGS += -Wall

#######################################################################
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
# Options for a debian 10 (buster) workstation with system libs
######################################################################
ifeq ($(MAKE_UNAME), debian10)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS += -O3 -Wno-unused-function
    CXXFLAGS += -Wall -g -Wno-unused-function
  else
    CXXFLAGS += -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  CXXFLAGS += -DCVODE3
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
################### ---  DEBIAN 10 (buster)  ---######################



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

  OPT += -DCVODE5
  MYLIB = $(CURDIR)/../../extra_libraries
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

######################################################################
# Options for an Ubuntu 16 workstation with locally installed libraries.
######################################################################
ifeq ($(MAKE_UNAME), ubuntu16)
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

  OPT += -DCVODE5
  MYLIB = $(CURDIR)/../../extra_libraries
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
################### --- UBUNTU 16 ---######################



###########################################
# Specific Options for OS X
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
  dataio_base.o  parameter_defs.o              \
  readparams.o dataio_text.o               \
  dataio_silo.o               \
  dataio_silo_utility.o       \
  dataio_silo_MPI.o           \
  uniform_grid_pllel.o        \
  comm_mpi.o                  \
  comm_files.o                \
  VectorOps.o                 \
  VectorOps_spherical.o       \
  stellar_wind_BC.o  stellar_wind_angle.o          \
  cell_interface.o            \
  file_status.o               \
  raytracer_SC_pllel.o        \
  raytracer_SC.o

OBJ_CHEM=integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
  hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
  microphysics_base.o mp_only_cooling.o \
  MPv3.o MPv5.o MPv6.o MPv7.o MPv10.o photo_xsections.o

OBJ_BCS = assign_update_bcs.o periodic_boundaries.o \
          outflow_boundaries.o inflow_boundaries.o \
          reflecting_boundaries.o oneway_out_boundaries.o \
          axisymmetric_boundaries.o \
          fixed_boundaries.o jet_boundaries.o \
          jetreflect_boundaries.o double_Mach_ref_boundaries.o \
          stellar_wind_boundaries.o MCMD_boundaries.o \
          periodic_boundaries_MPI.o assign_update_bcs_MPI.o \
          RT_MPI_boundaries.o

OBJ_NG_BCS = assign_update_bcs_NG_MPI.o assign_update_bcs_NG.o \
  NG_coarse_to_fine_boundaries.o \
  NG_fine_to_coarse_boundaries.o \
  NG_MPI_coarse_to_fine_boundaries.o \
  NG_MPI_fine_to_coarse_boundaries.o \
  NG_BC89flux.o \
  NG_MPI_BC89flux.o

OBJ_EQNS =  eqns_base.o eqns_hydro_adiabatic.o eqns_mhd_adiabatic.o 

OBJ_SOLVE = riemannMHD.o HLLD_MHD.o Roe_MHD_ConservedVar_solver.o \
            riemann.o findroot.o Riemann_FVS_hydro.o HLL_hydro.o \
            Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
            solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o



OBJECTS +=  $(OBJ_CHEM) $(OBJ_BCS) $(OBJ_NG_BCS) $(OBJ_EQNS) $(OBJ_SOLVE)

##########################
### -- EXECTUTABLES -- ###
##########################

silo2text: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o silo2text $(LDFLAGS)


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################
silo2text.o:	silo2text.cpp uniform_grid.h dataio_base.h ${GLOBAL}
global.o: global.cpp ${GLOBAL}

uniform_grid.o: uniform_grid.cpp uniform_grid.h ${GLOBAL} stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} comms.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h
dataio_base.o: dataio_base.cpp dataio_base.h ${GLOBAL} solver_eqn_base.h \
    readparams.h stellar_wind_BC.h file_status.h
dataio_silo.o:  dataio_silo.cpp dataio_silo.h dataio_base.h ${GLOBAL} \
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

MPv1.o:	MPv1.cpp MPv1.h microphysics_base.h \
   ${GLOBAL} cooling.h integrator.h

MPv0.o: MPv0.cpp MPv1.h \
   microphysics_base.h ${GLOBAL} cooling.h integrator.h

MPv9.o: MPv9.cpp MPv9.h \
   microphysics_base.h ${GLOBAL} integrator.h

hydrogen_mp.o: hydrogen_mp.cpp hydrogen_mp.h ${GLOBAL} hydrogen_recomb_Hummer94.h \
  hydrogen_photoion.h 

hydrogen_photoion.o: hydrogen_photoion.cpp hydrogen_photoion.h ${GLOBAL}

MPv2.o: MPv2.cpp MPv2.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

cvode_integrator.o: cvode_integrator.cpp cvode_integrator.h ${GLOBAL}

MPv3.o: MPv3.cpp MPv3.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

MPv10.o: MPv10.cpp MPv10.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h

MPv4.o: MPv4.cpp MPv4.h MPv3.h \
  ${GLOBAL} hydrogen_mp.h hydrogen_recomb_Hummer94.h hydrogen_photoion.h \
  cooling_SD93_cie.h microphysics_base.h

MPv5.o: MPv5.cpp MPv5.h \
  MPv3.h microphysics_base.h ${GLOBAL}

MPv6.o: MPv6.cpp MPv6.h \
  MPv3.h microphysics_base.h ${GLOBAL}

MPv7.o: MPv7.cpp MPv7.h \
  MPv3.h microphysics_base.h ${GLOBAL}

MPv8.o: MPv8.cpp \
  MPv8.h MPv3.h microphysics_base.h ${GLOBAL}

HHe_photoion.o: HHe_photoion.cpp HHe_photoion.h ${GLOBAL}
mpv9_HHe.o: mpv9_HHe.cpp mpv9_HHe.h microphysics_base.h \
  HHe_photoion.h cvode_integrator.h ${GLOBAL}

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

axisymmetric_boundaries.o: axisymmetric_boundaries.cpp \
  axisymmetric_boundaries.h \
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

assign_update_bcs.o: assign_update_bcs.cpp assign_update_bcs.h \
  periodic_boundaries.h outflow_boundaries.h \
  oneway_out_boundaries.h inflow_boundaries.h \
  reflecting_boundaries.h \
  boundaries.h ${GLOBAL} grid_base_class.h

assign_update_bcs_MPI.o: assign_update_bcs_MPI.cpp \
  assign_update_bcs_MPI.h periodic_boundaries_MPI.h \
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

#####################
# EQUATIONS CLASSES #
#####################

eqns_base.o: eqns_base.cpp eqns_base.h ${GLOBAL}

eqns_hydro_adiabatic.o: eqns_hydro_adiabatic.cpp eqns_hydro_adiabatic.h \
 eqns_base.h ${GLOBAL}

eqns_mhd_adiabatic.o: eqns_mhd_adiabatic.cpp eqns_mhd_adiabatic.h \
 eqns_base.h ${GLOBAL}

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

###########
# SOLVERS #
###########

solver_eqn_base.o: VectorOps.h solver_eqn_base.h  \
 eqns_base.h ${GLOBAL} cell_interface.h

solver_eqn_hydro_adi.o:	VectorOps.h VectorOps_spherical.h \
   solver_eqn_base.h  eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_hydro_adi.h   eqns_hydro_adiabatic.h \
   solver_eqn_hydro_adi.cpp grid_base_class.h

solver_eqn_mhd_adi.o: VectorOps.h grid_base_class.h \
   solver_eqn_base.h  eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_mhd_adi.h   eqns_mhd_adiabatic.h \
   solver_eqn_mhd_adi.cpp

################################################################



clean:
	rm -f *.o silo2text

