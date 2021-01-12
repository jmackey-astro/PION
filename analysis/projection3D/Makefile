# Makefile for sim_projection code, Jonathan Mackey
#
# 2010.12.13 JM: Changed aibn129 to standard workstation.
#   Updated Vpath to point to the tagged release used for MHD code.
# 2013.02.05 JM: updated to work with pion.
# 2013.10.15 JM: Updated to use microphysics to get temperature/n(H).
# 2015.07.03 JM: updated for pion_dev: uses MCMD, SimSetup, constants.h
#

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cpp .o

#######################
## COMPILATION RULES ##
#######################
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################

srcdir=.
bindir=.

SRCDIR=${CURDIR}/../../source
CODEDIR=${CURDIR}/../..

VPATH=../projection:../xray:${SRCDIR}:\
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
${SRCDIR}/contrib:\
${SRCDIR}/tools:\
${SRCDIR}/decomposition:\
${SRCDIR}/Riemann_solvers:\
${SRCDIR}/spatial_solvers:\
${SRCDIR}/boundaries


GLOBAL = functionality_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h rad_src_data.h


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

PION_OPTIMISE=HIGH
#PION_OPTIMISE=LOW

#MAKE_UNAME=standard
MAKE_UNAME=ubuntu18
#MAKE_UNAME=ubuntu16
#MAKE_UNAME=debian9
#MAKE_UNAME=debian10
#MAKE_UNAME=OSX
#MAKE_UNAME=SUPERMUC
#MAKE_UNAME=JUROPA

# All machines and options need readline and maths
LDFLAGS = -lm  
ifeq ($(PION_OPTIMISE), LOW)
  LDFLAGS += -lreadline
endif


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

  OPTION += -DCVODE2
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
    CXXFLAGS = -O3 -Wno-unused-function
#    CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif
  
  ifeq ($(IOMODE), USE_FILE_COMMS)
    CXX = g++
  endif
  ifeq ($(IOMODE), USE_MPI)
    CXX = mpicxx
  endif

  OPTION += -DCVODE3
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

  OPTION += -DCVODE5
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

  OPTION += -DCVODE5
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
################### --- UBUNTU 16 ---######################

###########################################
# Options for OSX 
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

  MYLIB = $(CURDIR)/../../extra_libraries
  INC += -I$(MYLIB)/include
  INC += -I/opt/local/include
  INC += -I${SRCDIR}

  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -L/opt/local/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5 
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
  CXXFLAGS = -O3 -Wall
  #CXXFLAGS += -g #-wd981 -wd383
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
  INC += -I${SRCDIR}
  INC += -I/usr/local/sundials/2.4.0/include
  INC += ${MKL_INC} ${MPI_INC}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += ${MKL_SHLIB} ${MPI_CXX_LIB}
  LDFLAGS += -lreadline
endif
##################### -- SuperMUC -- #################################


##################################################
##################################################

CXXFLAGS += $(OPTION)

OBJ_CHEM = integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
  hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
  microphysics_base.o mp_only_cooling.o \
  MPv3.o MPv5.o MPv6.o MPv7.o

OBJ_GEN = reporting.o mem_manage.o command_line_interface.o \
  timer.o interpolate.o constants.o sim_params.o \
  comm_mpi.o comm_files.o

OBJ_IO = parameter_defs.o file_status.o readparams.o \
  dataio_silo.o dataio_fits.o utility_fits_class.o \
  dataio_base.o dataio_text.o dataio_fits_MPI.o \
  dataio_silo_MPI.o dataio_silo_utility.o

OBJ_GRID = cell_interface.o uniform_grid.o setup_fixed_grid.o \
  stellar_wind_BC.o stellar_wind_angle.o MCMD_control.o \
  raytracer_SC.o setup_fixed_grid_MPI.o \
  uniform_grid_pllel.o raytracer_SC_pllel.o \
  setup_grid_NG_MPI.o setup_NG_grid.o

OBJ_CRD = VectorOps.o VectorOps_spherical.o 

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

OBJ_EQNS = eqns_base.o eqns_hydro_adiabatic.o eqns_mhd_adiabatic.o 

OBJ_SOLVE = riemannMHD.o HLLD_MHD.o Roe_MHD_ConservedVar_solver.o \
  riemann.o findroot.o Riemann_FVS_hydro.o HLL_hydro.o \
  Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
  solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o

OBJECTS = $(OBJ_CHEM) $(OBJ_GEN) $(OBJ_IO) $(OBJ_GRID) $(OBJ_CRD) \
  $(OBJ_BCS) $(OBJ_NG_BCS) $(OBJ_EQNS) $(OBJ_SOLVE) \
  point_quantities.o main_projection.o sim_projection.o image_io.o xray_emission.o

##########################
### -- EXECUTABLES -- ###
##########################

projection: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o projection $(LDFLAGS)


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################
reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
command_line_interface.o:  tools/command_line_interface.cpp ${GLOBAL}
timer.o: tools/timer.cpp tools/timer.h ${GLOBAL}
interpolate.o: tools/interpolate.cpp tools/interpolate.h ${GLOBAL}
constants.o: constants.cpp ${GLOBAL}
sim_params.o: sim_params.cpp ${GLOBAL}

main_NG_MPI.o: main_NG_MPI.cpp sim_control.h ${GLOBAL} comms.h
mainMPI.o: mainMPI.cpp sim_control.h ${GLOBAL} comms.h


icgen_NG_MPI.o: icgen_NG_MPI.cpp icgen.h get_sim_info.h \
  icgen_base.h grid_base_class.h \
  uniform_grid.h ${GLOBAL} microphysics_base.h dataio_base.h \
  dataio_silo.h dataio_fits.h setup_fixed_grid.h raytracer_base.h


#####################
# GRID CLASSES      #
#####################

cell_interface.o: cell_interface.cpp cell_interface.h ${GLOBAL}

uniform_grid.o: uniform_grid.cpp uniform_grid.h grid_base_class.h ${GLOBAL} \
    stellar_wind_angle.h stellar_wind_BC.h \
    VectorOps.h VectorOps_spherical.h

VectorOps.o:           VectorOps.cpp VectorOps.h ${GLOBAL} uniform_grid.h grid_base_class.h

VectorOps_spherical.o: VectorOps_spherical.cpp VectorOps_spherical.h VectorOps.h \
		       ${GLOBAL} uniform_grid.h grid_base_class.h

MCMD_control.o: MCMD_control.cpp MCMD_control.h ${GLOBAL}

setup_fixed_grid.o: setup_fixed_grid.cpp setup_fixed_grid.h ${GLOBAL} \
    grid_base_class.h uniform_grid.h assign_update_bcs.h

setup_NG_grid.o: setup_NG_grid.cpp setup_NG_grid.h \
  setup_fixed_grid.h ${GLOBAL} grid_base_class.h uniform_grid.h \
  assign_update_bcs.h assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h 

raytracer_SC.o: raytracer_SC.cpp raytracer_SC.h ${GLOBAL} \
    uniform_grid.h raytracer_base.h

raytracer_SC_pllel.o : raytracer_SC_pllel.cpp raytracer_SC.h \
    ${GLOBAL} uniform_grid.h comms.h raytracer_base.h

setup_fixed_grid_MPI.o: setup_fixed_grid_MPI.cpp setup_fixed_grid_MPI.h \
  setup_fixed_grid.h ${GLOBAL} uniform_grid_pllel.h \
  grid_base_class.h uniform_grid.h assign_update_bcs.h \
  periodic_boundaries_MPI.h assign_update_bcs_MPI.h \
  RT_MPI_boundaries.h VectorOps.h solver_eqn_base.h  eqns_base.h \
  solver_eqn_hydro_adi.h solver_eqn_hydro_adi.cpp

setup_grid_NG_MPI.o: setup_grid_NG_MPI.cpp setup_grid_NG_MPI.h \
  setup_fixed_grid_MPI.h setup_fixed_grid.h setup_NG_grid.h \
  ${GLOBAL} uniform_grid_pllel.h \
  grid_base_class.h uniform_grid.h assign_update_bcs.h \
  periodic_boundaries_MPI.h assign_update_bcs_MPI.h \
  RT_MPI_boundaries.h VectorOps.h solver_eqn_base.h  eqns_base.h \
  solver_eqn_hydro_adi.h solver_eqn_hydro_adi.cpp \
  assign_update_bcs_NG_MPI.h assign_update_bcs_MPI.h \
  assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.h \
  NG_MPI_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h periodic_boundaries_MPI.h \
  MCMD_boundaries.h boundaries.h \
  grid_base_class.h solver_eqn_base.h MCMD_control.h

uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid_pllel.h \
    uniform_grid.h ${GLOBAL} \
    comms.h stellar_wind_angle.h stellar_wind_BC.h \
    VectorOps.h VectorOps_spherical.h

#####################
# BOUNDARY CLASSES  #
#####################

periodic_boundaries.o: periodic_boundaries.cpp periodic_boundaries.h \
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
  boundaries.h ${GLOBAL} grid_base_class.h

NG_MPI_fine_to_coarse_boundaries.o: \
  NG_MPI_fine_to_coarse_boundaries.h NG_fine_to_coarse_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h

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
  
#####################
# CONTROL CLASSES   #
#####################

sim_control.o:  sim_control.cpp sim_control.h ${GLOBAL} solver_eqn_base.h \
  setup_fixed_grid.h  calc_timestep.h assign_update_bcs.h \
  parameter_defs.h dataio_silo.h dataio_fits.h microphysics_base.h \
  sim_init.h setup_fixed_grid.h grid_base_class.h assign_update_bcs.h \
  VectorOps.h

sim_init.o: sim_init.h sim_init.cpp ${GLOBAL} setup_fixed_grid.h \
 microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL} \
 grid_base_class.h assign_update_bcs.h 

sim_init_NG.o: sim_init_NG.h sim_init_NG.cpp sim_init.h \
 ${GLOBAL} setup_NG_grid.h setup_fixed_grid.h grid_base_class.h \
 microphysics_base.h solver_eqn_base.h raytracer_SC.h \
 assign_update_bcs_NG.h

time_integrator.o: time_integrator.cpp time_integrator.h sim_control.h \
 setup_fixed_grid.h calc_timestep.h assign_update_bcs.h grid_base_class.h \
 microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL} \
 assign_update_bcs.h

calc_timestep.o: calc_timestep.cpp calc_timestep.h \
 setup_fixed_grid.h microphysics_base.h grid_base_class.h \
 solver_eqn_base.h raytracer_SC.h ${GLOBAL} \
 assign_update_bcs.h

sim_control_NG.o : sim_control_NG.cpp sim_control_NG.h \
  sim_control.h calc_timestep.h setup_fixed_grid.h ${GLOBAL} \
  solver_eqn_base.h grid_base_class.h VectorOps.h \
  assign_update_bcs.h assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h 

sim_control_MPI.o: sim_control_MPI.cpp sim_control_MPI.h \
  sim_control.h MCMD_control.h ${GLOBAL} comms.h \
  periodic_boundaries_MPI.h assign_update_bcs_MPI.h \
  RT_MPI_boundaries.h setup_fixed_grid.h setup_fixed_grid_MPI.h \
  VectorOps.h

sim_control_NG_MPI.o: sim_control_NG_MPI.cpp sim_control_NG_MPI.h \
  sim_control_MPI.h sim_control_NG.h \
  sim_control.h MCMD_control.h ${GLOBAL} comms.h \
  periodic_boundaries_MPI.h assign_update_bcs_MPI.h \
  RT_MPI_boundaries.h setup_fixed_grid.h setup_fixed_grid_MPI.h \
  VectorOps.h  assign_update_bcs.h assign_update_bcs_NG.h \
  NG_coarse_to_fine_boundaries.h assign_update_bcs_NG_MPI.h \
  NG_fine_to_coarse_boundaries.h NG_MPI_coarse_to_fine_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.h 


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

######################
# MICROPHYSICS STUFF #
######################
microphysics_base.o: microphysics_base.cpp microphysics_base.h ${GLOBAL}

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

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h ${GLOBAL}

parameter_defs.o:  parameter_defs.cpp parameter_defs.h ${GLOBAL} \
    solver_eqn_base.h readparams.h \
    stellar_wind_angle.h stellar_wind_BC.h cell_interface.h

dataio_base.o:  dataio_base.cpp dataio_base.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h

dataio_text.o:  dataio_text.cpp dataio_text.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h grid_base_class.h dataio_base.h

dataio_silo.o:  dataio_silo.cpp dataio_silo.h ${GLOBAL} \
    parameter_defs.h  grid_base_class.h dataio_base.h \
    solver_eqn_base.h file_status.h

dataio_fits.o:  dataio_fits.cpp dataio_fits.h ${GLOBAL} \
    solver_eqn_base.h solver_eqn_hydro_adi.h  \
    parameter_defs.h  grid_base_class.h \
    file_status.h dataio_base.h

dataio_fits_MPI.o: dataio_fits_MPI.cpp dataio_fits_MPI.h \
    dataio_fits.h ${GLOBAL} dataio_base.h \
    solver_eqn_base.h solver_eqn_hydro_adi.h  \
    parameter_defs.h  grid_base_class.h \
    file_status.h uniform_grid_pllel.h

readparams.o: readparams.cpp readparams.h ${GLOBAL}

dataio_silo_MPI.o: dataio_silo_MPI.cpp dataio_silo.h dataio_base.h \
    parameter_defs.h  grid_base_class.h solver_eqn_base.h \
    ${GLOBAL} comms.h file_status.h uniform_grid_pllel.h

dataio_silo_utility.o: dataio_silo_utility.cpp dataio_silo_utility.h \
    dataio_silo.h dataio_base.h solver_eqn_base.h \
    parameter_defs.h  grid_base_class.h \
    ${GLOBAL} comms.h file_status.h uniform_grid_pllel.h

utility_fits_class.o: utility_fits_class.cpp dataio_fits.h dataio_base.h \
    parameter_defs.h  grid_base_class.h \
    ${GLOBAL} file_status.h uniform_grid_pllel.h


#######################
## MPI/files comms classes
########################
comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h ${GLOBAL}

comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}

###############################
# Stuff Specific to this code #
###############################

point_quantities.o: point_quantities.cpp point_quantities.h ${GLOBAL}

sim_projection.o: sim_projection.cpp sim_projection.h \
    point_quantities.h xray_emission.h \
    ${GLOBAL} uniform_grid.h dataIO/dataio_base.h

main_projection.o: main_projection.cpp sim_projection.h \
    point_quantities.h xray_emission.h \
    ${GLOBAL} uniform_grid.h dataIO/dataio_base.h
    
image_io.o: image_io.cpp image_io.h dataIO/dataio_base.h \
    dataIO/dataio_fits.h ${GLOBAL}

clean:
	rm -f *.o projection
