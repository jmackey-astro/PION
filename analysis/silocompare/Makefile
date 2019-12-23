# Makefile for silocompare code, Jonathan Mackey
#
# 2010.10.13 JM: Added new stellar_wind_BC.o object and new standard compile.
# 2010.12.07 JM: Added geometric grid support, and also cell_interface object.
# 2011.04.23 JM: now works on Juropa.
# 2013.09.05 JM: Updated for pion.
# 2015.03.26 JM: updated for pion v0.2.
# 2016.03.13 JM: updated for new pion version (.cc --> .cpp)

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
srcdir = .
SRCDIR=${CURDIR}/../../source
GLOBAL = functionality_flags.h testing_flags.h testing_flags.h reporting.h \
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
MAKE_UNAME=ubuntu18
#MAKE_UNAME=aibn129
#MAKE_UNAME=epona
#MAKE_UNAME=stokes

# All machines and options need readline and maths
#LDFLAGS = -lreadline
# display all warnings
CXXFLAGS = -O3 
CXXFLAGS += -g -Wall

######################################################################
# Options for a "standard" workstation.
# Need libsilo-dev libcfitsio-dev libsundials-serial-dev installed
# on debian/ubuntu systems.
######################################################################
ifeq ($(MAKE_UNAME), standard)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    #CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif
  CXX = mpicxx

  INC  = -I/usr/include
  INC += -I${SRCDIR}
  LDFLAGS += -lm
  # All libraries are in this path:
  LDFLAGS += -L/usr/lib
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
################### --- STANDARD WORKSTATION ---######################

######################################################################
# Options for an Ubuntu 18 workstation with locally installed libraries.
######################################################################
ifeq ($(MAKE_UNAME), ubuntu18)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    CXXFLAGS = -Wall -g -Wno-unused-function
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
../../source/tools:\
../../source/decomposition:\
../../source/boundaries


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
        sim_params.o \
        dataio_base.o parameter_defs.o \
	dataio_text.o dataio_silo.o \
        dataio_silo_MPI.o \
        dataio_silo_utility.o \
        cell_interface.o \
	stellar_wind_BC.o stellar_wind_angle.o \
        uniform_grid.o \
        uniform_grid_pllel.o \
        comm_mpi.o \
        comm_files.o \
        VectorOps.o \
        VectorOps_spherical.o \
	raytracer_SC_pllel.o \
        raytracer_SC.o \
  setup_grid_NG_MPI.o setup_NG_grid.o


OBJ_CHEM=integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
  hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
  microphysics_base.o mp_only_cooling.o \
  MPv3.o MPv5.o MPv6.o MPv7.o MPv10.o photo_xsections.o

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
  NG_MPI_fine_to_coarse_boundaries.o \
  NG_BC89flux.o \
  NG_MPI_BC89flux.o

OBJ_EQNS =  eqns_base.o eqns_hydro_adiabatic.o eqns_mhd_adiabatic.o 

OBJ_SOLVE = riemannMHD.o HLLD_MHD.o Roe_MHD_ConservedVar_solver.o \
            riemann.o findroot.o Riemann_FVS_hydro.o \
            Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
            solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o


OBJECTS += $(OBJ_BCS) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CHEM) $(OBJ_NG_BCS) 


##########################
### -- EXECTUTABLES -- ###
##########################

silocompare: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o silocompare $(LDFLAGS)



####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################

silocompare.o:	silocompare.cpp sim_params.h uniform_grid.h dataio_base.h ${GLOBAL}

reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
constants.o: constants.cpp ${GLOBAL}
timer.o: tools/timer.cpp tools/timer.h ${GLOBAL}
sim_params.o: 	  sim_params.cpp  ${GLOBAL}

comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h
comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}

uniform_grid.o: uniform_grid.cpp uniform_grid.h ${GLOBAL} stellar_wind_BC.h
uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid.h ${GLOBAL} comms.h
stellar_wind_BC.o: stellar_wind_BC.cpp stellar_wind_BC.h ${GLOBAL}

VectorOps.o:           VectorOps.cpp VectorOps.h uniform_grid.h ${GLOBAL}
VectorOps_spherical.o: VectorOps_spherical.cpp VectorOps_spherical.h VectorOps.h \
		       uniform_grid.h ${GLOBAL}

MCMD_control.o: MCMD_control.cpp MCMD_control.h ${GLOBAL}

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


########################
# Raytracing
########################
raytracer_SC.o: raytracer_SC.cpp raytracer_SC.h ${GLOBAL} \
    uniform_grid.h raytracer_base.h

raytracer_SC_pllel.o : raytracer_SC_pllel.cpp raytracer_SC.h \
    ${GLOBAL} uniform_grid.h comms.h raytracer_base.h

########################
# Data I/O
########################
file_status.o: file_status.cpp file_status.h

parameter_defs.o:  parameter_defs.cpp parameter_defs.h ${GLOBAL} \
    solver_eqn_base.h readparams.h \
    stellar_wind_angle.h stellar_wind_BC.h cell_interface.h

dataio_base.o: dataio_base.cpp dataio_base.h ${GLOBAL} solver_eqn_base.h \
    readparams.h stellar_wind_BC.h file_status.h

dataio_text.o:  dataio_text.cpp dataio_text.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h grid_base_class.h dataio_base.h

dataio_silo.o:  dataio_silo.cpp dataio_silo.h dataio_base.h ${GLOBAL} \
    solver_eqn_base.h file_status.h

dataio_silo_MPI.o:    dataio_silo_MPI.cpp dataio_silo.h dataio_base.h \
    ${GLOBAL} comms.h file_status.h

readparams.o:	readparams.cpp readparams.h ${GLOBAL}

dataio_silo_utility.o: dataio_silo_utility.cpp dataio_silo_utility.h\
    ${GLOBAL} file_status.h dataio_base.h dataio_silo.h comms.h

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
	rm -f *.o silocompare

