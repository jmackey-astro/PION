# Makefile for 'pion',
# Author: Jonathan Mackey
#
# 2011-2020

############################
### BASIC INFO AND PATHS ###
############################

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
bindir = ..
SRCDIR=${CURDIR}/../source
VPATH = \
../source:\
../source/ics:\
../source/Riemann_solvers:\
../source/defines:\
../source/microphysics:\
../source/raytracing:\
../source/microphysics/Harpreets_code:\
../source/grid:\
../source/dataIO:\
../source/comms:\
../source/coord_sys:\
../source/equations:\
../source/spatial_solvers:\
../source/sim_control:\
../source/contrib:\
../source/tools:\
../source/future:\
../source/decomposition\
../source/NG_grid\
../source/boundaries


GLOBAL = functionality_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h

ifneq ($(PION_OPTIMISE), HIGH)
  GLOBAL += command_line_interface.h
endif

ifneq (,$(findstring CODE_EXT_HHE,$(PION_OPTIONS)))
  #GLOBAL += constants.h
endif

#######################
## COMPILATION RULES ##
#######################
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################



######################################################################
# Options for a "standard" workstation, a fallback in case the OS was
# not detected by compile_code.sh.
######################################################################
ifeq ($(MAKE_UNAME), standard)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
    #CXXFLAGS = -O0 -Wall -g
  else
    CXXFLAGS = -O0 -Wall -g
  endif
  INC = -I/usr/include
  INC += -I${SRCDIR}
  LDFLAGS +=  -lm
  LDFLAGS += -L/usr/lib
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  #LDFLAGS += -lsilo
  LDFLAGS += -lgsl
  LDFLAGS +=-Wl,--rpath -Wl,/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################
######################################################################

######################################################################
# Options for a workstation with locally installed libraries.
######################################################################
ifeq ($(MAKE_UNAME), locallibs)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    #CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -Wl,-rpath,${MYLIB}/lib
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu/
  LDFLAGS += -lm
  LDFLAGS += -lcfitsio
  #LDFLAGS += -lsiloh5
  LDFLAGS += -lsilo
  LDFLAGS += -lgsl
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
################### --- local libs WORKSTATION ---######################

######################################################################
# Options for a debian 9 workstation, using system libraries
######################################################################
ifeq ($(MAKE_UNAME), debian9)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
    #CXXFLAGS = -O0 -Wall -g
  else
    CXXFLAGS = -O0 -Wall -g
  endif
  INC = -I/usr/include
  INC += -I${SRCDIR}
  LDFLAGS +=  -lm
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  LDFLAGS += -lgsl -lgslcblas
  LDFLAGS += -L/usr/lib
  LDFLAGS +=-Wl,--rpath -Wl,/usr/lib -lsundials_cvode -lsundials_nvecserial
  
  # For sundials 3/4/5 if needed
  #MYLIB = $(CURDIR)/../extra_libraries
  #INC  += -I${MYLIB}/include
  #LDFLAGS += -L${MYLIB}/lib
  #LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib
  #LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
######################################################################

######################################################################
# Options for a debian 10 workstation, using system libraries
######################################################################
ifeq ($(MAKE_UNAME), debian10)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
    #CXXFLAGS = -O0 -Wall -g
  else
    CXXFLAGS = -O0 -Wall -g
  endif
  INC = -I/usr/include
  INC += -I${SRCDIR}
  LDFLAGS += -lm
  LDFLAGS += -L/usr/lib/x86_64-linux-gnu
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  LDFLAGS += -lgsl -lgslcblas
  LDFLAGS += -L/usr/lib
  LDFLAGS += -Wl,--rpath -Wl,/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################

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

  MYLIB = $(CURDIR)/../extra_libraries
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
################### --- Ubuntu 18 WORKSTATION ---######################

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

  MYLIB = $(CURDIR)/../extra_libraries
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
################### --- Ubuntu 16 WORKSTATION ---######################


######################################################################
####  OS X with Macports and no self-installed libraries
######################################################################
ifeq ($(MAKE_UNAME), OSX-MP)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g
    LDFLAGS = -lreadline
  endif

  INC = -I/opt/local/include
  INC += -I${SRCDIR}
  LDFLAGS += -L/opt/local/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5
  LDFLAGS += -lgsl
  LDFLAGS += -lm  
  LDFLAGS += -Wl,-rpath,${MYLIB}/lib

endif
#######################  OSX-MP  #####################################

######################################################################
####  OS X with Homebrew libraries + locally compiled SILO
######################################################################
ifeq ($(MAKE_UNAME), OSX-HB)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
    #CXXFLAGS = -O0 -Wall -g
  else
    CXXFLAGS = -O0 -Wall -g
    LDFLAGS = -lreadline
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC = -I${MYLIB}/include
  #INC += -I/opt/local/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib
  #LDFLAGS += -L/opt/local/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  LDFLAGS += -lgsl
  LDFLAGS += -lm  
  LDFLAGS += -Wl,-rpath,${MYLIB}/lib

endif
#######################  OSX-HB  #####################################



#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################
CXXFLAGS += $(OPT) $(PION_OPTIONS)

OBJ_CHEM=integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
 hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
 microphysics_base.o mp_only_cooling.o \
 MPv3.o MPv5.o MPv6.o MPv7.o 

ifneq (,$(findstring LEGACY_CODE,$(PION_OPTIONS)))
  OBJ_CHEM +=  MPv0.o MPv1.o MPv2.o MPv4.o MPv8.o
endif

ifneq (,$(findstring HARPREETS_CODE_EXT,$(PION_OPTIONS)))
  OBJ_CHEM +=  MPv9.o HD_MetalFree_integrator.o HD_MetalFree_equations.o
endif

ifneq (,$(findstring CODE_EXT_HHE,$(PION_OPTIONS)))
  OBJ_CHEM += HHe_photoion.o mpv9_HHe.o
endif

OBJ_GEN = reporting.o mem_manage.o command_line_interface.o timer.o interpolate.o constants.o sim_params.o 

OBJ_IO =  parameter_defs.o file_status.o readparams.o \
          dataio_silo.o dataio_fits.o utility_fits_class.o \
          dataio_base.o dataio_text.o

OBJ_GRID =  cell_interface.o uniform_grid.o setup_fixed_grid.o \
            stellar_wind_BC.o stellar_wind_angle.o MCMD_control.o raytracer_SC.o

OBJ_EQNS =  eqns_base.o eqns_hydro_adiabatic.o eqns_mhd_adiabatic.o 

OBJ_SOLVE = riemannMHD.o HLLD_MHD.o Roe_MHD_ConservedVar_solver.o \
            riemann.o findroot.o Riemann_FVS_hydro.o HLL_hydro.o \
            Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
            solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o

OBJ_CRD = VectorOps.o VectorOps_spherical.o 

OBJ_ICS = shock_tube.o radiative_shock.o blast_wave.o basic_tests.o \
          shock_cloud.o jet.o laser_ablation.o spherical_clump.o \
          photoevaporating_clump.o photoevaporating_random_clumps.o \
          photoevaporating_multiclumps.o inside_sphere.o \
          icgen_base.o read_1Dto2D.o

OBJ_BCS = assign_update_bcs.o periodic_boundaries.o \
          outflow_boundaries.o inflow_boundaries.o \
          reflecting_boundaries.o oneway_out_boundaries.o \
          fixed_boundaries.o jet_boundaries.o \
          jetreflect_boundaries.o double_Mach_ref_boundaries.o \
          stellar_wind_boundaries.o axisymmetric_boundaries.o

OBJ_NG_BCS = assign_update_bcs_NG.o     \
          NG_coarse_to_fine_boundaries.o \
          NG_fine_to_coarse_boundaries.o NG_BC89flux.o


ifneq (,$(findstring CODE_EXT_SBII,$(PION_OPTIONS)))
  OBJ_ICS += StarBench_test.o
endif

ifneq (,$(findstring BBTURBULENCE_CODE_EXT,$(PION_OPTIONS)))
  OBJ_ICS += read_BBurkhart_data.o
endif

ifneq (,$(findstring HARPREETS_CODE_EXT,$(PION_OPTIONS)))
  OBJ_ICS += HD_2D_ShockCloud.o  # HD_3D_ShockCloud.o
endif

OBJ_SIM = sim_control.o time_integrator.o calc_timestep.o sim_init.o

OBJ_NG_ICS = setup_NG_grid.o sim_init.o

OBJ_NG_SIM = setup_NG_grid.o sim_control_NG.o

OBJECTS_PION = main.o $(OBJ_SIM) $(OBJ_CHEM) $(OBJ_GEN) $(OBJ_IO) \
  $(OBJ_GRID) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CRD) $(OBJ_BCS)

OBJECTS_ICGEN =  icgen.o get_sim_info.o $(OBJ_GEN) $(OBJ_GRID) $(OBJ_IO) \
  $(OBJ_CHEM) $(OBJ_ICS) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CRD) $(OBJ_BCS)

OBJECTS_NG_PION = main_NG.o $(OBJ_NG_SIM) $(OBJ_SIM) $(OBJ_CHEM) \
  $(OBJ_GEN) $(OBJ_IO) $(OBJ_GRID) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CRD) \
  $(OBJ_BCS) $(OBJ_NG_BCS)

OBJECTS_NG_ICGEN = icgen_NG.o get_sim_info.o $(OBJ_GEN) $(OBJ_GRID) \
  $(OBJ_IO) $(OBJ_CHEM) $(OBJ_ICS) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CRD) \
  $(OBJ_NG_ICS) $(OBJ_BCS)  $(OBJ_NG_BCS)

#####################################################
##### END OF LIST OF OBJECT FILES               #####
#####################################################

##########################
### -- EXECTUTABLES -- ###
##########################
all: ${bindir}/pion-ugs ${bindir}/icgen-ugs  \
     ${bindir}/pion-ngs ${bindir}/icgen-ngs

uni: ${bindir}/pion-ugs ${bindir}/icgen-ugs

NG: ${bindir}/icgen-ngs  ${bindir}/pion-ngs 

${bindir}/pion-ugs : $(OBJECTS_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_PION) \
        -o ${bindir}/pion-ugs $(LDFLAGS)

${bindir}/icgen-ugs: $(OBJECTS_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_ICGEN) \
        -o ${bindir}/icgen-ugs $(LDFLAGS)

${bindir}/pion-ngs : $(OBJECTS_NG_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_PION) \
        -o ${bindir}/pion-ngs $(LDFLAGS)

${bindir}/icgen-ngs: $(OBJECTS_NG_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_ICGEN) \
        -o ${bindir}/icgen-ngs $(LDFLAGS)

clean:
	rm -f *.o ${bindir}/pion-ugs ${bindir}/icgen-ugs \
        ${bindir}/pion-ngs ${bindir}/icgen-ngs


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

main.o:         main.cpp sim_control.h ${GLOBAL}
main_NG.o:         main_NG.cpp sim_control_NG.h sim_control.h ${GLOBAL}

icgen.o: icgen.cpp icgen.h get_sim_info.h icgen_base.h grid_base_class.h \
  uniform_grid.h ${GLOBAL} microphysics_base.h dataio_base.h \
  dataio_silo.h dataio_fits.h setup_fixed_grid.h raytracer_base.h

icgen_NG.o: icgen_NG.cpp icgen.h get_sim_info.h \
  icgen_base.h grid_base_class.h \
  uniform_grid.h ${GLOBAL} microphysics_base.h dataio_base.h \
  dataio_silo.h dataio_fits.h setup_fixed_grid.h raytracer_base.h \
  setup_NG_grid.h sim_init.h assign_update_bcs_NG.h

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
    grid_base_class.h uniform_grid.h

setup_NG_grid.o: setup_NG_grid.cpp setup_NG_grid.h \
  setup_fixed_grid.h ${GLOBAL} grid_base_class.h uniform_grid.h \
  assign_update_bcs_NG.h NG_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h boundaries.h NG_BC89flux.h

raytracer_SC.o: raytracer_SC.cpp raytracer_SC.h raytracer_base.h \
  ${GLOBAL} uniform_grid.h grid_base_class.h rad_src_data.h


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
  reflecting_boundaries.h axisymmetric_boundaries.h \
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

NG_BC89flux.o: NG_BC89flux.cpp NG_BC89flux.h ${GLOBAL} \
  grid_base_class.h uniform_grid.h

#####################
# CONTROL CLASSES   #
#####################
sim_control.o:  sim_control.cpp sim_control.h ${GLOBAL} solver_eqn_base.h \
  setup_fixed_grid.h  calc_timestep.h assign_update_bcs.h \
  parameter_defs.h dataio_silo.h dataio_fits.h microphysics_base.h \
  sim_init.h setup_fixed_grid.h grid_base_class.h

sim_init.o: sim_init.h sim_init.cpp ${GLOBAL} setup_fixed_grid.h \
 microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL} \
 grid_base_class.h assign_update_bcs.h

time_integrator.o: time_integrator.cpp time_integrator.h sim_control.h \
 setup_fixed_grid.h calc_timestep.h assign_update_bcs.h grid_base_class.h \
 microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL}

calc_timestep.o: calc_timestep.cpp calc_timestep.h \
 setup_fixed_grid.h microphysics_base.h grid_base_class.h \
 solver_eqn_base.h raytracer_SC.h ${GLOBAL}


sim_control_NG.o : sim_control_NG.cpp sim_control_NG.h \
  sim_control.h calc_timestep.h setup_fixed_grid.h ${GLOBAL} \
  solver_eqn_base.h grid_base_class.h assign_update_bcs_NG.h \
  NG_coarse_to_fine_boundaries.h \
  NG_fine_to_coarse_boundaries.h boundaries.h NG_BC89flux.h

#####################
# EQUATIONS CLASSES #
#####################

eqns_base.o:	          eqns_base.cpp eqns_base.h ${GLOBAL}

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

solver_eqn_base.o: VectorOps.h solver_eqn_base.h grid_base_class.h \
                   eqns_base.h ${GLOBAL} cell_interface.h

solver_eqn_hydro_adi.o:	VectorOps.h VectorOps_spherical.h \
   solver_eqn_base.h eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_hydro_adi.h eqns_hydro_adiabatic.h \
   solver_eqn_hydro_adi.cpp grid_base_class.h

solver_eqn_mhd_adi.o: VectorOps.h grid_base_class.h \
   solver_eqn_base.h  eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_mhd_adi.h   eqns_mhd_adiabatic.h \
   solver_eqn_mhd_adi.cpp


######################
# MICROPHYSICS STUFF #
######################

microphysics_base.o: microphysics_base.cpp microphysics_base.h ${GLOBAL}

integrator.o:	integrator.cpp integrator.h ${GLOBAL}
cooling.o:	cooling.cpp cooling.h tools/interpolate.h ${GLOBAL}

cooling_SD93_cie.o: cooling_SD93_cie.cpp cooling_SD93_cie.h \
 tools/interpolate.h ${GLOBAL}

hydrogen_recomb_Hummer94.o: hydrogen_recomb_Hummer94.h hydrogen_recomb_Hummer94.cpp \
   tools/interpolate.h  ${GLOBAL}

mp_only_cooling.o: mp_only_cooling.cpp mp_only_cooling.h ${GLOBAL} \
   cooling_SD93_cie.h microphysics_base.h hydrogen_recomb_Hummer94.h

MPv1.o:	MPv1.cpp MPv1.h microphysics_base.h \
   ${GLOBAL} cooling.h integrator.h tools/interpolate.h
MPv0.o: MPv0.cpp MPv1.h \
   microphysics_base.h ${GLOBAL} cooling.h integrator.h
MPv9.o: MPv9.cpp MPv9.h \
   microphysics_base.h ${GLOBAL} integrator.h

hydrogen_mp.o: hydrogen_mp.cpp hydrogen_mp.h ${GLOBAL} hydrogen_recomb_Hummer94.h \
  hydrogen_photoion.h tools/interpolate.h 

hydrogen_photoion.o: hydrogen_photoion.cpp hydrogen_photoion.h \
  tools/interpolate.h ${GLOBAL}

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

HHe_photoion.o: HHe_photoion.cpp HHe_photoion.h ${GLOBAL}
mpv9_HHe.o: mpv9_HHe.cpp mpv9_HHe.h microphysics_base.h \
  HHe_photoion.h cvode_integrator.h ${GLOBAL}

MPv8.o: MPv8.cpp \
  MPv8.h MPv3.h microphysics_base.h ${GLOBAL}


########################
# Data I/O
########################
file_status.o:  file_status.cpp file_status.h ${GLOBAL}
parameter_defs.o:       parameter_defs.cpp parameter_defs.h ${GLOBAL}\
    solver_eqn_base.h readparams.h \
    stellar_wind_angle.h stellar_wind_BC.h cell_interface.h
dataio_base.o:  dataio_base.cpp dataio_base.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h
dataio_text.o:  dataio_text.cpp dataio_text.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h grid_base_class.h
dataio_silo.o:  dataio_silo.cpp dataio_silo.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h grid_base_class.h
dataio_fits.o:  dataio_fits.cpp dataio_fits.h parameter_defs.h \
    ${GLOBAL} solver_eqn_base.h \
    solver_eqn_hydro_adi.h grid_base_class.h

readparams.o:   readparams.cpp readparams.h ${GLOBAL}
utility_fits_class.o: utility_fits_class.cpp dataio_fits.h parameter_defs.h ${GLOBAL} grid_base_class.h

stellar_wind_BC.o: stellar_wind_BC.cpp stellar_wind_BC.h \
   grid_base_class.h tools/interpolate.h ${GLOBAL}
stellar_wind_angle.o: stellar_wind_angle.cpp stellar_wind_angle.h stellar_wind_BC.h \
   grid_base_class.h tools/interpolate.h ${GLOBAL}

dataio_silo_NG.o : dataio_silo_NG.cpp \
  dataio_silo_NG.h dataio_silo.h parameter_defs.h ${GLOBAL} \
  solver_eqn_base.h grid_base_class.h

##############################
##### IC GENERATOR FILES #####
##############################
icgen_base.o: icgen_base.cpp icgen_base.h ${GLOBAL}
get_sim_info.o: get_sim_info.h get_sim_info.cpp ${GLOBAL} readparams.h stellar_wind_BC.h
inside_sphere.o: inside_sphere.cpp inside_sphere.h ${GLOBAL}
########################
shock_tube.o: shock_tube.cpp icgen_base.h icgen.h ${GLOBAL}
radiative_shock.o: radiative_shock.cpp icgen_base.h icgen.h ${GLOBAL}
blast_wave.o:      blast_wave.cpp icgen_base.h icgen.h ${GLOBAL} inside_sphere.h
photoevaporating_clump.o:  photoevaporating_clump.cpp icgen_base.h icgen.h ${GLOBAL} inside_sphere.h
basic_tests.o: basic_tests.cpp icgen_base.h icgen.h ${GLOBAL} 
photoevaporating_random_clumps.o: photoevaporating_random_clumps.cpp icgen_base.h icgen.h ${GLOBAL}
shock_cloud.o: shock_cloud.cpp icgen_base.h icgen.h ${GLOBAL} inside_sphere.h
photoevaporating_multiclumps.o: photoevaporating_multiclumps.cpp icgen_base.h icgen.h ${GLOBAL}
laser_ablation.o:   laser_ablation.cpp icgen_base.h icgen.h ${GLOBAL}
jet.o:		    jet.cpp icgen_base.h icgen.h ${GLOBAL}
spherical_clump.o: spherical_clump.cpp icgen_base.h icgen.h ${GLOBAL}
StarBench_test.o: StarBench_test.cpp icgen_base.h icgen.h ${GLOBAL}
read_1Dto2D.o: read_1Dto2D.cpp icgen_base.h icgen.h ${GLOBAL}
#######################

###############################################
## Contributed and/or project-specific code. ##
###############################################
read_BBurkhart_data.o:  read_BBurkhart_data.cpp icgen_base.h icgen.h ${GLOBAL} dataio_fits.h parameter_defs.h

HD_2D_ShockCloud.o:     HD_2D_ShockCloud.cpp icgen_base.h icgen.h ${GLOBAL}
HD_3D_ShockCloud.o:     HD_3D_ShockCloud.cpp icgen_base.h icgen.h ${GLOBAL}
HD_MetalFree_equations.o:  HD_MetalFree_equations.cpp HD_MetalFree.h ${GLOBAL}
HD_MetalFree_integrator.o: HD_MetalFree_integrator.cpp HD_MetalFree.h ${GLOBAL}
###############################################
###############################################


