# Makefile for 'pion',
# Author: Jonathan Mackey
#
# - 2010.06.24 JM: added section for aibn214 (2GB RAM, 2 cores)
# - 2010.07.21 JM: fixed aibn214 section to use system HDF5
# - 2010.07.29 JM: changed aibn214 for new computer aibn129 (8GB,4core,64-bit)
# - 2010-09-22 JM: Added Riemann_FVS_hydro files.
# - 2010-10-01 JM: Added spherical vector-ops file.
# - 2010-10-05 JM: Added stellar_wind_BC object file.
# - 2010-10-11 JM: Added 'standard' workstation
# - 2010.10.13 JM: Added MPv9 object file
# - 2010.10.13 JM: moved executables to ../exe/
# - 2010.11.09 JM: Added options for Phalanx
# - 2010.11.15 JM: Added cell_interface.o file.
# - 2010.12.22 JM: Added new files for solvers.
# - 2010.12.27 JM: added extra source files. Added GLOBAL
#                  variable for list of global header files.
# - 2010.12.28 JM: Added internal-energy integrating solver files.
# - 2011.01.11 JM: Added options for JUROPA at JSC Juelich.
# - 2011.01.17 JM: Added new microphysics stuff
# - 2011.02.17 JM: Added new raytracer_shielding.o object (with path).
# - 2011.03.17 JM: Added Harpreet's module files
# - 2011.04.12 JM: Added Hummer-recomb microphysics file.
# - 2011.04.14 JM: Added new MP interface, with CVODES library support.
#                  Simplified lib-paths.
# - 2011.05.12 JM: updated phalanx library paths
# - 2011.06.20 JM: updated iMac libraries to work on my MBP.
# - 2011.10.13 JM: added new microphysics files.
# - 2011.11.14 JM: Got rid of MAKE_UNAME Stuff (moved to compile_code.sh)
#                  Deleted some old machines from the list of computers.
# - 2011.12.20 JM: Merged code and icgenerator makefiles and moved most of the
#     logic to compile_code.sh.  Re-named to just "Makefile", and removed most
#     of the machine-specific options.
# - 2012.02.08 JM: Added new IC class for Harpreet's module.
# - 2012.02.22 JM/HD: Added section for dougal at UCL compile options.
# - 2012.07.15 JM: some new vpath entries and changed file names.
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.01.10 JM: Changed name of exe to pion_serial.
#                  Added StarBench test problems class to icgen.
# - 2013.01.17 JM: Got rid of readline/ncurses from link line in
#    production version of pion.
# - 2013.02.15 JM: Added NEW_METALLICITY flag for testing (with new
#    source files).
# - 2013.02.19 JM: Added file_status.o/cpp object file.
# - 2013.02.27 JM: Added extensions for contributed code.
# - 2013.08.19 JM: Changed filenames for Harpreet's module.
# - 2014.04.14 JM: Added options for JUDGE at JSC.
# - 2015.01.08 JM: Added tools/reporting.h.cpp and tools/mem_manage.h
# - 2015.01.12 JM: Added command_line_interface.h.cpp.
# - 2016.02.12 JM: renamed .cc files to .cpp
# - 2016.05.04 JM: Added FIONN to list of machines
# - 2016.08.29 JM: removed *Eint solvers.
# - 2017-07-29 JM: updated for latitude-dependent winds.
# - 2017-08-03 JM: updated standard workstation to use system libs.
# - 2018-[05-09] JM: added static mesh-refinement
#
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
# Options for a "standard" workstation, assumed to be like debian linux
# System CVODE (Sundials), Silo, and cfitsio libraries should be installed.
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
  LDFLAGS += -lm
  # All libraries are in this path:
  # Libraries for fits, silo, cvode
  LDFLAGS += -lcfitsio
  #LDFLAGS += -lsiloh5
  LDFLAGS += -lsilo
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
################### --- local libs WORKSTATION ---######################

######################## -- JUROPA -- ################################
# Options for the Juelich Supercomputer JUROPA
# Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), JUROPA)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g -wd981 -wd383
  endif
  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  #CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check

  MYLIB = ../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
#  INC += -I/usr/local/sundials/2.4.0/include
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
#  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lmkl_sequential -lmkl_core

endif
##################### -- JUROPA -- ###################################

######################## -- FIONN -- ################################
# Options for the ICHEC Supercomputer FIONN
# Intel Compilers
######################################################################
ifeq ($(MAKE_UNAME), FIONN)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g -wd981 -wd383
  endif
  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  #CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check

  MYLIB = ../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
##################### -- FIONN -- ###################################
######################################################################

######################## -- JUDGE -- ################################
# Options for the Juelich Supercomputer JUDGE.  Intel Compilers!
#####################################################################
ifeq ($(MAKE_UNAME), JUDGE)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g -wd981 -wd383
  endif
  MYLIB = ../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS  += -L${MYLIB}/lib -lcfitsio -lsilo
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
##################### -- JUDGE -- ###################################



######################## -- SuperMUC -- ##############################
# Options for the LRZ Supercomputer SuperMUC.  Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), SUPERMUC)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g -wd981 -wd383
  endif
  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  #CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check

  MYLIB = $(CURDIR)/../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  INC += -I/usr/local/sundials/2.4.0/include
  INC += ${MKL_INC}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += ${MKL_LIB}

endif
##################### -- SuperMUC -- ###################################
######################################################################


######################################################################
########## OPTIONS FOR PHALANX ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), phalanx)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -Wall -g
  endif
  MYLIB = ../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm
endif
#################### --- PHALANX --- #################################
######################################################################

######################################################################
########### OPTIONS FOR DOUGAL ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), dougal)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -Wall -g
  endif
  MYLIB = ../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm
endif
#################### --- DOUGAL --- #################################
######################################################################



######################################################################
####  OS X
######################################################################
ifeq ($(MAKE_UNAME), OSX)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
    #CXXFLAGS = -O0 -Wall -g
  else
    CXXFLAGS = -O0 -Wall -g
    LDFLAGS = -lreadline
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I/opt/local/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -L/opt/local/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lcfitsio
  #LDFLAGS += -lsiloh5
  LDFLAGS += -lsilo
  LDFLAGS += -lm  
  LDFLAGS += -Wl,-rpath,${MYLIB}/lib

endif
######################################################################
######################################################################



#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################
CXXFLAGS += $(OPT) $(PION_OPTIONS)

OBJ_CHEM=integrator.o cvode_integrator.o cooling.o cooling_SD93_cie.o \
 hydrogen_mp.o hydrogen_photoion.o hydrogen_recomb_Hummer94.o \
 microphysics_base.o mp_only_cooling.o \
 MPv3.o MPv5.o MPv6.o MPv7.o MPv10.o photo_xsections.o

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
            riemann.o findroot.o Riemann_FVS_hydro.o \
            Roe_Hydro_ConservedVar_solver.o Roe_Hydro_PrimitiveVar_solver.o \
            solver_eqn_base.o solver_eqn_hydro_adi.o solver_eqn_mhd_adi.o

OBJ_CRD = VectorOps.o VectorOps_spherical.o 

OBJ_ICS = shock_tube.o radiative_shock.o blast_wave.o basic_tests.o \
          shock_cloud.o jet.o laser_ablation.o spherical_clump.o \
          photoevaporating_clump.o photoevaporating_random_clumps.o \
          photoevaporating_multiclumps.o inside_sphere.o \
          icgen_base.o

OBJ_BCS = assign_update_bcs.o periodic_boundaries.o \
          outflow_boundaries.o inflow_boundaries.o \
          reflecting_boundaries.o oneway_out_boundaries.o \
          fixed_boundaries.o jet_boundaries.o \
          jetreflect_boundaries.o double_Mach_ref_boundaries.o \
          stellar_wind_boundaries.o

OBJ_NG_BCS = assign_update_bcs_NG.o     \
          NG_coarse_to_fine_boundaries.o \
          NG_fine_to_coarse_boundaries.o


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
all: ${bindir}/pion_serial ${bindir}/icgen_serial  \
     ${bindir}/pion_NG_serial ${bindir}/icgen_NG_serial

uni: ${bindir}/pion_serial ${bindir}/icgen_serial

NG: ${bindir}/icgen_NG_serial  ${bindir}/pion_NG_serial 

${bindir}/pion_serial : $(OBJECTS_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_PION) \
        -o ${bindir}/pion_serial $(LDFLAGS)

${bindir}/icgen_serial: $(OBJECTS_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_ICGEN) \
        -o ${bindir}/icgen_serial $(LDFLAGS)

${bindir}/pion_NG_serial : $(OBJECTS_NG_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_PION) \
        -o ${bindir}/pion_NG_serial $(LDFLAGS)

${bindir}/icgen_NG_serial: $(OBJECTS_NG_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_ICGEN) \
        -o ${bindir}/icgen_NG_serial $(LDFLAGS)

clean:
	rm -f *.o ${bindir}/pion_serial ${bindir}/icgen_serial \
        ${bindir}/pion_NG_serial ${bindir}/icgen_NG_serial


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
  NG_fine_to_coarse_boundaries.h boundaries.h

raytracer_SC.o: raytracer_SC.cpp raytracer_SC.h raytracer_base.h \
  ${GLOBAL} uniform_grid.h grid_base_class.h


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
  NG_fine_to_coarse_boundaries.h boundaries.h

#####################
# EQUATIONS CLASSES #
#####################

eqns_base.o:	          eqns_base.cpp eqns_base.h ${GLOBAL}

eqns_hydro_adiabatic.o: eqns_hydro_adiabatic.cpp eqns_hydro_adiabatic.h \
			eqns_base.h ${GLOBAL}

eqns_hydro_isothermal.o: eqns_hydro_isothermal.cpp \
                  eqns_hydro_isothermal.h eqns_base.h ${GLOBAL}

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

MPv10.o: MPv10.cpp MPv10.h ${GLOBAL} hydrogen_mp.h \
  hydrogen_recomb_Hummer94.h hydrogen_photoion.h cooling_SD93_cie.h \
  microphysics_base.h  photo_xsections.h

photo_xsections.o:  photo_xsections.cpp  photo_xsections.h ${GLOBAL}

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
    solver_eqn_hydro_adi.h solver_eqn_hydro_iso.h grid_base_class.h

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


