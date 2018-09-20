# Makefile for the Hydro/MHD Code "pion", parallel version.
# Author: Jonathan Mackey
#
# - 2009.12.18 JM: updated epona libs.
# - 2009.02.03 JM: updated stokes debug compile flags.
# - 2010.04.21 JM: added section for furfur (32GB RAM, 8 cores)
# - 2010.06.24 JM: added section for aibn214 (2GB RAM, 2 cores)
# - 2010.07.28 JM: added section for klaipedaX (64GB RAM, 16 cores)
# - 2010-09-22 JM: Added Riemann_FVS_hydro files.
# - 2010-10-11 JM: Added 'standard' workstation; added stellar_wind_BC file.
# - 2010.10.13 JM: Added MPv9 object file
# - 2010.10.13 JM: moved executables to ../exe/
# - 2010.11.09 JM: Added options for Phalanx at UCL
# - 2010.11.15 JM: Added cell_interface.o file.
# - 2010.12.22 JM: Added new files for solvers.
# - 2010.12.27 JM: added extra source files. Added GLOBAL
#                  variable for list of global header files.
# - 2010.12.28 JM: Added internal-energy integrating solver files.
# - 2011.01.10 JM: Added options for JUROPA at JSC Juelich.
# - 2011.01.17 JM: Added new microphysics stuff
# - 2011.02.17 JM: Added new raytracer_shielding.o object (with path).
# - 2011.03.22 JM: Added Harpreet's module files
# - 2011.04.12 JM: Added Hummer-recomb microphysics file.
# - 2011.04.22 JM: Added new MP interface, with CVODES library support.
#                  Simplified lib-paths.
# - 2011.05.12 JM: updated phalanx library paths
# - 2011.10.24 JM: updated RT and MP (removed shielding file, added MP files).
# - 2011.11.14 JM: Got rid of MAKE_UNAME Stuff (moved to compile_code.sh)
#                  Deleted some old machines from the list of computers.
# - 2011.12.20 JM: Merged code and icgenerator makefiles and moved most of the
#     logic to compile_code.sh.  Re-named to just "Makefile", and removed most
#     of the machine-specific options.
# - 2012.02.08 JM: Added new IC class for Harpreet's module.
# - 2012.02.22 JM/HD: Added section for compile options for dougal at UCL.
# - 2012.02.22 JM/HD: Added section for dougal at UCL compile options.
# - 2012.07.15 JM: some new vpath entries and changed file names.
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.01.10 JM: Changed name of exe to pion_parallel
# - 2013.01.14 JM: Added section for DIRAC/Complexity (it works now).
# - 2013.01.17 JM: Got rid of readline/ncurses from link line in
#    production version of pion.
# - 2013.02.27 JM: Added extensions for contributed code.
#    Added file_status.o/cpp object file. Added NEW_METALLICITY flag.
# - 2013.08.19 JM: Changed filenames for Harpreet's module.
# - 2014.04.14 JM: Section for Judge at JSC.
# - 2015.01.26-28 JM: Edited for new filenames and code structure.
# - 2016.03.14 JM: renamed .cc files to .cpp
# - 2016.05.04 JM: Added FIONN to list of machines
# - 2016.08.29 JM: removed *Eint solvers.
# - 2017-07-29 JM: updated for latitude-dependent winds.
# - 2017-08-03 JM: updated standard workstation to use system libs.
# - 2018-[07-09] JM: added static mesh-refinement
###

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
../source/grid:\
../source/dataIO:\
../source/comms:\
../source/coord_sys:\
../source/equations:\
../source/spatial_solvers:\
../source/sim_control:\
../source/time_integrators:\
../source/contrib:\
../source/tools:\
../source/decomposition:\
../source/NG_grid:\
../source/boundaries


GLOBAL = functionality_flags.h testing_flags.h reporting.h \
         sim_params.h sim_constants.h constants.h rad_src_data.h

#######################
## COMPILATION RULES ##
#######################
.cc.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC) -c $<
#######################


######################################################################
# Options for a "standard" workstation.
# Need libsilo-dev libcfitsio-dev libsundials-serial-dev installed
# on debian/ubuntu systems.
######################################################################
ifeq ($(MAKE_UNAME), standard)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif

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
# Options for a workstation with locally installled libraries.
######################################################################
ifeq ($(MAKE_UNAME), locallibs)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    CXXFLAGS = -Wall -g -Wno-unused-function
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
  LDFLAGS += -lsiloh5
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
################### --- local libs WORKSTATION ---######################

################### --- FIONN at ICHEC.IE ---######################
# Options for fionn.ichec.ie
######################################################################
ifeq ($(MAKE_UNAME), FIONN)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g
    # -wd981 -wd383
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
  LDFLAGS += -Wl,-L${MYLIB}/lib
  LDFLAGS += -lsilo 
  #LDFLAGS += -lcfitsio
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm

endif
################### --- FIONN at ICHEC.IE ---######################


######################## -- JUROPA -- ################################
# Options for the Juelich Supercomputer JUROPA.  Intel Compilers!
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
  LDFLAGS  += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
endif
##################### -- JUROPA -- ###################################



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
  INC += ${MKL_INC} ${MPI_INC}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += ${MKL_SHLIB} ${MPI_CXX_LIB}

endif
##################### -- SuperMUC -- #################################




######################################################################
########## OPTIONS FOR PHALANX ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), phalanx)
  # Special include for phalanx mpi libs:
  INC = -I/opt/sgi/mpt/mpt-1.26/include
  LDFLAGS += -lmpi
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
########## OPTIONS FOR dougal ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), dougal)
  # Special include for dougal mpi libs:
  INC = -I/usr/include/mpich2-x86_64
  LDFLAGS += -L/usr/lib64/mpich2/lib -lmpich -lmpichcxx -lopa  -lpthread   -lrt 
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -Wall -g
  endif
  MYLIB = ../extra_libraries
  INC += -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib
  LDFLAGS += -Wl,--rpath -Wl,/usr/lib64/mpich2/lib 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm
endif
#################### --- DOUGAL --- #################################

###########################################
# Specific Options for OSX
###########################################
ifeq ($(MAKE_UNAME), OSX)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3 -Wno-unused-function
    CXXFLAGS = -Wall -g -Wno-unused-function
  else
    CXXFLAGS = -Wall -g -Wno-unused-function
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I/opt/local/include
  INC += -I${SRCDIR}
  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -L/opt/local/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsiloh5 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  #LDFLAGS += -Wl,-no_compact_unwind

endif
#################### --- OS X Darwin --- #############################


######################################################################
#################### -- DIRAC-COMPLEXITY -- ##########################
# Options for the DIRAC-II Supercomputer Complexity. Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), DIRAC)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
#    CXXFLAGS = -O0 -Wall -g -wd981 -wd383
    CXXFLAGS = -O0 -Wall -g
  endif

  # -wd981 disables a warning about order of evaluation of operands:
  #  see http://software.intel.com/en-us/forums/showthread.php?t=62309
  #  and also http://software.intel.com/en-us/articles/cdiag981/
  # -wd383 also diasables a warning about passing a temporary object to a function.
  CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
  #CXXFLAGS += -fp-stack-check -fstack-security-check

  MYLIB = ../extra_libraries
  INC  = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm

endif
################## -- DIRAC-COMPLEXITY -- ############################



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

OBJ_GEN = reporting.o mem_manage.o command_line_interface.o \
          timer.o interpolate.o constants.o sim_params.o \
          comm_mpi.o comm_files.o

OBJ_IO =  parameter_defs.o file_status.o readparams.o \
          dataio_silo.o dataio_fits.o utility_fits_class.o \
          dataio_base.o dataio_text.o dataio_fits_MPI.o \
          dataio_silo_MPI.o dataio_silo_utility.o

OBJ_GRID =  cell_interface.o uniform_grid.o setup_fixed_grid.o \
            stellar_wind_BC.o stellar_wind_angle.o MCMD_control.o \
            raytracer_SC.o setup_fixed_grid_MPI.o \
            uniform_grid_pllel.o raytracer_SC_pllel.o

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
          stellar_wind_boundaries.o MCMD_boundaries.o \
          periodic_boundaries_MPI.o assign_update_bcs_MPI.o \
          RT_MPI_boundaries.o

OBJ_NG_BCS = assign_update_bcs_NG_MPI.o     \
          NG_coarse_to_fine_boundaries.o \
          NG_fine_to_coarse_boundaries.o \
          NG_MPI_coarse_to_fine_boundaries.o \
          NG_MPI_fine_to_coarse_boundaries.o \


ifneq (,$(findstring CODE_EXT_SBII,$(PION_OPTIONS)))
  OBJ_ICS += StarBench_test.o
endif

ifneq (,$(findstring BBTURBULENCE_CODE_EXT,$(PION_OPTIONS)))
  OBJ_ICS += read_BBurkhart_data.o
endif

ifneq (,$(findstring HARPREETS_CODE_EXT,$(PION_OPTIONS)))
  OBJ_ICS += HD_2D_ShockCloud.o  # HD_3D_ShockCloud.o
endif

OBJ_SIM = sim_control.o sim_control_MPI.o time_integrator.o \
  calc_timestep.o sim_init.o

OBJ_NG_SIM = sim_control_NG.o sim_control_NG_MPI.o \
  setup_grid_NG_MPI.o

OBJECTS_PION = mainMPI.o $(OBJ_SIM) $(OBJ_CHEM) $(OBJ_GEN) \
  $(OBJ_IO) $(OBJ_GRID) $(OBJ_EQNS) $(OBJ_SOLVE) $(OBJ_CRD) \
  $(OBJ_BCS)

OBJECTS_ICGEN = icgen_parallel.o get_sim_info.o $(OBJ_GEN) \
  $(OBJ_GRID) $(OBJ_IO) $(OBJ_CHEM) $(OBJ_ICS) $(OBJ_EQNS) \
  $(OBJ_CRD) $(OBJ_BCS) $(OBJ_SOLVE)

OBJECTS_NG_PION = main_NG_MPI.o $(OBJ_NG_SIM) $(OBJ_SIM) \
  $(OBJ_CHEM) $(OBJ_GEN) $(OBJ_IO) $(OBJ_GRID) $(OBJ_EQNS) \
  $(OBJ_SOLVE) $(OBJ_CRD) $(OBJ_BCS) $(OBJ_NG_BCS)

OBJECTS_NG_ICGEN = icgen_NG_MPI.o get_sim_info.o $(OBJ_GEN) \
  $(OBJ_GRID) $(OBJ_IO) $(OBJ_CHEM) $(OBJ_ICS) $(OBJ_EQNS) \
  $(OBJ_CRD) $(OBJ_NG_ICS) $(OBJ_SOLVE) $(OBJ_BCS)  $(OBJ_NG_BCS)


#####################################################
##### END OF LIST OF OBJECT FILES               #####
#####################################################


##########################
### -- EXECTUTABLES -- ###
##########################
all: ${bindir}/pion_parallel ${bindir}/icgen_parallel
uni: ${bindir}/pion_parallel ${bindir}/icgen_parallel
NG:  ${bindir}/pion_NG_parallel ${bindir}/icgen_NG_parallel

${bindir}/pion_parallel: $(OBJECTS_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_PION) -o ${bindir}/pion_parallel $(LDFLAGS)

${bindir}/icgen_parallel: $(OBJECTS_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_ICGEN) -o ${bindir}/icgen_parallel $(LDFLAGS)

${bindir}/pion_NG_parallel : $(OBJECTS_NG_PION)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_PION) -o ${bindir}/pion_NG_parallel $(LDFLAGS)

${bindir}/icgen_NG_parallel: $(OBJECTS_NG_ICGEN)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS_NG_ICGEN) -o ${bindir}/icgen_NG_parallel $(LDFLAGS)

clean:
	rm -f *.o ${bindir}/pion_parallel ${bindir}/icgen_parallel


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
  assign_update_bcs.h assign_update_bcs_NG.h

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


uniform_grid_pllel.o: uniform_grid_pllel.cpp uniform_grid_pllel.h \
    uniform_grid.h ${GLOBAL} \
    comms.h stellar_wind_angle.h stellar_wind_BC.h \
    VectorOps.h VectorOps_spherical.h

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

NG_MPI_coarse_to_fine_boundaries.o: \
  NG_MPI_coarse_to_fine_boundaries.h NG_coarse_to_fine_boundaries.h \
  NG_MPI_coarse_to_fine_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h

NG_MPI_fine_to_coarse_boundaries.o: \
  NG_MPI_fine_to_coarse_boundaries.h NG_fine_to_coarse_boundaries.h \
  NG_MPI_fine_to_coarse_boundaries.cpp solver_eqn_base.h \
  boundaries.h ${GLOBAL} grid_base_class.h

assign_update_bcs_MPI.o: assign_update_bcs_NG.cpp \
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
  assign_update_bcs.h assign_update_bcs_NG.h

sim_control_MPI.o: sim_control_MPI.cpp sim_control_MPI.h \
  sim_control.h MCMD_control.h ${GLOBAL} comms.h \
  periodic_boundaries_MPI.h assign_update_bcs_MPI.h \
  RT_MPI_boundaries.h setup_fixed_grid.h setup_fixed_grid_MPI.h \
  VectorOps.h

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
file_status.o: file_status.cpp file_status.h

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


#######################
## MPI/files comms classes
########################
comm_mpi.o: comm_mpi.cpp comm_mpi.h comms.h ${GLOBAL}
comm_files.o: comm_files.cpp comm_files.h comms.h ${GLOBAL}


##############################
##### IC GENERATOR FILES #####
##############################

icgen_parallel.o: icgen_parallel.cpp icgen_base.h icgen.h \
  ${GLOBAL} get_sim_info.h dataio_base.h \
  setup_fixed_grid_MPI.h setup_fixed_grid.h \
  uniform_grid_pllel.h grid_base_class.h uniform_grid.h


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
read_BBurkhart_data.o:  read_BBurkhart_data.cpp icgen_base.h icgen.h ${GLOBAL} dataio_fits.h dataio_base.h

HD_2D_ShockCloud.o:     HD_2D_ShockCloud.cpp icgen.h ${GLOBAL}
HD_3D_ShockCloud.o:     HD_3D_ShockCloud.cpp icgen.h ${GLOBAL}
HD_MetalFree_equations.o:  HD_MetalFree_equations.cpp HD_MetalFree.h ${GLOBAL}
HD_MetalFree_integrator.o: HD_MetalFree_integrator.cpp HD_MetalFree.h ${GLOBAL}

###############################################
###############################################

