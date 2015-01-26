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
# - 2010.10.13 JM: Added microphysics_lowZ object file
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
../source/flux_calc:\
../source/spatial_solvers:\
../source/time_integrators:\
../source/contrib:\
../source/tools:\
../source/future


GLOBAL = global.h functionality_flags.h testing_flags.h reporting.h sim_params.h sim_constants.h

ifneq ($(PION_OPTIMISE), HIGH)
  GLOBAL += command_line_interface.h
endif

ifneq (,$(findstring CODE_EXT_HHE,$(PION_OPTIONS)))
  GLOBAL += constants.h
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
# Options for a "standard" workstation.
# make sure cfitsio and silo are compiled in ../extra_libraries/
######################################################################
ifeq ($(MAKE_UNAME), standard)
  ifeq ($(PION_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC = -I${MYLIB}/include
  INC += -I${SRCDIR}
  LDFLAGS +=  -lm
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################
######################################################################


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
####  -- may need hard-coded path to work
######################################################################
ifeq ($(MAKE_UNAME), imac)
  CXXFLAGS = -O3
#CXXFLAGS = -O0
  CXXFLAGS += -Wall
  CXXFLAGS += -g    # enable debugging with gdb
#CXXFLAGS += -arch x86_64

# CURDIR may only work with GNU Make!  It gives absolute paths, so that
# the executable can always link to shared libs.
#MYLIB = /Users/jmackey/Library/JMs_stuff
  MYLIB = $(CURDIR)/../extra_libraries

  INC = -I${MYLIB}/include
    INC += -I${SRCDIR}
  LDFLAGS += -lm  
  LDFLAGS += -L${MYLIB}/lib
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo 
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif
######################################################################
######################################################################



#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################
CXXFLAGS += $(OPT) $(PION_OPTIONS)


OBJECTS = reporting.o mem_manage.o command_line_interface.o \
        main.o \
        uniform_grid.o cell_interface.o \
        sim_control.o global.o \
        time_integrator.o calc_timestep.o \
	VectorOps.o VectorOps_spherical.o \
	eqns_base.o solver_eqn_base.o \
	eqns_hydro_adiabatic.o riemann.o findroot.o \
	Riemann_FVS_hydro.o \
	Roe_Hydro_ConservedVar_solver.o \
	Roe_Hydro_PrimitiveVar_solver.o \
	Roe_MHD_ConservedVar_solver.o \
	flux_hydro_adiabatic.o solver_eqn_hydro_adi.o \
	eqns_hydro_isothermal.o flux_hydro_isothermal.o solver_eqn_hydro_iso.o \
	eqns_mhd_adiabatic.o riemannMHD.o flux_mhd_adiabatic.o \
	solver_eqn_mhd_adi.o \
	integrator.o cooling.o microphysics.o microphysics_v1.o raytracer_SC.o \
	dataio.o readparams.o dataio_silo.o dataio_fits.o utility_fits_class.o \
        file_status.o \
	stellar_wind_BC.o \
	flux_base.o \
	eqns_hydro_adi_Eint.o \
	flux_hydro_adi_Eint.o \
	solver_eqn_hydro_adi_Eint.o \
      	mp_only_cooling.o cooling_SD93_cie.o hydrogen_recomb_Hummer94.o \
        hydrogen_mp.o hydrogen_photoion.o \
        mp_explicit_H.o cvode_integrator.o mp_implicit_H.o mp_v2_aifa.o

IC_OBJ	= reporting.o mem_manage.o command_line_interface.o \
          icgen.o \
          uniform_grid.o cell_interface.o \
          global.o readparams.o file_status.o \
	dataio.o dataio_silo.o get_sim_info.o inside_sphere.o \
	integrator.o cooling.o microphysics.o microphysics_v1.o \
	shock_tube.o radiative_shock.o blast_wave.o \
	photoevaporating_clump.o basic_tests.o shock_cloud.o \
	photoevaporating_random_clumps.o dataio_fits.o utility_fits_class.o \
	VectorOps.o VectorOps_spherical.o photoevaporating_multiclumps.o \
	eqns_base.o eqns_hydro_adiabatic.o eqns_hydro_isothermal.o \
	eqns_mhd_adiabatic.o laser_ablation.o jet.o \
	stellar_wind_BC.o \
	spherical_clump.o \
        cooling_SD93_cie.o hydrogen_recomb_Hummer94.o \
        hydrogen_mp.o hydrogen_photoion.o \
        mp_explicit_H.o cvode_integrator.o mp_implicit_H.o mp_v2_aifa.o \
        StarBench_test.o

ifneq (,$(findstring HARPREETS_CODE_EXT,$(PION_OPTIONS)))
  OBJECTS +=  microphysics_lowZ.o HD_MetalFree_integrator.o \
              HD_MetalFree_equations.o 
  IC_OBJ +=   microphysics_lowZ.o HD_MetalFree_integrator.o \
              HD_MetalFree_equations.o \
              HD_2D_ShockCloud.o 
# HD_3D_ShockCloud.o
endif

ifneq (,$(findstring BBTURBULENCE_CODE_EXT,$(PION_OPTIONS)))
  IC_OBJ +=   read_BBurkhart_data.o
endif

ifneq (,$(findstring CODE_EXT_HHE,$(PION_OPTIONS)))
  OBJECTS += HHe_photoion.o mpv9_HHe.o
  IC_OBJ += HHe_photoion.o mpv9_HHe.o
endif

ifneq (,$(findstring CODE_EXT_SBII,$(PION_OPTIONS)))
  OBJECTS += mpv8_StarBench_heatcool.o
  IC_OBJ +=  mpv8_StarBench_heatcool.o
endif

# New chemistry classes.
  OBJECTS += mpv5_molecular.o mpv6_PureH.o mpv7_TwoTempIso.o
  IC_OBJ  += mpv5_molecular.o mpv6_PureH.o mpv7_TwoTempIso.o

#####################################################
##### END OF LIST OF OBJECT FILES               #####
#####################################################

##########################
### -- EXECTUTABLES -- ###
##########################
all: ${bindir}/pion_serial ${bindir}/icgen_serial

${bindir}/pion_serial : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o ${bindir}/pion_serial $(LDFLAGS)

${bindir}/icgen_serial: $(IC_OBJ)
	$(CXX) $(CXXFLAGS) $(IC_OBJ) -o ${bindir}/icgen_serial $(LDFLAGS)

clean:
	rm -f *.o ${bindir}/pion_serial ${bindir}/icgen_serial


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################

reporting.o:   tools/reporting.cpp ${GLOBAL}
mem_manage.o:  tools/mem_manage.cpp ${GLOBAL}
command_line_interface.o:  tools/command_line_interface.cpp ${GLOBAL}

main.o:         main.cpp sim_control.h ${GLOBAL}
global.o:       global.cc ${GLOBAL}

cell_interface.o: cell_interface.cc cell_interface.h ${GLOBAL}

uniform_grid.o:  uniform_grid.cc uniform_grid.h ${GLOBAL} stellar_wind_BC.h \
		VectorOps.h VectorOps_spherical.h

VectorOps.o:           VectorOps.cc VectorOps.h ${GLOBAL} uniform_grid.h
VectorOps_spherical.o: VectorOps_spherical.cc VectorOps_spherical.h VectorOps.h \
		       ${GLOBAL} uniform_grid.h

#####################
# EQUATIONS CLASSES #
#####################

eqns_base.o:	          eqns_base.cc eqns_base.h ${GLOBAL}

eqns_hydro_adiabatic.o: eqns_hydro_adiabatic.cc eqns_hydro_adiabatic.h \
			eqns_base.h ${GLOBAL}

eqns_hydro_isothermal.o: eqns_hydro_isothermal.cc \
                  eqns_hydro_isothermal.h eqns_base.h ${GLOBAL}

eqns_mhd_adiabatic.o: eqns_mhd_adiabatic.cc eqns_mhd_adiabatic.h \
		      eqns_base.h ${GLOBAL}

eqns_hydro_adi_Eint.o: eqns_hydro_adi_Eint.cc eqns_hydro_adi_Eint.h \
    eqns_hydro_adiabatic.h eqns_base.h ${GLOBAL}

###################
# RIEMANN SOLVERS #
###################

riemann.o: riemann.cc riemann.h findroot.h ${GLOBAL} \
	   eqns_hydro_adiabatic.h  eqns_base.h

findroot.o: findroot.cc findroot.h ${GLOBAL}

Riemann_FVS_hydro.o: Riemann_FVS_hydro.cc Riemann_FVS_hydro.h \
   ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

Roe_Hydro_ConservedVar_solver.o: Roe_Hydro_ConservedVar_solver.cc \
   Roe_Hydro_ConservedVar_solver.h \
   ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

Roe_Hydro_PrimitiveVar_solver.o: Roe_Hydro_PrimitiveVar_solver.cc \
   Roe_Hydro_PrimitiveVar_solver.h \
   ${GLOBAL} eqns_hydro_adiabatic.h eqns_base.h

riemannMHD.o: riemannMHD.cc riemannMHD.h ${GLOBAL} \
              eqns_mhd_adiabatic.h eqns_base.h ${GLOBAL} 

Roe_MHD_ConservedVar_solver.o : Roe_MHD_ConservedVar_solver.cc \
   Roe_MHD_ConservedVar_solver.h riemannMHD.h \
   ${GLOBAL} eqns_mhd_adiabatic.h eqns_base.h

####################
# FLUX CALCULATORS #
####################
flux_base.o : flux_base.cc flux_base.h eqns_base.h \
              ${GLOBAL} cell_interface.h

flux_hydro_adiabatic.o: flux_hydro_adiabatic.cc flux_hydro_adiabatic.h \
   eqns_hydro_adiabatic.h flux_base.h eqns_base.h riemann.h \
   Roe_Hydro_ConservedVar_solver.h Roe_Hydro_PrimitiveVar_solver.h \
   Riemann_FVS_hydro.h ${GLOBAL} cell_interface.h

flux_hydro_isothermal.o: flux_hydro_isothermal.cc flux_hydro_isothermal.h \
   eqns_hydro_isothermal.h flux_base.h eqns_base.h ${GLOBAL} cell_interface.h

flux_mhd_adiabatic.o: flux_mhd_adiabatic.cc flux_mhd_adiabatic.h \
   eqns_mhd_adiabatic.h riemannMHD.h flux_base.h eqns_base.h \
   ${GLOBAL} cell_interface.h

flux_hydro_adi_Eint.o: flux_hydro_adi_Eint.cc flux_hydro_adi_Eint.h \
   eqns_hydro_adi_Eint.h \
   flux_hydro_adiabatic.h eqns_hydro_adiabatic.h \
   flux_base.h eqns_base.h \
   riemann.h Roe_Hydro_PrimitiveVar_solver.h \
   ${GLOBAL} cell_interface.h

###########
# SOLVERS #
###########

solver_eqn_base.o: VectorOps.h solver_eqn_base.h flux_base.h \
                   eqns_base.h ${GLOBAL} cell_interface.h

solver_eqn_hydro_adi.o:	VectorOps.h VectorOps_spherical.h \
   solver_eqn_base.h flux_base.h eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_hydro_adi.h flux_hydro_adiabatic.h  eqns_hydro_adiabatic.h \
   solver_eqn_hydro_adi.cc

solver_eqn_hydro_iso.o: VectorOps.h solver_eqn_base.h eqns_base.h \
   flux_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_hydro_iso.h flux_hydro_isothermal.h eqns_hydro_isothermal.h \
   solver_eqn_hydro_iso.cc

solver_eqn_mhd_adi.o: VectorOps.h \
   solver_eqn_base.h flux_base.h eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_mhd_adi.h flux_mhd_adiabatic.h  eqns_mhd_adiabatic.h \
   solver_eqn_mhd_adi.cc

solver_eqn_hydro_adi_Eint.o:	VectorOps.h VectorOps_spherical.h \
   solver_eqn_hydro_adi_Eint.h solver_eqn_hydro_adi_Eint.cc \
   solver_eqn_base.h flux_base.h eqns_base.h ${GLOBAL} cell_interface.h \
   solver_eqn_hydro_adi.h flux_hydro_adiabatic.h  eqns_hydro_adiabatic.h \
   solver_eqn_hydro_adi_Eint.h solver_eqn_hydro_adi_Eint.cc


######################
# MICROPHYSICS STUFF #
######################

integrator.o:	integrator.cc integrator.h ${GLOBAL}
cooling.o:	cooling.cc cooling.h ${GLOBAL}

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
#header.h

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

HHe_photoion.o: HHe_photoion.cpp HHe_photoion.h ${GLOBAL}
mpv9_HHe.o: mpv9_HHe.cpp mpv9_HHe.h microphysics_base.h \
  HHe_photoion.h cvode_integrator.h ${GLOBAL}

mpv8_StarBench_heatcool.o: mpv8_StarBench_heatcool.cpp \
  mpv8_StarBench_heatcool.h mp_explicit_H.h microphysics_base.h ${GLOBAL}


########################
# Raytracing
########################
raytracer_SC.o: raytracer_SC.cc raytracer_SC.h ${GLOBAL} uniform_grid.h

sim_control.o:  sim_control.cpp sim_control.h ${GLOBAL} solver_eqn_base.h \
		microphysics.h microphysics_lowZ.h  mp_v2_aifa.h \
		dataio.h dataio_silo.h dataio_fits.h microphysics_base.h \
                mp_only_cooling.h mp_explicit_H.h mp_implicit_H.h

# functions for stepping forward in time
time_integrator.o: time_integrator.cpp sim_control.h microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL}
# functions for calculating the timestep
calc_timestep.o: calc_timestep.cpp sim_control.h microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL}

########################
# Data I/O
########################
file_status.o:  file_status.cpp file_status.h ${GLOBAL}
dataio.o:       dataio.cc dataio.h ${GLOBAL} solver_eqn_base.h readparams.h \
		stellar_wind_BC.h cell_interface.h
dataio_silo.o:  dataio_silo.cc dataio_silo.h dataio.h ${GLOBAL} solver_eqn_base.h
dataio_fits.o:  dataio_fits.cc dataio_fits.h dataio.h ${GLOBAL} solver_eqn_base.h \
        solver_eqn_hydro_adi.h solver_eqn_hydro_iso.h

readparams.o:   readparams.cc readparams.h ${GLOBAL}
utility_fits_class.o: utility_fits_class.cc dataio_fits.h dataio.h ${GLOBAL}

stellar_wind_BC.o: stellar_wind_BC.cc stellar_wind_BC.h ${GLOBAL}

##############################
##### IC GENERATOR FILES #####
##############################
icgen.o: icgen.cc icgen.h uniform_grid.h ${GLOBAL} microphysics_base.h  mp_v2_aifa.h \
         microphysics.h microphysics_lowZ.h
get_sim_info.o: get_sim_info.h get_sim_info.cc ${GLOBAL} readparams.h stellar_wind_BC.h
inside_sphere.o: inside_sphere.cc inside_sphere.h ${GLOBAL}
########################
shock_tube.o: shock_tube.cc icgen.h ${GLOBAL}
radiative_shock.o: radiative_shock.cc icgen.h ${GLOBAL}
blast_wave.o:      blast_wave.cc icgen.h ${GLOBAL} inside_sphere.h
photoevaporating_clump.o:  photoevaporating_clump.cc icgen.h ${GLOBAL} inside_sphere.h
basic_tests.o: basic_tests.cc icgen.h ${GLOBAL} 
photoevaporating_random_clumps.o: photoevaporating_random_clumps.cc icgen.h ${GLOBAL}
shock_cloud.o: shock_cloud.cc icgen.h ${GLOBAL} inside_sphere.h
photoevaporating_multiclumps.o: photoevaporating_multiclumps.cc icgen.h ${GLOBAL}
laser_ablation.o:   laser_ablation.cc icgen.h ${GLOBAL}
jet.o:		    jet.cc icgen.h ${GLOBAL}
spherical_clump.o: spherical_clump.cc icgen.h ${GLOBAL}
StarBench_test.o: StarBench_test.cpp icgen.h ${GLOBAL}
#######################

###############################################
## Contributed and/or project-specific code. ##
###############################################
read_BBurkhart_data.o:  read_BBurkhart_data.cpp icgen.h ${GLOBAL} dataio_fits.h dataio.h

HD_2D_ShockCloud.o:     HD_2D_ShockCloud.cpp icgen.h ${GLOBAL}
HD_3D_ShockCloud.o:     HD_3D_ShockCloud.cpp icgen.h ${GLOBAL}
HD_MetalFree_equations.o:  HD_MetalFree_equations.cpp HD_MetalFree.h ${GLOBAL}
HD_MetalFree_integrator.o: HD_MetalFree_integrator.cpp HD_MetalFree.h ${GLOBAL}
###############################################
###############################################


