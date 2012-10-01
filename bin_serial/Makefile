# Makefile for my Grid Hydro Code, Jonathan Mackey
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

############################
### BASIC INFO AND PATHS ###
############################

SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
.SUFFIXES: .cpp .o
srcdir = .
bindir = ../bin
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
../source/time_integrators


GLOBAL = global.h functionality_flags.h testing_flags.h

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
  ifeq ($(JMCODE_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -O0 -Wall -g
  endif

  MYLIB = $(CURDIR)/../extra_libraries
  INC = -I${MYLIB}/include
  LDFLAGS =  -lm  -lreadline -lncurses
  # All libraries are in this path:
  LDFLAGS += -L${MYLIB}/lib
  # Libraries for fits, silo
  LDFLAGS += -lcfitsio
  LDFLAGS += -lsilo
  # libs for CVODES
  LDFLAGS +=-Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
endif
######################################################################
# Options for a "standard" workstation.
######################################################################

######################## -- JUROPA -- ################################
# Options for the Juelich Supercomputer JUROPA
# Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), JUROPA)
  ifeq ($(JMCODE_OPTIMISE), HIGH)
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
  INC += -I/usr/local/sundials/2.4.0/include
  LDFLAGS  = -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lmkl_sequential -lmkl_core
  LDFLAGS += -lreadline

endif
##################### -- JUROPA -- ###################################

######################## -- SuperMUC -- ################################
# Options for the LRZ Supercomputer SuperMUC.  Intel Compilers!
######################################################################
ifeq ($(MAKE_UNAME), SUPERMUC)
  ifeq ($(JMCODE_OPTIMISE), HIGH)
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
  INC += -I/usr/local/sundials/2.4.0/include
  INC += ${MKL_INC}
  LDFLAGS  = -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -L/usr/local/sundials/2.4.0/lib
  LDFLAGS += -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += ${MKL_LIB}
  LDFLAGS += -lreadline

endif
##################### -- SuperMUC -- ###################################

######################################################################
########## OPTIONS FOR PHALANX ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), phalanx)
  ifeq ($(JMCODE_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -Wall -g
  endif
  MYLIB = ../extra_libraries
  INC = -I${MYLIB}/include
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm
endif
#################### --- PHALANX --- #################################

######################################################################
########## OPTIONS FOR DOUGAL ICC/ICPC ##############################
######################################################################
ifeq ($(MAKE_UNAME), dougal)
  ifeq ($(JMCODE_OPTIMISE), HIGH)
    CXXFLAGS = -O3
  else
    CXXFLAGS = -Wall -g
  endif
  MYLIB = ../extra_libraries
  INC = -I${MYLIB}/include
  LDFLAGS += -L${MYLIB}/lib -lcfitsio -lsilo 
  LDFLAGS += -Wl,--rpath -Wl,${MYLIB}/lib -lsundials_cvode -lsundials_nvecserial
  LDFLAGS += -lm
endif
#################### --- PHALANX --- #################################



###########################################
# Specific Options for my imac
###########################################
ifeq ($(MAKE_UNAME), imac)
CXXFLAGS = -O3
#CXXFLAGS = -O0
CXXFLAGS += -Wall
CXXFLAGS += -g    # enable debugging with gdb
#CXXFLAGS += -arch x86_64

# CURDIR may only work with GNU Make!  It gives absolute paths, so that
# the executable can always link to shared libs.
MYLIB = /Users/jmackey/Library/JMs_stuff
#MYLIB = $(CURDIR)/../extra_libraries

INC = -I${MYLIB}/include
LDFLAGS  =  -lm  
LDFLAGS +=-lreadline
LDFLAGS += -L${MYLIB}/lib
LDFLAGS += -lcfitsio
LDFLAGS += -lsilo 
LDFLAGS += -lsundials_cvode -lsundials_nvecserial

endif


#############################################
# Specific Options for Stokes - ICHEC cluster
#############################################
ifeq ($(MAKE_UNAME), stokes)
CXX = icpc
CXXFLAGS = -O3
MYLIB = /ichec/work/dsast002b/jm_code/extra_libs

INC  = -I$(MYLIB)/cfitsio/include
INC += -I$(MYLIB)/silo_intel_nohdf/include
LDFLAGS += -L$(MYLIB)/cfitsio/lib/ -lcfitsio
LDFLAGS += -L$(MYLIB)/silo_intel_nohdf/lib/ -lsilo
LDFLAGS += -lmkl_sequential -lmkl_core
endif

#############################################
# Specific Options for ICHEC Blue Gene /L
#############################################
# module load bgl
#############################################
ifeq ($(MAKE_UNAME), bgl)
CXX = xlc++

# The following option allows stuff to compile -- MPI complains without it.
#CXXFLAGS = -DMPICH_IGNORE_CXX_SEEK 
CXXFLAGS = -O3  # optimised code
#CXXFLAGS += -O0
#CXXFLAGS += -Wall
#CXXFLAGS += -g   # for debugging with ddt. 

MYLIB = /ichec/home/users/jmackey/jm_code/extra_libs

# Overwrite linking to get rid of readline, which is not on bg/l
LDFLAGS = -lnss_files -lnss_dns -lresolv -lc -lnss_files -lnss_dns -lresolv
LDFLAGS += -lm

INC  = -I$(MYLIB)/cfitsio_bgl_xlc/include
LDFLAGS += -L$(MYLIB)/cfitsio_bgl_xlc/lib/ -lcfitsio
INC += -I$(MYLIB)/silo_bgl_nohdf_xlc/include
LDFLAGS += -L$(MYLIB)/silo_bgl_nohdf_xlc/lib/ -lsilo
endif

#####################################################
##### LIST OF OBJECT FILES, DIRECTIVES FOR MAKE #####
#####################################################
CXXFLAGS += $(OPT) $(JMCODE_OPTIONS)


OBJECTS = main.o \
        uniform_grid.o cell_interface.o \
        gridMethods.o global.o \
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
	stellar_wind_BC.o microphysics_lowZ.o \
	flux_base.o \
	eqns_hydro_adi_Eint.o \
	flux_hydro_adi_Eint.o \
	solver_eqn_hydro_adi_Eint.o \
	test_equations.o test_integrator.o \
	mp_only_cooling.o cooling_SD93_cie.o hydrogen_recomb_Hummer94.o \
        hydrogen_mp.o hydrogen_photoion.o \
        mp_explicit_H.o cvode_integrator.o mp_implicit_H.o mp_v2_aifa.o

IC_OBJ	= icgen.o \
          uniform_grid.o cell_interface.o \
          global.o readparams.o \
	dataio.o dataio_silo.o get_sim_info.o inside_sphere.o \
	integrator.o cooling.o microphysics.o microphysics_v1.o \
	shock_tube.o radiative_shock.o blast_wave.o \
	photoevaporating_clump.o basic_tests.o shock_cloud.o \
	photoevaporating_random_clumps.o dataio_fits.o utility_fits_class.o \
	VectorOps.o VectorOps_spherical.o photoevaporating_multiclumps.o \
	eqns_base.o eqns_hydro_adiabatic.o eqns_hydro_isothermal.o \
	eqns_mhd_adiabatic.o laser_ablation.o jet.o \
	stellar_wind_BC.o microphysics_lowZ.o \
	test_equations.o test_integrator.o HD_2D_ShockCloud.o \
	spherical_clump.o \
        cooling_SD93_cie.o hydrogen_recomb_Hummer94.o \
        hydrogen_mp.o hydrogen_photoion.o \
        mp_explicit_H.o cvode_integrator.o mp_implicit_H.o mp_v2_aifa.o


##########################
### -- EXECTUTABLES -- ###
##########################
all: ${bindir}/main_serial ${bindir}/icgen_serial

${bindir}/main_serial : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) $(OBJECTS) -o ${bindir}/main_serial $(LDFLAGS)

${bindir}/icgen_serial: $(IC_OBJ)
	$(CXX) $(CXXFLAGS) $(IC_OBJ) -o ${bindir}/icgen_serial $(LDFLAGS)

clean:
	rm -f *.o ${bindir}/main_serial ${bindir}/icgen_serial


####################################
### DEPENDENCIES OF SOURCE FILES ###
####################################


main.o:         main.cc grid.h ${GLOBAL}
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
   microphysics_base.h ${GLOBAL} integrator.h header.h

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



########################
# Harpreet's files
########################

#RK5_functions.o:	header.h RK5_functions.cpp
#RKCK_function.o:	header.h RKCK_function.cpp
test_equations.o:	test_equations.cpp header.h
test_integrator.o: 	test_integrator.cpp header.h
HD_2D_ShockCloud.o:     HD_2D_ShockCloud.cpp icgen.h ${GLOBAL}
########################


raytracer_SC.o: raytracer_SC.cc raytracer_SC.h ${GLOBAL} uniform_grid.h

gridMethods.o:  gridMethods.cc grid.h ${GLOBAL} solver_eqn_base.h \
		microphysics.h microphysics_lowZ.h  mp_v2_aifa.h \
		dataio.h dataio_silo.h dataio_fits.h microphysics_base.h \
                mp_only_cooling.h mp_explicit_H.h mp_implicit_H.h

# functions for stepping forward in time
time_integrator.o: time_integrator.cpp grid.h microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL}
# functions for calculating the timestep
calc_timestep.o: calc_timestep.cpp grid.h microphysics_base.h solver_eqn_base.h raytracer_SC.h ${GLOBAL}


dataio.o:       dataio.cc dataio.h ${GLOBAL} solver_eqn_base.h readparams.h \
		stellar_wind_BC.h
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
#######################



