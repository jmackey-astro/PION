///
/// \file flux_mhd_adiabatic.cpp
/// \author Jonathan Mackey
/// Function Definitions of the adiabatic hydrodynamics flux solver class
///
/// History: 
///  - 2009-10-20 Started on the file
///  - 2009-10-24 Got it working.
///  - 2010-01-15 JM: Added in calculation of Pstar for Roe Flux solver
///  - 2010-01-18 JM: Fixed error in Roe Solver (Pstar not calculated when star state
///      is the left or right state).
///  - 2010-02-19 JM: Today and yesterday I changed some of the error reporting in
///     the Roe solver.  Now it should only complain about the first 1000 negative 
///     pressures and 1000 negative densities in the starred state.  Also the Roe solver()
///     will return 1 on negative pressure and 2 on negative density, so calculate_flux()
///     can then decide what to do about the solution.  Currently if an error is returned
///     I try Sam Falle's solver and use that solution regardless.
///  - 2010-02-19 JM: Removed pointless PtoU()->UtoP() calculation to get Pstar when
///     the resolved state is either the left or right state.
///  - 2010-02-20 JM: Set whether to use Ustar or meanP for the Pstar value to return to 
///     be in an ifdef.  Ustar is not great for stability, and I think meanP will be much
///     better, as well as quicker to calculate.  So I have it unset now.
///  - 2010-09-02 JM: Commented out unnecessary initialisation of
///     variables in the Roe solver.  Identified ways to speed up the
///     eigenvector calculation by avoiding doubling the calculation
///     for the left and right moving waves.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///    Renamed calculate_flux() to inviscid_flux() and moved AV
///    calculation to FV_solver_base class.
/// - 2010.12.07 JM: AVFalle() now uses the passed--in pointer to the
///   flux vector.  I need to get rid of the class variables -- they 
///   are just a recipe for disaster (AV was having no effect at all!)
/// - 2010.12.27 JM: Moved Roe flux solver to own class in Riemann_solvers/
///   Got rid of inherited class flux/left/right/pstar variables.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "flux_mhd_adiabatic.h"
using namespace std;


// ##################################################################
// ##################################################################



