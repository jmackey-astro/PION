///
/// \file flux_hydro_adiabatic.cc
/// \author Jonathan Mackey
/// Function Definitions of the adiabatic hydrodynamics flux solver class
///
/// ------------------------------------------------------------------
///
/// Wish-List:
///
/// - I think Lapidus viscosity could be in the base flux solver.  Or
///   maybe, since it is the only function here with grid cell pointers,
///   it should be in the solver class built on top of this flux class?
///   Or I can re-define it so that each direction has a div(v) for its
///   own interface, like the H-correction has an etamax?  That way I 
///   can avoid any grid-awareness.
///
/// ------------------------------------------------------------------
///
/// History: \n
///  - 2009-10-20 Started on the file
///  - 2010-01-13 JM: Added Resolved state calculation to Roe Solver (Toro,1999,eq.11.17).
///  - 2010-07-31 JM: Added Primitive variable linear solver with Roe average.
///  - 2010-09-21 JM: Shortened some long lines.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.11.12 JM: Changed ->col to use cell interface for 
///   extra_data.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///    Renamed calculate_flux() to inviscid_flux() and moved AV
///    calculation to FV_solver_base class.
///  - 2010.11.19 JM: Added H-correction to Roe conserved Var solver.
/// - 2010.11.21 JM: ifdef option to use the mean state for Pstar in
///   Roe conserved variable solver (4).  This is only relevant for
///   the FKJ98 viscosity function, and make sure that pstar has
///   positive definite pressure and density.  The ifdef statement is
///   at the top of this file.
/// - 2010.11.21 JM: Added symmetric version of the Roe conserved
///   variables flux solver, called Roe_flux_solver_symmetric().  It
///   is MUCH BETTER than the one-sided one, which I have renamed to
///   Roe_flux_solver_onesided().  For an axi-symmetric blast wave
///   with the H-correction, the one-sided solver developed spikes at
///   90 degree intervals, but the symmetric solver is perfectly
///   clean.  SO, DON'T USE ONE-SIDED ROE FLUX SOLVER ANYMORE!
/// - 2010.12.22 JM: Moved Roe PV and CV solvers to Riemann_solvers/
///   directory.  Got rid of the ROE_CV_MEANP ifdef (it is default).
/// - 2010.12.23 JM: Moved UtoP() etc. from solver to flux-solver.
///   Now all tracer stuff is dealt with here.  Added pstar[] array
///   to inviscid_flux() function since it is no longer inherited.
/// - 2010.12.27 JM: Enclosed Lapidus AV in an ifdef.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "flux_hydro_adiabatic.h"
using namespace std;


// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################




