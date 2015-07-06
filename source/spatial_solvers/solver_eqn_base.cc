/// \file solver_eqn_base.cc
/// 
/// \brief Class definition for various solvers implemented in my FV grid-code.
/// \author Jonathan Mackey
/// 
/// Modifications:
///  - 2007-07-10 Part Way through writing it.
///  - 2007-07-11 Still writing it and working out class heirarchy.
///  - 2007-07-12 Got the class heirarchy working (I think...).
///  - 2007-07-13 New Class structure implemented.
///  - 2007-07-16 Reverted to less complicated class hierarchy.  Runs fast.
///  - 2007-07-23 Started to add passive tracer variables
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08-01 cylindrical coordinates support.
///  - 2007-08-08 cylindrical coordinates hd/i-mhd/glm-mhd working.
///  - 2007-11-30 Worked on Toth's Field-CD method for keeping divB=0.
///  - 2009-10-20 New structure built on equations and flux classes...
///  - 2009-10-24 Cut out all the old solver classes.
///  - 2010.09.30 JM: Worked on Lapidus AV (added set_div_v() function for cells)
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.  Added Pre- and 
///   Post-flux viscosity functions for the H-correction and Lapidus
///   viscosity functions.
///   Added H-correction functions for calculating eta and getting
///   the maximum value of eta on a given H-stencil.
///   Made InterCellFlux general for all classes.
/// - 2010.11.19 JM: Debugged.
/// - 2010.12.08 JM: Got the H-correction to work for 1D grids.
/// - 2010.12.23 JM: Added new set_interface_tracer_flux(lp,rp,flux);
///   function.  Removed riemann_base calls.  Added pstar[] array
///   to intercell_flux() function since it is no longer inherited.
/// - 2010.12.27 JM: modified args to post_calc_viscous_terms to pass
///   in left, right, pstar (since these are no longer class vars).
/// - 2011.01.03 JM: Moved preprocess_data() and calc_Hcorrection from
///   gridMethods.cc
/// - 2011.02.25 JM: removed HCORR ifdef around new code
/// - 2011.04.06 JM: Added thermal-conduction to preprocess_data.
/// - 2012.08.05 JM: Added spacers between functions for readability.
///    Thermal Conduction module is unusable now, needs re-writing.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"

#include "solver_eqn_base.h"
using namespace std;






// ##################################################################
// ##################################################################


FV_solver_base::FV_solver_base(const int nv, ///< number of variables in state vector.
			       const int nd, ///< number of space dimensions in grid.
			       const double cflno,   ///< CFL number
			       const double delx,    ///< dx, cell size.
			       const double gam,     ///< gas eos gamma.
			       const double avcoeff, ///< Artificial Viscosity Parameter etav.
			       const int ntr         ///< Number of tracer variables.
			       )
  : eqns_base(nv),
    // riemann_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_gndim(nd), FV_cfl(cflno), FV_dx(delx)
{
  //
  // Can't hurt to set this everywhere...
  //
  eq_gamma = gam;
  //
  // Nothing else to do really...
  //
  return;
}

FV_solver_base::~FV_solver_base()
{
  return;
}




// ##################################################################
// ##################################################################


int FV_solver_base::get_LaxFriedrichs_flux(const double *l,
					   const double *r,
					   double *f,
					   const double
					   )
{
  //
  // This has to be in this solver because we need dx,dt,ndim
  //
  double u1[eq_nvar], u2[eq_nvar], f1[eq_nvar], f2[eq_nvar];
  PtoU(l, u1, eq_gamma);   PtoU(r, u2, eq_gamma);
  UtoFlux(u1, f1, eq_gamma);  UtoFlux(u2, f2, eq_gamma);
  // Then get inter-cell flux from this.
  for (int v=0;v<eq_nvar;v++)
    f[v] = 0.5*(f1[v]+f2[v] +FV_dx/FV_dt*(u1[v]-u2[v])/FV_gndim);

  //
  // Calculate tracer flux based on whether flow is to left or right.
  //
  if (FS_ntr>0) {
    if (f[eqRHO]>=0.) {
      for (int t=0;t<FS_ntr;t++)
	f[eqTR[t]] = l[eqTR[t]]*f[eqRHO];
    }
    else {
      for (int t=0;t<FS_ntr;t++)
	f[eqTR[t]] = r[eqTR[t]]*f[eqRHO];
    }
  }
  return 0;
}

void FV_solver_base::set_div_v(
      class cell *c,
      class GridBaseClass *grid
      )
{
  int indices[MAX_DIM];
  indices[0]=eqVX; indices[1]=eqVY; indices[2]=eqVZ;
  CI.set_DivV(c, Divergence(c,1,indices, grid));
  return;
}




// ##################################################################
// ##################################################################


///
/// Calculate Flux between a left and right state.
///
int FV_solver_base::InterCellFlux(
        class GridBaseClass *grid,
        const cell *Cl, ///< Left state cell pointer
        const cell *Cr, ///< Right state cell pointer
        double *lp, ///< Left Primitive State Vector.
        double *rp, ///< Right Primitive State Vector.
        double *f, ///< Flux Vector. (written to).
        const int solve_flag,    ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
        const int av_flag,    ///< Viscosity Flag (0=none,1=Falle's,2=Lapidus(broken),etc.)
        const double g, ///< gas EOS gamma.
        const double dx ///< Cell size dx.
        )
{
  //cout <<"FV_solver_base::InterCellFlux() gamma="<<eq_gamma<<" and passed in was g="<<g<<"\n";
  eq_gamma=g;
  double pstar[eq_nvar];

  //
  // Pre-calcualate anything needed for the viscosity (H-correction).
  // Data is stored in each cell, so no values returned.
  //
  pre_calc_viscous_terms(grid,Cl,Cr,av_flag);


  //
  // Get the flux from the flux solver:
  // I CAN GET RID OF CL,CR,G FROM THIS.
  //
  int err = inviscid_flux(Cl,Cr,lp, rp, f, pstar, solve_flag, g);

  //
  // Post-calculate anthing needed for the viscosity (calls the FKJ98
  // or the Lapidus viscosity functions which act after the flux has
  // been calculated).
  //
  post_calc_viscous_terms(Cl,Cr,lp,rp,pstar,f,av_flag);

  //
  // Calculate tracer flux based on whether flow is to left or right.
  //
  set_interface_tracer_flux(lp,rp,f);

  return err;
}




// ##################################################################
// ##################################################################


void FV_solver_base::pre_calc_viscous_terms(
        class GridBaseClass *grid,
        const cell *cl, ///< left-of-interface cell
        const cell *cr, ///< right-of-interface cell
        const int av_flag ///< what kind of AV?
        )
{
  //
  // Only the H-correction acts before the actual flux calculation:
  //
  switch (av_flag) {
  case AV_HCORRECTION: // H-correction
  case AV_HCORR_FKJ98: // H-correction +FKJ98 viscosity
    flux_solver_base::HC_etamax = select_Hcorr_eta(cl,cr, grid);
    break;
  default:
    // Just silently continue if not doing the H-correction.
    break;
  }

  return;
}




// ##################################################################
// ##################################################################



void FV_solver_base::post_calc_viscous_terms(const cell *cl, ///< left-of-interface cell
					     const cell *cr, ///< right-of-interface cell
					     const double *Pl,
					     const double *Pr,
					     const double *Pstar,
					     double *flux,   ///< flux vector
					     const int av_flag ///< what kind of AV?
					     )
{
  //  cout <<"etav="<<FS_etav<<"\t";  rep.printVec("flux",flux,eq_nvar);
  //  rep.printVec("flux",flux,eq_nvar);

  int err=0;
  switch (av_flag) {
  case AV_FKJ98_1D:    // FKJ98 Viscosity
  case AV_HCORR_FKJ98: // FKJ98+HCORR
    //cout <<"FKJ98 AV being applied!\n";
    err += AVFalle(Pl,Pr,Pstar,flux,FS_etav,eq_gamma);
    break;

#ifdef LAPIDUS_VISCOSITY_ENABLED
  case AV_LAPIDUS: // Broken Lapidus Viscosity
    err += AVLapidus(cl,cr,flux,FS_etav,eq_gamma);
    break;
#endif // LAPIDUS_VISCOSITY_ENABLED

  default:
    // Silently continue if av_flag==0 or 3.
    break;
  }

  //  rep.printVec("flux",flux,eq_nvar);
  return;
}




// ##################################################################
// ##################################################################


//
// Multi-dimensional calculations to be performed on every cell before
// the flux calculations.
//
int FV_solver_base::preprocess_data(
        const int csp,
        const int ctm,
        class GridBaseClass *grid
        )
{
  //  cout <<"\t\t\tpreprocess_data(): Starting: ndim = "<<SimPM.ndim<<"\n";
  int err=0;

#ifdef THERMAL_CONDUCTION
  //
  // If we are on the first half-step (or first order integration) we don't need
  // to calculate Edot, because it was already done in calc_dt().  But on the 
  // second half step we need to calculate it here.
  //
  if (ctm != OA1) {
    //cout <<"\tFV_solver_base::preprocess_data: setting Edot for 2nd step.\n";
    err = set_thermal_conduction_Edot();
    if (err) {
      rep.error("FV_solver_base::preprocess_data: set_thermal_conduction_Edot()",err);
    }
  }
  //
  // We need to multiply dU[ERG] by dt to convert from Edot to Delta-E
  // TODO: Wrap this in an if statement when I get the SimPM.EP.conduction
  //       flag integrated into the code.
  //
  class cell* c = grid->FirstPt();
  //cout <<"\tmultiplying conduction dU by dt.\n";
  do {
    c->dU[ERG] *= SimPM.dt;
    //if (c->id==39783) {
    //  cout <<"\tPPD cell id: "<<c->id;
    //  cout <<", dU[cond]="<<c->dU[ERG]<<", E="<<c->Ph[PG]/(SimPM.gamma-1.0)<<"\n";
    //}
  } while ( (c =grid->NextPt(c)) !=0);
#endif // THERMAL CONDUCTION

  //
  // If we have a Lapidus-type viscosity, then we need to save div(v)
  // for every cell.
  //
  if (SimPM.artviscosity==AV_LAPIDUS) {
#ifdef TEST_LAPIDUS
    //
    // TODO: NEED FIRST BOUNDARY CELL VALUES ALSO!
    //
    cout <<"\t\t Lapidus viscosity: Calculating div(v).\n";
    if (SimPM.artviscosity==AV_LAPIDUS) {
      class cell* c = grid->FirstPt();
      do {
	set_div_v(c);
      } while ( (c =grid->NextPt(c)) !=0);
    }
    //
#else
    cout <<"\t\t*** ASKED FOR LAPIDUS VISC, BUT IFDEF IS NOT ENABLED.\n";
#endif //TEST_LAPIDUS
  }
  //
  // For the H-correction we need a maximum speed in each direction
  // for every cell.  Note this calculation is very time-consuming
  // because all of the slopes and edge-states must be calculated.
  //
  else if (SimPM.artviscosity==AV_HCORRECTION ||
	   SimPM.artviscosity==AV_HCORR_FKJ98) {
    err += calc_Hcorrection(csp,ctm, grid);
  }
  return err;
}



// ##################################################################
// ##################################################################


int FV_solver_base::calc_Hcorrection(
        const int csp,
        const int ctm,
        class GridBaseClass *grid
        )
{
  //  cout <<"\t\t\tcalc_Hcorrection() ndim = "<<SimPM.ndim<<"\n";

  //
  // This function is quite similar to calc_dU() and dU_column()
  // because it steps through the grid in the same way.
  //

  // 
  // Allocate arrays
  //
  int err=0;
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP; posdirs[1] = YP; posdirs[2] = ZP;
  negdirs[0] = XN; negdirs[1] = YN; negdirs[2] = ZN;
  axis[0] = XX; axis[1] = YY; axis[2] = ZZ;
  //
  // Slope and edge state temporary arrays: This could be more
  // computationally efficient if these were cell members (i.e. if
  // each cell had a slope vector for each direction), but obviously
  // this would increase the memory overhead hugely.
  //
  double *slope_cpt=0, *slope_npt=0, *edgeR=0, *edgeL=0, *temp=0;
  slope_cpt = mem.myalloc(slope_cpt, SimPM.nvar);
  slope_npt = mem.myalloc(slope_npt, SimPM.nvar);
  edgeL     = mem.myalloc(edgeL,     SimPM.nvar);
  edgeR     = mem.myalloc(edgeR,     SimPM.nvar);

  //
  // Loop through each direction.
  //
  for (int idim=0;idim<SimPM.ndim;idim++) {
    //    cout <<"\t\t\tidim="<<idim<<"\n";
    SetDirection(axis[idim]);

    //
    // Want to start at first boundary cell eventually, but for now
    // stick to the grid cells.
    //
    class cell *start  = grid->FirstPt();
    class cell *marker = grid->FirstPt();

    //
    // Loop over z-planes (there must be at least one!)
    //
    bool zplanes_finished = false;
    bool xcolumns_finished = false;
    
    do {
      //
      // Loop over x-columns in the y-direction (at least one!)
      //
      xcolumns_finished = false;
      do {
	// --------------------------------------------------------
	// Calculate the H-correction coefficients for this column:
	// Start at the outermost boundary cell, and go to the end.
	// We need all interfaces between and including the
	// grid-boundary interface.
	// --------------------------------------------------------

	//
	// Set three cell pointers (2nd order slopes have a 3-point
	// stencil).
	//
	cell *cpt = start;
	while (grid->NextPt(cpt,negdirs[idim])) {
	  cpt = grid->NextPt(cpt,negdirs[idim]);
	}	
	cell *npt  = grid->NextPt(cpt,posdirs[idim]);
	cell *n2pt = grid->NextPt(npt,posdirs[idim]);
	if (npt==0 || n2pt==0)
	  rep.error("Couldn't find two real cells in column",0);

	// 
	// Need to get slopes and edge states if 2nd order (csp==OA2).
	//
	for (int v=0;v<SimPM.nvar;v++) { 
	  slope_cpt[v] = 0.;
	  edgeL[v]     = 0.;
	} // slope_npt[] and edgeR[] get initialised in next loop.
	
	//
	// Run through column, calculating slopes, edge-states, and
	// eta[] values as we go.
	//
	do {
	  err += SetEdgeState(cpt, posdirs[idim], SimPM.nvar, slope_cpt, edgeL, csp, grid);
	  err += SetSlope(npt, axis[idim], SimPM.nvar, slope_npt, csp, grid);
	  err += SetEdgeState(npt, negdirs[idim], SimPM.nvar, slope_npt, edgeR, csp, grid);
	  set_Hcorrection(cpt, axis[idim], edgeL, edgeR, SimPM.gamma);
	  //cout <<" Hcorr["<<axis[idim]<<"] = "<<CI.get_Hcorr(cpt,axis[idim])<<"\n";

	  //
	  // Set npt slope to cpt slope, set n2pt to npt, npt to cpt.
	  // Move to next cell.
	  //
	  temp = slope_cpt;
	  slope_cpt = slope_npt;
	  slope_npt = temp;
	  cpt = npt; npt = n2pt;
	}  while ( (n2pt=grid->NextPt(n2pt,posdirs[idim])) !=0);

	//
	// If 1st order, cpt is still a grid cell (2nd order we are done)
	//
	err += SetEdgeState(cpt, posdirs[idim], SimPM.nvar, slope_cpt, edgeL, csp, grid);
	for (int v=0;v<SimPM.nvar;v++) slope_npt[v] = 0.; // last cell must be 1st order.
	err += SetEdgeState(npt, negdirs[idim], SimPM.nvar, slope_npt, edgeR, csp, grid);
	set_Hcorrection(cpt, axis[idim], edgeL, edgeR, SimPM.gamma);
	//cout <<" Hcorr["<<axis[idim]<<"] = "<<CI.get_Hcorr(cpt,axis[idim])<<"\n";
	
	// --------------------------------------------------------
	// Finished H-correction calculation for the column.
	// --------------------------------------------------------

	//
	// Get next x-column, or if it doesn't exist set a flag to
	// indicate that we are finished.
	//
	if (SimPM.ndim==1) xcolumns_finished = true;
	else {
	  start = grid->NextPt(start,posdirs[(idim+1)%SimPM.ndim]);
	  if (!start || !start->isgd) xcolumns_finished = true;
	}
      } while (!xcolumns_finished);
	
      //
      // Get next z-plane, or if it doesn't exist set a flag to
      // indicate that we are finished.  Reset marker to the first
      // cell in the new plane.
      //
      if (SimPM.ndim<=2) zplanes_finished = true;
      else {
	start = grid->NextPt(marker,posdirs[(idim+2)%SimPM.ndim]);
	marker = start;
	if (!start || !start->isgd) zplanes_finished = true;
      }
    } while (!zplanes_finished);

  } // Loop over Ndim directions.
  SetDirection(axis[0]); // Reset fluxes to x-dir, (just to be safe!).

  slope_cpt = mem.myfree(slope_cpt);
  slope_npt = mem.myfree(slope_npt);
  edgeL     = mem.myfree(edgeL);
  edgeR     = mem.myfree(edgeR);
  return err;
} // calc_Hcorrection()




// ##################################################################
// ##################################################################


void FV_solver_base::set_Hcorrection(cell *c, ///< cell to operate on
				     const axes axis, ///< axis normal to interface.
				     const double *edgeL, ///< Left state
				     const double *edgeR, ///< right state
				     const double g       ///< gamma
				     )
{
  if (axis != GetDirection()) {
    cout <<GetDirection()<<"\t";
    rep.error("bad direction in FV_solver_base::set_Hcorrection()",axis);
  }
  //
  // Sanders, Morano, Druguet, (1998, JCP, 145, 511)  eq. 10
  //
  double eta = 0.5*(fabs(edgeR[eqVX]-edgeL[eqVX]) + 
		    fabs(maxspeed(edgeR,g)-maxspeed(edgeL,g)) );

  CI.set_Hcorr(c,axis,eta);
  return;
}




// ##################################################################
// ##################################################################


double FV_solver_base::select_Hcorr_eta(
        const cell *cl, ///< cell to left of interface
        const cell *cr, ///< cell to right of interface
        class GridBaseClass *grid
        )
{
  //
  // Based on Sanders, Morano, Druguet, (1998, JCP, 145, 511)
  // Fig. 9 and Eq. 16 (with obvious extension to 3D).
  //

  //
  // First get the current direction.
  //
  double eta = 0.0;
  enum axes axis = GetDirection();
  //
  // First the interface between left and right cells.
  //
  eta = CI.get_Hcorr(cl,axis);

  //
  // If we are in 1D this is the only interface, so we can return.
  // Otherwise we have another 4 or 8 more interfaces to calculate.
  //
  if (FV_gndim==1) return eta;

  //
  // Now the two in the positive direction of the 1st perpendicular
  // axis.  (i,j+1/2) and (i+1,j+1/2)
  //
  enum axes perp = static_cast<axes>((static_cast<int>(axis) +1)%FV_gndim);
  eta = max(eta, CI.get_Hcorr(cl,perp));
  eta = max(eta, CI.get_Hcorr(cr,perp));
  //
  // Positive direction of the 2nd perp. axis, if 3D grid:
  // (i,j,k+1/2) (i+1,j,k+1/2)
  //
  if (FV_gndim>2) {
    perp = static_cast<axes>((static_cast<int>(axis) +2)%FV_gndim);
    eta = max(eta, CI.get_Hcorr(cl,perp));
    eta = max(eta, CI.get_Hcorr(cr,perp));
  }
  //
  // Negative direction of the (up to) two perp. axes.  We need to
  // make sure the cells exist, and if they don't we just ignore the
  // non-existent interface.  It's up to the grid to have all the
  // cells it needs, with the correct connectivity.
  //
  enum direction negdir = NO;
  cell *cneg=0;
  for (int idim=1; idim < FV_gndim; idim++) {
    perp = static_cast<axes>((static_cast<int>(axis) +idim)%FV_gndim);
    negdir = static_cast<direction>(static_cast<int>(axis)*2);
    cneg = grid->NextPt(cl,negdir);
    if (cneg)
      eta = max(eta, CI.get_Hcorr(cneg,perp));
    cneg = grid->NextPt(cr,negdir);
    if (cneg)
      eta = max(eta, CI.get_Hcorr(cneg,perp));
  }
  
  //
  // Will want to comment this out later...
  //
  //  cout <<"cell id="<<cl->id<<" axis="<<axis<<", eta_max="<<eta<<"\n";
  
  return eta;
}




// ##################################################################
// ##################################################################


#ifdef THERMAL_CONDUCTION
int FV_solver_base::set_thermal_conduction_Edot()
{
  //
  // This function is quite similar to calc_dU() and dU_column()
  // because it steps through the grid in the same way.
  //
  //cout <<"IntUniformFV::calc_conduction_dt_and_Edot()\n\tcalculating Temperature.\n";
  //
  // First we need to calculate the temperature in every cell.  We do
  // this first for grid cells, and then we will calculate it for boundary
  // data later as they are encountered.  Temperature is stored in dU[RHO] --
  // We will reset it to zero at the end of this function.
  //
  if (!MP) {
    rep.error("Why do conductivity without having Microphysics?",MP);
  }
  cell *c = grid->FirstPt();
  do {
    //cout <<"dU[RHO]="<<c->dU[RHO];
    c->dU[RHO] = MP->Temperature(c->Ph,SimPM.gamma);
    //cout <<", replaced with T="<<c->dU[RHO]<<"\n";
  } while ((c=grid->NextPt(c)) !=0);
  
  //cout <<"\tT calculated, now calculating divQ.\n";
  // 
  // Allocate arrays
  //
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP; posdirs[1] = YP; posdirs[2] = ZP;
  negdirs[0] = XN; negdirs[1] = YN; negdirs[2] = ZN;
  axis[0] = XX; axis[1] = YY; axis[2] = ZZ;  
  double q_neg=0.0, q_pos=0.0, gradT=0.0, Qclassical=0.0, Qsaturated=0.0, T=0.0;
  double dx=grid->DX();

  //
  // Loop through each direction.
  //
  for (int idim=0;idim<SimPM.ndim;idim++) {
    //cout <<"\t\t\tidim="<<idim<<"\n";

    //
    // Start at first grid cells in each column, and we will move back
    // to the first boundary cell to get the correct edge-cell in/out
    // heat flux in each direction.
    //
    class cell *start  = grid->FirstPt();
    class cell *marker = grid->FirstPt();

    //
    // Loop over z-planes (there must be at least one!)
    //
    bool zplanes_finished = false;
    bool xcolumns_finished = false;
    
    do {
      //
      // Loop over x-columns in the y-direction (at least one!)
      //
      xcolumns_finished = false;
      do {
	// --------------------------------------------------------
	// Calculate the Heat fluxes coefficients for this column:
	// Start at the outermost boundary cell, and go to the end.
	// We need all interfaces between and including the
	// grid-boundary interface.
	// --------------------------------------------------------
	cell *cpt = start;
	while (grid->NextPt(cpt,negdirs[idim])) {
	  cpt = grid->NextPt(cpt,negdirs[idim]);
	}	
	cell *npt  = grid->NextPt(cpt,posdirs[idim]);
	if (npt==0)
	  rep.error("Couldn't find two cells in column",0);

        q_neg = 0.0; // no flux coming in from non-existent boundary data.
        q_pos = 0.0;

	//
	// Run through column, calculating slopes, edge-states, and
	// eta[] values as we go.
	//
	do {
          //
          // Check if we have boundary data, b/c we need to set T for these cells.
          //
          if (!cpt->isgd) {
            //cout <<"\tBoundary cpt dU[RHO]="<<cpt->dU[RHO];
            cpt->dU[RHO] = MP->Temperature(cpt->Ph,SimPM.gamma);
            //cout <<", replaced with T="<<cpt->dU[RHO]<<"\n";
          }
          if (!npt->isgd) {
            //cout <<"\tBoundary npt dU[RHO]="<<npt->dU[RHO];
            npt->dU[RHO] = MP->Temperature(npt->Ph,SimPM.gamma);
            //cout <<", replaced with T="<<npt->dU[RHO]<<"\n";
          }
          //
          // Now use the Slavin & Cox (1992) formula for conduction to get
          // the conductive heat flux from cpt to npt in direction posdir[idim].
          // We use the more opaque function idifference_cell2cell() since it gives
          // the correct distance between centres-of-volume of cells on curvilinear
          // grids.
          //
          gradT = (npt->dU[RHO]-cpt->dU[RHO])/(grid->idifference_cell2cell(cpt,npt,axis[idim])*CI.phys_per_int());
          //cout <<"\tT2="<<npt->dU[RHO]<<", T1="<<cpt->dU[RHO]<<", grad(T)="<<gradT<<"\n";
          //
          // If flow is from npt to cpt, we use npt's values for calculating Q.
          // Else we use cpt's values. (note if gradT>0, then T2>T1, flow from 2->1
          // in the *negative* direction).
          //
          if (gradT>0.0) c=npt;
          else           c=cpt;
          //
          // First we get ln(Lambda) and then the classical and saturated fluxes.
          // For ln(Lambda) the formula is only valid for T>4.2e5.  The value of
          // 4.2735e23 is (1.4*m_p)^{-1}.
          //
          T = c->dU[RHO];
          if (T<4.2e5) Qclassical = 29.7;
          else         Qclassical = 29.7+log(T/(1.0e6*sqrt(c->Ph[RO]*4.2735e23)));
          Qclassical = -1.84e-5*pow(T,2.5)*gradT/Qclassical;
          //
          // For saturated Q we follow S&C(1992) and use phi_s=0.3
          //
          Qsaturated = -1.5*pow(c->Ph[PG],1.5)/sqrt(c->Ph[RO]);
          if (gradT<0.0) Qsaturated *= -1.0;
          //
          // now Q = Qs*(1-exp(-Qc/Qs)).   (Qs>>Qc)=>(Q->Qc). (Qc>>Qs)=>(Q->Qs).
          //
          q_pos = Qsaturated*(1.0-exp(-Qclassical/Qsaturated));
          //
          // Finally cpt needs an updated -div(q) value from the current direction.
          // This is a hack, because there is no VectorOps function to do this for me.
          // I should write a function in VectorOps which I can call to do this... I
          // should also make the base grid derive from base-VectorOps, so that the 
          // functions are accessible!
          //
          if      (SimPM.coord_sys==COORD_CYL && axis[idim]==Rcyl) {
            double rp = CI.get_dpos(c,Rcyl)+0.5*dx; double rn=rp-dx;
            cpt->dU[ERG] += 2.0*(rn*q_neg-rp*q_pos)/(rp*rp-rn*rn);
          }
          else if (SimPM.coord_sys==COORD_SPH && axis[idim]==Rsph) {
            double rc = CI.get_dpos(c,Rsph);
            double rp = rc+0.5*dx; double rn=rp-dx;
            rc = (pow(rp,3.0) -pow(rn,3.0))/3.0;
            cpt->dU[ERG] += (rn*rn*q_neg-rp*rp*q_pos)/rc;
          }
          else {
            cpt->dU[ERG] += (q_neg-q_pos)/dx;
          }
          //cout <<"\tQc="<<Qclassical<<", Qs="<<Qsaturated<<", Q="<<q_pos<<", Edot="<<cpt->dU[ERG]<<"\n";

 	  //
	  // Set npt to cpt, set current q_pos to q_neg for next cell.
	  // Move to next cell.
	  //
          q_neg = q_pos;
	  cpt = npt;
	}  while ( (npt=grid->NextPt(npt,posdirs[idim])) !=0);

	// --------------------------------------------------------
	// Finished Heat conduction calculation for the column.
	// --------------------------------------------------------

	//
	// Get next x-column, or if it doesn't exist set a flag to
	// indicate that we are finished.
	//
	if (SimPM.ndim==1) xcolumns_finished = true;
	else {
	  start = grid->NextPt(start,posdirs[(idim+1)%SimPM.ndim]);
	  if (!start || !start->isgd) xcolumns_finished = true;
	}
      } while (!xcolumns_finished);
	
      //
      // Get next z-plane, or if it doesn't exist set a flag to
      // indicate that we are finished.  Reset marker to the first
      // cell in the new plane.
      //
      if (SimPM.ndim<=2) zplanes_finished = true;
      else {
	start = grid->NextPt(marker,posdirs[(idim+2)%SimPM.ndim]);
	marker = start;
	if (!start || !start->isgd) zplanes_finished = true;
      }
    } while (!zplanes_finished);

  } // Loop over Ndim directions.

  //
  // Now reset the temporary storage of Temperature in dU[RHO] to zero.
  //
  c = grid->FirstPt();
  do {
    c->dU[RHO]=0.0;
  } while ((c=grid->NextPt(c)) !=0);

  return 0;
}
#endif // THERMAL_CONDUCTION




// ##################################################################
// ##################################################################



