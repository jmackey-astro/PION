///
/// \file stellar_wind_BC.cc
/// \author Jonathan Mackey
/// \date 2010.10.05
///
/// Stellar wind boundary condition class.
/// Moved from global.cc b/c it got too big.
///
/// Modifications:
///
/// - 2010.10.05 JM: Added spherical coordinate possibility.  Also
///    changed pressure calculation to calculate it at the stellar
///    radius and scale outwards assuming adiabatic expansion.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  So VX,VY,VZ are now set by using calls to
///   distance() functions in the grid class.  It doesn't seem to help
///   much.
/// - 2011.01.07 JM: I debugged the geometric grid functions, and now
///   it works very well!  I have a nice spherical expansion.
/// - 2011.01.18 JM: Added #def which sets pressure so that Tmin=10K
/// - 2011.02.14 JM: Added stellar_wind_evolution class for winds with 
///    evolving properties, determined by a stellar evolution model, and
///    fitted with spline functions.
///    02.15 JM: Debugged. 02.16 JM: Debugged
/// - 2011.04.29 JM: Now in add_cell(), the c->isbd bool is set to true to
///    indicate that it has become boundary data (microphysics updates 
///    will skip it in this case and, more importantly, microphysics
///    timescales calculations.
/// - 2011.06.20 JM: Got rid of non-ANSI-C exp10 functions
/// - 2011.11.22 JM: Added t_scalefactor parameter for stellar winds.
/// - 2011.12.01 JM: Switched from spline to linear interpolation for
///    winds.
/// - 2012.12.07/10 JM: Changed min. ion frac. in wind from 0 to 1e-7.
/// - 2013.04.15 JM: removed lots of comments (or put in TESTING def)
/// - 2013.04.16 JM: Fixed bug where Set_Temp() was called when 
///    tracer variables were still (potentially) unset in wind cells.
/// - 2013.08.19 JM: got rid of cm_per_km() function.
/// - 2013.09.06 JM: removed the slowly-expanding switch-on wind
///    because it didn't help with anything, ever.
///    Removed the integer positions because they created potential
///    errors in the wind properties from rounding errors.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.10.19 JM: Fixed wind-tracer to always use pion_flt.
/// - 2017.07.26 JM: cleaned up code a bit.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/interpolate.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "grid/grid_base_class.h"
#include "grid/stellar_wind_BC.h"
#include "microphysics/microphysics_base.h"
#include <sstream>
using namespace std;

//#define TESTING


// ##################################################################
// ##################################################################



stellar_wind::stellar_wind(
      const int nd, ///< ndim
      const int nv, ///< nvar
      const int nt, ///< ntracer
      const int ft, ///< ftr
      const int cs, ///< coord_sys
      const int eq, ///< eqn_type
      const double mt ///< Minimum temperature allowed on grid
      )
: ndim(nd), nvar(nv), ntracer(nt), ftr(ft), coordsys(cs), eqntype(eq), Tmin(mt)
{
  nsrc=0;
}


// ##################################################################
// ##################################################################




stellar_wind::~stellar_wind()
{
  //
  // Need to delete the wind_cell structs in each wind_source struct,
  // and then delete the wind_source structs.
  //
  struct wind_source *ws;

  for (int n=0; n<nsrc; n++) {
    ws = wlist[n];
    for (int m=0; m<ws->ncell; m++) {
      ws->wcells[m]->p = mem.myfree(ws->wcells[m]->p);
      ws->wcells[m] = mem.myfree(ws->wcells[m]);
    }
    ws->wcells.clear();
    ws->tracers = mem.myfree(ws->tracers);
    ws = mem.myfree(ws);
  }
  return;      
}



// ##################################################################
// ##################################################################


// Function to replace pow(a, b) - exp(b*log(a)) is twice as fast
double stellar_wind::pow_fast(
	double a,
	double b
	)
{
	return exp(b*log(a));
}


// ##################################################################
// ##################################################################



int stellar_wind::add_source(
      const double *pos, ///< position (cm, w.r.t. grid origin)
      const double  rad, ///< radius (cm)
      const int    type, ///< type (0=constant, only option here)
      const double mdot, ///< Mdot (Msun/yr)
      const double vinf, ///< Vinf (km/s)
      const double temp, ///< Wind Temperature (K) (actually p_g.m_p/(rho.k_b))
      const double Rstar, ///< Radius at which T=T* (cm).
      const pion_flt *trv  ///< Tracer values of wind (if any)
      )
{
  struct wind_source *ws = 0;
  ws = mem.myalloc(ws,1);
  ws->id = wlist.size();
  ws->ncell = 0;
  ws->type = type;
  switch (type) {
  case WINDTYPE_CONSTANT: case WINDTYPE_EVOLVING:
    cout <<"\tAdding constant wind source as id="<<ws->id<<"\n";
    break;
  default:
    rep.error("What type of source is this?  add a new type?",type);
    break;
  }

  for (int v=0;v<ndim;v++)
    ws->dpos[v] = pos[v];
  rep.printVec("ws->dpos",ws->dpos,ndim);

  for (int v=ndim;v<MAX_DIM;v++)
    ws->dpos[v] = VERY_LARGE_VALUE;

  ws->radius = rad;

  //
  // Mdot and Vinf are passed to the function in Msun/yr and km/s,
  // but are stored internally in cgs units.
  //
  ws->Mdot  = mdot *pconst.Msun()/pconst.year();
  ws->Vinf  = vinf *1.0e5;

  ws->Tw    = temp;
  ws->Rstar = Rstar;

  ws->tracers=0;
  ws->tracers = mem.myalloc(ws->tracers,ntracer);
  for (int v=0;v<ntracer; v++) {
    ws->tracers[v] = trv[v];
    cout <<"ws->tracers[v] = "<<ws->tracers[v]<<"\n";
  }
  // if using microphysics, find H+ tracer variable, if it exists.
  int hplus=-1;
  if (MP) {
    hplus = MP->Tr("H1+");
  }
  ws->Hplus = hplus;

  ws->cells_added = false;
  if (!ws->wcells.empty())
    rep.error("wind_source: wcells not empty!",ws->wcells.size());

  //
  // Make sure the source position is compatible with the geometry:
  //
  if (coordsys==COORD_SPH) {
    if (!pconst.equalD(ws->dpos[Rsph],0.0))
      rep.error("Spherical symmetry but source not at origin!",
                ws->dpos[Rsph]);
  }
  if (coordsys==COORD_CYL && ndim==2) {
    if (!pconst.equalD(ws->dpos[Rcyl],0.0))
      rep.error("Axisymmetry but source not at R=0!",ws->dpos[Rcyl]);
  }

  wlist.push_back(ws);
  nsrc++;
#ifdef TESTING
  cout <<"\tAdded wind source id="<<nsrc-1<<" to list of ";
  cout <<nsrc<<" elements.\n";
#endif // TESTING
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind::Nsources()
{
  return nsrc;
}



// ##################################################################
// ##################################################################



int stellar_wind::add_cell(
      class GridBaseClass *grid,
      const int id, ///< src id
      cell *c       ///< cell to add to list.
      )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);
  struct wind_source *WS = wlist[id];

  //
  // Setup a wind_cell struct
  //
  struct wind_cell *wc = 0;
  wc = mem.myalloc(wc,1);
  //
  // GEOMETRY: this is a distance from a vertex to a cell--centre, so
  // we call the appropriate function here.
  //
  wc->dist = grid->distance_vertex2cell(WS->dpos, c);
  if (wc->dist > WS->radius) {
    rep.warning("stellar_wind::add_cell() cell is outside radius",
                WS->radius,wc->dist);
    CI.print_cell(c);
    return 1;
  }
  //
  // Now set wc cell pointer to this one.  Also set c->isbd to indicate
  // that it is now boundary data (while also grid data).
  //
  c->isbd=true;
  c->isdomain=false;
  wc->c = c;

  //
  // Calculate the polar angle theta
  //

  // Set theta to 0 if 1D - no angle dependent wind in this case (should add exit if angle + 1D)
  if (ndim == 1) {
	wc->theta = 0;
  }

  // Polar angle in 2D
  else if (ndim == 2) {
    // Opposite and adjacent of cell angle
    double opp = grid->difference_vertex2cell(WS->dpos, c, Rcyl);
    double adj = grid->difference_vertex2cell(WS->dpos, c, Zcyl);
    wc->theta = atan(fabs(opp/adj));
  }

  // Polar angle in 3D
  else if (ndim == 3) {
    // Opposite and adjacent in X-Y plane
    double opp1 = grid->difference_vertex2cell(WS->dpos, c, XX);
    double adj1 = grid->difference_vertex2cell(WS->dpos, c, YY);
    // Opposite and adjacent in Y-Z plane
    double opp2 = grid->difference_vertex2cell(WS->dpos, c, ZZ);
    double adj2 = sqrt(pow_fast(opp1, 2.0) + pow_fast(adj1, 2.0));
    wc->theta = atan(fabs(adj2/opp2));
  }

  //
  // Allocate memory for wind_cell reference state vector.
  //
  wc->p = 0;
  wc->p = mem.myalloc(wc->p, nvar);
  
  //
  // Now assign values to the state vector:
  //
  /// NB: This function has the EOS Gamma hardcoded to 5/3
  set_wind_cell_reference_state(grid, wc, WS, 5./3.);

  WS->wcells.push_back(wc);
  WS->ncell += 1;

#ifdef TESTING
  cout <<"*** dist="<<wc->dist<<". "; 
  rep.printVec("Wind BC cell pos",wc->c->pos,ndim);
  rep.printVec("Wind BC cell values", wc->p, nvar);
  CI.print_cell(c);
  cout <<" Added cell: array size="<<WS->wcells.size()<<"\tncell="<<WS->ncell<<"\n";
#endif
  return 0;
}


// ##################################################################
// ##################################################################




void stellar_wind::set_wind_cell_reference_state(
      class GridBaseClass *grid,
      struct wind_cell *wc,
      const struct wind_source *WS,
      const double gamma ///< EOS gamma
      )
{
  //
  // In this function we set the density, pressure, velocity, and tracer values
  // for the reference state of the cell.  Every timestep the cell-values will
  // be reset to this reference state.
  //
  double pp[ndim];
  CI.get_dpos(wc->c,pp);
  //
  // Density at cell position: rho = Mdot/(4.pi.R^2.v_inf) (for 3D)
  // or in 2D (slab-symmetry) rho = Mdot/(2.pi.R.v_inf)
  //
  if (ndim==2 && coordsys==COORD_CRT) {
    //
    // 2D slab symmetry --> 1/r force laws and density profile
    //
    wc->p[RO] = WS->Mdot/(WS->Vinf*2.0*M_PI*wc->dist);
    //
    // Set pressure based on wind density/temperature at the stellar radius:
    // rho_star = Mdot/(2.pi.R_star.v_inf), p_star = rho_star.k.T_star/(mu.m_p)
    // and then p(r) = p_star (rho(r)/rho_star)^gamma
    // ******************************************************************************
    // *********** WARNING MU=1 HERE, PROBABLY SHOULD BE O.6 (IONISED) 1.3 (NEUTRAL).
    // ******************************************************************************
    //
    wc->p[PG] = pconst.kB()*WS->Tw/pconst.m_p();
    wc->p[PG]*= exp((gamma-1.0)*log(2.0*M_PI*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((gamma)*log(wc->p[RO]));
  }

  else {
    //
    // 3D geometry, so either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
    // rho = Mdot/(4.pi.R^2.v_inf) 
    //
    wc->p[RO] = 1.0/(wc->dist);
    wc->p[RO] *= wc->p[RO];
    wc->p[RO] *= WS->Mdot/(WS->Vinf*4.0*M_PI);
    //
    // Set pressure based on wind density/temperature at the stellar radius,
    // assuming adiabatic expansion outside Rstar, and that we don't care what
    // the temperature is inside Rstar (because this function will make it 
    // hotter than Teff, which is not realistic):
    //
    // rho_star = Mdot/(4.pi.R_star^2.v_inf),
    //   p_star = rho_star.k.T_star/(mu.m_p)
    // So then p(r) = p_star (rho(r)/rho_star)^gamma
    // ******************************************************************************
    // *********** WARNING MU=1 HERE, PROBABLY SHOULD BE O.6 (IONISED) 1.3 (NEUTRAL).
    // ******************************************************************************
    //
    wc->p[PG] = pconst.kB()*WS->Tw/pconst.m_p();
    wc->p[PG]*= exp((gamma-1.0)*log(4.0*M_PI*WS->Rstar*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((gamma)*log(wc->p[RO]));
  }


  //
  // VELOCITIES: These should be cell-average values, which are
  // the values at the centre-of-volume, so we call the geometry-aware
  // grid functions.
  //
  cell *c = wc->c;
  wc->p[VX] = WS->Vinf*grid->difference_vertex2cell(WS->dpos,c,XX)/wc->dist;
  if (ndim>1)
    wc->p[VY] = WS->Vinf*grid->difference_vertex2cell(WS->dpos,c,YY)/wc->dist;
  else
    wc->p[VY] = 0.0;
  if (ndim>2)
    wc->p[VZ] = WS->Vinf*grid->difference_vertex2cell(WS->dpos,c,ZZ)/wc->dist;
  else
    wc->p[VZ] = 0.0;

  // Add in statement for magnetic field of the stellar wind (B=100G, R=10Ro)....
  if (eqntype==EQMHD || eqntype==EQGLM) {
    double R=695508e5;  // R_sun in cm
    double x = grid->difference_vertex2cell(WS->dpos,c,XX);
    wc->p[BX] = (100.0/sqrt(4.0*M_PI))*pow(10.0*R/wc->dist,2)*
                fabs(x)/wc->dist;
    if (ndim>1) {
      wc->p[BY] = (100.0/sqrt(4.0*M_PI))*pow(10.0*R/wc->dist,2)*
                grid->difference_vertex2cell(WS->dpos,c,YY)/wc->dist;
      wc->p[BY] = (x>0.0) ? wc->p[BY] : -1.0*wc->p[BY];
    }
    else
      wc->p[BY] = 0.0;
    if (ndim>2) {
      wc->p[BZ] = (100.0/sqrt(4.0*M_PI))*pow(10.0*R/wc->dist,2)*
                grid->difference_vertex2cell(WS->dpos,c,ZZ)/wc->dist;
      wc->p[BZ] = (x>0.0) ? wc->p[BZ] : -1.0*wc->p[BZ];
    }
    else
      wc->p[BZ] = 0.0;
  }
  if (eqntype==EQGLM) {
    wc->p[SI] = 0.0;
  }
    
  //if (eqntype!=EQEUL && eqntype!=EQEUL_EINT)
    //rep.error("Need to code B into winds model!",eqntype);

  // update tracers
  for (int v=0;v<ntracer;v++)
    wc->p[ftr+v] = WS->tracers[v];

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // Set a minimum temperature in the wind
  //
  if (MP) {
    if (MP->Temperature(wc->p,gamma) <WS->Tw) {
      MP->Set_Temp(wc->p,WS->Tw,gamma);
    }
  }
  else {
    // appropriate for a neutral medium, He+M mass fraction 0.285.
    wc->p[PG] = max(static_cast<double>(wc->p[PG]), 
      Tmin*wc->p[RO]*pconst.kB()*0.78625/pconst.m_p());
  }
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return;
}



// ##################################################################
// ##################################################################



int stellar_wind::set_num_cells(
      const int id, ///< src id
      const int nc  ///< number of cells.
      )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);

  if (nc != wlist[id]->ncell) {
    rep.warning("stellar_wind::set_num_cells() COUNTING PROBLEM!!",
                nc,wlist[id]->ncell);
    return 1;
  }
  return 0;
}



// ##################################################################
// ##################################################################




int stellar_wind::get_num_cells(
      const int id ///< src id
      )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);
  return wlist[id]->ncell;
}



// ##################################################################
// ##################################################################



int stellar_wind::set_cell_values(
      class GridBaseClass *,
      const int id,  ///< src id
      const double t ///< simulation time
      )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);

  //
  // go through every cell in one go and update them all
  //
  struct wind_source *WS = wlist[id];
#ifdef TESTING
  cout <<"updating source "<<id<<" which has "<<WS->ncell<<" cells.\n";
#endif
  for (int i=0; i<WS->ncell; i++) {
    //cout <<"i="<<i<<", density = "<<WS->wcells[i]->p[RO]<<"\n";
    for (int v=0;v<nvar;v++)
      WS->wcells[i]->c->P[v]  = WS->wcells[i]->p[v];
    for (int v=0;v<nvar;v++)
      WS->wcells[i]->c->Ph[v] = WS->wcells[i]->p[v];
  }
    
  return 0;
}



// ##################################################################
// ##################################################################



void stellar_wind::get_src_posn(
      const int id, ///< src id
      double *x     ///< position vector (output)
      )
{
  for (int v=0;v<ndim;v++) x[v] = wlist[id]->dpos[v];
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_Mdot(
      const int id, ///< src id
      double *x     ///< mdot (output)
      )
{
  *x = wlist[id]->Mdot *pconst.year()/pconst.Msun();
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_drad(
      const int id, ///< src id
      double *x     ///< radius (output)
      )
{
  *x = wlist[id]->radius;
  return;
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_Vinf(
      const int id, ///< src id
      double *x   ///< Vinf (output)
      )
{
  *x = wlist[id]->Vinf/1.0e5;
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_Tw(
      const int id, ///< src id
      double *x   ///< Stellar Radius (output)
      )
{
  *x = wlist[id]->Tw;
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_Rstar(
      const int id, ///< src id
      double *x   ///< Stellar radius (output)
      )
{
  *x = wlist[id]->Rstar;
}


// ##################################################################
// ##################################################################



void stellar_wind::get_src_trcr(
      const int id, ///< src id
      pion_flt *x   ///< tracers (output)
      )
{
  for (int v=0;v<ntracer; v++)
    x[v] = wlist[id]->tracers[v];
}



// ##################################################################
// ##################################################################



void stellar_wind::get_src_type(
      const int id, ///< src id
      int *x   ///< type of wind (=0 for now) (output)
      )
{
  *x=wlist[id]->type;
}
//------------------------------------------------






// ##################################################################
// ##################################################################



double stellar_wind::beta(const double Teff)
{
  //
  // Eldridge et al. (2006, MN, 367, 186).
  // Beta = Zeta^2
  //
  double beta;
  double rsg=0.125; // Eldridge value
  //double rsg=0.04;  // Mackey+2012 Betelgeuse value

  if (Teff <= 3600.0)
    beta = rsg;
  else if (Teff >= 22000.0)
    beta = 2.6;
  else {
    //
    // Linear interpolation for beta from Eldridge et al. Table 1.
    //
    double b0, b1, T0, T1;
    if      (Teff<6000.0) {
      T0 = 3600.0; b0 = rsg;
      T1 = 6000.0; b1 = 0.5;
    }
    else if (Teff <8000.0) {
      T0 = 6000.0; b0 = 0.5;
      T1 = 8000.0; b1 = 0.7;
    }
    else if (Teff <10000.0) {
      T0 = 8000.0; b0 = 0.7;
      T1 = 10000.0; b1 = 1.3;
    }
    else if (Teff <20000.0) {
      T0 = 10000.0; b0=1.3;
      T1 = 20000.0; b1=1.3;
    }
    else {
      T0 = 20000.0; b0 = 1.3;
      T1 = 22000.0; b1 = 2.6;
    }
    beta = b0 + (Teff-T0)*(b1-b0)/(T1-T0);
  }

  return beta;
};




// ##################################################################
// ##################################################################



// ------------------------------------------------------------------
// ----------  STELLAR WIND WITH STELLAR EVOLUTION ------------------
// ------------------------------------------------------------------

stellar_wind_evolution::stellar_wind_evolution(
      const int nd, ///< ndim
      const int nv, ///< nvar
      const int nt, ///< ntracer
      const int ft, ///< ftr
      const int cs, ///< coord_sys
      const int eq, ///< eqn_type
      const double mt, ///< Minimum temperature allowed on grid
      const double ss, ///< Simulation start time.
      const double sf  ///< Simulation finish time.
      )
: stellar_wind(nd,nv,nt,ft,cs,eq,mt), sim_start(ss), sim_finish(sf)
{
#ifdef TESTING
  cout <<"Stellar wind with time evolution, constructor.\n";
#endif
}



// ##################################################################
// ##################################################################



stellar_wind_evolution::~stellar_wind_evolution()
{
#ifdef TESTING
  cout <<"Stellar wind with time evolution, destructor.\n";
#endif

  //
  // Delete arrays
  //
  // delete each wdata_evol[] element, making sure to delete the
  // arrays in each element first (arrays, pos, tracers).
  //
  struct evolving_wind_data *wd=0;
  while (wdata_evol.size() >0) {
    wd = wdata_evol.back();
    wdata_evol.pop_back();
    wd = mem.myfree(wd);
  }
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::add_source(
        const double *pos, ///< position (physical units)
        const double  rad, ///< radius (physical units)
        const int    type, ///< type (0=fixed in time,1=slow switch on)
        const double mdot, ///< Mdot (Msun/yr)
        const double vinf, ///< Vinf (km/s)
        const double Twnd, ///< Wind Temperature (actually p_g.m_p/(rho.k_b))
        const double Rstar, ///< Stellar radius (for T*-->gas pres.).
        const pion_flt *trv  ///< Tracer values of wind (if any)
        )
{
  //
  // This function sets up an evolving_wind_data struct to hold the 
  // constant wind.  It sets all the evolving wind stuff to zero, and
  // then calls the stellar_wind:: version to set up the wind source.
  //
  struct evolving_wind_data *temp=0;
  temp = mem.myalloc(temp,1);
  //
  // First set all the stuff we don't need for a constant wind.
  //
  temp->Npt = 0;
  temp->offset = 0.0;
  temp->tstart = sim_start -1.0e10; // so it has already started.
  temp->tfinish= 2.0*sim_finish;    // so it keeps going to end of sim.
  temp->update_freq = 1.0e99;             // so it never updates.
  temp->t_next_update = 1.0e99;           // so it never updates.
  
  //
  // Set source to be active
  //
  temp->is_active = true;
  //
  // Now add source using constant wind version.
  //
  stellar_wind::add_source(pos,rad,type,mdot,vinf,Twnd,Rstar,trv);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);

  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::add_evolving_source(
  const double *pos,        ///< position (physical units).
  const double  rad,        ///< radius (physical units).
  const int    type,        ///< type (must be 3, for variable wind).
  const double Rstar,       ///< Radius at which to get gas pressure from Teff
  const pion_flt *trv,        ///< Any (constant) wind tracer values.
  const string infile,      ///< file name to read data from.
  const int ,   ///< enhance mdot based on rotation (0=no,1=yes).
  const double time_offset, ///< time offset = [t(sim)-t(wind_file)] in years
  const double t_now,       ///< current simulation time, to see if src is active.
  const double update_freq, ///< frequency with which to update wind properties.
  const double t_scalefactor ///< wind evolves this factor times faster than normal
  )
{
  if (type != WINDTYPE_EVOLVING) {
    rep.error("Bad wind type for evolving stellar wind!",type);
  }
  //
  // First we will read the file, and see when the source should
  // switch on in the simulation (it may not be needed for a while).
  //
#ifdef TESTING
  cout <<"\t\tsw-evo: adding source from file "<<infile<<"\n";
#endif

  //
  // Wind source struct, to be added to class vector wdata_evol
  //
  struct evolving_wind_data *temp=0;
  temp = mem.myalloc(temp,1);
  //temp->Npt = Npt;
  //temp->t = 0;
  //temp->mdot = 0;
  //temp->vinf = 0;
  //temp->Teff = 0;

  //
  // Read in stellar evolution data
  // Format: time	M	L	Teff	Mdot	vrot   vcrit
  //
  FILE *wf = 0;
  wf = fopen(infile.c_str(), "r");
  if (!wf) rep.error("can't open wind file, stellar_wind_angle",wf);
  // Skip first two lines
  char line[512];
  char *rval=0;
  rval = fgets(line,512,wf);
  if (!rval) rep.error("stwind_angle: failed to get line 1",line);
  //printf("%s",line);
  rval = fgets(line,512,wf);
  if (!rval) rep.error("stwind_angle: failed to get line 2",line);
  //printf("%s",line);

  // Temp. variables for column values
  double t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0, t6=0.0, t7=0.0;
  while ((rval = fgets(line,512,wf))  != 0) {
    sscanf(line, "   %lE   %lE %lE %lE %lE %lE %lE",
                      &t1, &t2, &t3, &t4, &t5, &t6, &t7);
    //cout.precision(16);
#ifdef TESTING
    cout <<t1 <<"  "<<t2  <<"  "<< t3  <<"  "<< t4 <<"  ";
    cout << t5 <<"  "<< t6 <<"\n";
#endif
    // Set vector value
    temp->time_evo.push_back(t1);
    temp->M_evo.push_back(t2);
    temp->L_evo.push_back(t3);
    temp->Teff_evo.push_back(t4);

    // Stellar radius
    t6 = sqrt( t3/ (4.0*pconst.pi()*pconst.StefanBoltzmannConst()*
                                                pow_fast(t4, 4.0)));
    temp->R_evo.push_back(t6);
    
    // Hydrogen mass fraction (should make this a sim parameter?) 
    double H_X = 0.7;

    // Eddington luminosity (taking the opacity as the electron
    // scattering cross section)
    double L_edd = (4.0*pconst.pi()*pconst.c()*pconst.G()*t2)/
                                                    (0.2*(1 + H_X));

    // Escape velocity
    temp->vinf_evo.push_back(sqrt(beta(t4)) 
                        *sqrt(2.0*pconst.G()*t2*(1 - t3/L_edd)/t6));

    
    // Mdot: 
    temp->Mdot_evo.push_back(t5);
  }
  fclose(wf);

  // Column length
  size_t Npt = temp->time_evo.size();
  temp->Npt = Npt;

  //
  // Next we set up the interpolation, first modifying the
  // time array using the time-offset so that it has the same zero
  // offset as the simulation time.
  //
  for (size_t i=0; i<Npt; i++) {
    //
    // times in the file are measured in seconds, so offset should be
    // in years.
    //
    temp->time_evo[i] += time_offset;
    // scale times by scale factor.
    temp->time_evo[i] /= t_scalefactor;
    //cout <<"t="<<temp->time_evo[i]<<"\n";
  }


  //
  // Offset is not used in the code past here.  It's just here for I/O b/c a
  // restart will need to read the data-file again.
  // We reset all variables to be in seconds here.  So the interpolation uses
  // years, but everything else is in seconds.  Note that the global SWP struct
  // still has times in years though.
  //
  temp->offset = time_offset/t_scalefactor;    // in seconds
  temp->tstart = temp->time_evo[0];            // in seconds (already scaled)
  temp->tfinish= temp->time_evo[Npt-1];        // in seconds (already scaled)
  temp->update_freq = update_freq/t_scalefactor; // in seconds
  temp->t_next_update = max(temp->tstart,t_now);
#ifdef TESTING
  cout <<"\t\t tstart="<<temp->tstart;
  cout <<", next update="<<temp->t_next_update;
  cout <<", and tfinish="<<temp->tfinish<<"\n";
#endif

  //
  // We need to decide if the wind src is active yet.  If it is, then
  // we also set up a constant wind source for updating its
  // properties.  We set it to be active if the current time is
  // within update_freq of tstart.
  //
  double mdot=0.0, vinf=0.0, Twind=0.0, rstar=0.0;
  if ( ((t_now+temp->update_freq)>temp->tstart ||
        pconst.equalD(temp->tstart, t_now))
       && t_now<temp->tfinish) {
    temp->is_active = true;
    //
    // Get the current values for mdot, vinf, Teff, and setup a wind
    // source using the constant-wind function.
    //
    interpolate.root_find_linear_vec(temp->time_evo, temp->Teff_evo, t_now, Twind);
    interpolate.root_find_linear_vec(temp->time_evo, temp->Mdot_evo, t_now, mdot);
    interpolate.root_find_linear_vec(temp->time_evo, temp->vinf_evo, t_now, vinf);
    interpolate.root_find_linear_vec(temp->time_evo, temp->R_evo, t_now, rstar);
#ifdef TESTING
    cout <<"Source is Active\n";
    cout <<"T = "<<Twind<<",  mdot="<<mdot<<",  vinf="<<vinf;
    cout <<",  rstar="<<rstar<<"\n";
    cout.flush();
#endif
  }
  else {
    cout <<"WARNING: Source is not yet active: tnow="<<t_now;
    cout <<", tstart=";
    cout <<temp->tstart<<". Setting wind source to INACTIVE.\n";
    temp->is_active = false;
    mdot=-100.0; vinf=-100.0; Twind=-100.0;
  }


  //
  // Now add source using constant wind version.
  //
  stellar_wind::add_source(pos,rad,type,mdot,vinf,Twind,Rstar,trv);
  temp->ws = wlist.back();

  //
  // So now we have all of the properties of the wind source in the
  // struct evolving_wind_data 'temp', and we have added it to the
  // list of constant wind sources, wlist[].  We can now add temp to
  // the list of evolving wind sources and return (so that wlist[i]
  // and wdata_evol[i] point to the same source).
  //
  wdata_evol.push_back(temp);
  //NSRC_TOTAL++;

  return 0;
}



// ##################################################################
// ##################################################################



void stellar_wind_evolution::update_source(
        class GridBaseClass *grid,
        struct evolving_wind_data *wd,
        const double t_now,
        const double gamma
        )
{
  //
  // We have a source that needs updating.  If it is not active, and
  // needs activating then we set that.
  //
  if (!wd->is_active) {
    cout <<"stellar_wind_evo::update_source activating source id=";
    cout << wd->ws->id <<" at Simulation time t="<<t_now<<"\n";
    rep.printVec("Source position",wd->ws->dpos,ndim);
    wd->is_active=true;
  }

  if (t_now < wd->tstart) {
    rep.error("Updating source, not yet active!",wd->tstart-t_now);
  }

  wd->t_next_update = t_now; // (We update every timestep now)
  wd->t_next_update = min(wd->t_next_update, wd->tfinish);
  
  //
  // Now we update Mdot, Vinf, Teff by linear interpolation.
  //
  double mdot=0.0, vinf=0.0, Twind=0.0, rstar=0.0;
  interpolate.root_find_linear_vec(wd->time_evo, wd->Teff_evo, t_now, Twind);
  interpolate.root_find_linear_vec(wd->time_evo, wd->Mdot_evo, t_now, mdot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vinf_evo, t_now, vinf);
  interpolate.root_find_linear_vec(wd->time_evo, wd->R_evo, t_now, rstar);
  //
  // Assign new values to wd->ws (the wind source struct), all in CGS
  // because the arrays are in CGS.
  //
#ifdef TESTING
  cout <<"updating wind: old [Mdot,Vinf,Tstar,Rstar] = [";
  cout <<wd->ws->Mdot<<", "<<wd->ws->Vinf<<", "<<wd->ws->Tw<<", "<<wd->ws->Rstar<<"]\n";
  cout <<"               new [Mdot,Vinf,Tstar,Rstar] = [";
  cout <<mdot<<", "<<vinf<<", "<<Twind<<", "<<rstar<<"]\n";
#endif
  wd->ws->Mdot = mdot;  // already cgs.
  wd->ws->Vinf = vinf;  // this is in cm/s already.
  wd->ws->Tw   = Twind; // This is in K.
  wd->ws->Rstar = rstar;

  //
  // Now re-assign state vector of each wind-boundary-cell with
  // updated values.
  //
  for (int i=0; i<wd->ws->ncell; i++) {
    set_wind_cell_reference_state(grid,wd->ws->wcells[i],wd->ws,gamma);
  }

  //
  // Now the source is updated, and the reference states are all set
  // for the new values, and the next update time has been set.  So
  // we can return.
  //
  return;
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::set_cell_values(
        class GridBaseClass *grid,
        const int id,  ///< src id
        const double t_now ///< simulation time
        )
{
  int err=0;
  if (id<0 || id>=nsrc) {
    rep.error("bad src id",id);
  }
  
  struct evolving_wind_data *wd = wdata_evol[id];

  //
  // Check if the source needs to be updated, and if so, update it
  //
  if (t_now >= wd->t_next_update) {
    /// NB: This function has hardcoded EOS Gamma to 5/3
    update_source(grid, wd, t_now, 5./3.);
  }

  //
  // Now we have an up-to-date source.  We check if it is active,
  // and if so we update the wind cell values using the constant wind
  // functions.  If not, then we ignore it and return.
  //
  if (wd->is_active) {
    err += stellar_wind::set_cell_values(grid, id,t_now);
  }

  return err;
}



// ##################################################################
// ##################################################################





