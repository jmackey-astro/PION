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
      const std::string *tr,  ///< List of tracer variable names.
      const int cs, ///< coord_sys
      const int eq, ///< eqn_type
      const double mt ///< Minimum temperature allowed on grid
      )
: ndim(nd), nvar(nv), ntracer(nt), ftr(ft), coordsys(cs), eqntype(eq), Tmin(mt)
{
  nsrc=0;

  for (int v=0;v<nt;v++) tracers.push_back(tr[v]);
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
      const double vrot, ///< Vrot (km/s)
      const double temp, ///< Wind Temperature (K)
      const double Rstar, ///< Radius of star (cm).
      const double Bstar, ///< Surface Magnetic field of star (Gauss).
      pion_flt *trv  ///< Tracer values of wind (if any)
      )
{
  struct wind_source *ws = 0;
  ws = mem.myalloc(ws,1);
  ws->id = wlist.size();
  ws->ncell = 0;
  ws->type = type;
  switch (type) {
  case WINDTYPE_CONSTANT: case WINDTYPE_EVOLVING:
    //cout <<"\tAdding wind source as id="<<ws->id<<"\n";
    break;
  default:
    rep.error("What type of source is this?  add a new type?",type);
    break;
  }

  for (int v=0;v<ndim;v++)
    ws->dpos[v] = pos[v];
  //rep.printVec("ws->dpos",ws->dpos,ndim);

  for (int v=ndim;v<MAX_DIM;v++)
    ws->dpos[v] = VERY_LARGE_VALUE;

  ws->radius = rad;

  //
  // Mdot and Vinf are passed to the function in Msun/yr and km/s,
  // but are stored internally in cgs units.
  //
  ws->Mdot  = mdot *pconst.Msun()/pconst.year();
  ws->Vinf  = vinf *1.0e5;
  ws->v_rot  = vrot *1.0e5;

  // temperature, radius and surface B-field already in CGS
  ws->Tw    = temp;
  ws->Rstar = Rstar;
  ws->Bstar = Bstar;

  ws->tracers=0;
  ws->tracers = mem.myalloc(ws->tracers,ntracer);
  for (int v=0;v<ntracer; v++) {
    ws->tracers[v] = trv[v];
    //cout <<"ws->tracers[v] = "<<ws->tracers[v]<<"\n";
  }

  // if using microphysics, find H+ tracer variable, if it exists.
  ws->Hplus = -1;
  ws->iHplus = -1;
  int hplus=-1;
  if (MP) {
    hplus = MP->Tr("H1+");
  }
  ws->Hplus = hplus;
  if (hplus >=0) ws->iHplus = hplus - nvar + ntracer; 

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
  if (wc-dist < 0.75*WC->radius) {
    wc->p[RO] = 1.0e-31;
    wc->p[PG] = 1.0e-31;
    wc->p[VX] = 0.0;
    wc->p[VY] = 0.0;
    wc->p[VZ] = 0.0;
    if (eqntype==EQMHD || eqntype==EQGLM) {
      wc->p[BX] = 0.0;
      wc->p[BY] = 0.0;
      wc->p[BZ] = 0.0;
    }
    if (eqntype==EQGLM) {
      wc->p[SI] = 0.0;
    }
  }


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
    //
    wc->p[PG] = pconst.kB()*WS->Tw/pconst.m_p();
    wc->p[PG]*= exp((gamma-1.0)*log(4.0*M_PI*WS->Rstar*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((gamma)*log(wc->p[RO]));
  }

  // Velocities and magnetic fields: get coordinates relative to star
  // for calculating sin/cos angles in theta and phi.
  cell *c = wc->c;
  double x,y,z;
  switch (ndim) {
  case 1:
    x = grid->difference_vertex2cell(WS->dpos,c,XX);
    y = 0.0;
    z = 0.0;
    break;
  case 2:
    x = grid->difference_vertex2cell(WS->dpos,c,XX);
    y = grid->difference_vertex2cell(WS->dpos,c,YY);
    z = 0.0;
    break;
  case 3:
    x = grid->difference_vertex2cell(WS->dpos,c,XX);
    y = grid->difference_vertex2cell(WS->dpos,c,YY);
    z = grid->difference_vertex2cell(WS->dpos,c,ZZ);
    break;
  default:
    rep.error("bad ndim in set_wind_cell_reference_state",ndim);
    break;
  }

  // Velocities: cell-average values, i.e. values at the
  // centre-of-volume.
  // TODO: for general J vector in 3D, what is rotational component.
  switch (ndim) {
  case 1:
    wc->p[VX] = WS->Vinf * x / wc->dist;
    wc->p[VY] = 0.0;
    wc->p[VZ] = 0.0;
    break;

  case 2:
    wc->p[VX] = WS->Vinf * x / wc->dist;
    wc->p[VY] = WS->Vinf * y / wc->dist;
    // J is hardcoded to be parallel to positive z-axis
    wc->p[VZ] = WS->v_rot *  WS->Rstar * y / pow_fast(wc->dist,2);
    break;

  case 3:
    wc->p[VX] = WS->Vinf * x / wc->dist;
    wc->p[VY] = WS->Vinf * y / wc->dist;
    wc->p[VZ] = WS->Vinf * z / wc->dist;

    // add non-radial component to x/y-dir from rotation.
    // J is hardcoded to be parallel to positive z-axis
    wc->p[VX] += -WS->v_rot * WS->Rstar * y / pow_fast(wc->dist,2);
    wc->p[VY] +=  WS->v_rot * WS->Rstar * x / pow_fast(wc->dist,2);
    break;

  default:
    rep.error("bad ndim in set_wind_cell_reference_state",ndim);
    break;
  }


  // Magnetic field: cell-average values, i.e. values at the
  // centre-of-volume.
  // Use a split monopole plus a rotational term adding toroidal
  // component.
  // TODO: for general J vector, what is rotational component.
  if (eqntype==EQMHD || eqntype==EQGLM) {
    //double t=0.0;
    double B_s = WS->Bstar/sqrt(4.0*M_PI); // code units for B_surf
    double D_s = WS->Rstar/wc->dist;     // 1/d in stellar radii
    double D_2 = D_s*D_s;                // 1/d^2 in stellar radii
    // this multiplies the toroidal component:
    double beta_B_sint = (WS->v_rot /  WS->Vinf) * B_s * D_s;

    switch (ndim) {
    case 1:
      rep.error("1D spherical but MHD?",ndim);
      break;
    case 2:
      // split monopole
      wc->p[BX] = B_s * D_2 * fabs(x)/wc->dist;
      wc->p[BY] = B_s * D_2 / wc->dist;
      wc->p[BY] = (x>0.0) ? y * wc->p[BY] : -y * wc->p[BY];
      // toroidal component
      beta_B_sint = beta_B_sint * y / wc->dist;
      wc->p[BZ] = (x>0.0) ? -beta_B_sint : beta_B_sint;
      break;

    case 3:
      // split monopole along z-axis, parallel to J
      wc->p[BX] = B_s * D_2 / wc->dist;
      wc->p[BX] = (z>0.0) ? x * wc->p[BX] : -x * wc->p[BX];

      wc->p[BY] = B_s * D_2 / wc->dist;
      wc->p[BY] = (z>0.0) ? y * wc->p[BY] : -y * wc->p[BY];

      wc->p[BZ] = B_s * D_2 * fabs(z) / wc->dist;

      // toroidal component in x-y plane from rotation, such that
      // we have a Parker spiral, inward winding for z<0 and 
      // outwards for z>0.
      beta_B_sint *= sqrt(x*x+y*y)/wc->dist;
      beta_B_sint = (z>0.0) ? -beta_B_sint : beta_B_sint;

      // modulate strength near the equator by linearly reducing 
      // torodial component for |theta|<1 degree from equator
      // See Pogorelov et al (2006,ApJ,644,1299).  This is for
      // testing the code.
      //t = fabs(z)/wc->dist * 180.0 / M_PI; // angle in degrees.
      //if (t < 2.0) beta_B_sint *= 0.5*t;

      wc->p[BX] += - beta_B_sint * y / wc->dist;
      wc->p[BY] +=   beta_B_sint * x / wc->dist;
      break;

    default:
      rep.error("bad ndim in set_wind_cell_reference_state",ndim);
      break;
    }
  }
  if (eqntype==EQGLM) {
    wc->p[SI] = 0.0;
  }
    
  // update tracers: should be set already.
  for (int v=0;v<ntracer;v++)
    wc->p[ftr+v] = WS->tracers[v];

  if (WS->Hplus>=0) {
    //cout <<"WS->Hplus = "<<WS->Hplus;
    //cout <<" itr = "<<WS->iHplus<<endl;
    if      (WS->Tw < 1.0e4) WS->tracers[WS->iHplus] = 1.0e-10;
    else if (WS->Tw > 1.5e4) WS->tracers[WS->iHplus] = 1.0;
    else    WS->tracers[WS->iHplus] = 1.0e-10 + (WS->Tw-1.0e4)*(1.0-1.0e-10)/0.5e4;
  }


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
  // Their eq. 2 has a typo, because the hot star values from Vink+
  // (2001) are for v_inf/v_esc, not for (v_inf/v_esc)^2.
  // So the beta value scales the escape velocity itself, not the
  // square.
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
      const std::string *tr,  ///< List of tracer variable names.
      const int cs, ///< coord_sys
      const int eq, ///< eqn_type
      const double mt, ///< Minimum temperature allowed on grid
      const double ss, ///< Simulation start time.
      const double sf  ///< Simulation finish time.
      )
: stellar_wind(nd,nv,nt,ft,tr,cs,eq,mt), sim_start(ss), sim_finish(sf)
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
      const double vrot, ///< Vrot (km/s)
      const double Twnd, ///< Wind Temperature (K)
      const double Rstar, ///< Stellar radius.
      const double Bstar, ///< Surface Magnetic field of star (Gauss).
      pion_flt *trv  ///< Tracer values of wind (if any)
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

  // find indices of elements and dust in tracers list.
  set_element_indices(temp);

  // Now add source using constant wind version.
  stellar_wind::add_source(pos,rad,type,mdot,vinf,vrot,Twnd,Rstar,Bstar,trv);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);

  return 0;
}



// ##################################################################
// ##################################################################



void stellar_wind_evolution::set_element_indices(
      struct evolving_wind_data *ewd
      )
{
  ewd->i_XH  = -1;
  ewd->i_XHe = -1;
  ewd->i_XC  = -1;
  ewd->i_XN  = -1;
  ewd->i_XO  = -1;
  ewd->i_XZ  = -1;
  ewd->i_XD  = -1;

  // Check for elements if using microphysics
  if (MP) {
    for (int v=0; v<ntracer;v++) {
      if      (tracers[v] == "X_H")  ewd->i_XH  = v;
      else if (tracers[v] == "X_He") ewd->i_XHe = v;
      else if (tracers[v] == "X_C")  ewd->i_XC = v;
      else if (tracers[v] == "X_N")  ewd->i_XN = v;
      else if (tracers[v] == "X_O")  ewd->i_XO = v;
      else if (tracers[v] == "X_Z")  ewd->i_XZ = v;
      else if (tracers[v] == "X_D")  ewd->i_XD = v;
    }
  }
  return;
}



// ##################################################################
// ##################################################################

int stellar_wind_evolution::read_evolution_file(
      const string infile,      ///< file name to read data from.
      struct evolving_wind_data *data   ///< where to put the data
      )
{

  //
  // Read in stellar evolution data
  // Format: time M L Teff Mdot vrot vcrit vinf X_H X_He X_C X_N X_O X_Z X_D
  //
  FILE *wf = 0;
  wf = fopen(infile.c_str(), "r");
  if (!wf) rep.error("can't open wind file, stellar_wind_angle",wf);

  // Skip first two lines
  char line[512];
  char *rval=0;
  rval = fgets(line,512,wf);
  if (!rval) rep.error("stwind_angle: failed to get line 1",line);
  //printf("Star Calculation Source: %s",line);
  cout <<"Star Calculation Source: "<<string(line)<<"\n";
  rval = fgets(line,512,wf);
  if (!rval) rep.error("stwind_angle: failed to get line 2",line);
  //printf("%s",line);

  // read file line by line and add to struct vectors.
  // Everthing must be in CGS units already.
  double time=0.0, mass=0.0, lumi=0.0, teff=0.0, radi=0.0, mdot=0.0,
    vrot=0.0, vcrt=0.0, vinf=0.0;
  double xh=0.0, xhe=0.0, xc=0.0, xn=0.0, xo=0.0, xz=0.0, xd=0.0;
  while ((rval = fgets(line,512,wf))  != 0) {
    sscanf(line, "   %lE   %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
            &time, &mass, &lumi, &teff, &mdot, &vrot, &vcrt, &vinf,
            &xh, &xhe, &xc, &xn, &xo, &xz, &xd);
    //cout.precision(16);
#ifdef TESTING
    cout <<time <<"  "<<mass  <<"  "<< lumi<<"  "<< teff <<"  ";
    cout << mdot <<"  "<< vrot <<"  "<<vcrt<<"  "<<vinf<<"\n";
#endif
    // Set vector value
    data->time_evo.push_back(time);
    data->M_evo.push_back(mass);
    data->L_evo.push_back(lumi);
    data->Teff_evo.push_back(teff);

    // Stellar radius
    radi = sqrt( lumi/ (4.0*pconst.pi()*pconst.StefanBoltzmannConst()*
                                                pow_fast(teff, 4.0)));
    data->R_evo.push_back(radi);
    
    data->Mdot_evo.push_back(mdot);
    data->vrot_evo.push_back(vrot);
    data->vcrt_evo.push_back(vcrt);
    data->vinf_evo.push_back(vinf);
    data->X_H_evo.push_back(xh);
    data->X_He_evo.push_back(xhe);
    data->X_C_evo.push_back(xc);
    data->X_N_evo.push_back(xn);
    data->X_O_evo.push_back(xo);
    data->X_Z_evo.push_back(xz);
    data->X_D_evo.push_back(xd);
  }
  fclose(wf);

  // Column length
  size_t Npt = data->time_evo.size();
  data->Npt = Npt;

  return 0;
}



// ##################################################################
// ##################################################################




int stellar_wind_evolution::add_evolving_source(
  const double *pos,        ///< position (physical units).
  const double  rad,        ///< radius (physical units).
  const int    type,        ///< type (must be 3, for variable wind).
  pion_flt *trv,        ///< Any (constant) wind tracer values.
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
  int err = read_evolution_file(infile,temp);
  if (err) rep.error("couldn't read wind evolution file",infile);

  //
  // Optional time offset between simulation time and evolutionary
  // time.  Also optional scaling.
  //
  for (int i=0; i<temp->Npt; i++) {
    temp->time_evo[i] += time_offset;
    temp->time_evo[i] /= t_scalefactor;
    //cout <<"t="<<temp->time_evo[i]<<"\n";
  }

  //
  // Offset is not used in the code past here.  All times are in
  // seconds.
  //
  temp->offset = time_offset/t_scalefactor;
  temp->tstart = temp->time_evo[0];       
  temp->tfinish= temp->time_evo[temp->Npt-1];  
  temp->update_freq = update_freq/t_scalefactor;
  temp->t_next_update = max(temp->tstart,t_now);
#ifdef TESTING
  cout <<"\t\t tstart="<<temp->tstart;
  cout <<", next update="<<temp->t_next_update;
  cout <<", and tfinish="<<temp->tfinish<<"\n";
#endif

  //
  // Decide if the wind src is active yet.  If it is, then
  // set up a constant wind source for updating its properties.
  //
  double mdot=0.0, vinf=0.0, vrot=0.0, Twind=0.0, rstar=0.0;
  double xh=0.0, xhe=0.0, xc=0.0, xn=0.0, xo=0.0, xz=0.0, xd=0.0;

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
    interpolate.root_find_linear_vec(temp->time_evo, temp->vrot_evo, t_now, vrot);
    interpolate.root_find_linear_vec(temp->time_evo, temp->R_evo,    t_now, rstar);

    // get tracer values for elements.
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_H_evo, t_now, xh);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_He_evo, t_now, xhe);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_C_evo, t_now, xc);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_N_evo, t_now, xn);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_O_evo, t_now, xo);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_Z_evo, t_now, xz);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_D_evo, t_now, xd);

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

  // set tracer values for elements
  set_element_indices(temp);
  if (temp->i_XH>=0)  trv[temp->i_XH] = xh;
  if (temp->i_XHe>=0) trv[temp->i_XHe]= xhe;
  if (temp->i_XC>=0)  trv[temp->i_XC] = xc;
  if (temp->i_XN>=0)  trv[temp->i_XN] = xn;
  if (temp->i_XO>=0)  trv[temp->i_XO] = xo;
  if (temp->i_XZ>=0)  trv[temp->i_XZ] = xz;
  if (temp->i_XD>=0)  trv[temp->i_XD] = xd;

  // Set B-field of star
  // TODO: Decide how to set this better!  For now pick B=10G at
  //       radius 10 R_sun, and scale with R^-2 for constant flux.
  //
  double Bstar= 10.0*pow_fast(10.0*pconst.Rsun()/rstar,2.0);

  // Now add source using constant wind version.
  stellar_wind::add_source(pos,rad,type,mdot,vinf,vrot,Twind,rstar,Bstar,trv);
  temp->ws = wlist.back();

  // Add evolutionary data to list of wind sources and return
  wdata_evol.push_back(temp);
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
  double mdot=0.0, vinf=0.0, vrot=0.0, Twind=0.0, rstar=0.0;
  double xh=0.0, xhe=0.0, xc=0.0, xn=0.0, xo=0.0, xz=0.0, xd=0.0;

  interpolate.root_find_linear_vec(wd->time_evo, wd->Teff_evo, t_now, Twind);
  interpolate.root_find_linear_vec(wd->time_evo, wd->Mdot_evo, t_now, mdot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vrot_evo, t_now, vrot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vinf_evo, t_now, vinf);
  interpolate.root_find_linear_vec(wd->time_evo, wd->R_evo, t_now, rstar);

  wd->ws->Mdot = mdot;  // already cgs.
  wd->ws->Vinf = vinf;  // this is in cm/s already.
  wd->ws->v_rot = vrot;  // this is in cm/s already.
  wd->ws->Tw   = Twind; // This is in K.
  wd->ws->Rstar = rstar;

  // get tracer values for elements.
  if (wd->i_XH>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_H_evo, t_now, xh);
    wd->ws->tracers[wd->i_XH] = xh;
  }
  if (wd->i_XHe>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_He_evo, t_now, xhe);
    wd->ws->tracers[wd->i_XHe]= xhe;
  }
  if (wd->i_XC>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_C_evo, t_now, xc);
    wd->ws->tracers[wd->i_XC] = xc;
  }
  if (wd->i_XN>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_N_evo, t_now, xn);
    wd->ws->tracers[wd->i_XN] = xn;
  }
  if (wd->i_XO>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_O_evo, t_now, xo);
    wd->ws->tracers[wd->i_XO] = xo;
  }
  if (wd->i_XZ>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_Z_evo, t_now, xz);
    wd->ws->tracers[wd->i_XZ] = xz;
  }
  if (wd->i_XD>=0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_D_evo, t_now, xd);
    wd->ws->tracers[wd->i_XD] = xd;
  }


  // Set B-field of star
  // TODO: Decide how to set this better!  For now pick B=10G at
  //       radius 10 R_sun, and scale with R^-2 for constant flux.
  //
  wd->ws->Bstar= 10.0*pow_fast(10.0*pconst.Rsun()/rstar,2.0);
  
  //
  // Now re-assign state vector of each wind-boundary-cell with
  // updated values.
  //
  for (int i=0; i<wd->ws->ncell; i++) {
    set_wind_cell_reference_state(grid,wd->ws->wcells[i],wd->ws,gamma);
  }

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





