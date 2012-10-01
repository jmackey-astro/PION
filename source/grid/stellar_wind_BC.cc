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
///
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  So VX,VY,VZ are now set by using calls to
///   distance() functions in the grid class.  It doesn't seem to help
///   much.
/// - 2011.01.07 JM: I debugged the geometric grid functions, and now
///   it works very well!  I have a nice spherical expansion.
///
/// - 2011.01.18 JM: Added #def which sets pressure so that Tmin=10K
///
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
/// - 2011.12.01 JM: Switched from spline to linear interpolation for winds.
//------------------------------------------------
//------------ STELLAR WIND CLASS ----------------
//------------------------------------------------
#include "stellar_wind_BC.h"
#include <sstream>
using namespace std;

stellar_wind::stellar_wind()
{
  nsrc=0;
}

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

int stellar_wind::add_source(const double *pos, ///< position (physical units)
           const double  rad, ///< radius (physical units)
           const int    type, ///< type (0=fixed in time,1=slow switch on,2 switch on in radius)
           const double mdot, ///< Mdot (Msun/yr)
           const double vinf, ///< Vinf (km/s)
           const double temp, ///< Wind Temperature (actually p_g.m_p/(rho.k_b))
           const double Rstar, ///< Stellar radius (for T*-->gas pres.).
           const double *trv  ///< Tracer values of wind (if any)
           )
{
  struct wind_source *ws = 0;
  ws = mem.myalloc(ws,1);
  for (int v=0;v<MAX_DIM;v++) 
    ws->ipos[v] = 0;
  ws->id = wlist.size();
  ws->ncell = 0;
  ws->type = type;
  switch (type) {
  case 0: case 3:
    cout <<"Adding constant wind source as id="<<ws->id<<"\n";
    break;
  case 1:
    cout <<"Adding wind source with gradual switch-on over 100kyr, as id="<<ws->id<<"\n";
    break;
  case 2:
    cout <<"Adding wind source with switch on based on r<f*Vinf*t, as id="<<ws->id<<"\n";
    break;
  default:
    rep.error("What type of source is this?  add a new type?",type);
    break;
  }

  for (int v=0;v<SimPM.ndim;v++)       ws->dpos[v] = pos[v];
  for (int v=SimPM.ndim;v<MAX_DIM;v++) ws->dpos[v] = VERY_LARGE_VALUE;

  ws->radius = rad;
  //
  // Mdot and Vinf are passed to the function in Msun/yr and km/s, but
  // are stored internally in cgs units.
  //
  ws->Mdot  = mdot *GS.Msun()/GS.s_per_yr();
  ws->Vinf  = vinf *GS.cm_per_km();

  ws->Tw    = temp;
  ws->Rstar = Rstar;

  ws->tracers=0;
  ws->tracers = mem.myalloc(ws->tracers,SimPM.ntracer);
  for (int v=0;v<SimPM.ntracer; v++) ws->tracers[v] = trv[v];
    
  ws->ipos_set    = false;
  ws->cells_added = false;

  if (!ws->wcells.empty())
    rep.error("wind_source: wcells not empty!",ws->wcells.size());

  //
  // Make sure the source position is compatible with the geometry!
  //
  if (SimPM.coord_sys==COORD_SPH) {
    if (!GS.equalD(ws->dpos[Rsph],0.0))
      rep.error("Spherical symmetry but source not at origin!",ws->dpos[Rsph]);
  }
  if (SimPM.coord_sys==COORD_CYL && SimPM.ndim==2) {
    //
    // Axisymmetry
    //
    if (!GS.equalD(ws->dpos[Rcyl],0.0))
      rep.error("Axi-symmetry but source not at R=0!",ws->dpos[Rsph]);
  }

  wlist.push_back(ws);
  //nsrc++;
  cout <<" Added wind source id="<<nsrc<<" to list of "<<++nsrc<<" elements.\n";
  return ws->id;
}



int stellar_wind::Nsources()
{
  return nsrc;
}



void stellar_wind::get_src_ipos(
        const int id, ///< src id
        int *pos      ///< position vector (output)
        )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);

  struct wind_source *ws = wlist[id];
  if (!(ws->ipos_set)) set_integer_positions(ws);

  for (int v=0;v<SimPM.ndim;v++) pos[v] = ws->ipos[v];
  return;
}



void stellar_wind::set_integer_positions(struct wind_source *ws)
{
  //
  // Need to set position, and convert radius to internal units too.
  //
  if (!grid) rep.error("don't call me until grid is set up!",grid);
  CI.get_ipos_vec(ws->dpos, ws->ipos);
  cout <<"<<<<set_integer_positions: source ";rep.printVec("ipos",ws->ipos,SimPM.ndim);
  cout <<"<<<<set_integer_positions: source ";rep.printVec("dpos",ws->dpos,SimPM.ndim);
  ws->radius /= CI.phys_per_int();
  ws->ipos_set = true;
  return;
}



void stellar_wind::get_src_irad(
        const int id, ///< src id
        double *rad   ///< radius (output)
        )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);
  //
  // Make sure we've converted radius to internal units.
  //
  if (!(wlist[id]->ipos_set))
    set_integer_positions(wlist[id]);

  *rad = wlist[id]->radius;
  return;
}



int stellar_wind::add_cell(
          const int id, ///< src id
           cell *c       ///< cell to add to list.
           )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);
  struct wind_source *WS = wlist[id];
  //
  // Make sure we've converted radius/position to internal units.
  //
  if (!(WS->ipos_set)) set_integer_positions(WS);

  //
  // Setup a wind_cell struct
  //
  struct wind_cell *wc = 0;
  wc = mem.myalloc(wc,1);
  //
  // GEOMETRY: this is a distance from a vertex to a cell--centre, so
  // we call the appropriate function here.
  //
#ifdef GEOMETRIC_GRID
  wc->dist = grid->idistance_vertex2cell(WS->ipos, c);
#else  // GEOMETRIC_GRID
  wc->dist = GS.idistance(c->pos, WS->ipos, SimPM.ndim);
#endif // GEOMETRIC_GRID
  if (wc->dist > WS->radius) {
    rep.warning("stellar_wind::add_cell() cell is outside radius",WS->radius,wc->dist);
    CI.print_cell(c);
    return 1;
  }
  //
  // Now set wc cell pointer to this one.  Also set c->isbd to indicate
  // that it is now boundary data (while also grid data).
  //
  c->isbd=true;
  wc->c = c;

  //
  // Allocate memory for wind_cell reference state vector.
  //
  wc->p = 0; wc->p = mem.myalloc(wc->p, SimPM.nvar);
  
  //
  // Now assign values to the state vector:
  //
  set_wind_cell_reference_state(wc,WS);

  WS->wcells.push_back(wc);
  WS->ncell += 1;
  //cout <<"*** dist="<<wc->dist<<". "; 
  //rep.printVec("Wind BC cell pos",wc->c->pos,SimPM.ndim);
  //rep.printVec("Wind BC cell values", wc->p, SimPM.nvar);
  //cout <<" Added cell: array size="<<WS->wcells.size()<<"\tncell="<<WS->ncell<<"\n";
  return 0;
}

void stellar_wind::set_wind_cell_reference_state(struct wind_cell *wc,
                                                const struct wind_source *WS
                                                )
{
  //
  // In this function we set the density, pressure, velocity, and tracer values
  // for the reference state of the cell.  Every timestep the cell-values will
  // be reset to this reference state.
  //

  //
  // Density at cell position: rho = Mdot/(4.pi.R^2.v_inf) (for 3D)
  // or in 2D (slab-symmetry) rho = Mdot/(2.pi.R.v_inf)
  //
  if (SimPM.ndim==2 && SimPM.coord_sys==COORD_CRT) {
    //
    // 2D slab symmetry --> 1/r force laws and density profile
    //
    wc->p[RO] = WS->Mdot/(WS->Vinf*2.0*M_PI*wc->dist*CI.phys_per_int());
    //
    // Set pressure based on wind density/temperature at the stellar radius:
    // rho_star = Mdot/(2.pi.R_star.v_inf), p_star = rho_star.k.T_star/(mu.m_p)
    // and then p(r) = p_star (rho(r)/rho_star)^gamma
    // ******************************************************************************
    // *********** WARNING MU=1 HERE, PROBABLY SHOULD BE O.6 (IONISED) 1.3 (NEUTRAL).
    // ******************************************************************************
    //
    wc->p[PG] = GS.kB()*WS->Tw/GS.m_p();
    wc->p[PG]*= exp((SimPM.gamma-1.0)*log(2.0*M_PI*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((SimPM.gamma)*log(wc->p[RO]));
#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
    //
    // Set the minimum temperature to be 10K in the wind...
    //
    if (MP) {
      if (MP->Temperature(wc->p,SimPM.gamma) <SimPM.EP.MinTemperature) {
        MP->Set_Temp(wc->p,SimPM.EP.MinTemperature,SimPM.gamma);
      }
    }
    else {
      wc->p[PG] = max(wc->p[PG], 10.0*wc->p[RO]*GS.kB()/GS.m_p());
    }
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  }


  else {
    //
    // 3D geometry, so either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
    // rho = Mdot/(4.pi.R^2.v_inf) 
    //
    wc->p[RO] = 1.0/(wc->dist*CI.phys_per_int());
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
    wc->p[PG] = GS.kB()*WS->Tw/GS.m_p();
    wc->p[PG]*= exp((SimPM.gamma-1.0)*log(4.0*M_PI*WS->Rstar*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((SimPM.gamma)*log(wc->p[RO]));
#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
    //
    // Set the minimum temperature to be 10K in the wind...
    //
    if (MP) {
      if (MP->Temperature(wc->p,SimPM.gamma) <SimPM.EP.MinTemperature) {
        MP->Set_Temp(wc->p,SimPM.EP.MinTemperature,SimPM.gamma);
      }
    }
    else {
      wc->p[PG] = max(wc->p[PG], 12.0*wc->p[RO]*GS.kB()/GS.m_p());
    }
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  }


  //
  // NOW VELOCITIES: These should be cell-average values, which are
  // the values at the centre-of-volume, so we call the geometry-aware
  // grid functions rather than GS.idistance().
  //
  cell *c = wc->c;
#ifdef GEOMETRIC_GRID
  wc->p[VX] = WS->Vinf*grid->idifference_vertex2cell(WS->ipos,c,XX)/wc->dist;
  if (SimPM.ndim>1)
    wc->p[VY] = WS->Vinf*grid->idifference_vertex2cell(WS->ipos,c,YY)/wc->dist;
  else
    wc->p[VY] = 0.0;
  if (SimPM.ndim>2)
    wc->p[VZ] = WS->Vinf*grid->idifference_vertex2cell(WS->ipos,c,ZZ)/wc->dist;
  else
    wc->p[VZ] = 0.0;

#else  // not GEOMETRIC_GRID

  wc->p[VX] = WS->Vinf*(c->pos[XX]-WS->ipos[XX])/wc->dist;
  if (SimPM.ndim>1)
    wc->p[VY] = WS->Vinf*(c->pos[YY]-WS->ipos[YY])/wc->dist;
  else
    wc->p[VY] = 0.0;
  if (SimPM.ndim>2)
    wc->p[VZ] = WS->Vinf*(c->pos[ZZ]-WS->ipos[ZZ])/wc->dist;
  else
    wc->p[VZ] = 0.0;
#endif // GEOMETRIC_GRID

  if (SimPM.eqntype!=EQEUL && SimPM.eqntype!=EQEUL_EINT)
    rep.error("Need to code B into winds model!",SimPM.eqntype);


  // update tracers
  for (int v=0;v<SimPM.ntracer;v++)
    wc->p[SimPM.ftr+v] = WS->tracers[v];
  //
  // HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK!
  // Assume the first tracer variable is the H+ ion fraction, and set it so
  // that it goes from y=1 at T>12500K to y=0 at T<10000, with linear
  // interpolation.
  //
#ifdef HACK_WARNING
#error "REMOVE HACK setting ion fraction of H+ in winds"
#endif
  if (SimPM.ntracer>0) {
    if      (WS->Tw >1.25e4)
      wc->p[SimPM.ftr] = 1.0;
    else if (WS->Tw <1.00e4)
      wc->p[SimPM.ftr] = 0.0;
    else
      wc->p[SimPM.ftr] = (WS->Tw-1.0e4)*4e-4;
  }
  return;
}


int stellar_wind::set_num_cells(const int id, ///< src id
        const int nc  ///< number of cells.
        )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);

  if (nc != wlist[id]->ncell) {
    rep.warning("stellar_wind::set_num_cells() COUNTING INCONSISTENCY!!",nc,wlist[id]->ncell);
    return 1;
  }
  return 0;
}

int stellar_wind::get_num_cells(const int id ///< src id
        )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);
  return wlist[id]->ncell;
}





int stellar_wind::set_cell_values(const int id,  ///< src id
          const double t ///< simulation time
          )
{
  if (id<0 || id>=nsrc)
    rep.error("bad src id",id);

  //
  // Since we have the list here we may as well go through
  // every cell in one go and update them all.  This is hard-coded
  // for non-varying winds, so if I want to model that I'll need to
  // add quite a bit here.
  //
  struct wind_source *WS = wlist[id];
  //cout <<"updating source "<<id<<" which has "<<WS->ncell<<" cells.\n";
  if      (WS->type==WINDTYPE_CONSTANT ||
           WS->type==WINDTYPE_EVOLVING) {
    //
    // Constant wind (type==0 or type==3)
    //
    //cout <<"UPDATING WIND -- CONSTANT ID=0\n";
    for (int i=0; i<WS->ncell; i++) {
      for (int v=0;v<SimPM.nvar;v++)
        WS->wcells[i]->c->P[v]  = WS->wcells[i]->p[v];
      for (int v=0;v<SimPM.nvar;v++)
        WS->wcells[i]->c->Ph[v] = WS->wcells[i]->p[v];
      //cout <<"TW="<<WS->Tw<<",  ftr val="<<WS->wcells[i]->c->Ph[SimPM.ftr]<<"\n";
    }
  }

  else if (WS->type==WINDTYPE_SWITCHON) {
    //
    // Slow switch on over 10^5 years (type==1)
    //
    //cout <<"UPDATING WIND -- SLOW SWITCH ON OVER TIME ID=1\n";
    rep.error("Don't use switch-on wind -- no advantage",WS->type);
    cell *c=0;
    double tfrac = t/3.16e12;
    bool slow=false; if (tfrac<1) slow=true;

    for (int i=0; i<WS->ncell; i++) {
      c = WS->wcells[i]->c;

      for (int v=0;v<SimPM.nvar;v++)
        c->P[v]  = WS->wcells[i]->p[v];
      //
      // Slow down velocity if t<10^5 yrs
      //
      if (slow) {
        //cout <<"tfrac="<<tfrac<<" max-vel="<<c->P[VX]*tfrac<<"\n";
        c->P[VX] *= tfrac;
        c->P[VY] *= tfrac;
        c->P[VZ] *= tfrac;
      }
      //
      // Set Ph[]=P[]
      //
      for (int v=0;v<SimPM.nvar;v++)
        c->Ph[v] = c->P[v];
    }
  }
  else if (WS->type==WINDTYPE_EXPAND) {
    //
    // Only apply boundary condition if r<f*v*t (where f<=1)
    // EXPERIMENTAL CODE!!!
    //
    //    cout <<"UPDATING WIND -- SLOW SWITCH ON SPATIALLY ID=2\n";
    rep.error("Don't use expanding wind -- no advantage.",WS->type);
    cell *c=0;
    double vt= 0.25*WS->Vinf*t/CI.phys_per_int();
    double d2= max(2.5,min(0.5*vt,0.25*WS->radius));
    vt = max(vt,d2);
    //cout <<"vt="<<vt<<"\tradius="<<WS->radius<<"\n";
    for (int i=0; i<WS->ncell; i++) {
      c = WS->wcells[i]->c;
      if (WS->wcells[i]->dist < vt) {
  //if (WS->wcells[i]->dist < d2) {
  c->P[RO] = WS->wcells[i]->p[RO];
  c->P[PG] = WS->wcells[i]->p[PG];
  //}
  c->P[VX] = WS->wcells[i]->p[VX];
  c->P[VY] = WS->wcells[i]->p[VY];
  c->P[VZ] = WS->wcells[i]->p[VZ];
  for (int v=SimPM.ftr;v<SimPM.nvar;v++)
    c->P[v]  = WS->wcells[i]->p[v];
  //
  // Set Ph[]=P[]
  //
  for (int v=0;v<SimPM.nvar;v++)
    c->Ph[v] = c->P[v];
      }
      // else {
      //   c->P[VX] = WS->wcells[i]->p[VX]*vt/WS->wcells[i]->dist;
      //   c->P[VY] = WS->wcells[i]->p[VY]*vt/WS->wcells[i]->dist;
      //   c->P[VZ] = WS->wcells[i]->p[VZ]*vt/WS->wcells[i]->dist;
      // }
    }
  }
  else {
    rep.error("set_cell_values(): What type of source is this?  add a new type?",WS->type);
  }
    
  return 0;
}

void stellar_wind::get_src_posn(const int id, ///< src id
        double *x   ///< position vector (output)
        )
{
  for (int v=0;v<SimPM.ndim;v++) x[v] = wlist[id]->dpos[v];
}
void stellar_wind::get_src_Mdot(const int id, ///< src id
        double *x   ///< mdot (output)
        )
{
  *x = wlist[id]->Mdot *GS.s_per_yr()/GS.Msun();
}
void stellar_wind::get_src_drad(const int id, ///< src id
             double *x   ///< radius (output)
             )
{
  if (!(wlist[id]->ipos_set)) 
    *x = wlist[id]->radius;
  else 
    *x = wlist[id]->radius *CI.phys_per_int();
  return;
}
void stellar_wind::get_src_Vinf(const int id, ///< src id
        double *x   ///< Vinf (output)
        )
{
  *x = wlist[id]->Vinf/GS.cm_per_km();
}
void stellar_wind::get_src_Tw(const int id, ///< src id
            double *x   ///< Stellar Radius (output)
            )
{
  *x = wlist[id]->Tw;
}
void stellar_wind::get_src_Rstar(const int id, ///< src id
         double *x   ///< Stellar radius (output)
         )
{
  *x = wlist[id]->Rstar;
}
void stellar_wind::get_src_trcr(const int id, ///< src id
        double *x   ///< tracers (output)
        )
{
  for (int v=0;v<SimPM.ntracer; v++)
    x[v] = wlist[id]->tracers[v];
}
void stellar_wind::get_src_type(const int id, ///< src id
        int *x   ///< type of wind (=0 for now) (output)
        )
{
  *x=wlist[id]->type;
}
//------------------------------------------------


// ------------------------------------------------------------------
// ----------  STELLAR WIND WITH STELLAR EVOLUTION ------------------
// ------------------------------------------------------------------

stellar_wind_evolution::stellar_wind_evolution()
{
  cout <<"Stellar wind with time evolution, constructor.\n";
}

stellar_wind_evolution::~stellar_wind_evolution()
{
  cout <<"Stellar wind with time evolution, destructor.\n";
  //
  // Delete arrays
  //
  // delete each wdata_evol[] element, making sure to delete the arrays in each
  // element first (spline arrays, pos, tracers).
  //
  struct evolving_wind_data *wd=0;
  while (wdata_evol.size() >0) {
    wd = wdata_evol.back();
    wd->t = mem.myfree(wd->t);
    wd->mdot  = mem.myfree(wd->mdot);
    wd->mdot2 = mem.myfree(wd->mdot2);
    wd->vinf  = mem.myfree(wd->vinf);
    wd->vinf2 = mem.myfree(wd->vinf2);
    wd->Teff  = mem.myfree(wd->Teff);
    wd->Teff2 = mem.myfree(wd->Teff2);
    wdata_evol.pop_back();
    wd = mem.myfree(wd);
  }
}

int stellar_wind_evolution::add_source(const double *pos, ///< position (physical units)
           const double  rad, ///< radius (physical units)
           const int    type, ///< type (0=fixed in time,1=slow switch on,2 switch on in radius)
           const double mdot, ///< Mdot (Msun/yr)
           const double vinf, ///< Vinf (km/s)
           const double Twnd, ///< Wind Temperature (actually p_g.m_p/(rho.k_b))
           const double Rstar, ///< Stellar radius (for T*-->gas pres.).
           const double *trv  ///< Tracer values of wind (if any)
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
  temp->Nspl = 0;
  temp->t = 0;
  temp->mdot  = 0;
  temp->mdot2 = 0;
  temp->vinf  = 0;
  temp->vinf2 = 0;
  temp->Teff  = 0;
  temp->Teff2 = 0;
  temp->offset = 0.0;
  temp->tstart = SimPM.starttime -1.0e10; // so it has already started.
  temp->tfinish= 2.0*SimPM.finishtime;    // so it keeps going to end of sim.
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



int stellar_wind_evolution::add_evolving_source(
  const double *pos,        ///< position (physical units).
  const double  rad,        ///< radius (physical units).
  const int    type,        ///< type (must be 3, for variable wind).
  const double Rstar,       ///< Radius at which to get gas pressure from Teff
  const double *trv,        ///< Any (constant) wind tracer values.
  const string infile,      ///< file name to read data from.
  const double time_offset, ///< time offset = [t(sim)-t(wind_file)]
  const double t_now,       ///< current simulation time, to see if src is active.
  const double update_freq, ///< frequency with which to update wind properties.
  const double t_scalefactor ///< wind evolves this factor times faster than normal
  )
{
  if (type != WINDTYPE_EVOLVING) {
    rep.error("Bad wind type for evolving stellar wind! If one is, all must be.",type);
  }
  //
  // First we will read the file, and see when the source should switch on in 
  // the simulation (it may not be needed for a while).
  //
  cout <<"\t\tsw-evo: adding source from file "<<infile<<"\n";
  ifstream wf;
  wf.open(infile.c_str());
  if (!wf.is_open()) {cerr<<"Error opening file.\n";return 1;}

  string line;
  getline(wf, line);
  while (line.empty() == true || line.substr(0,1) == "#") {
    getline(wf, line);
  }
  //
  // Now first non-comment line should tell me Nspline
  //
  istringstream iss2(line);
  string junk;
  int Nspl=0;
  iss2 >> junk >> Nspl;
  cout <<"\t\tgetting Nspl:: "<<junk<<": "<<Nspl<<"\n";
  if (!isfinite(Nspl) || Nspl<2) {
    rep.error("Bad Nspline in stellar_wind_evolution",Nspl);
  }
  //
  // From here on in we read in Nspl lines, and assign the values to arrays,
  // which we now need to set up.  We need time[], vinf[], vinf2[], mdot[],
  // mdot2[], Teff, Teff2[].
  //
  double *t=0, *vi=0, *vi2=0, *md=0, *md2=0, *Tf=0, *Tf2=0;
  t   = mem.myalloc(t  , Nspl);
  vi  = mem.myalloc(vi , Nspl);
  vi2 = mem.myalloc(vi2, Nspl);
  md  = mem.myalloc(md , Nspl);
  md2 = mem.myalloc(md2, Nspl);
  Tf  = mem.myalloc(Tf , Nspl);
  Tf2 = mem.myalloc(Tf2, Nspl);
  int i=0; // iterator for which element of spline we are at.

  while (!wf.eof()) {
    getline(wf, line);
    istringstream iss(line);
    if (!line.empty()) {
      //cout <<"line["<<i<<"]: "<<line<<"\n";
      //cout <<"iss ["<<i<<"]: "<<iss.str()<<"\n";
      // File format:
      // Time/yr log10(L/Lsun) log10(Teff/K) log10(Mdot/Msun/yr) log10(vinf/cm/s)
      //
      if (i>= Nspl) {
        cout <<"ERROR: line: "<<line<<"\n";
        rep.error("Too many lines in file for given Nspl", Nspl);
      }
      iss >> t[i] >> junk >> Tf[i] >> md[i] >> vi[i];
      //double temp;
      //iss >> temp; cout <<"j="<<temp<<"\n";
      //cout <<iss.str()<<"\n";
      //cout <<t[i] <<"  "<< junk <<"  "<< Tf[i] <<"  "<< md[i] <<"  "<< vi[i]<<"\n";
      i++;
    }
  }
  if (i!=Nspl) {
    rep.error("Too few lines in file!",Nspl-i);
  }

  //
  // Next we set up the spline interpolation, first modifying the time array
  // using the time-offset so that it has the same zero offset as the simulation
  // time.
  //
  for (i=0; i<Nspl; i++) {
    //cout <<"i="<<i<<"  t[i]="<<t[i];
    // times in the file are measured in years, so offset should be in years.
    t[i] += time_offset;
    // scale times by scale factor.
    t[i] /= t_scalefactor;
    //cout <<" t[i]="<<t[i]<<"\n";
  }
  //GS.spline(t,md,Nspl,1.e99,1.e99,md2);
  //GS.spline(t,vi,Nspl,1.e99,1.e99,vi2);
  //GS.spline(t,Tf,Nspl,1.e99,1.e99,Tf2);


  //
  // Some properties of the wind source are specific to this module, such as the
  // spline interpolation array, the number of points in the spline array, and
  // the timings (offset, update frequency).  They are stored in local data.
  struct evolving_wind_data *temp=0;
  temp = mem.myalloc(temp,1);
  temp->Nspl = Nspl;
  temp->t = t;  // in years
  temp->mdot  = md;
  temp->mdot2 = md2;
  temp->vinf  = vi;
  temp->vinf2 = vi2;
  temp->Teff  = Tf;
  temp->Teff2 = Tf2;
  //
  // Offset is not used in the code past here.  It's just here for I/O b/c a
  // restart will need to read the data-file again.
  // We reset all variables to be in seconds here.  So the spline/splint uses
  // years, but everything else is in seconds.  Note that the global SWP struct
  // still has times in years though.
  //
  temp->offset = time_offset*GS.s_per_yr()/t_scalefactor; // now in seconds
  temp->tstart = t[0]       *GS.s_per_yr(); // now in seconds (already scaled)
  temp->tfinish= t[Nspl-1]  *GS.s_per_yr(); // now in seconds (already scaled)
  temp->update_freq = update_freq*GS.s_per_yr()/t_scalefactor; // now in seconds
  temp->t_next_update = max(temp->tstart,t_now);
  cout <<"\t\t tstart="<<temp->tstart<<", next update="<<temp->t_next_update<<", and tfinish="<<temp->tfinish<<"\n";

  //
  // We need to decide if the wind src is active yet.  If it is, then we also set
  // up a constant wind source for updating its properties.  We set it to be active
  // if the current time is within update_freq of tstart.
  //
  double mdot=0.0, vinf=0.0, Twind=0.0;
  if (((t_now+temp->update_freq)>temp->tstart || GS.equalD(temp->tstart, t_now))
      && t_now<temp->tfinish) {
    temp->is_active = true;
    //
    // Get the current values for mdot, vinf, Teff, and setup a wind source
    // using the constant-wind function.
    //
    //GS.splint(t, Tf, Tf2, Nspl, t_now/GS.s_per_yr(), &Twind);
    //GS.splint(t, md, md2, Nspl, t_now/GS.s_per_yr(), &mdot);
    //GS.splint(t, vi, vi2, Nspl, t_now/GS.s_per_yr(), &vinf);
    GS.root_find_linear(t, Tf, Nspl, t_now/GS.s_per_yr(), &Twind);
    GS.root_find_linear(t, md, Nspl, t_now/GS.s_per_yr(), &mdot);
    GS.root_find_linear(t, vi, Nspl, t_now/GS.s_per_yr(), &vinf);
    cout <<"Source is Active\n";
  }
  else {
    cout <<"WARNING: Source is not yet active: tnow="<<t_now<<", tstart=";
    cout <<temp->tstart<<". Setting wind source to INACTIVE.\n";
    temp->is_active = false;
    mdot=-100.0; vinf=-100.0; Twind=-100.0;
  }
  //
  // Need to convert units.  add_source expects mdot in msun/yr, vinf in km/s, 
  // and Twind in K.
  //
  mdot = exp(GS.ln10()*mdot);
  vinf = exp(GS.ln10()*vinf) /GS.cm_per_km();
  Twind = exp(GS.ln10()*Twind);
  //
  // Now add source using constant wind version.
  //
  stellar_wind::add_source(pos,rad,type,mdot,vinf,Twind,Rstar,trv);
  temp->ws = wlist.back();

  //
  // So now we have all of the properties of the wind source in the struct
  // evolving_wind_data 'temp', and we have added it to the list of constant
  // wind sources, wlist[].  We can now add temp to the list of evolving
  // wind sources and return (so that wlist[i] and wdata_evol[i] point to the
  // same source).
  //
  wdata_evol.push_back(temp);
  //NSRC_TOTAL++;

  return 0;
}

void stellar_wind_evolution::update_source(struct evolving_wind_data *wd,
                                           const double t_now
                                           )
{
  //
  // We have a source that needs updating.  If it is not active, and needs
  // activating then we set that.
  //
  if (!wd->is_active) {
    cout <<"stellar_wind_evolution::update_source() activating source id=";
    cout << wd->ws->id <<" at Simulation time t="<<t_now<<"\n";
    rep.printVec("Source position",wd->ws->dpos,SimPM.ndim);
    wd->is_active=true;
  }

  if (t_now < wd->tstart) {
    rep.error("Requested updating source, but it hasn't switched on yet!",wd->tstart-t_now);
  }

  //wd->t_next_update = t_now + wd->update_freq;
  wd->t_next_update = t_now+SimPM.dt;
  wd->t_next_update = min(wd->t_next_update, wd->tfinish);
  
  //cout <<"updating source: current [mdot,vinf,Tw]=[";
  //cout <<wd->ws->Mdot<<", "<<wd->ws->Vinf<<", "<<wd->ws->Tw<<"]\n";
  //
  // Now we update Mdot, Vinf, Teff by spline interpolation.
  //
  double mdot=0.0, vinf=0.0, Twind=0.0;
  //GS.splint(wd->t, wd->Teff, wd->Teff2, wd->Nspl, t_now/GS.s_per_yr(), &Twind);
  //GS.splint(wd->t, wd->mdot, wd->mdot2, wd->Nspl, t_now/GS.s_per_yr(), &mdot);
  //GS.splint(wd->t, wd->vinf, wd->vinf2, wd->Nspl, t_now/GS.s_per_yr(), &vinf);
  GS.root_find_linear(wd->t, wd->Teff, wd->Nspl, t_now/GS.s_per_yr(), &Twind);
  GS.root_find_linear(wd->t, wd->mdot, wd->Nspl, t_now/GS.s_per_yr(), &mdot);
  GS.root_find_linear(wd->t, wd->vinf, wd->Nspl, t_now/GS.s_per_yr(), &vinf);
  //
  // Assign new values to wd->ws (the wind source struct), converting from log10
  // to actual values, and also unit conversions.
  //
  wd->ws->Mdot = exp(GS.ln10()*mdot) *GS.Msun()/GS.s_per_yr();  // this was log10(msun/yr)
  wd->ws->Vinf = exp(GS.ln10()*vinf);  // this is in log10(cm/s) already.
  wd->ws->Tw   = exp(GS.ln10()*Twind); // This is in log10(K).

  //cout <<"updating source: updated [mdot,vinf,Tw]=[";
  //cout <<wd->ws->Mdot<<", "<<wd->ws->Vinf<<", "<<wd->ws->Tw<<"]\n";

  //
  // Now re-assign state vector of each wind-boundary-cell with updated values.
  //
  for (int i=0; i<wd->ws->ncell; i++) {
    set_wind_cell_reference_state(wd->ws->wcells[i],wd->ws);
  }

  //
  // Now the source is updated, and the reference states are all set for the new
  // values, and the next update time has been set.  So we can return.
  //
  return;
}

int stellar_wind_evolution::set_cell_values(
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
  // Check if the source has reached its last time, because we want to stop  
  // the simulation in this case!
  //
  if (t_now >= wd->tfinish) {
    cout <<"stellar_wind_evolution::set_cell_values() reached final time.";
    cout <<"\t t(now)="<<t_now<<" and t(final)="<<wd->tfinish<<"\n";
    SimPM.maxtime=true;
  }
  //
  // Check if the source needs to be updated, and if so, update it (but only if 
  // we haven't reached the end of the sim, b/c we could run off the end of a 
  // spline interpolation, which can have very bad results).
  // Not too worried about testing for equality b/c we will be at most one
  // timestep wrong, which is not going to be serious.
  //
  else if (t_now >= wd->t_next_update) {
    update_source(wd, t_now);
  }

  //
  // So now we have an up-to-date source.  We check if it is active, and if so we
  // update the wind cell values using the constant wind functions.  If not, then
  // we ignore it and return.
  //
  if (wd->is_active) {
    err += stellar_wind::set_cell_values(id,t_now);
  }

  return err;
}



