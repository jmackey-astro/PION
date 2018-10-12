/// \file uniform_grid_pllel.cpp
/// \author Jonathan Mackey
/// 
/// Definitions for parallel uniform grid.
/// 
/// Modifications:\n
/// - 2007-11-06 started to get it going.
/// - 2010-01-22 JM: Changed RT source finding if we have non-cell-centred sources.
/// - 2010-01-26 JM: Debugging 3D RT for vertex-centred sources.
/// - 2010-02-03 JM: changed variable names in destructor ('i' was defined twice!)
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; Added oneway-outflow BC
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.  Switched 'endl' statements for '\n'.
/// - 2010.11.19 JM: Got rid of testing myalloc() myfree() functions.
/// - 2010.12.07 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  Will probably keep it since it is the way
///   things will go eventually.  The new grid classes have extra
///   function for the distance between two points, two cells, and
///   between a vertex and a cell.
/// - 2011.02.24 JM: Worked on generalising multi-source radiative transfer.
///    Still some way to go to getting it functional.  Wraps old setup() 
///    function for point sources into a new function.  The setup function now
///    calls either this new function, or another new one for sources at infinity.
/// - 2011.02.25 JM: removed HCORR ifdef around new code. 
///     Added possibility of multiple sources in send/recv.
/// - 2011.02.28 JM: Got rid of RSP radiation-sources-parameters class.
/// - 2011.03.02 JM: Fixed bugs in diffuse radiation boundary comms.
/// - 2011.03.21 JM: Got rid of zero-ing of column densities.  It is done in
///    the cell constructor, so we don't need to worry about it.
/// - 2011.03.22 JM: small bugfixes.
/// - 2011.04.22 JM: Fixed setup_recv_boundaries() functions so that new boundary cells
///    are only created if they don't already exist.  It now works for multiple point
///    sources, and I don't think they need to be at the same place.
/// - 2013.09.05 JM: Debugged for new get/set col functions.
/// - 2015.01.28 JM: Removed GEOMETRIC_GRID (it is default now), and
///    updated for the new code structure.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
/// - 2016.03.21 JM: Worked on simplifying RT boundaries (first hack!
///    compiles but buggy...) 03.24:fixed some bugs, redefined dir_XP
/// - 2017.11.07-22 JM: updating boundary setup.
/// - 2018.08.09 JM: moved all boundary stuff to new classes.
/// - 2018.10.12 JM: added flux_update functions for BC89 refinement.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"


#include "grid/uniform_grid_pllel.h"
#include "microphysics/microphysics_base.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef PARALLEL



//-------------------------------------------------------------
//------------------- CARTESIAN GRID START --------------------
//-------------------------------------------------------------



// ##################################################################
// ##################################################################



UniformGridParallel::UniformGridParallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn, ///< array of minimum values of x,y,z for this grid.
      double *xp, ///< array of maximum values of x,y,z for this grid.
      int *nc,    ///< array of number of cells in x,y,z directions.
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{

#ifdef TESTING
  cout <<"UniformGridParallel constructor.\n";
  rep.printVec("Local Xmin",xn,nd);
  rep.printVec("Local Xmax",xp,nd);
  rep.printVec("Local Npt ",nc,nd);

  rep.printVec("SIM Xmin ", Sim_xmin, G_ndim);
  rep.printVec("SIM Xmax ", Sim_xmax, G_ndim);
  rep.printVec("SIM Range", Sim_range,G_ndim);

  rep.printVec("SIM iXmin ", Sim_ixmin, G_ndim);
  rep.printVec("SIM iXmax ", Sim_ixmax, G_ndim);
  rep.printVec("SIM iRange", Sim_irange,G_ndim);

  cout <<"UniformGridParallel constructor done.\n";
#endif

  return;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::setup_flux_recv(
      class SimParams &par, ///< simulation params (including BCs)
      const int l           ///< level to receive from
      )
{
//#ifdef DEBUG_NG
  cout <<"UniformGridParallel::setup_flux_recv() recv from level=";
  cout <<l<<"\n";
//#endif

  if (par.levels[0].MCMD.get_nproc()==1) {
    cout <<"setup_flux_recv(): nproc=1, calling serial version.\n";
    return UniformGrid::setup_flux_recv(par,l);
  }
  else return 0;

  //
  // Get size of interface region and number of cells.
  //
  size_t nc  = 1; // number of cells in each interface
  int ixmin[MAX_DIM], ixmax[MAX_DIM], ncell[MAX_DIM]; // interface
  int lxmin[MAX_DIM], lxmax[MAX_DIM]; // finer grid
  bool recv[2*G_ndim];  // whether to get data in this direction
  size_t nel[2*G_ndim]; // number of interfaces in each direction
  struct flux_interface *fi = 0;
  CI.get_ipos_vec(par.levels[l].Xmin, lxmin);
  CI.get_ipos_vec(par.levels[l].Xmax, lxmax);

  // Modify the following so that there can be up to 2^Ndim child
  // grids, all of which can have an external level boundary.  Should
  // be straightforward...

  // define interface region of fine and coarse grids, and whether
  // each direction is to be included or not.  Note that this allows
  // for a fine grid that is not completely encompassed by the coarse
  // grid.
  for (int v=0;v<G_ndim;v++) {
    ixmin[v] = std::max(G_ixmin[v], lxmin[v]);
    if (G_ixmin[v] < lxmin[v]) recv[2*v] = true;
    else                       recv[2*v] = false;
    
    ixmax[v] = std::min(G_ixmax[v], lxmax[v]);
    if (G_ixmax[v] > lxmax[v]) recv[2*v+1] = true;
    else                       recv[2*v+1] = false;

    ncell[v] = (ixmax[v]-ixmin[v])/G_idx;
    if ( (ixmax[v]-ixmin[v]) % G_idx !=0) {
      rep.error("interface region not divisible!",ixmax[v]-ixmin[v]);
    }
  }
  for (int v=0;v<2*G_ndim;v++) nel[v]=0;

  // different number of interfaces depending on dimensionality.
  switch (G_ndim) {
  case 1:
    if (recv[XN]) nel[XN] = 1;
    if (recv[XP]) nel[XP] = 1;
    break;
  case 2:
    if (recv[XN]) nel[XN] = ncell[YY];
    if (recv[XP]) nel[XP] = ncell[YY];
    if (recv[YN]) nel[YN] = ncell[XX];
    if (recv[YP]) nel[YP] = ncell[XX];
    break;
  case 3:
    if (recv[XN]) nel[XN] = ncell[YY]*ncell[ZZ];
    if (recv[XP]) nel[XP] = ncell[YY]*ncell[ZZ];
    if (recv[YN]) nel[YN] = ncell[XX]*ncell[ZZ];
    if (recv[YP]) nel[YP] = ncell[XX]*ncell[ZZ];
    if (recv[ZN]) nel[ZN] = ncell[XX]*ncell[YY];
    if (recv[ZP]) nel[ZP] = ncell[XX]*ncell[YY];
    break;
  default:
    rep.error("bad ndim in setup_flux_recv",G_ndim);
    break;
  }

  // initialize arrays
  flux_update_recv.resize(2*G_ndim);
  for (int v=0; v<2*G_ndim; v++) {
    if (recv[v] == true) {
      flux_update_recv[v].resize(nel[v]);
    }
    else {
      flux_update_recv[v].resize(1);
      flux_update_recv[v][0] = 0;
    }
    for (size_t i=0; i<nel[v]; i++) {
      flux_update_recv[v][i] = mem.myalloc(flux_update_recv[v][i],1);
      fi = flux_update_recv[v][i];
      fi->Ncells = nc;
      fi->c.resize(nc);
      fi->area.resize(nc);
      fi->flux = mem.myalloc(fi->flux,G_nvar);
      for (int v=0;v<G_nvar;v++) fi->flux[v]=0.0;
    }
  }

  // For each interface, find the cell that is outside the fine grid
  // and that includes the interface.
  for (int v=0; v<2*G_ndim; v++) {
    if (recv[v]) {
      add_cells_to_face(static_cast<enum direction>(v),ixmin,ixmax,ncell,1,flux_update_recv[v]);
    }
  }

  return 0;
}




// ##################################################################
// ##################################################################



int UniformGridParallel::setup_flux_send(
      class SimParams &par, ///< simulation params (including BCs)
      const int l           ///< level to send to
      )
{
//#ifdef DEBUG_NG
  cout <<" UniformGridParallel::setup_flux_send() send to level=";
  cout <<l<<"\n";
//#endif

  if (par.levels[0].MCMD.get_nproc()==1) {
    cout <<"setup_flux_recv(): nproc=1, calling serial version.\n";
    return UniformGrid::setup_flux_send(par,l);
  }
  else return 0;

}


// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//------------------- CARTESIAN GRID END ----------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID START ------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_cyl_parallel::uniform_grid_cyl_parallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn,
      double *xp,
      int *nc,
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  uniform_grid_cyl(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel:: cyl. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim>2)
    rep.error("Need to write code for 3 dimensions",G_ndim);

#ifdef TESTING
  cout <<"cylindrical grid: dR="<<G_dx<<"\n";
#endif 
  return;
}



// ##################################################################
// ##################################################################


uniform_grid_cyl_parallel::~uniform_grid_cyl_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel destructor. Present and correct!\n";
#endif 
}



// ##################################################################
// ##################################################################


double uniform_grid_cyl_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  //cout <<" Cell radius: "<< R_com(c,G_dx)/CI.phys_per_int() +G_ixmin[Rcyl];
  //rep.printVec("  cell centre",c->pos,G_ndim);
  return (R_com(c,G_dx)-Sim_xmin[Rcyl])/CI.phys_per_int() +Sim_ixmin[Rcyl];
}



// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID END --------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- SPHERICAL GRID START --------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_sph_parallel::uniform_grid_sph_parallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn,
      double *xp,
      int *nc,
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd),
  uniform_grid_sph(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel:: sph. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim!=1)
    rep.error("Need to write code for >1 dimension",G_ndim);

#ifdef TESTING
  cout <<"spherical grid: dr="<<G_dx<<"\n";
#endif 
  return;
}


// ##################################################################
// ##################################################################



uniform_grid_sph_parallel::~uniform_grid_sph_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel destructor. Present and correct!\n";
#endif 
}


// ##################################################################
// ##################################################################



double uniform_grid_sph_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  return (R_com(c,G_dx)-Sim_xmin[Rsph])/CI.phys_per_int() +Sim_ixmin[Rsph];
}


// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------
#endif //PARALLEL
