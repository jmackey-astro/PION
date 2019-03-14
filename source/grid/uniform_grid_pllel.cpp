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


//#define TEST_BC89FLUX


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
      double *lev_xn, // level xmin
      double *lev_xp, // level xmax
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp)
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
      const int lp1         ///< level to receive from
      )
{
#ifdef TEST_BC89FLUX
  cout <<"UniformGridParallel::setup_flux_recv() recv from level=";
  cout <<lp1<<"\n";
#endif
  int l = lp1-1;  // my level

  // if only one MPI process, then no send/recv necessargy and we
  // call the serial version:
  if (par.levels[0].MCMD.get_nproc()==1) {
#ifdef TEST_BC89FLUX
    cout <<"setup_flux_recv(): nproc=1, calling serial version.\n";
#endif
    int err = UniformGrid::setup_flux_recv(par,lp1);
    rep.errorTest("UniformGrid flux setup",0,err);
    err = flux_update_recv.size();
    for (int i=0;i<err;i++) {
      flux_update_recv[i].rank.push_back(0);
      flux_update_recv[i].dir = i;
    }
    return 0;
  }

  //
  // Get size of interface region and number of cells.
  //
  size_t nc  = 1; // number of cells in each interface
  int ixmin[MAX_DIM], ixmax[MAX_DIM], ncell[MAX_DIM]; // interface region
  int lxmin[MAX_DIM], lxmax[MAX_DIM]; // finer grid
  int dxmin[MAX_DIM], dxmax[MAX_DIM]; // full domain
  struct flux_interface *fi = 0;
  CI.get_ipos_vec(par.levels[lp1].Xmin, lxmin);
  CI.get_ipos_vec(par.levels[lp1].Xmax, lxmax);
  CI.get_ipos_vec(par.levels[0].Xmin, dxmin);
  CI.get_ipos_vec(par.levels[0].Xmax, dxmax);

  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int nchild = MCMD->child_procs.size();
//#ifdef TEST_BC89FLUX
  cout <<"UniformGridParallel::setup_flux_recv: "<<nchild<<" child grids\n";
//#endif

  // Two cases: if there are children, then part or all of grid is
  // within l+1 level.  If there are no children, then at most one
  // face of my grid could be an outer face of the l+1 grid.
  if (nchild >0) {
    bool recv[nchild*2*G_ndim];  // whether to get data in this dir
    size_t nel[nchild*2*G_ndim]; // number of interfaces in each dir

    for (int ic=0;ic<nchild;ic++) {
      int off = ic*2*G_ndim;
      CI.get_ipos_vec(MCMD->child_procs[ic].Xmin,ixmin);
      CI.get_ipos_vec(MCMD->child_procs[ic].Xmax,ixmax);

      // define interface region of fine and coarse grids, and if 
      // each direction is to be included or not.
      for (int ax=0;ax<G_ndim;ax++) {
        if (lxmin[ax] == ixmin[ax] &&
            lxmin[ax] != dxmin[ax]) recv[off +2*ax] = true;
        else                        recv[off +2*ax] = false;
        
        if (lxmax[ax] == ixmax[ax] &&
            lxmax[ax] != dxmax[ax]) recv[off +2*ax+1] = true;
        else                        recv[off +2*ax+1] = false;

        ncell[ax] = (ixmax[ax]-ixmin[ax])/G_idx;
        if ( (ixmax[ax]-ixmin[ax]) % G_idx !=0) {
          rep.error("interface region not divisible!",
                    ixmax[ax]-ixmin[ax]);
        }
      } // all dimensions
      for (int d=0;d<2*G_ndim;d++) nel[off +d]=0;
    } // all child grids.

    for (int ic=0;ic<nchild;ic++) {
      int off = ic*2*G_ndim;
      // different number of interfaces depending on dimensionality.
      switch (G_ndim) {
      case 1:
        if (recv[off +XN]) nel[off +XN] = 1;
        if (recv[off +XP]) nel[off +XP] = 1;
        break;
      case 2:
        if (recv[off +XN]) nel[off +XN] = ncell[YY];
        if (recv[off +XP]) nel[off +XP] = ncell[YY];
        if (recv[off +YN]) nel[off +YN] = ncell[XX];
        if (recv[off +YP]) nel[off +YP] = ncell[XX];
        break;
      case 3:
        if (recv[off +XN]) nel[off +XN] = ncell[YY]*ncell[ZZ];
        if (recv[off +XP]) nel[off +XP] = ncell[YY]*ncell[ZZ];
        if (recv[off +YN]) nel[off +YN] = ncell[XX]*ncell[ZZ];
        if (recv[off +YP]) nel[off +YP] = ncell[XX]*ncell[ZZ];
        if (recv[off +ZN]) nel[off +ZN] = ncell[XX]*ncell[YY];
        if (recv[off +ZP]) nel[off +ZP] = ncell[XX]*ncell[YY];
        break;
      default:
        rep.error("bad ndim in setup_flux_recv",G_ndim);
        break;
      } // dims
    } // all child grids

    // initialize arrays
    flux_update_recv.resize(nchild*2*G_ndim);
    for (int ic=0;ic<nchild;ic++) {
#ifdef TEST_BC89FLUX
      cout <<"UniformGridParallel::setup_flux_recv: ";
      cout <<ic<<" has "<< 2*G_ndim <<"boundaries,\n";
#endif
      int off = ic*2*G_ndim;
      for (int d=0; d<2*G_ndim; d++) {
        int el = off+d;
#ifdef TEST_BC89FLUX
        cout <<"UniformGridParallel::setup_flux_recv: ";
        cout << ic <<" d="<< d <<", recv[el]="<<recv[el]<<", nel=";
        cout << nel[el] <<"\n";
#endif
        if (recv[el] == true) {
          flux_update_recv[el].fi.resize(nel[el]);
          flux_update_recv[el].dir = d;
          flux_update_recv[el].ax  = d/2;
          flux_update_recv[el].Ncells = nc;
          flux_update_recv[el].rank.push_back(
                                          MCMD->child_procs[ic].rank);
          for (size_t i=0; i<nel[el]; i++) {
            flux_update_recv[el].fi[i] = 
                            mem.myalloc(flux_update_recv[el].fi[i],1);
            fi = flux_update_recv[el].fi[i];
            fi->c.resize(nc);
            fi->area.resize(nc);
            fi->flux = mem.myalloc(fi->flux,G_nvar);
            for (int v=0;v<G_nvar;v++) fi->flux[v]=0.0;
          }
        }
        else {
          flux_update_recv[el].fi.resize(1);
          flux_update_recv[el].fi[0] = 0;
          flux_update_recv[el].dir = d;
          flux_update_recv[el].ax  = d/2;
          flux_update_recv[el].Ncells = nc;
        }
      } // loop over dims
      CI.get_ipos_vec(MCMD->child_procs[ic].Xmin,ixmin);
      CI.get_ipos_vec(MCMD->child_procs[ic].Xmax,ixmax);


      // For each interface, find the cell that is off the fine grid
      // and that includes the interface.
      for (int d=0; d<2*G_ndim; d++) {
        if (recv[off+d]) {
#ifdef TEST_BC89FLUX
          cout <<"parallel setup_flux_recv(): adding cells.\n";
#endif
          add_cells_to_face(static_cast<enum direction>(d),
                    ixmin,ixmax,ncell,1,flux_update_recv[off+d]);
        }
      } // loop over dims.
    }  // loop over child grids.
  } // if there are child grids.

  else  {
    // There are no children, but my grid might have a boundary in 
    // common with the l+1 level's outer boundary.
    // First try to exclude this:
#ifdef TEST_BC89FLUX
    cout <<"pllel setup_flux_recv(): checking for external faces\n";
#endif
    int edge=-1, axis=-1;
    for (int ax=0;ax<G_ndim;ax++) {
      if (G_ixmin[ax] == lxmax[ax]) {
        edge=2*ax;
        axis=ax;
      }
      if (G_ixmax[ax] == lxmin[ax]) {
        edge=2*ax+1;
        axis=ax;
      }
    }
    if (edge<0) {
#ifdef TEST_BC89FLUX
      cout <<"no edges adjacent to l+1 level\n";
#endif
      flux_update_recv.resize(1);
      flux_update_recv[0].fi.resize(1);
      flux_update_recv[0].fi[0] = 0;
      return 0;
    }

    size_t nel = 0;
    double pos[G_ndim];
    int ch = -1;
    // If we get to here, then one edge borders the l+1 domain, so we
    // see which edge it is, see how many grids on the l+1 domain are
    // facing onto this grid, and add these faces to the list.
    if (G_ndim==1) {
      // this is easy, max one child grid, with one face.
      flux_update_recv.resize(1);
      nel = 1;
      rep.error("Write 1D flux recv setup code",0);
      if (edge==XN) {
        pos[0] = G_xmin[0]-G_dx;
      }
      else {
        pos[0] = G_xmax[0]+G_dx;
      }
      ch = MCMD->get_grid_rank(par,pos,l+1);
      flux_update_recv[0].rank.push_back(ch);
      // dir is the outward normal of the grid edge at level l+1
      flux_update_recv[0].dir = 
                        OppDir(static_cast<enum direction>(edge));
      flux_update_recv[0].ax  = axis;
      flux_update_recv[0].fi.resize(1);
      flux_update_recv[0].fi[0] = 
                      mem.myalloc(flux_update_recv[0].fi[0],1);
      fi = flux_update_recv[0].fi[0];
      fi->c.resize(nc);
      fi->area.resize(nc);
      fi->flux = mem.myalloc(fi->flux,G_nvar);
      for (int v=0;v<G_nvar;v++) fi->flux[v]=0.0;
      for (int v=0;v<G_ndim;v++) ncell[v] = G_ng[v];
      // change G_ixmin,G_ixmax to min/max of neighbouring child.
      add_cells_to_face(OppDir(static_cast<enum direction>(edge)),
                      G_ixmin,G_ixmax,ncell,1,flux_update_recv[0]);
    }
    else if (G_ndim==2) {
      // up to 2 grids on l+1, so go through them one by one, see
      // if they exist, and then add the cells.
      flux_update_recv.resize(2);
      int perp = (axis+1+G_ndim) % G_ndim;
      nel = G_ng[perp]/2;  // child covers half of the grid.
      pos[perp] = G_xmin[perp]+0.25*G_range[perp];
      double xmin[2], xmax[2];
      for (int ic=0;ic<2;ic++) {
        // find xmin/xmax of boundary region, and rank of child
        // procs.
        xmin[perp] = G_xmin[perp] + ic*0.5*G_range[perp];
        xmax[perp] =   xmin[perp]   + 0.5*G_range[perp];
        pos[perp] += ic*0.5*G_range[perp];
        if (edge%2==0) {
          // negative direction
          pos[axis] = G_xmin[axis] - G_dx; // just off grid.
          xmin[axis] = G_xmin[axis] - 0.5*G_range[axis];
          xmax[axis] = G_xmin[axis];
        }
        else {
          // positive direction
          pos[axis] = G_xmax[axis] + G_dx; // just off grid.
          xmin[axis] = G_xmax[axis];
          xmax[axis] = G_xmax[axis] + 0.5*G_range[axis];
        }
#ifdef TEST_BC89FLUX
        cout <<"ic="<<ic<<", perp="<<perp<<", ax="<<axis<<"\n";
        rep.printVec("xmin",xmin,G_ndim);
        rep.printVec("xmax",xmax,G_ndim);
#endif
        
        ch = MCMD->get_grid_rank(par,pos,l+1);
        if (ch>=0) {
          flux_update_recv[ic].rank.push_back(ch);
          // dir is the outward normal of the grid edge at level l+1
          flux_update_recv[ic].dir = 
                          OppDir(static_cast<enum direction>(edge));
          flux_update_recv[ic].ax  = axis;
          flux_update_recv[ic].fi.resize(nel);
          for (size_t i=0; i<nel; i++) {
            flux_update_recv[ic].fi[i] = 
                          mem.myalloc(flux_update_recv[ic].fi[i],1);
            fi = flux_update_recv[ic].fi[i];
            fi->c.resize(nc);
            fi->area.resize(nc);
            fi->flux = mem.myalloc(fi->flux,G_nvar);
            for (int v=0;v<G_nvar;v++) fi->flux[v]=0.0;
          }
          for (int v=0;v<G_ndim;v++) ncell[v] = G_ng[v]/2;
          CI.get_ipos_vec(xmin,ixmin);
          CI.get_ipos_vec(xmax,ixmax);
          
#ifdef TEST_BC89FLUX
          cout <<"FLUX: adding "<<nel<<" cells to recv boundary."<<endl;
#endif
          add_cells_to_face(OppDir(static_cast<enum direction>(edge)),
                    ixmin,ixmax,ncell,1,flux_update_recv[ic]);
        }
        else {
          // no child here, so just create one null element
          flux_update_recv[ic].fi.resize(1);
          flux_update_recv[ic].fi[0] = 0;
        }
      } // loop over child grids (up to 2)
    } // 2D

    else {
      // up to 4 grids on l+1, so go through them one by one, see if
      // they exist, and add the cells.
      flux_update_recv.resize(4);
      int perp[2];
      perp[0] = (axis+1+G_ndim) % G_ndim;
      perp[1] = (axis+2+G_ndim) % G_ndim;
      nel = G_ng[perp[0]]*G_ng[perp[1]]/4;
      //pos[perp[0]] = G_xmin[perp[0]]+0.25*G_range[perp[0]];
      //pos[perp[1]] = G_xmin[perp[1]]+0.25*G_range[perp[1]];
      
      double xmin[3], xmax[3];
      for (int ic=0;ic<4;ic++) {
        // find xmin/xmax of boundary region, for each of the four
        // potential child grids.
        // First the in-plane directions:
        switch (ic) {
          case 0: // (-,-)
          for (int d=0;d<2;d++) {
            xmin[perp[d]] = G_xmin[perp[d]];
            xmax[perp[d]] =   xmin[perp[d]] + 0.5 *G_range[perp[d]];
            pos[perp[d]]  = G_xmin[perp[d]] + 0.25*G_range[perp[d]];
          }
          break;

          case 1: // (+,-)
          xmin[perp[0]] = G_xmin[perp[0]] + 0.5*G_range[perp[0]];
          xmin[perp[1]] = G_xmin[perp[1]];
          for (int d=0;d<2;d++) {
            xmax[perp[d]] = xmin[perp[d]] + 0.5*G_range[perp[d]];
          }
          pos[perp[0]] = G_xmin[perp[0]]+0.75*G_range[perp[0]];
          pos[perp[1]] = G_xmin[perp[1]]+0.25*G_range[perp[1]];
          break;

          case 2: // (-,+)
          xmin[perp[0]] = G_xmin[perp[0]];
          xmin[perp[1]] = G_xmin[perp[1]] + 0.5*G_range[perp[1]];
          for (int d=0;d<2;d++) {
            xmax[perp[d]] = xmin[perp[d]] + 0.5*G_range[perp[d]];
          }
          pos[perp[0]] = G_xmin[perp[0]]+0.25*G_range[perp[0]];
          pos[perp[1]] = G_xmin[perp[1]]+0.75*G_range[perp[1]];
          break;

          case 3: // (+,+)
          for (int d=0;d<2;d++) {
            xmin[perp[d]] = G_xmin[perp[d]] + 0.5 *G_range[perp[d]];
            xmax[perp[d]] =   xmin[perp[d]] + 0.5 *G_range[perp[d]];
            pos[perp[d]]  = G_xmin[perp[d]] + 0.75*G_range[perp[d]];
          }
          break;

          default:
          rep.error("loopy error",ic);
        }
        // Then the normal direction
        if (edge%2==0) {
          pos[axis] = G_xmin[axis] - G_dx; // just off grid.
          xmin[axis] = G_xmin[axis] - 0.5*G_range[axis];
          xmax[axis] = G_xmin[axis];
        }
        else {
          // positive direction
          pos[axis] = G_xmax[axis] + G_dx; // just off grid.
          xmin[axis] = G_xmax[axis];
          xmax[axis] = G_xmax[axis] + 0.5*G_range[axis];
        }
#ifdef TEST_BC89FLUX
        cout <<"ic="<<ic<<", perp=["<<perp[0]<<","<<perp[1]<<"] ";
        cout <<"ax="<<axis<<"\n";
        rep.printVec("xmin",xmin,G_ndim);
        rep.printVec("xmax",xmax,G_ndim);
        rep.printVec(" pos", pos,G_ndim);
#endif

        ch = MCMD->get_grid_rank(par,pos,l+1);
        if (ch>=0) {
#ifdef TEST_BC89FLUX
          cout <<"ic="<<ic<<": from rank "<<ch<<"...\n";
#endif
          flux_update_recv[ic].rank.push_back(ch);
          // dir is the outward normal of the grid edge at level l+1
          flux_update_recv[ic].dir = 
                          OppDir(static_cast<enum direction>(edge));
          flux_update_recv[ic].ax  = axis;
          flux_update_recv[ic].fi.resize(nel);
          for (size_t i=0; i<nel; i++) {
            flux_update_recv[ic].fi[i] = 
                          mem.myalloc(flux_update_recv[ic].fi[i],1);
            fi = flux_update_recv[ic].fi[i];
            fi->c.resize(nc);
            fi->area.resize(nc);
            fi->flux = mem.myalloc(fi->flux,G_nvar);
            for (int v=0;v<G_nvar;v++) fi->flux[v]=0.0;
          }
          for (int v=0;v<G_ndim;v++) ncell[v] = G_ng[v]/2;
          CI.get_ipos_vec(xmin,ixmin);
          CI.get_ipos_vec(xmax,ixmax);
          
#ifdef TEST_BC89FLUX
          cout <<"ic="<<ic<<"FLUX: adding "<<nel<<" cells to recv boundary."<<endl;
#endif
          add_cells_to_face(OppDir(static_cast<enum direction>(edge)),
                    ixmin,ixmax,ncell,1,flux_update_recv[ic]);
        }
        else {
          // no child here, so just create one null element
          flux_update_recv[ic].fi.resize(1);
          flux_update_recv[ic].fi[0] = 0;
        }
      } // loop over child grids (up to 4)
    } // 3D
  } // child grids?

  return 0;
}




// ##################################################################
// ##################################################################



int UniformGridParallel::setup_flux_send(
      class SimParams &par, ///< simulation params (including BCs)
      const int lm1         ///< level to send to
      )
{
#ifdef TEST_BC89FLUX
  cout <<" UniformGridParallel::setup_flux_send() send to level=";
  cout <<lm1<<" from MY LEVEL l="<<lm1+1<<"\n";
#endif

  int err = UniformGrid::setup_flux_send(par,lm1);
  rep.errorTest("UniformGrid::setup_flux_send",0,err);

  // Add ranks for each send, based on parent rank.
  int l = lm1+1; // my level
  int ns = flux_update_send.size();
  if (ns != 2*G_ndim) rep.error("bad flux send size",ns);

  if (par.levels[l].MCMD.get_nproc()==1)  {
    for (int d=0;d<ns;d++) {
      flux_update_send[d].rank.push_back(0);
      flux_update_send[d].dir = d;
    }
  }
  else {
    // always send to parent, sometimes send to neighbour of parent
    // if boundary between parent and neighbour sits on the edge of
    // my level.
    class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
    int pproc = MCMD->parent_proc;
    for (int d=0;d<ns;d++) {
      // if boundary element is not empty, send data to parent.
      if (flux_update_send[d].fi[0] !=0) {
        flux_update_send[d].rank.push_back(pproc);
      }
    }
  
    for (int ax=0;ax<G_ndim;ax++) {
      // check if parent boundary is also my boundary, in which
      // case we need to send data to parent's neighbour
      double pos[G_ndim];

      // First negative direction along this axis
      int d = 2*ax, p1=-1, p2=-1;
      if (flux_update_send[d].fi[0] !=0) {
        for (int i=0;i<G_ndim;i++) pos[i] = G_xmin[i] + G_dx;
        //pos[ax] += G_dx;
        p1 = MCMD->get_grid_rank(par,pos,lm1);
        if (p1!=pproc) rep.error("BC89 finding parents 1",p1-pproc);
        pos[ax] -= 2*G_dx;
        p2 = MCMD->get_grid_rank(par,pos,lm1);
#ifdef TEST_BC89FLUX
        cout <<"ax="<<ax<<", d="<<d<<", parent="<<pproc<<", p1="<<p1;
        cout<<", and ngb="<<p2<<"\n";
#endif
        if (p2!=pproc) {
#ifdef TEST_BC89FLUX
          cout <<"Adding 2nd parent to dir="<<d<<", "<<p2<<"\n";
#endif
          flux_update_send[d].rank.push_back(p2);
        }
      }
      flux_update_send[d].dir = d;
      flux_update_send[d].ax  = ax;

      // now positive direction
      d = 2*ax+1;
      p1=-1;
      p2=-1;
      if (flux_update_send[d].fi[0] !=0) {
        for (int i=0;i<G_ndim;i++) pos[i] = G_xmax[i] - G_dx;
        //pos[ax] -= G_dx;
        rep.printVec("p1 pos",pos,G_ndim);
        p1 = MCMD->get_grid_rank(par,pos,lm1);
        if (p1!=pproc) rep.error("BC89 finding parents 2",p1-pproc);
        pos[ax] += 2*G_dx;
        rep.printVec("p2 pos",pos,G_ndim);
        p2 = MCMD->get_grid_rank(par,pos,lm1);
#ifdef TEST_BC89FLUX
        cout <<"ax="<<ax<<", d="<<d<<", parent="<<pproc<<", p1="<<p1;
        cout<<", and p2="<<p2<<endl;
#endif
        if (p2!=pproc) {
#ifdef TEST_BC89FLUX
          cout <<"Adding 2nd parent to dir="<<d<<", "<<p2<<endl;
#endif
          flux_update_send[d].rank.push_back(p2);
        }
      }
      flux_update_send[d].dir = d;
      flux_update_send[d].ax  = ax;
    } // loop over axes
  } // if nproc>1

  return 0;
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
      double *lev_xn, // level xmin
      double *lev_xp, // level xmax
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  uniform_grid_cyl(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp)
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
      double *lev_xn, // level xmin
      double *lev_xp, // level xmax
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp  ///< array of max. x/y/z for full simulation.
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd),
  uniform_grid_sph(nd,nv,eqt,nbc,xn,xp,nc,lev_xn,lev_xp,sim_xn,sim_xp)
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
