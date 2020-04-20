///
/// \file uniform_grid.cpp
/// 
/// \author Jonathan Mackey
/// 
/// Function definitions for UniformGrid class:
/// - setup boundaries.
/// - allocate/delete grid data and setup linked lists.
/// - update boundaries.
/// 
/// Modifications:\n
///  - 2007-08-07 minor fixes to boundary conditions.
///  - 2007-11-06 changed BC setup structure.  tested and it works
///     identically. This is for ease of use when adding parallel
///     boundaries.
/// - 2010-01-26 JM: Added get_iXmin/iXmax/iRange() functions to
///     grid class to give integer positions.
/// - 2010-01-26 JM: Took out unused cstep,maxstep vars in some 
///     update_boundary functions.
///  - 2010-02-03 JM: fixed compiler warnings about re-used and 
///     unused variables.
/// - 2010-03-12 JM: Started work on one-way boundaries which allow 
///     inflow/outflow but not both. For outflow-only, this will work 
///     by checking the normal velocity.  If it is off the grid then 
///     we'll use zero gradient to set the boundary value; if it's 
///     onto the grid then we'll use reflecting condition to set the 
///     boundary value (but what about field??? Maybe just setting 
///     the velocity to be either outflowing or zero in the boundary 
///     cells is a better method...).
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; 
///     added oneway-outflow boundary.
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010-07-24 JM: Added stellar wind class, and boundary setup/ 
///     update functions.
/// - 2010-07-28 JM: Slightly changed boundary update logic -- we
///    check for internal and external BCs in both functions and just
///    do nothing for the ones we're not interested in.  Then if there
///    is really an unhandled boundary we report a warning.
/// - 2010-09-27 JM: fixed one comment.
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.10.04 JM: Added last-point-set flag.
/// - 2010.10.05 JM: Added spherical coordinates to the 
///    BC_assign_STWIND() function, since external data needs to be 
///    considered.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///    GEOMETRIC_GRID ifdef.  Will probably keep it since it is the
///    way things will go eventually.  The new grid classes have extra
///    function for the distance between two points, two cells, and
///    between a vertex and a cell.
/// - 2011.01.06 JM: New stellar wind interface.
/// - 2011.01.07 JM: I debugged the geometric grid functions, and now
///    it works very well!  I have a nice spherical expansion for
///    stellar winds.
/// - 2011.02.15 JM: Added support for time-varying stellar winds.  I
///   think it is working now, but it needs testing.
/// - 2011.02.16 JM: Fixed a bug in spherical grid iR_cov() function, 
///    for where the grid does not begin at the origin.
/// - 2011.02.25 JM: Set column densities in boundary data to zero in 
///    boundary setup function.  Removed HCORR ifdef.
/// - 2011.03.21 JM: Got rid of zero-ing of column-densities -- done
///    in the cell constructor.
/// - 2011.04.06 JM: Added idifference_cell2cell() function for 
///    thermal conduction calculation.
/// - 2011.11.22 JM: Added t_scalefactor parameter for stellar winds.
/// - 2012.01.20 JM: Check that we have a stellar-wind source before 
///    trying to set up the wind boundary!  (avoids seg.faults).
///
/// - 2013.01.11 JM: Fiddled with radiative shock boundary to try to
///    make it better. (also 2013.01.23, and 2013.01.20).
/// - 2013.01.14 JM: Added GRIDV2 ifdef to eventually retire this...
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.04.15 JM: Tidied up stdIO, fixed bug in iR_cov() (cyl).
/// - 2013.06.13 JM: Added STARBENCH1 boundary to put in a wall that
///    shields the ISM from ionising radiation, to study the thermal
///    physics and recombination in the shadowed region.
/// - 2013.06.17 JM: Changed STARBENCH1 from internal to external BC
///    and fixed some bugs with it.  It works well now.
/// - 2013.08.20 JM: Changed cell_interface for radiative transfer
///    variables.
/// - 2013.09.06 JM: Added difference_vertex2cell() functions.
///    Set stellar wind boundary conditions to use physical distances
///    not integer distances, to avoid rounding errors.
/// - 2015.01.10 JM: New include statements for new file structure.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2016.02.11 JM: Worked on Grid_v2 update (full boundaries).
/// - 2016.02.12 JM: including source code from static_grid.cc, and
///    renamed to uniform_grid.cpp, worked on new grid structure.
/// - 2016.02.19 JM: new grid structure finished, compiles and runs
///    the DMR test.
/// - 2016.02.22 JM: bugfixes for periodic boundaries.
/// - 2016.03.08 JM: bugfixes for outflow boundaries, grid setup.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
/// - 2017.11.07-22 JM: updating boundary setup.
/// - 2017.12.09 JM: got rid of SimPM references.
/// - 2018.05.09 JM: set cell positions using physical pos, not
///    integer pos (b/c NG grids don't all have cell-size=2
///    units) and removed STARBENCH1 boundary.
/// - 2018.05.10 JM: Moved boundary assignment to setup_fixed_grid,
///    and boundary updates to update_boundaries.cpp

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "grid/uniform_grid.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include <fstream>
#include <iostream>
using namespace std;



//#define TEST_BC89FLUX

// ##################################################################
// ##################################################################



UniformGrid::UniformGrid(
    int nd,
    int nv,
    int eqt,
    int Nbc,       ///< Number of boundary cells to use.
    double *g_xn,  // this grid xmin
    double *g_xp,  // this grid xmax
    int *g_nc,
    double *lev_xn, // level xmin
    double *lev_xp, // level xmax
    double *sim_xn, // sim xmin
    double *sim_xp  // sim xmax
    )
  : 
  VectorOps_Cart(nd),
  G_ndim(nd),
  G_nvar(nv),
  G_eqntype(eqt),
  BC_nbc(Nbc)
{

#ifdef TESTING
  cout <<"Setting up UniformGrid with G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
  rep.printVec("\tXmin",g_xn,nd);
  rep.printVec("\tXmax",g_xp,nd);
  rep.printVec("\tNpt ",g_nc,nd);
#endif
  G_coordsys = COORD_CRT;  // Cartesian Coordinate system

  //
  // Allocate arrays for dimensions of grid.
  //
  G_ng=0;
  G_ng_all=0;
  G_xmin=0;
  G_xmax=0;
  G_range=0;
  
  G_ng     = mem.myalloc(G_ng,    MAX_DIM);
  G_ng_all = mem.myalloc(G_ng_all,MAX_DIM);
  G_xmin   = mem.myalloc(G_xmin,  G_ndim);
  G_xmax   = mem.myalloc(G_xmax,  G_ndim);
  G_range  = mem.myalloc(G_range, G_ndim);
  G_xmin_all  = mem.myalloc(G_xmin_all, G_ndim);
  G_xmax_all  = mem.myalloc(G_xmax_all, G_ndim);
  G_range_all = mem.myalloc(G_range_all,G_ndim);
  G_ixmin  = mem.myalloc(G_ixmin, G_ndim);
  G_ixmax  = mem.myalloc(G_ixmax, G_ndim);
  G_irange = mem.myalloc(G_irange,G_ndim);
  G_ixmin_all  = mem.myalloc(G_ixmin_all, G_ndim);
  G_ixmax_all  = mem.myalloc(G_ixmax_all, G_ndim);
  G_irange_all = mem.myalloc(G_irange_all,G_ndim);

  //
  // Useful for NG/parallel grids, where simulation domain
  // may be larger than this grid.
  //
  Sim_xmin   = mem.myalloc(Sim_xmin,  G_ndim);
  Sim_xmax   = mem.myalloc(Sim_xmax,  G_ndim);
  Sim_range  = mem.myalloc(Sim_range, G_ndim);
  Sim_ixmin  = mem.myalloc(Sim_ixmin, G_ndim);
  Sim_ixmax  = mem.myalloc(Sim_ixmax, G_ndim);
  Sim_irange = mem.myalloc(Sim_irange,G_ndim);

  L_xmin   = mem.myalloc(L_xmin,  G_ndim);
  L_xmax   = mem.myalloc(L_xmax,  G_ndim);
  L_range  = mem.myalloc(L_range, G_ndim);
  L_ixmin  = mem.myalloc(L_ixmin, G_ndim);
  L_ixmax  = mem.myalloc(L_ixmax, G_ndim);
  L_irange = mem.myalloc(L_irange,G_ndim);

  //
  // Assign values to grid dimensions and set the cell-size
  //
  G_ncell = 1;
  G_ncell_all = 1;
  // initialise to 1 for routines that loop over unused dimensions.
  for (int i=0;i<MAX_DIM;i++) {
    G_ng_all[i] = 1;
    G_ng[i] = 1;
  }

  for (int i=0;i<G_ndim;i++) {
    G_ng[i] = g_nc[i];
    G_ng_all[i] = g_nc[i] + 2*BC_nbc;
    G_ncell *= G_ng[i]; // Get total number of cells.
    G_ncell_all *= G_ng[i] +2*BC_nbc; // all cells including boundary.
    G_xmin[i] = g_xn[i];
    G_xmax[i] = g_xp[i];
    G_range[i] = G_xmax[i]-G_xmin[i];
    if(G_range[i]<0.) rep.error("Negative range in direction i",i);

    //
    // For single static grid, Sim_xmin/xmax/range are the same as
    // for the grid (becuase there is only one grid).  But we set
    // them here anyway so that it is correctly placed in a multigrid
    // setup.
    //
    Sim_xmin[i] = sim_xn[i];
    Sim_xmax[i] = sim_xp[i];
    Sim_range[i] = sim_xp[i]-sim_xn[i];

    L_xmin[i] = lev_xn[i];
    L_xmax[i] = lev_xp[i];
    L_range[i] = lev_xp[i]-lev_xn[i];

  }

#ifdef TESTING
  cout <<"MIN.MAX for x = "<<G_xmin[XX]<<"\t"<<G_xmax[XX]<<"\n";
  cout <<"Setting cell size...\n";
#endif
  //
  // Checks grid dimensions and discretisation is reasonable,
  // and that it gives cells with equal length in each dimension.
  //
  set_cell_size();

  // grid dimensions including boundary data
  for (int v=0;v<G_ndim;v++)
    G_xmin_all[v] = G_xmin[v] - BC_nbc*G_dx;
  for (int v=0;v<G_ndim;v++)
    G_xmax_all[v] = G_xmax[v] + BC_nbc*G_dx;
  for (int v=0;v<G_ndim;v++)
    G_range_all[v] = G_xmax_all[v] - G_xmin_all[v];

  //
  // Now create the first cell, and then allocate data from there.
  // Safe to assume we have at least one cell...
  //
#ifdef TESTING
  cout <<" done.\n Initialising first cell...\n";
#endif // TESTING
  G_fpt_all = CI.new_cell();
  G_fpt_all->id = 0;
#ifdef TESTING
  cout <<" done.\n";
#endif // TESTING
  if(G_fpt_all==0) {
    rep.error("Couldn't assign memory to first cell in grid.",
              G_fpt_all);
  }
  
  //
  // assign memory for all cells, including boundary cells.
  //
#ifdef TESTING
  cout <<"\t allocating memory for grid.\n";
#endif // TESTING
  int err = allocate_grid_data();
  if(err!=0)
    rep.error("Error setting up grid, allocate_grid_data",err);

  //
  // assign grid structure on cells, setting positions and ngb
  // pointers.
  //
#ifdef TESTING
  cout <<"\t assigning pointers to neighbours.\n";
#endif // TESTING
  err += assign_grid_structure();
  if(err!=0)
    rep.error("Error setting up grid, assign_grid_structure",err);
  
#ifdef TESTING
  rep.printVec("\tFirst Pt. integer position",FirstPt()->pos,nd);
  rep.printVec("\tLast  Pt. integer position", LastPt()->pos,nd);
#endif // TESTING

  //
  // Leave boundaries uninitialised.
  //
  BC_nbd=-1;

  //
  // Set integer dimensions/location of grid.
  //
  CI.get_ipos_vec(G_xmin, G_ixmin );
  CI.get_ipos_vec(G_xmax, G_ixmax );
  for (int v=0;v<G_ndim;v++)
    G_irange[v] = G_ixmax[v]-G_ixmin[v];
  G_idx = G_irange[XX]/G_ng[XX];
  for (int v=0;v<G_ndim;v++)
    G_ixmin_all[v] = G_ixmin[v] - BC_nbc*G_idx;
  for (int v=0;v<G_ndim;v++)
    G_ixmax_all[v] = G_ixmax[v] + BC_nbc*G_idx;
  for (int v=0;v<G_ndim;v++)
    G_irange_all[v] = G_ixmax_all[v] - G_ixmin_all[v];

  CI.get_ipos_vec(Sim_xmin, Sim_ixmin);
  CI.get_ipos_vec(Sim_xmax, Sim_ixmax);
  for (int v=0;v<G_ndim;v++) {
    Sim_irange[v] = Sim_ixmax[v]-Sim_ixmin[v];
  }
  CI.get_ipos_vec(L_xmin, L_ixmin);
  CI.get_ipos_vec(L_xmax, L_ixmax);
  for (int v=0;v<G_ndim;v++) {
    L_irange[v] = L_ixmax[v]-L_ixmin[v];
  }

  // set external boundary cells to be not part of the simulation
  // domain if they are outside of Sim_Xmin and Sim_Xmax.
  class cell *c = FirstPt_All();
  bool dom=true;
  do {
    c->isdomain=true;
    dom=true;
    for (int v=0;v<G_ndim;v++) {
      if (c->pos[v]<Sim_ixmin[v] || c->pos[v]>Sim_ixmax[v]) dom=false;
    }
    if (!dom) c->isdomain = false;
    c->isleaf=true; // assume all cells are leaves unless changed.
    if (dom) c->timestep = true;
  } while ( (c=NextPt_All(c))!=0);


#ifdef TESTING
  rep.printVec("grid ixmin ", G_ixmin, G_ndim);
  rep.printVec("grid ixmax ", G_ixmax, G_ndim);
  rep.printVec("grid irange", G_irange,G_ndim);
  rep.printVec("grid ixmin_all ", G_ixmin_all, G_ndim);
  rep.printVec("grid ixmax_all ", G_ixmax_all, G_ndim);
  rep.printVec("grid irange_all", G_irange_all,G_ndim);
  rep.printVec("grid xmin ", G_xmin, G_ndim);
  rep.printVec("grid xmax ", G_xmax, G_ndim);
  rep.printVec("grid range", G_range,G_ndim);
  rep.printVec("grid xmin_all ", G_xmin_all, G_ndim);
  rep.printVec("grid xmax_all ", G_xmax_all, G_ndim);
  rep.printVec("grid range_all", G_range_all,G_ndim);
  rep.printVec("Sim xmin ", Sim_xmin, G_ndim);
  rep.printVec("Sim xmax ", Sim_xmax, G_ndim);
  rep.printVec("Sim range", Sim_range,G_ndim);
#endif


#ifdef TESTING
  cout <<"Cartesian grid: dr="<<G_dx<<"\n";
#endif // TESTING
  RT=0;

#ifdef TESTING
  cout <<"UniformGrid Constructor done.\n";
#endif
} //UniformGrid Constructor



// ##################################################################
// ##################################################################


UniformGrid::~UniformGrid()
{
  //
  // Delete the grid data.
  //
  cell *cpt=FirstPt_All();
  cell *npt=NextPt_All(cpt);
  do {
    //cout <<"deleting cell id: "<<cpt->id<<"\n";
    CI.delete_cell(cpt);
    cpt = npt;
  } while ( (npt=NextPt_All(cpt))!=0 );
  //cout <<"deleting cell id: "<<cpt->id<<"\n";
  CI.delete_cell(cpt);

  G_ng    = mem.myfree(G_ng);
  G_xmin  = mem.myfree(G_xmin);
  G_xmax  = mem.myfree(G_xmax);
  G_range = mem.myfree(G_range);
  G_ixmin  = mem.myfree(G_ixmin);
  G_ixmax  = mem.myfree(G_ixmax);
  G_irange = mem.myfree(G_irange);

  Sim_xmin  = mem.myfree(Sim_xmin);
  Sim_xmax  = mem.myfree(Sim_xmax);
  Sim_range = mem.myfree(Sim_range);
  Sim_ixmin  = mem.myfree(Sim_ixmin);
  Sim_ixmax  = mem.myfree(Sim_ixmax);
  Sim_irange = mem.myfree(Sim_irange);

  BC_deleteBoundaryData();
  
#ifdef TESTING
  cout <<"UniformGrid Destructor:\tdone.\n";
#endif
} // Destructor



// ##################################################################
// ##################################################################



int UniformGrid::assign_grid_structure()
{
#ifdef TESTING
  cout<<"\tAssignGridStructure.\n";
#endif // TESTING

  /// \section Structure
  /// There is a base grid, Nx,Ny,Nz elements, which is 
  /// allocated at the start of the simulation.  These are
  /// ordered so that element (ix,iy,iz) has 
  /// \f$ \mbox{ID} = \mbox{firstPtID} + i_x +N_x*i_y + N_x*N_y*i_z \f$.
  /// NextPt(cpt) will return the next
  /// cell in this set, and null if it is at the last one. 
  /// NextPt(cpt,direction) will return the next cell in a given
  /// direction, and null if it has reached an edge cell. Each cell 
  /// has pointers to its neighbours.
  /// 

  // ----------------------------------------------------------------
  //
  // The cell integer positions are based off Sim_xmin (not G_xmin)
  // for reasons that become clear in the parallel version -- I want
  // positions to be globally applicable.  So we set the offset based
  // on (Sim_xmin-G_xmin).
  // First cell is at [1,1,1], so the offset has 1 added to it
  // (because the cell centre is 1 unit from xmin).
  //
  //int ipos[MAX_DIM], offset[MAX_DIM];
  //for (int i=0; i<MAX_DIM; i++) ipos[i] = offset[i] = 0;
  //for (int i=0; i<G_ndim; i++) {
  //  offset[i] = 1+ 2*static_cast<int>(ONE_PLUS_EPS*(G_xmin[i]-Sim_xmin[i])/G_dx);
#ifdef TESTING
  //  cout <<"************OFFSET["<<i<<"] = "<<offset[i]<<"\n";
#endif
  //}

  // ----------------------------------------------------------------
  // ---------------------- SET CELL POSITIONS ----------------------
  // Go through full grid, setting positions, where the first cell is
  // BC_nbc cells negative of G_xmin in all directions.
  //
  int    ix[MAX_DIM];  // ix[] is a counter for cell number.
  for (int i=0;i<MAX_DIM;i++) ix[i]=-BC_nbc;

  double dpos[MAX_DIM];
  for (int i=0;i<MAX_DIM;i++) dpos[i] = 0.0;

  class cell *c = FirstPt_All();
  do {
#ifdef TESTING
    //cout <<"Cell positions: id = "<<c->id<<"\n";
#endif // TESTING

    //
    // Assign positions, for integer positions the first on-grid cell
    // is at [1,1,1], and a cell is 2 units across, so the second
    // cell is at [3,1,1] etc.
    // Boundary cells can have negative positions.
    //
    // Here we use the physical position as input, and the
    // cell-interface class assigns an integer position based on
    // this physical position.
    //
    for (int i=0; i<G_ndim; i++) dpos[i] = G_xmin[i] + G_dx*(ix[i]+0.5);
    CI.set_pos(c,dpos);
#ifdef TESTING
    //rep.printVec("    pos", c->pos, G_ndim);
#endif // TESTING

    //
    // Initialise the cell data to zero.
    //
    for(int v=0;v<G_nvar;v++) c->P[v] = 0.0;
    if (!CI.query_minimal_cells()) {
      for(int v=0;v<G_nvar;v++) c->Ph[v] = c->dU[v] =0.;
    }

    //
    // See if cell is on grid or not, and set isgd/isbd accordingly.
    //
    CI.get_dpos(c,dpos);
    bool on_grid=true;
    for (int v=0;v<G_ndim;v++)
      if (dpos[v]<G_xmin[v] || dpos[v]>G_xmax[v])
        on_grid=false;
#ifdef TESTING
    //rep.printVec("    dpos",dpos,G_ndim);
    //rep.printVec("    xmax",G_xmax,G_ndim);
#endif // TESTING



    if (on_grid) {
      c->isgd = true;
      c->isbd = false;
#ifdef TESTING
      //cout <<"    cell is on grid";
#endif // TESTING
    }
    else {
      c->isgd = false;
      c->isbd = true;
#ifdef TESTING
      //cout <<"    cell NOT on grid";
#endif // TESTING
    }

    ///
    /// \section Edges
    /// Edge cells can be at corners in multi-dimensions, and there
    /// are 27 possible configurations ranging from not-an-edge to a
    /// corner with three boundary cells.  I track this with an 
    /// integer flag 'isedge', defined so that:
    /// - \f$ i \bmod 3=0 \f$ Not an X-edge.
    /// - \f$ i \bmod 3=1 \f$ Neg. X-edge.
    /// - \f$ i \bmod 3=2 \f$ Pos. X-edge.
    /// - \f$ i \bmod 9<3 \f$ Not a Y-edge.
    /// - \f$ 3 \leq i\bmod 9 <6 \f$ Neg. Y-edge.
    /// - \f$ 6 \leq i\bmod 9 \f$ Pos. Y-edge.
    /// - \f$ i <9 \f$ Not a Z-edge.
    /// - \f$ 9 \leq i <18 \f$ Neg. Z-edge.
    /// - \f$ 18\leq i <27 \f$ Pos. Z-edge.
    /// - \f$ i \geq 27 \f$ Out of Range error.
    ///
    c->isedge =0;
    if      (ix[XX]==0)            c->isedge +=1;
    else if (ix[XX]==G_ng[XX]-1)   c->isedge +=2;
    if (G_ndim>1){
      if      (ix[YY]==0)          c->isedge += 1*3;
      else if (ix[YY]==G_ng[YY]-1) c->isedge += 2*3;
    }
    if (G_ndim>2) {
      if      (ix[ZZ]==0)          c->isedge += 1*3*3;
      else if (ix[ZZ]==G_ng[ZZ]-1) c->isedge += 2*3*3;
    }

    //
    // Boundary values for isedge:  if the cell is off-grid, then set
    // isedge to equal the number of cells it is from the grid.
    //
    if      (ix[XX]<0)         {
      c->isedge = ix[XX];
      //cout <<"ix[x]="<<ix[XX]<<", ";
      //rep.printVec("pos",c->pos,G_ndim);
    }
    else if (ix[XX]>=G_ng[XX]) {
      c->isedge = G_ng[XX]-1-ix[XX];
      //cout <<"ix[x]="<<ix[XX]<<", ";
      //rep.printVec("pos",c->pos,G_ndim);
    }
    if (G_ndim>1) {
      if      (ix[YY]<0)         {
        c->isedge = ix[YY];
        //cout <<"ix[y]="<<ix[YY]<<", ";
        //rep.printVec("pos",c->pos,G_ndim);
      }
      else if (ix[YY]>=G_ng[YY]) {
        c->isedge = G_ng[YY]-1-ix[YY];
        //cout <<"ix[y]="<<ix[YY]<<", ";
        //rep.printVec("pos",c->pos,G_ndim);
      }
    }
    if (G_ndim>2) {
      if      (ix[ZZ]<0)         c->isedge = ix[ZZ];
      else if (ix[ZZ]>=G_ng[ZZ]) c->isedge = G_ng[ZZ]-1-ix[ZZ];
    }

    
    //
    // Increment counters
    //
    ix[XX]++;
    //
    // if we got to the end of a row, then loop back and increment iy
    //
    if (ix[XX] == G_ng[XX]+BC_nbc) {
      ix[XX] = -BC_nbc;
      if (G_ndim > 1) {
        ix[YY]++;
        //
        // if we got to the end of a y-column, then loop back and
        // increment z.
        //
        if (ix[YY] == G_ng[YY]+BC_nbc) {
          ix[YY] = -BC_nbc;
          if (G_ndim > 2) {
            ix[ZZ]++;
            if (ix[ZZ]>G_ng[ZZ]+BC_nbc) {
              cerr <<"\tWe should be done by now.\n";
              return(ix[ZZ]);
            }
          } // If 3D.
        } // If at end of YY row.
      } // If 2D.
    } // If at end of XX row.
  } while ( (c=NextPt_All(c))!=0);
  // ---------------------- SET CELL POSITIONS ----------------------
  // ----------------------------------------------------------------
#ifdef TESTING
  //c = FirstPt_All();
  //do {
  //  cout <<"Cell id = "<<c->id<<"\n";
  //  rep.printVec("cell pos",c->pos,G_ndim);
  //} while ( (c=NextPt_All(c))!=0);
#endif // TESTING

  
  // ----------------------------------------------------------------
  // ----------------------  SET npt POINTERS  ----------------------
  // Now run through grid and set npt pointer for on-grid cells.
  // Also set fpt for the first on-grid point, and lpt for the last.
  //
  bool set_fpt=false, set_lpt=false;
  cell *ctemp=0;
  c=FirstPt_All();
  do {
    //
    // If cell is on-grid, set npt to be the next cell in the list
    // that is also on-grid.
    //
    if (c->isgd) {
      //
      // set pointer to first on-grid cell.
      //
      if (!set_fpt) {
        G_fpt = c;
        set_fpt=true;
      }
      //
      // if next cell in list is also on-grid, this is npt.
      //
      if (NextPt_All(c)!=0 && NextPt_All(c)->isgd) {
        c->npt = NextPt_All(c);
      }
      //
      // Otherwise, go through the list until we get to the end (so
      // we are at the last on-grid cell) or we find another on-grid
      // cell, which we set npt to point to.
      //
      else {
        ctemp = c;
        do {
          ctemp = NextPt_All(ctemp);
        } while (ctemp!=0 && !ctemp->isgd);

        if (ctemp) {
          // there is another on-grid point.
          if (!ctemp->isgd) rep.error("Setting npt error",ctemp->id);
          c->npt = ctemp;
        }
        else {
          // must be last point on grid.
          c->npt = 0;
          G_lpt = c;
          set_lpt = true;
        }
      }
    } // if c->isgd
  } while ( (c=NextPt_All(c))!=0);
  if (!set_lpt) rep.error("failed to find last on-grid point",1);
  if (!set_fpt) rep.error("failed to find first on-grid point",1);
  // ----------------------  SET npt POINTERS  ----------------------
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // ----------------------  SET ngb POINTERS  ----------------------
  //
  // Set up pointers to neighbours in x-direction
  //
  c = FirstPt_All();
  cell *c_prev=0, *c_next = c;
  do {
    c_next = NextPt_All(c_next);
    //cout <<" Starting X ngb loop, c="<<c->id;
    //rep.printVec("pos",c->pos,G_ndim);
    c->ngb[XN] = c_prev;

    if (!c_next) {
      //cout <<"*** c_next=0, cell id ="<<c->id<<", pos=";
      //rep.printVec("pos",c->pos,G_ndim);
      c->ngb[XP] = 0;
    }
    else if (c_next->pos[XX] > c->pos[XX]) {
      c->ngb[XP] = c_next;
      c_prev = c;
      //cout <<"*** c_next[X]>c[X], cell id ="<<c->id<<", pos=";
      //rep.printVec("pos",c->pos,G_ndim);
    }
    else {
      //cout <<"*** c_next[X]<=c[X], cell id ="<<c->id<<"\n";
      //rep.printVec("this pos",c->pos,G_ndim);
      //rep.printVec("next pos",c_next->pos,G_ndim);
      c->ngb[XP] = 0;
      c_prev = 0;
    }

    c = c_next;
  }
  while (c != 0);

  //
  // Pointers to neighbours in the Y-Direction, if it exists.
  // This is not a very smart algorithm and is probably quite slow,
  // but it does the job.
  //
  if (G_ndim>1) {
    c = FirstPt_All();
    do {
      c_next=c;
      do {
        c_next = NextPt_All(c_next);
      } while ( (c_next!=0) && 
                (c_next->pos[YY] >= c->pos[YY]) &&
                (c_next->pos[XX] != c->pos[XX]) );
      //
      // So now maybe c_next is zero, in which case ngb[YP]=0
      //
      if      (!c_next) {
        c->ngb[YP] = 0;
      }
      //
      // Or else we have looped onto the next z-plane, so ngb[YP]=0
      //
      else if (c_next->pos[YY] < c->pos[YY]) {
        c->ngb[YP] = 0;
      }
      //
      // Or else cells have the same pos[XX], so they are neighbours.
      //
      else {
        c->ngb[YP] = c_next;
        c_next->ngb[YN] = c;
      }
    } while ((c=NextPt_All(c)) !=0);
  } // if G_ndim>1

  //
  // Pointers to neighbours in the Z-Direction, if it exists.
  //
  if (G_ndim>2) {
    //
    // First get two cells, one above the other:
    //
    c = FirstPt_All();
    c_next = NextPt_All(c);
    //CI.print_cell(c);
    //CI.print_cell(c_next);
    while (c_next->pos[ZZ] == c->pos[ZZ]) {
      c_next=NextPt_All(c_next);
      //cout <<"c_next, id="<<c_next->id<<"\n";
      // Let c_next loop through x-y plane until i get to above the
      // first point.
    }
    //
    // Now c_next should be c's ZP neighbour.
    //
    do {
      //rep.printVec("C-",c->pos,G_ndim);
      //rep.printVec("C+",c_next->pos,G_ndim);
      c->ngb[ZP] = c_next;
      c_next->ngb[ZN] = c;
      c = NextPt_All(c);
      c_next = NextPt_All(c_next);
    } while (c_next!=0);
  }
  // ----------------------  SET ngb POINTERS  ----------------------
  // ----------------------------------------------------------------

  //  cout<<"\tAssignGridStructure Finished.\n";
  return 0;
} // assign_grid_structure()



// ##################################################################
// ##################################################################



int UniformGrid::allocate_grid_data()
{
#ifdef TESTING
  cout <<"\tAllocating grid data... G_ncell="<<G_ncell<<"\n";
#endif
  //
  // First generate G_ncell_all points and add a npt_all counter so
  // that they are in a singly linked list.
  //
  cell *c = G_fpt_all;
  size_t count=0;
  while (c->id < static_cast<long int>(G_ncell_all-1)) {
    c->npt_all = CI.new_cell();
    c = c->npt_all;
    count++;
    c->id = count;
  }
  c->npt_all = 0;
  G_lpt_all = c; 
  
#ifdef TESTING
  cout <<"\tFinished Allocating Data.\n";
#endif // TESTING
  return 0;
} // allocate_grid_data



// ##################################################################
// ##################################################################



int UniformGrid::set_cell_size()
{
  //
  // Uniform Cartesian grid, with cells that have the same length in
  // each direction, so this is very easy...
  //
#ifdef TESTING
  cout <<"\t Setting G_dx=constant for all cells.\n";
#endif

  G_dx = G_range[0]/(G_ng[0]);


  if(G_ndim>1) {
    if (!pconst.equalD(G_range[1]/G_dx, static_cast<double>(G_ng[1])))
      rep.error("Cells must be same length in each direction! Set the range and number of points appropriately.", G_range[1]/G_dx/G_ng[1]);
  }

  if (G_ndim>2) {
    if (!pconst.equalD(G_range[2]/G_dx, static_cast<double>(G_ng[2])))
      rep.error("Cells must be same length in each direction! Set the range and number of points appropriately.", G_range[2]/G_dx/G_ng[2]);
  }

  //
  // Surface area of interface: It is assumed extra dimensions are
  // per unit length.
  //
  if (G_ndim==1) G_dA = 1.; 
  else if (G_ndim==2) G_dA = G_dx;
  else G_dA = G_dx*G_dx;

  //
  // Volume of cell.
  //
  if (G_ndim==1) G_dV = G_dx;
  else if (G_ndim==2) G_dV = G_dx*G_dx;
  else  G_dV = G_dx*G_dx*G_dx;
  
  return 0;
} // set_cell_size



// ##################################################################
// ##################################################################



enum direction UniformGrid::OppDir(enum direction dir)
{
  // This should be in VectorOps, or somewhere like that...
  if (dir==XP) return(XN);
  else if (dir==XN) return(XP);
  else if (dir==YP) return(YN);
  else if (dir==YN) return(YP);
  else if (dir==ZP) return(ZN);
  else if (dir==ZN) return(ZP);
  else {rep.error("Bad direction given to OppDir",dir); return(NO);}
}



// ##################################################################
// ##################################################################



class cell* UniformGrid::FirstPt()
{
  ///
  /// First point is always the at the negative corner of the grid in
  /// all directions.
  ///
  return(G_fpt);
} // FirstPt



// ##################################################################
// ##################################################################



class cell* UniformGrid::FirstPt_All()
{
  ///
  /// First point is always the at the negative corner of the grid in
  /// all directions.
  ///
  return(G_fpt_all);
} // FirstPt



// ##################################################################
// ##################################################################



class cell* UniformGrid::LastPt()
{
//#ifdef TESTING
//  cout <<"Last Point is :"<<G_lpt; CI.print_cell(G_lpt);
//#endif
  return(G_lpt);
} // LastPt



// ##################################################################
// ##################################################################



class cell* UniformGrid::LastPt_All()
{
  return(G_lpt_all);
} // LastPt




// ##################################################################
// ##################################################################



class cell* UniformGrid::PrevPt(
      const class cell* p,
      enum direction dir
      )
{
  // Returns previous cell, including virtual boundary cells.
  //
  /// \section direction
  /// There is a good argument here for having the direction enums
  /// with negative values, so that XN=-XP.  But I am using the same
  /// enum for addressing elements of arrays, so I'm not sure what
  /// the way around this is.  The best thing is probably not to use
  /// this function, but to call NextPt in the reverse direction.
  ///
  enum direction opp = OppDir(dir);
#ifdef TESTING
  // This is going to be very inefficient...
  cout <<"This function is very inefficient and probably shouldn't be used.\n";
#endif // TESTING
  return(p->ngb[opp]);
}



// ##################################################################
// ##################################################################



int UniformGrid::SetupBCs(
      class SimParams &par  ///< List of simulation params (including BCs)
      )
{
  //
  // Set ntracer and ftr, bsed on SimParams
  //
  G_ntracer = par.ntracer;
  G_ftr     = G_nvar - par.ntracer;
  
  ///
  /// \section Ordering of Cells in Boundary data
  /// Cells are always ordered by increasing cell id.  So the cell
  /// with the most negative ZYZ coords comes first, then increasing
  /// in X, then in Y, and finally in Z.  i.e. the most slowly
  /// increasing index is Z.
  ///

  //
  // This function loops through all points on the grid, and if it is
  // an edge point, it finds which direction(s) are off the grid and
  // adds the boundary cells to boundary data.
  //
  class cell *c, *cy, *cz;

  // ----------------------------------------------------------------
  //
  // Set up X-boundaries, these are all cells adjacent to the grid 
  // in the XN and XP directions.  First we go up the XN boundary
  // along the YP and ZP directions (if multi-D).
  // - c is a temporary cell pointer to move through the boundary
  // data.
  // - cy is a pointer to the cells in increasing y
  // - cz is a pointer to the cells in increasing z
  //
  cz = FirstPt();
  // loop in ZP direction
  do {
    cy = cz;
    // loop in YP direction
    do {
      c=cy;
      // add boundary cells beside the grid.
      for (int v=0; v<BC_nbc; v++) c = NextPt(c, XN);
      if (!c) rep.error("Got lost on grid! XN",cy->id);
      for (int v=0; v<BC_nbc; v++) {
        BC_bd[XN]->data.push_back(c);
#ifdef TESTING
        cout << " Adding cell "<<c->id<<" to XN boundary.\n";
#endif // TESTING
        c = NextPt(c, XP);
      }
      if (G_ndim>1) cy=NextPt(cy,YP);
    } while (G_ndim>1 && cy!=0 && cy->isgd);
    if (G_ndim>2) cz=NextPt(cz,ZP);
  } while (G_ndim>2 && cz!=0 && cz->isgd);
#ifdef TESTING
  cout <<"** Setup XN boundary, got "<<BC_bd[XN]->data.size();
  cout <<" grid cells.\n";
#endif // TESTING
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // Now do the same for XP, beginning at the most negative point
  // at the XP boundary and moving to the most positive.
  //
  cz = FirstPt();
  while (NextPt(cz,XP)->isgd) cz=NextPt(cz,XP);
  // loop in ZP direction
  do {
    cy = cz;
    // loop in YP direction
    do {
      c=cy;
      // add boundary cells beside the grid.
      for (int v=0; v<BC_nbc; v++) {
        c = NextPt(c, XP);
        if (!c) {
          CI.print_cell(cy);
          rep.error("Got lost on grid! XP",cy->id);
        }
        BC_bd[XP]->data.push_back(c);
#ifdef TESTING
        cout << " Adding cell "<<c->id<<" to XP boundary.\n";
#endif // TESTING
      }
      if (G_ndim>1) cy=NextPt(cy,YP);
    } while (G_ndim>1 && cy!=0 && cy->isgd);
    if (G_ndim>2) cz=NextPt(cz,ZP);
  } while (G_ndim>2 && cz!=0 && cz->isgd);
#ifdef TESTING
  cout <<"** Setup XP boundary, got "<<BC_bd[XP]->data.size();
  cout <<" grid cells.\n";
#endif // TESTING
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // The YN/YP boundaries are different because they have the
  // corner data too.
  //
  if (G_ndim>1) {
    // --------------------------------------------------------------
    //
    // First YN boundary.
    // Get to most negative cell, and add cells, progressing to
    // the next one.
    //
    cz = FirstPt();
    while (NextPt(cz,XN) != 0) cz = NextPt(cz,XN);
    while (NextPt(cz,YN) != 0) cz = NextPt(cz,YN);
    //
    // loop in ZP-direction, at least once because there must be
    // at least one plane of boundary cells.
    //
    do {
      cy = cz;
      //
      // loop in the X-Y plane while we are in the YN boundary
      // region, and add all the cells.
      //
      do {
        BC_bd[YN]->data.push_back(cy);
        //cout << " Adding cell "<<cy->id<<" to YN boundary.\n";
        cy = NextPt_All(cy);
      } while (cy->pos[YY] < G_ixmin[YY]);

      if (G_ndim>2) cz = NextPt(cz, ZP);
    } while (G_ndim>2 && cz!=0 && cz->pos[ZZ]<G_ixmax[ZZ]);
#ifdef TESTING
    cout <<"** Setup YN boundary, got "<<BC_bd[YN]->data.size();
    cout <<" grid cells.\n";
#endif // TESTING
    // --------------------------------------------------------------

    // --------------------------------------------------------------
    //
    // Then YP boundary.
    // Get to the first row of YP boundary cells, move to the 
    // first cell in the XY plane, and go from there.
    //
    cz = FirstPt();
    while (cz->isgd) cz = NextPt(cz,YP);
    while (NextPt(cz,XN) != 0) cz = NextPt(cz,XN);

    //
    // loop in ZP-direction, at least once because there must be
    // at least one plane of boundary cells.
    //
    do {
      cy = cz;
      //
      // loop in the X-Y plane while we are in the YP boundary
      // region, and add all the cells.
      //
      do {
        BC_bd[YP]->data.push_back(cy);
        //cout << " Adding cell "<<cy->id<<" to YP boundary.\n";
        cy = NextPt_All(cy);
      } while ((cy !=0) && (cy->pos[YY] > G_ixmax[YY]));
      
      if (G_ndim>2) cz = NextPt(cz, ZP);
    } while (G_ndim>2 && cz!=0 && cz->pos[ZZ]<G_ixmax[ZZ]);
#ifdef TESTING
    cout <<"** Setup YP boundary, got "<<BC_bd[YP]->data.size();
    cout <<" grid cells.\n";
#endif // TESTING
    // --------------------------------------------------------------
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // Now the Z-boundary
  //
  if (G_ndim>2) {
    //
    // ZN is easy... all points until pos[Z] > xmin[Z]
    //
    cz = FirstPt_All();
    do {
      BC_bd[ZN]->data.push_back(cz);
      //cout << " Adding cell "<<cz->id<<" to ZN boundary.\n";
      cz = NextPt_All(cz);
    } while (cz->pos[ZZ] < G_ixmin[ZZ]);
#ifdef TESTING
    cout <<"** Setup ZN boundary, got "<<BC_bd[ZN]->data.size();
    cout <<" grid cells.\n";
#endif // TESTING
    //
    // ZP is also easy... all points with pos[Z] > xmax[Z]
    //
    cz = LastPt();
    while (cz->pos[ZZ] < G_ixmax[ZZ]) cz = NextPt_All(cz);
    do {
      BC_bd[ZP]->data.push_back(cz);
      //cout << " Adding cell "<<cz->id<<" to ZP boundary.\n";
      cz = NextPt_All(cz);
    } while (cz!=0);
#ifdef TESTING
    cout <<"** Setup ZP boundary, got "<<BC_bd[ZP]->data.size();
    cout <<" grid cells.\n";
#endif // TESTING
  }
  // ----------------------------------------------------------------
  return 0;
}


// ##################################################################
// ##################################################################



int UniformGrid::BC_printBCdata(boundary_data *b)
{
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {  
    CI.print_cell(*c);
  }
  return 0;
}



// ##################################################################
// ##################################################################



void UniformGrid::BC_deleteBoundaryData(
      boundary_data *b
      )
{
  if (b->refval !=0) {
    b->refval = mem.myfree(b->refval);
  }

  list<cell *>::iterator i=b->data.begin();
  if (b->data.empty()) {
#ifdef TESTING
    cout <<"BC destructor: No boundary cells to delete.\n";
#endif
  }
  else {
    do {
      b->data.erase(i);
      i=b->data.begin();
    }  while(i!=b->data.end());
  }

  if (b->send_data.size()>0) {
    i=b->send_data.begin();
    do {
      b->send_data.erase(i);
      i=b->send_data.begin();
    }  while(i!=b->send_data.end());
  }

  for (unsigned int j=0; j<b->NGrecvF2C.size(); j++) {
    i=b->NGrecvF2C[j].begin();
    do {
      b->NGrecvF2C[j].erase(i);
      i=b->NGrecvF2C[j].begin();
    }  while(i!=b->NGrecvF2C[j].end());
  }
  b->NGrecvF2C.clear();

  for (unsigned int j=0; j<b->NGrecvC2F.size(); j++) {
    i=b->NGrecvC2F[j].begin();
    do {
      b->NGrecvC2F[j].erase(i);
      i=b->NGrecvC2F[j].begin();
    }  while(i!=b->NGrecvC2F[j].end());
  }
  b->NGrecvC2F.clear();

  for (unsigned int j=0; j<b->NGsendC2F.size(); j++) {
    b->NGsendC2F[j]->c.clear();
    struct c2f *t = b->NGsendC2F[j];
    delete t;
    //b->NGsendC2F[j] = mem.myfree( b->NGsendC2F[j]);
  }
  b->NGsendC2F.clear();

  for (unsigned int j=0; j<b->avg.size(); j++) {
    b->avg[j].avg_state = mem.myfree( b->avg[j].avg_state);
    b->avg[j].c.clear();
  }
  b->avg.clear();

  return;
}


// ##################################################################
// ##################################################################



void UniformGrid::BC_deleteBoundaryData()
{
#ifdef TESTING
  cout <<"BC destructor: deleting Boundary data...\n";
#endif
  struct boundary_data *b;
  for (unsigned int ibd=0; ibd<BC_bd.size(); ibd++) {
    b = BC_bd[ibd];
    BC_deleteBoundaryData(b);
  } // loop over all boundaries.
  BC_bd.clear();
  return;
}
  

// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in physical units, and outputs in 
// physical units.
//
double UniformGrid::distance(
    const double *p1, ///< position 1 (physical)
    const double *p2  ///< position 2 (physical)
    )
{
  double temp=0.0;
  for (int i=0;i<G_ndim;i++)
    temp += (p1[i]-p2[i])*(p1[i]-p2[i]);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in code integer units, and outputs in
// integer units (but obviously the answer is not an integer).
//
double UniformGrid::idistance(
    const int *p1, ///< position 1 (integer)
    const int *p2  ///< position 2 (integer)
    )
{
  double temp=0.0;
  for (int i=0;i<G_ndim;i++)
    temp += static_cast<double>((p1[i]-p2[i])*(p1[i]-p2[i]));
  return sqrt(temp);
}
   


// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in physical units (e.g. centimetres).
//
double UniformGrid::distance_cell2cell(
    const cell *c1, ///< cell 1
    const cell *c2  ///< cell 2
    )
{
  double temp = 0.0;
  for (int i=0;i<G_ndim;i++)
    temp += pow(CI.get_dpos(c1,i)-CI.get_dpos(c2,i), 2.0);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in grid--integer units (one cell has a diameter
// two units).
//
double UniformGrid::idistance_cell2cell(
    const cell *c1, ///< cell 1
    const cell *c2  ///< cell 2
    )
{
  double temp = 0.0;
  for (int i=0;i<G_ndim;i++)
    temp += 
      pow(static_cast<double>(CI.get_ipos(c1,i)-CI.get_ipos(c2,i)),
    2.0);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are physical units.
//
double UniformGrid::distance_vertex2cell(
    const double *v, ///< vertex (physical)
    const cell *c    ///< cell
    )
{
  double temp = 0.0;
  for (int i=0;i<G_ndim;i++)
    temp += pow(v[i]-CI.get_dpos(c,i), 2.0);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// As distance_vertex2cell(double[],cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double UniformGrid::difference_vertex2cell(
    const double *v,  ///< vertex (double)
    const cell *c, ///< cell
    const axes a   ///< Axis to calculate.
    )
{
  return (CI.get_dpos(c,a)-v[a]);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are code-integer units.
//
double UniformGrid::idistance_vertex2cell(
    const int *v, ///< vertex (integer)
    const cell *c ///< cell
    )
{
  double temp = 0.0;
  for (int i=0;i<G_ndim;i++)
    temp += pow(static_cast<double>(v[i]-CI.get_ipos(c,i)), 2.0);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// As idistance_vertex2cell(int,cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double UniformGrid::idifference_vertex2cell(
    const int *v,  ///< vertex (integer)
    const cell *c, ///< cell
    const axes a   ///< Axis to calculate.
    )
{
  return (CI.get_ipos(c,a)-v[a]);
}



// ##################################################################
// ##################################################################



//
// As idifference_vertex2cell(int,cell,axis) but for the coordinate
// difference between two cell positions along a given axis.
// It returns *cell2* coordinate minus *cell1* coordinate.
//
double UniformGrid::idifference_cell2cell(
    const cell *c1, ///< cell 1
    const cell *c2, ///< cell 2
    const axes a    ///< Axis.
    )
{
  return (CI.get_ipos(c2,a)-CI.get_ipos(c1,a));
}



// ##################################################################
// ##################################################################



bool UniformGrid::point_on_grid(
    const double *pos ///< position
    )
{
  bool on=true;
  for (int v=0;v<G_ndim;v++) {
    if (pos[v] < G_xmin[v]  ||  pos[v] > G_xmax[v]) on=false;
  }
  return on;
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


//
// Constructor
//
uniform_grid_cyl::uniform_grid_cyl(
    int nd,         ///< ndim, length of position vector.
    int nv,         ///< nvar, length of state vectors.
    int eqt,        ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,        ///< Number of boundary cells to use.
    double *g_xn,   ///< array of minimum values of x,y,z for this grid.
    double *g_xp,   ///< array of maximum values of x,y,z for this grid.
    int *g_nc,      ///< array of number of cells in x,y,z directions.
    double *lev_xn, // level xmin
    double *lev_xp, // level xmax
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,lev_xn,lev_xp,sim_xn,sim_xp),
  VectorOps_Cyl(nd)
{
#ifdef TESTING
  cout <<"Setting up cylindrical uniform grid with";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif
  if (G_ndim!=2)
    rep.error("Need to write code for !=2 dimensions",G_ndim);
  G_coordsys = COORD_CYL;  // Cylindrical Coordinate system

#ifdef TESTING
  cout <<"cylindrical grid: dr="<<G_dx<<"\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



uniform_grid_cyl::~uniform_grid_cyl()
{
#ifdef TESTING
  cout <<"uniform_grid_cyl destructor. Present and correct!\n";
#endif
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
#ifdef PARALLEL
  rep.error("redefine iR_cov() for parallel grid","please");
#endif
  //cout <<" Cell radius: "<< R_com(c)/CI.phys_per_int() +G_ixmin[Rcyl];
  //rep.printVec("  cell centre",c->pos,G_ndim);
  return (R_com(c, G_dx)-G_xmin[Rcyl])/CI.phys_per_int() +G_ixmin[Rcyl];
}
  

// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in physical units, and outputs in 
// physical units.
//
double uniform_grid_cyl::distance(
          const double *p1, ///< position 1 (physical)
          const double *p2  ///< position 2 (physical)
          )
{
  //
  // This is the same as for cartesian as long as there is no 3rd
  // dimension.
  //
  double temp=0.0;
  for (int i=0;i<G_ndim;i++)
    temp += (p1[i]-p2[i])*(p1[i]-p2[i]);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in code integer units, and outputs in
// integer units (but obviously the answer is not an integer).
//
double uniform_grid_cyl::idistance(
          const int *p1, ///< position 1 (integer)
          const int *p2  ///< position 2 (integer)
          )
{
  //
  // This is the same as for cartesian as long as there is no 3rd
  // dimension.
  //
  double temp=0.0;
  for (int i=0;i<G_ndim;i++)
    temp += static_cast<double>((p1[i]-p2[i])*(p1[i]-p2[i]));
  return sqrt(temp);
}
   


// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in physical units (e.g. centimetres).
//
double uniform_grid_cyl::distance_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2  ///< cell 2
      )
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d=0.0, temp=0.0;
  temp = CI.get_dpos(c1,Zcyl)-CI.get_dpos(c2,Zcyl);
  d += temp*temp;
  temp = R_com(c1,G_dx) - R_com(c2,G_dx);
  d += temp*temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in grid--integer units (one cell has a diameter
// two units).
//
double uniform_grid_cyl::idistance_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2  ///< cell 2
      )
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d=0.0, temp;
    //r1,r2;
  temp = static_cast<double>(CI.get_ipos(c1,Zcyl)-CI.get_ipos(c2,Zcyl));
  d += temp*temp;
  temp = iR_cov(c1) - iR_cov(c2);
  d += temp*temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are physical units.
//
double uniform_grid_cyl::distance_vertex2cell(
          const double *v, ///< vertex (physical)
          const cell *c    ///< cell
          )
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d=0.0, temp=0.0;
  temp = v[Zcyl] -CI.get_dpos(c,Zcyl);
  d += temp*temp;
  temp = v[Rcyl] - R_com(c,G_dx);
  d += temp*temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



//
// As distance_vertex2cell(double[],cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double uniform_grid_cyl::difference_vertex2cell(
      const double *v,  ///< vertex (double)
      const cell *c, ///< cell
      const axes a   ///< Axis to calculate.
      )
{
  if      (a==Zcyl) {
    return (CI.get_dpos(c,a)-v[a]);
  }
  else if (a==Rcyl) {
    return R_com(c,G_dx) - v[Rcyl];
  }
  else {
    cerr <<" Requested cylindrical distance in theta dir.\n";
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################


//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are code-integer units.
//
double uniform_grid_cyl::idistance_vertex2cell(
          const int *v, ///< vertex (integer)
          const cell *c ///< cell
          )
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  //cout <<"cylindrical vertex2cell distance...\n";
  //rep.printVec("idist vertex",v,G_ndim);
  //rep.printVec("idist cell  ",c->pos,G_ndim);
  double d=0.0, temp=0.0;
  temp = static_cast<double>(v[Zcyl]  -CI.get_ipos(c,Zcyl));
  d += temp*temp;
  //cout <<"Z-dist = "<<temp;
  // iR_cov() gives coordinates of centre-of-volume radius in integer
  // units.
  temp = static_cast<double>(v[Rcyl]) -iR_cov(c);
  d += temp*temp;
  //cout <<",  idist : R-dist = "<<temp<< " total dist = "<<sqrt(d)<<"\n";
  return sqrt(d);

}



// ##################################################################
// ##################################################################


//
// As idistance_vertex2cell(int,cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double uniform_grid_cyl::idifference_vertex2cell(const int *v,  ///< vertex (integer)
             const cell *c, ///< cell
             const axes a   ///< Axis to calculate.
             )
{
  if      (a==Zcyl)
    return (CI.get_ipos(c,a)-v[a]);
  else if (a==Rcyl) {
    //cout <<"cov="<<iR_cov(c)<<", v="<<v[a]<<"\n";
    return (iR_cov(c) -v[a]);
  }
  else 
    rep.error("Bad axis for uniform_grid_cyl::idifference_vertex2cell()",a);

  return -1.0;
}



// ##################################################################
// ##################################################################



//
// As idifference_vertex2cell(int,cell,axis) but for the coordinate
// difference between two cell positions along a given axis.
// It returns *cell2* coordinate minus *cell1* coordinate.
//
double uniform_grid_cyl::idifference_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2, ///< cell 2
      const axes a    ///< Axis.
      )
{
  if      (a==Zcyl)
    return (CI.get_ipos(c2,a)-CI.get_ipos(c1,a));
  else if (a==Rcyl)
    return (iR_cov(c2) -iR_cov(c1));
  else 
    rep.error("Bad axis for uniform_grid_cyl::idifference_cell2cell()",a);
  return -1.0;
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



uniform_grid_sph::uniform_grid_sph(
    int nd,         ///< ndim, length of position vector.
    int nv,         ///< nvar, length of state vectors.
    int eqt,        ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,        ///< Number of boundary cells to use.
    double *g_xn,   ///< array of minimum values of x,y,z for this grid.
    double *g_xp,   ///< array of maximum values of x,y,z for this grid.
    int *g_nc,      ///< array of number of cells in x,y,z directions.
    double *lev_xn, // level xmin
    double *lev_xp, // level xmax
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,lev_xn,lev_xp,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd)
{
#ifdef TESTING
  cout <<"Setting up spherical uniform grid with";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif
  if (G_ndim!=1)
    rep.error("Need to write code for >1 dimension",G_ndim);
  G_coordsys = COORD_SPH;  // Spherical Coordinate system

#ifdef TESTING
  cout <<"spherical grid: dr="<<G_dx<<"\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



uniform_grid_sph::~uniform_grid_sph()
{
#ifdef TESTING
  cout <<"uniform_grid_sph destructor. Present and correct!\n";
#endif
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
#ifdef PARALLEL
  rep.error("redefine iR_cov() for parallel grid","please");
#endif
  return ((R_com(c,G_dx)-G_xmin[Rsph])/CI.phys_per_int() +G_ixmin[Rsph]);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in physical units, and outputs in 
// physical units.
//
double uniform_grid_sph::distance(const double *p1, ///< position 1 (physical)
      const double *p2  ///< position 2 (physical)
      )
{
  return fabs(p1[Rsph]-p2[Rsph]);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two points, where the two position
// are interpreted in the appropriate geometry.
// This function takes input in code integer units, and outputs in
// integer units (but obviously the answer is not an integer).
//
double uniform_grid_sph::idistance(const int *p1, ///< position 1 (integer)
      const int *p2  ///< position 2 (integer)
      )
{
  return fabs(static_cast<double>(p1[Rsph]-p2[Rsph]));
}
   


// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in physical units (e.g. centimetres).
//
double uniform_grid_sph::distance_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2  ///< cell 2
      )
{
  return fabs(R_com(c1,G_dx)-R_com(c2,G_dx));
}



// ##################################################################
// ##################################################################



//
// Calculate distance between two cell--centres (will be between
// centre-of-volume of cells if non-cartesian geometry).
// Result returned in grid--integer units (one cell has a diameter
// two units).
//
double uniform_grid_sph::idistance_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2  ///< cell 2
      )
{
  return fabs(R_com(c1,G_dx)-R_com(c2,G_dx))/CI.phys_per_int();
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are physical units.
//
double uniform_grid_sph::distance_vertex2cell(
      const double *v, ///< vertex (physical)
      const cell *c    ///< cell
      )
{
  return fabs(v[Rsph] - R_com(c,G_dx));
}




// ##################################################################
// ##################################################################



//
// As distance_vertex2cell(double[],cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double uniform_grid_sph::difference_vertex2cell(
      const double *v,  ///< vertex (double)
      const cell *c, ///< cell
      const axes a   ///< Axis to calculate.
      )
{
  if      (a==Rsph) {
    return R_com(c,G_dx) - v[Rsph];
  }
  else {
    cerr <<" Requested spherical distance in not-radial dir.\n";
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are code-integer units.
//
double uniform_grid_sph::idistance_vertex2cell(
      const int *v, ///< vertex (integer)
      const cell *c ///< cell
      )
{
  //cout <<"idist_v2c: iV[0]="<<v[Rsph]<<", iR(c)="<<iR_cov(c);
  //cout <<", idist="<<fabs(static_cast<double>(v[Rsph]) -iR_cov(c))<<"\n";
  return
    fabs(static_cast<double>(v[Rsph]) -iR_cov(c));
}



// ##################################################################
// ##################################################################



//
// As idistance_vertex2cell(int,cell) but for a single component
// of the position vector, and not the absolute value.  It returns
// the *cell* coordinate minus the *vertex* coordinate.
//
double uniform_grid_sph::idifference_vertex2cell(
      const int *v,  ///< vertex (integer)
      const cell *c, ///< cell
      const axes a   ///< Axis to calculate.
      )
{
  if (a==Rsph)
    return (iR_cov(c) -static_cast<double>(v[a]));
  else 
    rep.error("Bad axis for uniform_grid_sph::idifference_vertex2cell()",a);

  return -1.0;
}



// ##################################################################
// ##################################################################



//
// As idifference_vertex2cell(int,cell,axis) but for the coordinate
// difference between two cell positions along a given axis.
// It returns *cell2* coordinate minus *cell1* coordinate.
//
double uniform_grid_sph::idifference_cell2cell(
      const cell *c1, ///< cell 1
      const cell *c2, ///< cell 2
      const axes a    ///< Axis.
      )
{
  if (a==Rsph)
    return (iR_cov(c2) -iR_cov(c1));
  else 
    rep.error("Bad axis for uniform_grid_sph::idifference_cell2cell()",a);
  return -1.0;
}




// ##################################################################
// ##################################################################



//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------

