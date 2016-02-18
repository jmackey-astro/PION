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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "grid/uniform_grid.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include <fstream>
#include <iostream>
using namespace std;


//#define GLM_ZERO_BOUNDARY ///< Set this flag to make Psi=0 on boundary cells.
#define GLM_NEGATIVE_BOUNDARY ///< Set this flag for Psi[boundary cell]=-Psi[edge cell]



// ##################################################################
// ##################################################################



UniformGrid::UniformGrid(
    int nd,
    int nv,
    int eqt,
    int Nbc,       ///< Number of boundary cells to use.
    double *g_xn,
    double *g_xp,
    int *g_nc,
    double *sim_xn,
    double *sim_xp
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

  //
  // Allocate arrays for dimensions of grid.
  //
  G_ng=0;
  G_xmin=0;
  G_xmax=0;
  G_range=0;
  
  G_ng     = mem.myalloc(G_ng,    G_ndim);

  G_xmin   = mem.myalloc(G_xmin,  G_ndim);
  G_xmax   = mem.myalloc(G_xmax,  G_ndim);
  G_range  = mem.myalloc(G_range, G_ndim);

  G_ixmin  = mem.myalloc(G_ixmin, G_ndim);
  G_ixmax  = mem.myalloc(G_ixmax, G_ndim);
  G_irange = mem.myalloc(G_irange,G_ndim);

  //
  // For serial code the Sim_XXX arrays are the same as G_XXX arrays,
  // but for parallel code with many grids this is useful.
  //
  Sim_xmin   = mem.myalloc(Sim_xmin,  G_ndim);
  Sim_xmax   = mem.myalloc(Sim_xmax,  G_ndim);
  Sim_range  = mem.myalloc(Sim_range, G_ndim);
  Sim_ixmin  = mem.myalloc(Sim_ixmin, G_ndim);
  Sim_ixmax  = mem.myalloc(Sim_ixmax, G_ndim);
  Sim_irange = mem.myalloc(Sim_irange,G_ndim);


  //
  // Assign values to grid dimensions and set the cell-size
  //
  G_ncell = 1;
  G_ncell_all = 1;

  for (int i=0;i<G_ndim;i++) {
    G_ng[i] = g_nc[i];
    G_ncell *= G_ng[i]; // Get total number of cells.
    G_ncell_all *= G_ng[i] +2*BC_nbc; // all cells including boundary.
    G_xmin[i] = g_xn[i];
    G_xmax[i] = g_xp[i];
    G_range[i] = G_xmax[i]-G_xmin[i];
    if(G_range[i]<0.) rep.error("Negative range in direction i",i);

    //
    // For single static grid, Sim_xmin/xmax/range are the same as
    // for the grid (becuase there is only one grid).
    //
    Sim_xmin[i] = G_xmin[i];
    Sim_xmax[i] = G_xmax[i];
    Sim_range[i] = G_range[i];
  }

#ifdef TESTING
  cout <<"MIN.MAX for x = "<<G_xmin[XX]<<"\t"<<G_xmax[XX]<<"\n";
  cout <<"Setting cell size...\n";
#endif
  //
  // Checks grid dimensions and discretisation is reasonable,
  // and that it gives cubic cells.
  //
  set_cell_size();

  //
  // Now create the first cell, and then allocate data from there.
  // Safe to assume we have at least one cell...
  //
  cout <<" done.\n Initialising first cell...\n";
  G_fpt_all = CI.new_cell();
  G_fpt_all->id = 0;
  cout <<" done.\n";
  if(G_fpt_all==0) {
    rep.error("Couldn't assign memory to first cell in grid.",
              G_fpt_all);
  }
  
  //
  // assign memory for all cells, including boundary cells.
  //
  cout <<"\t allocating memory for grid.\n";
  int err = allocate_grid_data();
  if(err!=0)
    rep.error("Error setting up grid, allocate_grid_data",err);

  //
  // assign grid structure on cells, setting positions and ngb
  // pointers.
  //
  cout <<"\t assigning pointers to neighbours.\n";
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
  BC_bd=0; BC_nbd=-1;
  Wind =0;

  //
  // Set integer dimensions/location of grid.
  //
  CI.get_ipos_vec(G_xmin, G_ixmin );
  CI.get_ipos_vec(G_xmax, G_ixmax );
  for (int v=0;v<G_ndim;v++) G_irange[v] = G_ixmax[v]-G_ixmin[v];
#ifdef TESTING
  rep.printVec("grid ixmin ", G_ixmin, G_ndim);
  rep.printVec("grid ixmax ", G_ixmax, G_ndim);
  rep.printVec("grid irange", G_irange,G_ndim);
#endif

  //
  // Set integer dimensions/location of Simulation (same as grid)
  //
  for (int v=0;v<G_ndim;v++) {
    Sim_irange[v] = G_irange[v];
    Sim_ixmax[v]  = G_ixmax[v];
    Sim_ixmin[v]  = G_ixmin[v];
  }

  cout <<"Cartesian grid: dr="<<G_dx<<"\n";
  set_dx(G_dx);

#ifdef TESTING
  cout <<"UniformGrid Constructor done.\n";
#endif
} //UniformGrid Constructor



// ##################################################################
// ##################################################################


UniformGrid::~UniformGrid()
{

  //  cout <<"\tUniformGrid Destructor: delete boundary data if it exists...\n";
  if (BC_bd !=0) BC_deleteBoundaryData();
  //  cout <<"\tdone.\n";   
  // Delete the list of data I have just created...
  //  cout <<"\tDeleting grid data...\n";

  cell *cpt=FirstPt(), *npt=NextPt(FirstPt());
  do {
    CI.delete_cell(cpt);
  } while ( (npt=NextPt(cpt=npt))!=0 );
  CI.delete_cell(cpt);

  if (Wind) {delete Wind; Wind=0;}

  G_ng    = mem.myfree(G_ng);
  G_xmin  = mem.myfree(G_xmin);
  G_xmax  = mem.myfree(G_xmax);
  G_range = mem.myfree(G_range);
  G_ixmin  = mem.myfree(G_ixmin);
  G_ixmax  = mem.myfree(G_ixmax);
  G_irange = mem.myfree(G_irange);

  Sim_ixmin  = mem.myfree(Sim_ixmin);
  Sim_ixmax  = mem.myfree(Sim_ixmax);
  Sim_irange = mem.myfree(Sim_irange);


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
  int ipos[MAX_DIM], offset[MAX_DIM];
  for (int i=0; i<MAX_DIM; i++) ipos[i] = offset[i] = 0;
  for (int i=0; i<G_ndim; i++) {
    offset[i] = 1+ 2*static_cast<int>(ONE_PLUS_EPS*(G_xmin[i]-Sim_xmin[i])/G_dx);
#ifdef TESTING
    cout <<"************OFFSET["<<i<<"] = "<<offset[i]<<"\n";
#endif
  }

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
    cout <<"Cell positions: id = "<<c->id<<"\n";
#endif // TESTING

    //
    // Assign positions, for integer positions the first on-grid cell
    // is at [1,1,1], and a cell is 2 units across, so the second
    // cell is at [3,1,1] etc.
    // Boundary cells can have negative positions.
    //
    for (int i=0; i<G_ndim; i++) ipos[i] = offset[i] + 2*ix[i];
    CI.set_pos(c,ipos);
#ifdef TESTING
    rep.printVec("    pos", c->pos, G_ndim);
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
    rep.printVec("    dpos",dpos,G_ndim);
    rep.printVec("    xmax",G_xmax,G_ndim);
#endif // TESTING



    if (on_grid) {
      c->isgd = true;
      c->isbd = false;
#ifdef TESTING
      cout <<"    cell is on grid";
#endif // TESTING
    }
    else {
      c->isgd = false;
      c->isbd = true;
#ifdef TESTING
      cout <<"    cell NOT on grid";
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
  cell *c_prev=0, *c_next = NextPt_All(c);
  do {
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
    c_next = NextPt_All(c_next);
  }
  while (c_next != 0);

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
    while (CI.get_ipos(c_next,ZZ) == CI.get_ipos(c,ZZ)) {
      c_next=NextPt(c_next);
      // Let c_next loop through x-y plane until i get to above the
      // first point.
    }
    //
    // Now c_next should be c's ZP neighbour.
    //
    do {
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
  cout <<"\tAllocating grid data... G_ncell="<<G_ncell<<"\n";
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
  // Uniform Cartesian grid, with cubic cells, so this is easy, and
  // every cell has the same dimensions.
  //
#ifdef TESTING
  cout <<"\t Setting G_dx=constant for all cells.\n";
#endif

  G_dx = G_range[0]/(G_ng[0]);

  //
  // Set Cell dx in cell interface class, and also xmin.
  //
  CI.set_dx(G_dx);
  CI.set_xmin(Sim_xmin);

  if(G_ndim>1) {
    if (!pconst.equalD(G_range[1]/G_dx, static_cast<double>(G_ng[1])))
      rep.error("Cells are not cubic! Set the range and number of \
        points appropriately.", G_range[1]/G_dx/G_ng[1]);
  }

  if (G_ndim>2) {
    if (!pconst.equalD(G_range[2]/G_dx, static_cast<double>(G_ng[2])))
      rep.error("Cells are not cubic! Set the range and number of \
        points appropriately.", G_range[2]/G_dx/G_ng[2]);
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
#ifdef TESTING
  cout <<"Last Point is :"<<G_lpt; CI.print_cell(G_lpt);
#endif
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
      const int Nbc,
      string typeofbc
      )
{
  if (BC_setBCtypes(typeofbc) !=0)
    rep.error("UniformGrid::setupBCs() Failed to set types of BC",1);
  
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
  int err=0;

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
        BC_bd[XN].data.push_back(c);
        cout << " Adding cell "<<c->id<<" to XN boundary.\n";
        c = NextPt(c, XP);
      }
      cy=NextPt(cy,YP);
    } while (G_ndim>1 && cy!=0 && cy->isgd);
    cz=NextPt(cz,ZP);
  } while (G_ndim>2 && cz!=0 && cz->isgd);
  cout <<"** Setup XN boundary, got "<<BC_bd[XN].data.size();
  cout <<" grid cells.\n";
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
        BC_bd[XP].data.push_back(c);
        cout << " Adding cell "<<c->id<<" to XP boundary.\n";
      }
      cy=NextPt(cy,YP);
    } while (G_ndim>1 && cy!=0 && cy->isgd);
    cz=NextPt(cz,ZP);
  } while (G_ndim>2 && cz!=0 && cz->isgd);
  cout <<"** Setup XP boundary, got "<<BC_bd[XP].data.size();
  cout <<" grid cells.\n";
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
        BC_bd[YN].data.push_back(cy);
        cout << " Adding cell "<<cy->id<<" to YN boundary.\n";
        cy = NextPt_All(cy);
      } while (cy->pos[YY] < G_xmin[YY]);

      cz = NextPt(cz, ZP);
    } while (G_ndim>2 && cz!=0 && cz->isgd);
    cout <<"** Setup YN boundary, got "<<BC_bd[YN].data.size();
    cout <<" grid cells.\n";
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
        BC_bd[YP].data.push_back(cy);
        cout << " Adding cell "<<cy->id<<" to YP boundary.\n";
        cy = NextPt_All(cy);
      } while (cy !=0 && cy->pos[YY] > G_xmax[YY]);

      cz = NextPt(cz, ZP);
    } while (G_ndim>2 && cz!=0 && cz->isgd);
    cout <<"** Setup YP boundary, got "<<BC_bd[YP].data.size();
    cout <<" grid cells.\n";
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
      BC_bd[ZN].data.push_back(cz);
      cout << " Adding cell "<<cz->id<<" to ZN boundary.\n";
      cz = NextPt_All(cz);
    } while (cz->pos[ZZ] < G_xmin[ZZ]);
    //
    // ZP is also easy... all points with pos[Z] > xmax[Z]
    //
    cz = LastPt();
    while (cz->pos[ZZ] < G_xmax[ZZ]) cz = NextPt(cz);
    do {
      BC_bd[ZP].data.push_back(cz);
      cout << " Adding cell "<<cz->id<<" to ZP boundary.\n";
      cz = NextPt_All(cz);
    } while (cz!=0);
  }
  // ----------------------------------------------------------------


  }
  // ----------------------------------------------------------------

  //
  // Loop through all boundaries, and assign data to them.
  //
  for (int i=0; i<BC_nbd; i++) {
    switch (BC_bd[i].itype) {
     case PERIODIC:   err += BC_assign_PERIODIC(  &BC_bd[i]); break;
     case OUTFLOW:    err += BC_assign_OUTFLOW(   &BC_bd[i]); break;
     case ONEWAY_OUT: err += BC_assign_ONEWAY_OUT(&BC_bd[i]); break;
     case INFLOW:     err += BC_assign_INFLOW(    &BC_bd[i]); break;
     case REFLECTING: err += BC_assign_REFLECTING(&BC_bd[i]); break;
     case FIXED:      err += BC_assign_FIXED(     &BC_bd[i]); break;
     case JETBC:      err += BC_assign_JETBC(     &BC_bd[i]); break;
     case JETREFLECT: err += BC_assign_JETREFLECT(&BC_bd[i]); break;
     case DMACH:      err += BC_assign_DMACH(     &BC_bd[i]); break;
     case DMACH2:     err += BC_assign_DMACH2(    &BC_bd[i]); break;
     case RADSHOCK:   err += BC_assign_RADSHOCK(  &BC_bd[i]); break;
     case RADSH2:     err += BC_assign_RADSH2(    &BC_bd[i]); break;
     case STWIND:     err += BC_assign_STWIND(    &BC_bd[i]); break;
     case STARBENCH1: err += BC_assign_STARBENCH1(&BC_bd[i]); break;
     default:
      rep.warning("Unhandled BC",BC_bd[i].itype,-1); err+=1; break;
    }

  }
  return(err);
}



// ##################################################################
// ##################################################################


int UniformGrid::BC_setBCtypes(
  std::string bctype
  )
{
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif
  if(bctype=="FIXED" || bctype=="PERIODIC" || bctype=="ABSORBING") {
#ifdef TESTING
    cout <<"using old-style boundary condition specifier: "<<bctype;
    cout <<" on all sides.\n";
    cout <<"Converting to new style specifier.\n";
#endif

    if (bctype=="FIXED") {
      bctype = "XNfix_XPfix_";
      if (G_ndim>1) bctype += "YNfix_YPfix_";
      if (G_ndim>2) bctype += "ZNfix_ZPfix_";
    }
    else if (bctype=="PERIODIC") {
      bctype = "XNper_XPper_";
      if (G_ndim>1) bctype += "YNper_YPper_";
      if (G_ndim>2) bctype += "ZNper_ZPper_";
    }
    else if (bctype=="ABSORBING") {
      bctype = "XNout_XPout_";
      if (G_ndim>1) bctype += "YNout_YPout_";
      if (G_ndim>2) bctype += "ZNout_ZPout_";
    }
#ifdef TESTING
    cout <<"New bctype string = "<<bctype<<"\n";
#endif
  }    

  // Set up boundary data struct.  Assumes bctype is in format "XPper XNper " etc.,
  // so that each boundary is defined by 6 characters, and number of boundaries is
  // given by length/6.
  int len = bctype.length(); len = (len+5)/6;
  if (len < 2*G_ndim) rep.error("Need boundaries on all sides!",len);
#ifdef TESTING
  cout <<"Got "<<len<<" boundaries to set up.\n";
#endif
  BC_bd=0;
  BC_bd = mem.myalloc(BC_bd, len);
  UniformGrid::BC_nbd = len;
  
  int i=0;
  string::size_type pos;
  string d[6] = {"XN","XP","YN","YP","ZN","ZP"};
  for (i=0; i<2*G_ndim; i++) {
    BC_bd[i].dir = static_cast<direction>(i); //XN=0,XP=1,YN=2,YP=3,ZN=4,ZP=5
    BC_bd[i].baxis = static_cast<axes>(i/2);
    //
    // odd values of i are positive boundaries, others are negative.
    //
    if ((i+2)%2 !=0) {
      BC_bd[i].bloc  = G_xmax[BC_bd[i].baxis];
      BC_bd[i].bpos  = true;
    }
    else {
      BC_bd[i].bloc  = G_xmin[BC_bd[i].baxis];
      BC_bd[i].bpos  = false;
    }
    //
    // find boundary condition specified:
    //
    if ( (pos=bctype.find(d[i])) == string::npos)
      rep.error("Couldn't find boundary condition for ",d[i]);
    BC_bd[i].type = bctype.substr(pos+2,3);

    if      (BC_bd[i].type=="per") {
      BC_bd[i].itype=PERIODIC;
      BC_bd[i].type="PERIODIC";
    }
    else if (BC_bd[i].type=="out" || BC_bd[i].type=="abs") {
      BC_bd[i].itype=OUTFLOW;
      BC_bd[i].type="OUTFLOW";
    }
    else if (BC_bd[i].type=="owo") {
      BC_bd[i].itype=ONEWAY_OUT;
      BC_bd[i].type="ONEWAY_OUT";
    }
    else if (BC_bd[i].type=="inf") {
      BC_bd[i].itype=INFLOW ;
      BC_bd[i].type="INFLOW";
    }
    else if (BC_bd[i].type=="ref") {
      BC_bd[i].itype=REFLECTING;
      BC_bd[i].type="REFLECTING";
    }
    else if (BC_bd[i].type=="jrf") {
      BC_bd[i].itype=JETREFLECT;
      BC_bd[i].type="JETREFLECT";
    }
    else if (BC_bd[i].type=="fix") {
      BC_bd[i].itype=FIXED;
      BC_bd[i].type="FIXED";
    }
    else if (BC_bd[i].type=="dmr") {
      BC_bd[i].itype=DMACH;
      BC_bd[i].type="DMACH";
    }
    else if (BC_bd[i].type=="sb1") {
      BC_bd[i].itype=STARBENCH1;
      BC_bd[i].type="STARBENCH1";  // Wall for Tremblin mixing test.
    }
    else {
      rep.error("Don't know this BC type",BC_bd[i].type);
    }

    if(!BC_bd[i].data.empty())
      rep.error("Boundary data not empty in constructor!",BC_bd[i].data.size());
    BC_bd[i].refval=0;
#ifdef TESTING
    cout <<"\tBoundary type "<<i<<" is "<<BC_bd[i].type<<"\n";
#endif
  }
  if (i<BC_nbd) {
#ifdef TESTING
    cout <<"Got "<<i<<" boundaries, but have "<<BC_nbd<<" boundaries.\n";
    cout <<"Must have extra BCs... checking if I know what they are.\n";
#endif
    do {
      BC_bd[i].dir = NO;
      if ( (pos=bctype.find("IN",i*6)) ==string::npos) {
        for (int ii=i; ii<BC_nbd; ii++)
          BC_bd[ii].refval=0;
          cout <<"\tEncountered an unrecognised boundary condition.\n";
          cout <<"\tPerhaps boundaries are specified for Y or Z";
          cout <<" directions in a lower dimensional simulation.\n";
          cout <<"\tPlease check the parameter BC in the parameter ";
          cout <<"file used for setting up this simulation.\n";
          rep.error("Couldn't find boundary condition for extra boundary condition","IN");
      }
      BC_bd[i].type = bctype.substr(pos+2,3);
      if      (BC_bd[i].type=="jet") {BC_bd[i].itype=JETBC;    BC_bd[i].type="JETBC";}
      else if (BC_bd[i].type=="dm2") {BC_bd[i].itype=DMACH2;   BC_bd[i].type="DMACH2";}
      else if (BC_bd[i].type=="rsh") {BC_bd[i].itype=RADSHOCK; BC_bd[i].type="RADSHOCK";}
      else if (BC_bd[i].type=="rs2") {BC_bd[i].itype=RADSH2;   BC_bd[i].type="RADSH2";}
      else if (BC_bd[i].type=="wnd") {BC_bd[i].itype=STWIND;   BC_bd[i].type="STWIND";}
      else rep.error("Don't know this BC type",BC_bd[i].type);
      if(!BC_bd[i].data.empty())
        rep.error("Boundary data not empty in constructor!",BC_bd[i].data.size());
      BC_bd[i].refval=0;
//#ifdef TESTING
      cout <<"\tBoundary type "<<i<<" is "<<BC_bd[i].type<<"\n";
//#endif
      i++;
    } while (i<BC_nbd);
  }
#ifdef TESTING
  cout <<"BC types and data set up.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_PERIODIC(  boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);

  if (b->data.empty())
    rep.error("BC_assign_PERIODIC: empty boundary data",b->itype);

  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  do{
    temp=(*bpt);
    while (NextPt(temp,ondir) && NextPt(temp,ondir)->isgd)
      temp=NextPt(temp,ondir);

    if(!temp->isgd)
      rep.error("BC_assign_PERIODIC: Got lost assigning periodic BCs",
                temp->id);
    //
    // Now temp is last grid point before i hit the opposite boundary.
    // So copy this cells value into bpt, and the next one too if 2nd order.
    //
    for (int i=0; i<BC_nbc; i++) {
      for (int v=0;v<G_nvar;v++) (*bpt)->P[v] = (*bpt)->Ph[v] = temp->P[v];
      (*bpt)->npt = temp;
      ct++;
      ++bpt;
      temp = NextPt(temp,offdir);
      if (!temp) rep.error("no grid point!",temp);
    }
  } while (bpt !=b->data.end());

  if (ct != b->data.size())
    rep.error("BC_assign_PERIODIC: missed some cells!",ct-b->data.size());

  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_OUTFLOW(   boundary_data *b)
{
  //
  // Zero order extrapolation, if edge cell is index k, and boundary cells
  // k+1 and k+2, then P_{k+1} = P_{k+2} = P_{k}
  // First order extrapolation would be P_{k+1} = 2P_{k} - P_{k-1} but is
  // said to have stability problems sometimes (LeVeque, S.7.2.1, p131-132)
  //
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);

  if (b->data.empty()) {
    rep.error("BC_assign_OUTFLOW: empty boundary data",b->itype);
  }
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  //
  // loop through all boundary cells, set npt to point to the grid
  // cell where we get the data.
  //
  do{
    //
    // Find the cell to point to.  This should be the first cell on
    // the other side of the boundary, located at x[baxis]=bloc.
    // --------END OF 18.02.2016--------
    temp = NextPt((*bpt),ondir);
    if(temp==0) {
      rep.error("Got lost assigning Absorbing bcs.",temp->id);
    }

    for (int i=0; i<BC_nbc; i++) {
      for (int v=0;v<G_nvar;v++) (*bpt)->P[v]  = temp->P[v];
      for (int v=0;v<G_nvar;v++) (*bpt)->Ph[v] = temp->P[v];
      (*bpt)->npt = temp;

#ifdef GLM_ZERO_BOUNDARY
      if (G_eqntype==EQGLM)
        (*bpt)->P[SI]=(*bpt)->Ph[SI]=0.0;
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
      if (G_eqntype==EQGLM) {
        if (G_ndim==1) {
          rep.error("Psi outflow boundary condition doesn't work for 1D!",99);
        }
        if ((*bpt)->id==-1)
          (*bpt)->P[SI]=(*bpt)->Ph[SI]=-temp->P[SI];
        else if ((*bpt)->id==-2) {
          // need to set a pointer to cell [edge-1] for Psi
          // Have to be a bit sneaky and use a null ngb pointer.
          // Can't use ondir/offdir b/c we rely on getting 0 from
          // NextPt() at the edge of the grid.  Use isedge index
          // to tell us which ngb[] is to the [edge-1] cell.
          //cout <<"cell isedge: "<<(*bpt)->isedge<<"\n";
          int n=0; (*bpt)->isedge =0;
          do {
            if (n!=ondir && n!=offdir && (*bpt)->ngb[n]==0) {
              (*bpt)->isedge = n;
              (*bpt)->ngb[n] = temp->ngb[ondir];
            }
            n++;
          } while (n<2*G_ndim && (*bpt)->isedge==0);    
          (*bpt)->P[SI]=(*bpt)->Ph[SI]=-temp->ngb[ondir]->P[SI];
        }
        else rep.error("only know 1st/2nd order bcs",(*bpt)->id);
      }
#endif // GLM_NEGATIVE_BOUNDARY
      ct++;
      ++bpt;
    }
  } while (bpt !=b->data.end());

  if (ct != b->data.size()) {
    rep.error("BC_assign_OUTFLOW: missed some cells!",
              ct-b->data.size());
  }
  return 0;
}



// ##################################################################
// ##################################################################


int UniformGrid::BC_assign_ONEWAY_OUT(boundary_data *b)
{
  //
  // The setup for this is identical to outflow, so just call
  // outflow() setup function.
  //
  int err=BC_assign_OUTFLOW(b);
  return err;
}



// ##################################################################
// ##################################################################


int UniformGrid::BC_assign_INFLOW(    boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);
  if (b->data.empty()) rep.error("BC_assign_INFLOW: empty boundary data",b->itype);
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  do{
    temp = NextPt((*bpt),ondir); if(temp==0) rep.error("Got lost assigning Absorbing bcs.",temp->id);
    for (int i=0; i<BC_nbc; i++) {
      for (int v=0;v<G_nvar;v++) {
  (*bpt)->P[v] = (*bpt)->Ph[v] = temp->P[v];
  (*bpt)->dU[v] = 0.;
      }
      ct++;
      ++bpt;
     }
  } while (bpt !=b->data.end());
  if (ct != b->data.size()) rep.error("BC_assign_INFLOW: missed some cells!",ct-b->data.size());
  
  return 0;
}



// ##################################################################
// ##################################################################


int UniformGrid::BC_assign_REFLECTING(boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);
  if (b->data.empty()) rep.error("BC_assign_: empty boundary data",b->itype);
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, G_nvar);
    for (int v=0;v<G_nvar;v++) b->refval[v] = 1;
    switch (offdir) {
     case XN: case XP:
      b->refval[VX] = -1; break;
     case  YN: case YP:
      b->refval[VY] = -1; break;
     case  ZN: case ZP:
      b->refval[VZ] = -1; break;
     default:
      rep.error("BAD DIRECTION",offdir);
    } // Set Normal velocity direction.
    if (G_eqntype==EQMHD || G_eqntype==EQGLM || G_eqntype==EQFCD) {
      switch (offdir) {
       case XN: case XP:
  b->refval[BX] = -1; break;
       case  YN: case YP:
  b->refval[BY] = -1; break;
       case  ZN: case ZP:
  b->refval[BZ] = -1; break;
       default:
  rep.error("BAD DIRECTION",offdir);
      } // Set normal b-field direction.
    } // Setting up reference value.
  } // if we needed to set up refval.
  // Now go through each column of boundary points and assign values to them.
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  do{
    temp = NextPt((*bpt),ondir);
    if(!temp) rep.error("Got lost assigning reflecting bcs.",temp->id);
    for (int i=0; i<BC_nbc; i++) {
      for (int v=0;v<G_nvar;v++) {
  (*bpt)->P[v] = (*bpt)->Ph[v] = temp->P[v]*b->refval[v];
  (*bpt)->dU[v] = 0.;
      }
      (*bpt)->npt = temp;
      ++bpt; ct++; temp=NextPt(temp,ondir);
      if(!temp) rep.error("Got lost assigning reflecting bcs.",temp->id);
    }
  } while (bpt !=b->data.end());
  if (ct != b->data.size()) rep.error("BC_assign_: missed some cells!",ct-b->data.size());
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_FIXED(     boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);
  if (b->data.empty()) rep.error("BC_assign_FIXED: empty boundary data",b->itype);
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  temp = NextPt((*bpt),ondir);
  if(!temp) rep.error("Got lost assigning FIXED bcs.",temp->id);
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, G_nvar);
    for (int v=0;v<G_nvar;v++) b->refval[v] = temp->P[v];
  }
  // Initialise all the values to be the fixed value.
  do{
    for (int v=0;v<G_nvar;v++) {
      (*bpt)->P[v] = (*bpt)->Ph[v] = b->refval[v];
      (*bpt)->dU[v] = 0.;
    }
    ++bpt; ct++;
  } while (bpt !=b->data.end());
  if (ct != b->data.size())
    rep.error("BC_assign_FIXED: missed some cells!",ct-b->data.size());
  
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_JETBC(     boundary_data *b)
{
  if (!JP.jetic) rep.error("BC_assign_JETBC: not a jet simulation!",JP.jetic);
  if (b->dir != NO) rep.error("BC_assign_JETBC: boundary is not an internal one!",b->dir);
    cell *c = FirstPt(); cell *temp=0, *cy=0;
  int ct=0; int ctot=0; int maxnv=0;

  double jr = JP.jetradius*G_dx; // Physical radius of jet.
  //cout <<"jetrad="<<JP.jetradius<<" dx="<<G_dx<<"\n";
  //
  // Assign reference values, containing Jet parameters:
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, G_nvar);
  }
  if (G_eqntype==EQEUL || G_eqntype==EQEUL_ISO || G_eqntype==EQEUL_EINT ||
      G_eqntype==EQMHD || G_eqntype==EQGLM || G_eqntype==EQFCD) {
    rep.printVec("JetState",JP.jetstate,MAX_NVAR);
    b->refval[RO] = JP.jetstate[RO];
    b->refval[PG] = JP.jetstate[PG];
    b->refval[VX] = JP.jetstate[VX];
    b->refval[VY] = JP.jetstate[VY];
    b->refval[VZ] = JP.jetstate[VZ];
    maxnv=5;
  } else rep.error("BC_assign_JETBC: bad equation type",G_eqntype);
  if (G_eqntype==EQMHD || G_eqntype==EQGLM || G_eqntype==EQFCD) {
    b->refval[BX] = JP.jetstate[BX];
    b->refval[BY] = JP.jetstate[BY];
    b->refval[BZ] = JP.jetstate[BZ];
    maxnv=8;
  }
  if (G_eqntype==EQGLM) {
    b->refval[SI] = JP.jetstate[SI];
    maxnv=9;
  }
  for (int v=maxnv; v<G_nvar; v++) b->refval[v] = JP.jetstate[v];
  
  if (!b->data.empty()) rep.error("BC_assign_JETBC: boundary data exists!",b->itype);


  if (JP.jetic==1) { // Simplest jet -- zero opening angle.
    //
    // Axi-symmetry first -- this is relatively easy to set up.
    //
    if (G_ndim==2 && SimPM.coord_sys==COORD_CYL) {
      do {
  temp=c;
  while ( (temp=NextPt(temp,XN))!=0 ) {
    for (int v=0;v<G_nvar;v++) temp->P[v] = temp->Ph[v] = b->refval[v];
# ifdef SOFTJET
    //
    // jetradius is in number of cells, jr is in integer grid units (dx=2).
    // Jet centre is along Z_cyl axis, centred on origin.
    //
    temp->P[VX] = temp->Ph[VX] = b->refval[VX] *min(1., 4.0-4.0*CI.get_dpos(temp,YY)/jr);
    //cout <<"Incoming VX = "<< temp->P[VX]<<" pos="<<CI.get_dpos(temp,YY)<<" jr="<<jr<<"\n";
# endif //SOFTJET
    b->data.push_back(temp); ctot++;
    if (temp->isgd) rep.error("Looking for Boundary cells! setupjet",temp);
  }
  ct++;
      } while ( (c=NextPt(c,YP)) && ct<JP.jetradius);
      if (ct!=JP.jetradius) rep.error("Not enough cells for jet",ct);
      cout <<"Got "<<ctot<<" Cells in total for jet boundary.\n";
    } // 2D Axial Symmetry

    //
    // 3D now, more difficult.
    //
    else if (G_ndim==3 && SimPM.coord_sys==COORD_CRT) {
      double dist=0.;
      // 3D, so we need to convert the jet radius to a real length.
      // Also, the jet will come in at the centre of the XN boundary,
      // which must be the origin.
      do { // loop over ZZ axis
  cy=c; do { // loop over YY axis
    if ((dist=sqrt(CI.get_dpos(cy,YY)*CI.get_dpos(cy,YY) +CI.get_dpos(cy,ZZ)*CI.get_dpos(cy,ZZ))) <= jr) {
      temp = cy;
      while ( (temp=NextPt(temp,XN))!=0 ) {
        for (int v=0;v<G_nvar;v++) temp->P[v] = temp->Ph[v] = b->refval[v];
# ifdef SOFTJET
        temp->P[VX] = temp->Ph[VX] = b->refval[VX] *min(1., 4.0-4.0*dist/jr);
        //    cout <<"Incoming VX = "<< temp->P[VX]<<"\n";
# endif //SOFTJET
        b->data.push_back(temp); ctot++;
        if (temp->isgd) rep.error("Looking for Boundary cells! setupjet",temp);
      }
    } // if within jet radius
  } while ( (cy=NextPt(cy,YP)) ); // scroll through cells on YY axis.
      } while ( (c=NextPt(c,ZP)) );     // scroll through cells on ZZ axis.
      //      BC_printBCdata(b);
    } // 3D Cartesian

    else rep.error("Only know how to set up jet in 2Dcyl or 3Dcart",G_ndim);

  } // jetic==1
  else rep.error("Only know simple jet with jetic=1",JP.jetic);

  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_JETREFLECT(boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = OppDir(offdir);
  if (b->data.empty()) rep.error("BC_assign_: empty boundary data",b->itype);

  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, G_nvar);
    for (int v=0;v<G_nvar;v++) b->refval[v] = 1;
    switch (offdir) {
     case XN: case XP:
      b->refval[VX] = -1; break;
     case  YN: case YP:
      b->refval[VY] = -1; break;
     case  ZN: case ZP:
      b->refval[VZ] = -1; break;
     default:
      rep.error("BAD DIRECTION",offdir);
    } // Set Normal velocity direction.
    if (G_eqntype==EQMHD || G_eqntype==EQGLM || G_eqntype==EQFCD) {
      switch (offdir) {
       case XN: case XP:
  b->refval[BY] = b->refval[BZ] = -1; break;
       case  YN: case YP:
  b->refval[BX] = b->refval[BZ] = -1; break;
       case  ZN: case ZP:
  b->refval[BX] = b->refval[BY] = -1; break;
       default:
  rep.error("BAD DIRECTION",offdir);
      } // Set normal b-field direction.
    } // if B-field exists
  } // if we needed to set up refval.

  // Now go through each column of boundary points and assign values to them.
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  do{
    temp = NextPt((*bpt),ondir);
    if(!temp) rep.error("Got lost assigning jet-reflecting bcs.",temp->id);
    for (int i=0; i<BC_nbc; i++) {
      for (int v=0;v<G_nvar;v++) {
  (*bpt)->P[v] = (*bpt)->Ph[v] = temp->P[v]*b->refval[v];
  (*bpt)->dU[v] = 0.;
      }
      (*bpt)->npt = temp;
      ++bpt; ct++; temp=NextPt(temp,ondir);
      if(!temp) rep.error("Got lost assigning jet-reflecting bcs.",temp->id);
    }
  } while (bpt !=b->data.end());
  if (ct != b->data.size()) rep.error("BC_assign_: missed some cells!",ct-b->data.size());

  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_DMACH(     boundary_data *b)
{
#ifdef TESTING
  cout <<"Setting up DMACH boundary... starting.\n";
#endif // TESTING

  if (b->data.empty()) rep.error("BC_assign_DMACH: empty boundary data",b->itype);
  // Set reference value to be downstream values.
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, G_nvar);
    b->refval[RO] = 1.4; b->refval[PG] = 1.0; b->refval[VX] = b->refval[VY] = b->refval[VZ] = 0.0;
    //cout <<"SimPM.ftr="<<SimPM.ftr<<"\tG_nvar="<<G_nvar<<"\tSimPM.nvar="<<SimPM.nvar<<"\n";
    for (int v=SimPM.ftr; v<G_nvar; v++) b->refval[v] = -1.0;
  }
  // Run through all boundary cells, and give them either upstream or downstream
  // value, depending on their position.
  list<cell*>::iterator bpt=b->data.begin();
  unsigned int ct=0;
  do {
    if (CI.get_dpos(*bpt,XX) <= (10.0*SimPM.simtime/sin(M_PI/3.) +1./6. +CI.get_dpos(*bpt,YY)/tan(M_PI/3.))) {
      (*bpt)->P[RO] = (*bpt)->Ph[RO] = 8.0;
      (*bpt)->P[PG] = (*bpt)->Ph[PG] = 116.5;
      (*bpt)->P[VX] = (*bpt)->Ph[VX] = 7.14470958;
      (*bpt)->P[VY] = (*bpt)->Ph[VY] = -4.125;
      (*bpt)->P[VZ] = (*bpt)->Ph[VZ] = 0.0;
      for (int v=SimPM.ftr; v<G_nvar; v++) (*bpt)->P[v] = (*bpt)->Ph[v] = 1.0;
    }
    else for (int v=0;v<G_nvar;v++) (*bpt)->P[v] = (*bpt)->Ph[v] = b->refval[v];
    ++bpt; ct++;
  } while (bpt !=b->data.end());
  if (ct != b->data.size()) rep.error("BC_assign_: missed some cells!",ct-b->data.size());

#ifdef TESTING
  cout <<"Setting up DMACH boundary... finished.\n";
#endif // TESTING
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_DMACH2(    boundary_data *b) 
{
#ifdef TESTING
  cout <<"Setting up DMACH2 boundary... starting.\n";
#endif // TESTING

  if (b->dir != NO) rep.error("DMACH2 not internal boundary!",b->dir);
  cout <<"DMACH2 boundary, from x=0 to x=1/6 at y=0, fixed bd.\n";
  if (b->refval) rep.error("Already initialised memory in DMACH2 boundary refval",b->refval);

  b->refval = mem.myalloc(b->refval, G_nvar);

  b->refval[RO] = 8.0;
  b->refval[PG] = 116.5;
  b->refval[VX] = 7.14470958;
  b->refval[VY] = -4.125;
  b->refval[VZ] = 0.0;
  for (int v=SimPM.ftr; v<G_nvar; v++) b->refval[v] = 1.0;
  // Now have to go from first point onto boundary and across to x<=1/6
  if (!b->data.empty()) rep.error("BC_assign_DMACH2: Not empty boundary data",b->itype);
  cell *c = FirstPt();
  cell *temp=0;
  do {
    // check if we are <1/6 in case we are in a parallel grid where the domain is outside
    // the DMR region.  This saves having to rewrite the function in the parallel uniform grid.
    if (CI.get_dpos(c,XX)<=1./6.) {
      temp=c;
      while ( (temp=NextPt(temp,YN))!=0 ) {
        for (int v=0;v<G_nvar;v++) temp->P[v] = temp->Ph[v] = b->refval[v];
        b->data.push_back(temp);
        if (temp->isgd) rep.error("BC_assign_DMACH2: Looking for Boundary cells!",temp);
      }
    }
  } while ( (c=NextPt(c,XP)) && (CI.get_dpos(c,XX)<=1./6.));

#ifdef TESTING
  cout <<"Setting up DMACH2 boundary... finished.\n";
#endif // TESTING
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_RADSHOCK(  boundary_data *b)
{
  if (b->dir != NO) rep.error("RADSHOCK not internal boundary!",b->dir);
  cout <<"Assigning data to RADSHOCK boundary at XN boundary.\n";
  if (G_ndim !=2 && G_ndim !=1) rep.error("RADSHOCK must be 1d or 2d",G_ndim);
  if (b->refval) rep.error("Already initialised memory in RADSHOCK boundary refval",b->refval);
  b->refval = mem.myalloc(b->refval, G_nvar);
  if (!b->data.empty()) rep.error("BC_assign_RADSHOCK: Not empty boundary data",b->itype);
  cell *c=FirstPt();
  for (int v=0;v<G_nvar;v++) b->refval[v]=0.0;
  //
  // Try assuming that outflow will be at approx 3000 K, and shock
  // will have cooled (with a fudge factor for lower velocity).
  // so rho_1 = (v/(5km/s))^2 rho_0, and
  //    v_1   = v_0 /(v/(5km/s))^2
  //
  double ratio = c->P[VX]*c->P[VX]/25.0e10;
  if (fabs(c->P[VX])<55.0e5) ratio *= 0.5*fabs(c->P[VX])/50.0e5;
  b->refval[RO] = c->P[RO]*ratio; // This boundary requires density to be M^2 times initial density.
  b->refval[VX] = c->P[VX]/ratio; // mass flux out of domain is same as into domain.
  
  int ct=0;
  if (G_ndim==1) {
    do {
      c->P[RO] = c->Ph[RO] = min(c->P[RO], b->refval[RO]);
      b->data.push_back(c);
      ct++;
    } while (CI.get_dpos((c=NextPt(c)),XX) <= (G_xmin[XX] + (G_range[XX]*0.02)));
  }
  else if (G_ndim==2) {
    do {
      b->data.push_back(c);
      ct++;
      c->P[RO] = c->Ph[RO] = min(c->P[RO], b->refval[RO]);
    } while ( (NextPt(c,YP)!=0) && (c=NextPt(c,YP))->isgd );
    if (ct != G_ng[YY]) rep.error("Didn't get all cells in 2d RADSHOCK boundary",ct-G_ng[YY]);
  }
  cout <<"******************************************Added "<<ct<<" cells to RADSHOCK internal boundary.\n";
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_RADSH2(     boundary_data *b)
{
  if (b->dir != NO) rep.error("RADSH2 not external boundary!",b->dir);
  cout <<"Assigning data to RADSH2 boundary at XN boundary.\n";
  if (G_ndim !=2 && G_ndim !=1) rep.error("RADSH2 must be 1d or 2d",G_ndim);
  if (b->refval) 
    rep.error("Already initialised memory in RADSH2 boundary refval",b->refval);
  b->refval = mem.myalloc(b->refval, G_nvar);
  if (!b->data.empty())
    rep.error("BC_assign_RADSH2: Not empty boundary data",b->itype);
  cell *c=FirstPt();
  for (int v=0;v<G_nvar;v++) b->refval[v]=0.0;
  cell *temp =c;
  do {temp=NextPt(temp,XP);} while (NextPt(temp,XP)->isgd);
  cout <<"endpoint: "<<temp->P[VX]<<"  firstpt: "<<c->P[VX]<<"  ";
  b->refval[VX] = temp->P[VX]*temp->P[RO]/c->P[RO]; // This boundary requires fixed velocity.
  cout <<"Refval[VX] = "<< b->refval[VX] <<"\n"; //c->P[VX]<<"\n";
  
  int ct=0;
  if (G_ndim==1) {
    temp=c; 
    while (CI.get_dpos(NextPt(temp),XX) <= G_xmax[XX]/128.)
      temp=NextPt(temp);
    do {
      temp->P[VX] = temp->Ph[VX] = b->refval[VX]; // /c->P[VX]; //max(b->refval[VX],temp->P[VX]); // XN boundary, and want to limit outflow speed.
      b->data.push_back(temp);
      ct++;
    } while ( (temp=NextPt(temp,XN)) !=0);
  }
  else if (G_ndim==2) {
    do {
      while (CI.get_dpos(NextPt(c,XP),XX) <= (G_xmax[XX]-G_xmin[XX])/128.)
  c=NextPt(c,XP);
      temp=c;
      do {
  b->data.push_back(temp);
  ct++;
  temp->P[VX] = temp->Ph[VX] = b->refval[VX]; // /c->P[VX]; //max(b->refval[VX],temp->P[VX]);
      } while ( (temp = NextPt(temp,XN)) !=0);
    } while ( (NextPt(c,YP)!=0) && (c=NextPt(c,YP))->isgd );
//    cout <<"Got "<<ct<<" cells in 2d RADSH2 boundary.\n";
//    if (ct != G_ng[YY]) rep.error("Didn't get all cells in 2d RADSH2 boundary",ct-G_ng[YY]);
  }
  cout <<"******************************************Added "<<ct<<" cells to RADSH2 internal boundary.\n";
  return 0;
}



// ##################################################################
// ##################################################################



///
/// Add internal stellar wind boundaries -- these are (possibly
/// time-varying) winds defined by a mass-loss-rate and a terminal
/// velocity.  A region within the domain is given fixed values
/// corresponding to a freely expanding wind from a
/// cell-vertex-located source.
///
int UniformGrid::BC_assign_STWIND(boundary_data *b)
{
  //
  // Check that we have an internal boundary struct, and that we have a stellar
  // wind source to set up.
  //
  if (b->dir != NO)
    rep.error("STWIND not external boundary!",b->dir);
#ifdef TESTING
  cout <<"Assigning data to STWIND boundary. Nsrc="<<SWP.Nsources<<"\n";
#endif
  if (SWP.Nsources<1) {
    rep.error("UniformGrid::BC_assign_STWIND() No Wind Sources!",SWP.Nsources);
  }

  //
  // Setup reference state vector and initialise to zero.
  //
  if (b->refval) 
    rep.error("Already initialised memory in STWIND boundary refval",
        b->refval);
  b->refval = mem.myalloc(b->refval, G_nvar);
  if (!b->data.empty())
    rep.error("BC_assign_STWIND: Not empty boundary data",b->itype);
  for (int v=0;v<G_nvar;v++)
    b->refval[v]=0.0;

  //
  // New structure: we need to initialise the stellar wind class with
  // all of the wind sources in the global parameter list (this was
  // formerly done in DataIOBase::read_simulation_parameters()).
  //
  // The type of class we set up is determined first.  If any of the 
  // wind sources have type==3==evolving, we set up stellar_wind_evolving(),
  // otherwise we set up stellar_wind().
  //
  int err=0;
  int Ns = SWP.Nsources;
  for (int isw=0; isw<Ns; isw++) {
    if (SWP.params[isw]->type ==3) err+=1;
  }
  Wind = 0;
  if (Ns>0) {
//#ifdef TESTING
    cout <<"\n------------------ SETTING UP STELLAR WIND CLASS -----------------\n";
//#endif
    if (!err) {
      Wind = new stellar_wind ();
    }
    else {
      Wind = new stellar_wind_evolution();
      err=0;
    }
  }

  //
  // Run through sources and add sources.
  //
  for (int isw=0; isw<Ns; isw++) {
//#ifdef TESTING
    cout <<"\tUniformGrid::BC_assign_STWIND: Adding source "<<isw<<"\n";
//#endif
    if (SWP.params[isw]->type==3) {
      err = Wind->add_evolving_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->Rstar,
        SWP.params[isw]->tr,
        SWP.params[isw]->evolving_wind_file,
        SWP.params[isw]->time_offset,
        SimPM.simtime,
        SWP.params[isw]->update_freq,
        SWP.params[isw]->t_scalefactor
        );
    }
    else {
      err = Wind->add_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->Mdot,
        SWP.params[isw]->Vinf,
        SWP.params[isw]->Tstar,
        SWP.params[isw]->Rstar,
        SWP.params[isw]->tr
        );
    }
    if (err) rep.error("Error adding wind source",isw);
  }

  //
  // loop over sources, adding cells to boundary data list in order.
  //
  for (int id=0;id<Ns;id++) {
//#ifdef TESTING
    cout <<"\tUniformGrid::BC_assign_STWIND: Adding cells to source "<<id<<"\n";
//#endif
    BC_assign_STWIND_add_cells2src(id,b);
  }
  //
  // Now we should have set everything up, so we assign the boundary
  // cells with their boundary values.
  //
  err += BC_update_STWIND(b,0,0);
//#ifdef TESTING
  cout <<"\tFinished setting up wind parameters\n";
  cout <<"------------- DONE SETTING UP STELLAR WIND CLASS -----------------\n\n";
//#endif
  return err;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_STWIND_add_cells2src(
        const int id, ///< source id
        struct boundary_data *b
        )
{
  //
  // this is for cartesian geometry, with cubic cells, so things are
  // quite simple.  We run through each cell, and if it is within the 
  // source's radius of influence, then we add it to the lists.
  //
  int err=0;
  int ncell=0;
  double srcpos[MAX_DIM];
  Wind->get_src_posn(id,srcpos);
  double srcrad;
  Wind->get_src_drad(id,&srcrad);
#ifdef TESTING
  cout <<"*** srcrad="<<srcrad<<"\n";
  rep.printVec("src", srcpos, G_ndim);
#endif

  cell *c = FirstPt();
  do {
    //
    // GEOMETRY: This should be to centre--of--volume of cell!
    // It makes no difference for Cartesian grids b/c the centre--of--
    // volume coincides the midpoint.
    //
    if (distance_vertex2cell(srcpos,c) <= srcrad) {
        ncell++;
        //b->data.push_back(c); // don't need b to have a lit too.
        err += Wind->add_cell(this, id,c);
        //cout <<"CART: adding cell "<<c->id<<" to list. d=";
        //rep.printVec("src", srcpos, G_ndim);
        //rep.printVec("pos", c->pos, G_ndim);
      }
  } while ((c=NextPt(c))!=0);
  
  err += Wind->set_num_cells(id,ncell);
//#ifdef TESTING
  cout <<"UniformGrid: Added "<<ncell<<" cells to wind boundary for WS "<<id<<"\n";
//#endif
  return err;
}



// ##################################################################
// ##################################################################



///
/// Update internal stellar wind boundaries -- these are (possibly time-varying)
/// winds defined by a mass-loss-rate and a terminal velocity.  If fixed in time
/// the wind is updated with b->refval, otherwise with a (slower) call to the 
/// stellar wind class SW
///
int UniformGrid::BC_update_STWIND(
        boundary_data *b, ///< Boundary to update.
        const int ,  ///< current fractional step being taken.
        const int    ///< final step (not needed b/c fixed BC).
        )
{
  //
  // The stellar_wind class already has a list of cells to update
  // for each source, together with pre-calculated state vectors,
  // so we just call the set_cell_values() function.
  //
  int err=0;
  for (int id=0;id<Wind->Nsources();id++) {
    //cout <<" updating source "<<id<<"\n";
    err += Wind->set_cell_values(this, id,SimPM.simtime);
    //cout <<" finished source "<<id<<"\n";
  }

  return err;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_assign_STARBENCH1(boundary_data *b)
{
  //
  // First set up outflow boundaries.
  //
  int err=BC_assign_OUTFLOW(b);
  if (err) rep.error("BC_assign_STARBENCH1 error in OUTFLOW",err);

  //
  // Now do some checks, and then set the column density in the 
  // cells that we want to shadow.
  //
  if (b->dir != XN) rep.error("RADSH2 not XN boundary!",b->dir);
  //enum direction offdir = b->dir;
  //enum direction ondir  = OppDir(offdir);
  if (b->data.empty())
    rep.error("BC_assign_STARBENCH1: empty boundary data",b->itype);

  //
  // Set the column density (in g/cm2) to some large values in the
  // range 1.4pc<y<2.6pc, and to zero for the others.
  //
  list<cell*>::iterator bpt=b->data.begin();
  unsigned int ct=0;
  double dpos[G_ndim];
  double column = 1.67; // g/cm^2 so, in number this is 1.0e24.
  do{
    //
    // Set cell optical depth to be something big in the 
    // relevant region for the first radiation source.
    //
    CI.get_dpos((*bpt),dpos);
    if (dpos[YY]>4.3204e18 && dpos[YY]<8.0235e18 && dpos[XX]<SimPM.Xmin[XX]+G_dx) {
      column = 1.67; // g/cm^2 so, in number this is 1.0e24.
      for (int s=0; s<SimPM.RS.Nsources; s++)
        CI.set_col((*bpt), s, &column);
    }
    else {
      column = 0.00; // zero optical depth outside the barrier
      for (int s=0; s<SimPM.RS.Nsources; s++)
        CI.set_col((*bpt), s, &column);
    }

    
    ++bpt; ct++;
  } while (bpt !=b->data.end());
  if (ct != b->data.size())
    rep.error("BC_assign_STARBENCH1: missed some cells!",ct-b->data.size());
  
  return 0;
}

// ##################################################################
// ##################################################################

int UniformGrid::BC_update_STARBENCH1(
        struct boundary_data *b,
        const int cstep,
        const int maxstep
        )
{
  //
  // Outflow or Absorbing BCs; boundary cells are same as edge cells.
  // This is zeroth order outflow bcs.
  //
  // For the one-way boundary we set the normal velocity to be zero if
  // the flow wants to be onto the domain.
  //

  //
  // First get the normal velocity component, and whether the offdir is
  // positive or negative.
  //
  enum direction offdir = b->dir;
  enum primitive Vnorm = RO;
  int norm_sign = 0;
  switch (offdir) {
  case XN: case XP:
    Vnorm = VX; break;
  case YN: case YP:
    Vnorm = VY; break;
  case ZN: case ZP:
    Vnorm = VZ; break;
  default:
    rep.error("bad dir",offdir);
  }
  if (Vnorm == RO)
    rep.error("Failed to set normal velocity component",Vnorm);

  switch (offdir) {
  case XN: case YN: case ZN:
    norm_sign = -1; break;
  case XP: case YP: case ZP:
    norm_sign =  1; break;
  default:
    rep.error("bad dir",offdir);
  }

  //
  // Now run through all cells in the boundary
  //
  double dpos[G_ndim];
  double column = 1.67; // g/cm^2 so, in number this is 1.0e24.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //
    //exactly same routine as for periodic and zero-gradient.
    //
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v];
      (*c)->dU[v] = 0.;
    }
    //
    // ONEWAY_OUT: overwrite the normal velocity if it is inflow:
    //
    (*c)->Ph[Vnorm] = norm_sign*max(static_cast<pion_flt>(0.0),
                                    (*c)->Ph[Vnorm]*norm_sign);


    //
    // Set cell optical depth to be something big in the 
    // relevant region for the first radiation source.
    //
    CI.get_dpos((*c),dpos);
    if (dpos[YY]>4.3204e18 && dpos[YY]<8.0235e18 && dpos[XX]<SimPM.Xmin[XX]+G_dx) {
      column = 1.67; // g/cm^2 so, in number this is 1.0e24.
      for (int s=0; s<SimPM.RS.Nsources; s++)
        CI.set_col((*c), s, &column);
    }
    else {
      column = 0.00; // zero optical depth outside the barrier
      for (int s=0; s<SimPM.RS.Nsources; s++)
        CI.set_col((*c), s, &column);
    }


    if (cstep==maxstep)
      for (int v=0;v<G_nvar;v++)
        (*c)->P[v] = (*c)->npt->P[v];

#ifdef GLM_ZERO_BOUNDARY
    if (G_eqntype==EQGLM) {
      (*c)->P[SI]=(*c)->Ph[SI]=0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (G_eqntype==EQGLM) {
      if ((*c)->id==-1) {
        (*c)->P[SI] =-(*c)->npt->P[SI];
        (*c)->Ph[SI]=-(*c)->npt->Ph[SI];
      }
      else if ((*c)->id==-2) {
        (*c)->P[SI] =-(*c)->ngb[(*c)->isedge]->P[SI];
        (*c)->Ph[SI]=-(*c)->ngb[(*c)->isedge]->Ph[SI];
      }
      else rep.error("only know 1st/2nd order bcs",(*c)->id);
    }
#endif // GLM_NEGATIVE_BOUNDARY

  } // all cells.
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



int UniformGrid::TimeUpdateInternalBCs(const int cstep, const int maxstep)
{
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<BC_nbd;i++) {
    b = &BC_bd[i];
    switch (b->itype) {
     case RADSHOCK:   err += BC_update_RADSHOCK(   b, cstep, maxstep); break;
     case RADSH2:     err += BC_update_RADSH2(     b, cstep, maxstep); break;
     case STWIND:     err += BC_update_STWIND(     b, cstep, maxstep); break;
    case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
    case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2: case BCMPI:
    case STARBENCH1:
      //
      // External BCs updated elsewhere
      //     
      break;
      
    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.warning("Unhandled BC: serial update internal",b->itype,-1); err+=1; break;
      break;
    }
  }
  return(0);
}
  


// ##################################################################
// ##################################################################



int UniformGrid::TimeUpdateExternalBCs(const int cstep, const int maxstep)
{
  // TEMP_FIX
  //BC_printBCdata(&BC_bd[0]);
  // TEMP_FIX
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<BC_nbd;i++) {
    b = &BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
    case PERIODIC:   err += BC_update_PERIODIC(   b, cstep, maxstep); break;
    case OUTFLOW:    err += BC_update_OUTFLOW(    b, cstep, maxstep); break;
    case ONEWAY_OUT: err += BC_update_ONEWAY_OUT( b, cstep, maxstep); break;
    case INFLOW:     err += BC_update_INFLOW(     b, cstep, maxstep); break;
    case REFLECTING: err += BC_update_REFLECTING( b, cstep, maxstep); break;
    case FIXED:      err += BC_update_FIXED(      b, cstep, maxstep); break;
    case JETBC:      err += BC_update_JETBC(      b, cstep, maxstep); break;
    case JETREFLECT: err += BC_update_JETREFLECT( b, cstep, maxstep); break;
    case DMACH:      err += BC_update_DMACH(      b, cstep, maxstep); break;
    case DMACH2:     err += BC_update_DMACH2(     b, cstep, maxstep); break;
    case STARBENCH1: err += BC_update_STARBENCH1( b, cstep, maxstep); break;
    case RADSHOCK: case RADSH2: case STWIND: case BCMPI:
      //
      // internal bcs updated separately
      //
      break;
    default:
      //      cout <<"do i have a bc to update?? no.\n";
      rep.warning("Unhandled BC: serial update external",b->itype,-1); err+=1; break;
      break;
    }
  }
  return(0);
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_PERIODIC(   struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v];
      (*c)->dU[v] = 0.;
    }
    if (cstep==maxstep) {
      for (int v=0;v<G_nvar;v++) {
        (*c)->P[v] = (*c)->npt->P[v];
      }
    }
  } // all cells.
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_OUTFLOW(    struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  // Outflow or Absorbing BCs; boundary cells are same as edge cells.
  // This is zeroth order outflow bcs.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
#ifdef TESTING
    //    if ((*c)->npt->id==5887) {
    //cout <<"****************Boundary Cells********************\n";
    //cout <<"**** bc before:";
    //CI.print_cell((*c));
    //cout <<"**** gc before:";
    //CI.print_cell((*c)->npt);
    //}
#endif //TESTING
    //exactly same routine as for periodic.
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v];
      (*c)->dU[v] = 0.;
    }
    if (cstep==maxstep) for (int v=0;v<G_nvar;v++)  (*c)->P[v] = (*c)->npt->P[v];
#ifdef GLM_ZERO_BOUNDARY
    if (G_eqntype==EQGLM) {
      (*c)->P[SI]=(*c)->Ph[SI]=0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (G_eqntype==EQGLM) {
      //(*c)->P[SI] =-(*c)->npt->P[SI];
      //(*c)->Ph[SI]=-(*c)->npt->Ph[SI];
      //}
      if ((*c)->id==-1) {
  (*c)->P[SI] =-(*c)->npt->P[SI];
  (*c)->Ph[SI]=-(*c)->npt->Ph[SI];
      }
      else if ((*c)->id==-2) {
  //cout <<"bc x="<<(*c)->x[XX]/G_dx<<" and npt->x="<<(*c)->ngb[(*c)->isedge]->x[XX]/G_dx<<"\n";
  (*c)->P[SI] =-(*c)->ngb[(*c)->isedge]->P[SI];
  (*c)->Ph[SI]=-(*c)->ngb[(*c)->isedge]->Ph[SI];
      }
      else rep.error("only know 1st/2nd order bcs",(*c)->id);
    }
#endif // GLM_NEGATIVE_BOUNDARY

  } // all cells.
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_ONEWAY_OUT( struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  //
  // Outflow or Absorbing BCs; boundary cells are same as edge cells.
  // This is zeroth order outflow bcs.
  //
  // For the one-way boundary we set the normal velocity to be zero if
  // the flow wants to be onto the domain.
  //

  //
  // First get the normal velocity component, and whether the offdir is
  // positive or negative.
  //
  enum direction offdir = b->dir;
  enum primitive Vnorm = RO;
  int norm_sign = 0;
  switch (offdir) {
  case XN: case XP:
    Vnorm = VX; break;
  case YN: case YP:
    Vnorm = VY; break;
  case ZN: case ZP:
    Vnorm = VZ; break;
  default:
    rep.error("bad dir",offdir);
  }
  if (Vnorm == RO)
    rep.error("Failed to set normal velocity component",Vnorm);

  switch (offdir) {
  case XN: case YN: case ZN:
    norm_sign = -1; break;
  case XP: case YP: case ZP:
    norm_sign =  1; break;
  default:
    rep.error("bad dir",offdir);
  }

  //
  // Now run through all cells in the boundary
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //
    //exactly same routine as for periodic.
    //
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v];
      (*c)->dU[v] = 0.;
    }
    //
    // ONEWAY_OUT: overwrite the normal velocity if it is inflow:
    //
    (*c)->Ph[Vnorm] = norm_sign*max(static_cast<pion_flt>(0.0),
                                    (*c)->Ph[Vnorm]*norm_sign);

    if (cstep==maxstep)
      for (int v=0;v<G_nvar;v++)
  (*c)->P[v] = (*c)->npt->P[v];

#ifdef GLM_ZERO_BOUNDARY
    if (G_eqntype==EQGLM) {
      (*c)->P[SI]=(*c)->Ph[SI]=0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (G_eqntype==EQGLM) {
      //(*c)->P[SI] =-(*c)->npt->P[SI];
      //(*c)->Ph[SI]=-(*c)->npt->Ph[SI];
      //}
      if ((*c)->id==-1) {
  (*c)->P[SI] =-(*c)->npt->P[SI];
  (*c)->Ph[SI]=-(*c)->npt->Ph[SI];
      }
      else if ((*c)->id==-2) {
  //cout <<"bc x="<<(*c)->x[XX]/G_dx<<" and npt->x="<<(*c)->ngb[(*c)->isedge]->x[XX]/G_dx<<"\n";
  (*c)->P[SI] =-(*c)->ngb[(*c)->isedge]->P[SI];
  (*c)->Ph[SI]=-(*c)->ngb[(*c)->isedge]->Ph[SI];
      }
      else rep.error("only know 1st/2nd order bcs",(*c)->id);
    }
#endif // GLM_NEGATIVE_BOUNDARY

  } // all cells.
  return 0;
  
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_INFLOW(     struct boundary_data *b,
               const int , /// current ooa step e.g. 1 (not used here)
               const int   /// overall ooa      e.g. 2 (not used here)
               )
{
  // Inflow means BC is constant at the initial value of the neighbouring edge cell.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) (*c)->dU[v]=0.;
  } // all cells.
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_REFLECTING( struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  //exactly same routine as for outflow, except multiply v_n,B_n by -1.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //  cout <<"refval: "; for (int i=0;i<G_nvar;i++) cout <<b->refval[i]<<" "; cout <<"\n";
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v] *b->refval[v];
      (*c)->dU[v] = 0.;
    }
    if (cstep==maxstep) for (int v=0;v<G_nvar;v++)  (*c)->P[v] = (*c)->npt->P[v];
  } // all cells.
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_FIXED(      struct boundary_data *b,
               const int , /// current ooa step e.g. 1 (not used here)
               const int   /// overall ooa      e.g. 2 (not used here)
               )
{
  // Fixed means all boundary points have same fixed value, stored in refval.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) {
      (*c)->dU[v]=0.;
      (*c)->P[v] = (*c)->Ph[v] = b->refval[v];
    }
  }    // all cells.   
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_JETREFLECT( struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  //exactly same routine as for reflecting, except the normal B is unchanged, but the tangential is reversed.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //  cout <<"refval: "; for (int i=0;i<G_nvar;i++) cout <<b->refval[i]<<" "; cout <<"\n";
    for (int v=0;v<G_nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v] *b->refval[v];
      (*c)->dU[v] = 0.;
    }
    if (cstep==maxstep) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->npt->P[v];
    }
  } // all cells.
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_DMACH(      struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    if (CI.get_dpos(*c,XX) <= (10.0*SimPM.simtime/sin(M_PI/3.) +1./6. +CI.get_dpos(*c,YY)/tan(M_PI/3.))) {
      (*c)->Ph[RO] = 8.0;
      (*c)->Ph[PG] = 116.5;
      (*c)->Ph[VX] = 7.14470958;
      (*c)->Ph[VY] = -4.125;
      (*c)->Ph[VZ] = 0.0;
      //      cout <<"SimPM.ftr = "<<SimPM.ftr<<" and nvar = "<<G_nvar<<"\n";
      for (int v=SimPM.ftr; v<G_nvar; v++) (*c)->Ph[v] = 1.0;
    }
    else for (int v=0;v<G_nvar;v++) (*c)->Ph[v] = b->refval[v];
    for (int v=0;v<G_nvar;v++) (*c)->dU[v] = 0.0;
    if (cstep==maxstep) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    }
  }
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_DMACH2(     struct boundary_data *b,
               const int,
               const int
               )
{
  // Fixed at all times, so no difference between full and half step.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) {
      (*c)->dU[v]=0.;
      (*c)->P[v] = (*c)->Ph[v] = b->refval[v];
    }
  }    // all cells.   
  return 0;
}


// ##################################################################
// ##################################################################



int UniformGrid::BC_update_RADSHOCK(   struct boundary_data *b,
               const int cstep,
               const int maxstep
               )
{
  // We set the density adjacent to the wall to the minimum of 
  // the current density and some reference density.
  list<cell*>::iterator c=b->data.begin();
  //double temp=0.0;
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //temp = (*c)->Ph[RO];
    (*c)->Ph[RO] = min((*c)->Ph[RO], b->refval[RO]);
    (*c)->Ph[VX] = max((*c)->Ph[VX], b->refval[VX]);
    if (cstep==maxstep)
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    
  }
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_RADSH2(   struct boundary_data *b,
             const int cstep,
             const int maxstep
             )
{
  // We set the velocity at the boundary to be the initial outflow value.
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    (*c)->Ph[VX] = b->refval[VX]; //max(b->refval[VX],(*c)->Ph[VX]);
    if (cstep==maxstep)
      (*c)->P[VX] = (*c)->Ph[VX];
  }
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_update_JETBC(      struct boundary_data *b,
               const int,
               const int
               )
{
# ifdef SOFTJET
  double dist=0.0, jr=JP.jetradius*G_dx;
# endif //SOFTJET
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) {
      (*c)->dU[v]=0.;
      (*c)->P[v] = (*c)->Ph[v] = b->refval[v];
    }
# ifdef SOFTJET
    dist =0.0;
    if (G_ndim==2) dist = CI.get_dpos(*c,YY);
    else if (G_ndim==3) dist = sqrt(CI.get_dpos(*c,YY)*CI.get_dpos(*c,YY)+CI.get_dpos(*c,ZZ)*CI.get_dpos(*c,ZZ));
    else rep.error("Jet BC, but not 2d or 3d!!!",G_ndim);      
    (*c)->P[VX] = (*c)->Ph[VX] = b->refval[VX] *min(1., 4-4.0*dist/jr);
# endif //SOFTJET
  }    // all cells.   
  //  if (SimPM.timestep%10==0)
  //    BC_printBCdata(b);
  return 0;
}



// ##################################################################
// ##################################################################



void UniformGrid::BC_deleteBoundaryData()
{
#ifdef TESTING
  cout <<"BC destructor: deleting Boundary data...\n";
#endif
  struct boundary_data *b;
  for (int ibd=0; ibd<BC_nbd; ibd++) {
    b = &BC_bd[ibd];
    if (b->refval !=0) {
      b->refval = mem.myfree(b->refval);
    }
    if (ibd<2*G_ndim) {
      if (b->data.empty()) {
#ifdef TESTING
        cout <<"BC destructor: No boundary cells to delete.\n";
#endif
      }
      else {
        list<cell *>::iterator i=b->data.begin();
        do {
          CI.delete_cell(*i);  // This seems to work in terms of actually freeing the memory.
          b->data.erase(i);
          i=b->data.begin();
        }  while(i!=b->data.end());
        if (b->data.empty()) {
#ifdef TESTING
          cout <<"\t done.\n";
#endif
        }
        else {
#ifdef TESTING
          cout <<"\t not empty list! FIX ME!!!\n";
#endif
        }
      }
    }
    else {
      if (b->data.empty()) {
#ifdef TESTING
        cout <<"BC destructor: No boundary cells to delete.\n";
#endif
      }
      else {
        list<cell *>::iterator i=b->data.begin();
        do {
          // Don't need to delete the cell here because the extra boundaries just have
          // pointers to cells initialised and listed elsewhere.
          b->data.erase(i);
          i=b->data.begin();
        }  while(i!=b->data.end());
        if(b->data.empty()) {
#ifdef TESTING
          cout <<"\t done.\n";
#endif
        }
        else {
#ifdef TESTING
          cout <<"\t not empty list! FIX ME!!!\n";
#endif
        }
      }
    }
  } // loop over all boundaries.
  BC_bd = mem.myfree(BC_bd);
  return;
}
  


// ##################################################################
// ##################################################################



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in physical units, and outputs in 
/// physical units.
///
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



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in code integer units, and outputs in
/// integer units (but obviously the answer is not an integer).
///
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



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in physical units (e.g. centimetres).
///
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



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in grid--integer units (one cell has a diameter
/// two units).
///
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



///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are physical units.
///
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



///
/// As distance_vertex2cell(double[],cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
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



///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are code-integer units.
///
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



///
/// As idistance_vertex2cell(int,cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
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



///
/// As idifference_vertex2cell(int,cell,axis) but for the coordinate
/// difference between two cell positions along a given axis.
/// It returns *cell2* coordinate minus *cell1* coordinate.
///
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
uniform_grid_cyl::uniform_grid_cyl(
    int nd,         ///< ndim, length of position vector.
    int nv,         ///< nvar, length of state vectors.
    int eqt,        ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,        ///< Number of boundary cells to use.
    double *g_xn,   ///< array of minimum values of x,y,z for this grid.
    double *g_xp,   ///< array of maximum values of x,y,z for this grid.
    int *g_nc,      ///< array of number of cells in x,y,z directions.
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd)
{
#ifdef TESTING
  cout <<"Setting up cylindrical uniform grid with";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif
  if (G_ndim!=2)
    rep.error("Need to write code for !=2 dimensions",G_ndim);

#ifdef TESTING
  cout <<"cylindrical grid: dr="<<G_dx<<"\n";
#endif
  set_dx(G_dx);
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
  return (R_com(c)-G_xmin[Rcyl])/CI.phys_per_int() +G_ixmin[Rcyl];
}
  


// ##################################################################
// ##################################################################



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in physical units, and outputs in 
/// physical units.
///
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



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in code integer units, and outputs in
/// integer units (but obviously the answer is not an integer).
///
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



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in physical units (e.g. centimetres).
///
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
  temp = R_com(c1) - R_com(c2);
  d += temp*temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in grid--integer units (one cell has a diameter
/// two units).
///
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



///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are physical units.
///
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
  temp = v[Rcyl] - R_com(c);
  d += temp*temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



///
/// As distance_vertex2cell(double[],cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
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
    return R_com(c) - v[Rcyl];
  }
  else {
    cerr <<" Requested cylindrical distance in theta dir.\n";
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################


///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are code-integer units.
///
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


///
/// As idistance_vertex2cell(int,cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
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



///
/// As idifference_vertex2cell(int,cell,axis) but for the coordinate
/// difference between two cell positions along a given axis.
/// It returns *cell2* coordinate minus *cell1* coordinate.
///
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



int uniform_grid_cyl::BC_assign_STWIND_add_cells2src(
        const int id, ///< source id
        struct boundary_data *b
        )
{
  //
  // this is for cylindrical geometry, with cubic cells, so things are
  // quite simple.  We run through all cells, and if any are within the 
  // source's radius of influence, then we add them to the lists.
  //
  int err=0;
  int ncell=0;
  double srcpos[MAX_DIM];
  Wind->get_src_posn(id,srcpos);
  double srcrad;
  Wind->get_src_drad(id,&srcrad);
#ifdef TESTING
  cout <<"*** srcrad="<<srcrad<<"\n";
  rep.printVec("src", srcpos, G_ndim);
#endif

  cell *c = FirstPt();
  do {
    //
    // GEOMETRY: This is to centre--of--volume of cell!
    //
    //cout <<"CYL: "<<idistance_vertex2cell(srcpos,c)<<"\n";
    //rep.printVec("***pos", c->pos, G_ndim);
    //rep.printVec("***src", srcpos, G_ndim);
    //
    if (distance_vertex2cell(srcpos,c) <= srcrad) {

      ncell++;
      //b->data.push_back(c); // don't need b to have a lit too.
      err += Wind->add_cell(this, id,c);
#ifdef TESTING
      cout <<"CYL adding cell "<<c->id<<" to list.\n";
      rep.printVec("***src", srcpos, G_ndim);
      rep.printVec("***pos", c->pos, G_ndim);
      cout <<"*** distance="<<distance_vertex2cell(srcpos,c)<<"\n";
#endif
      //rep.printVec("src", srcpos, G_ndim);
      //rep.printVec("pos", c->pos, G_ndim);
    }
  } while ((c=NextPt(c))!=0);
  
  err += Wind->set_num_cells(id,ncell);

  return err;
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
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd)
{
#ifdef TESTING
  cout <<"Setting up spherical uniform grid with";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif
  if (G_ndim!=1)
    rep.error("Need to write code for >1 dimension",G_ndim);

#ifdef TESTING
  cout <<"spherical grid: dr="<<G_dx<<"\n";
#endif
  set_dx(G_dx);
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
  //cout <<"R_com(c)="<<R_com(c)<<", ppi="<<CI.phys_per_int()<<", ixmin="<<G_ixmin[Rsph]<<"\n";
  //cout <<"So iR_cov(c)="<<(R_com(c)-G_xmin[Rsph])/CI.phys_per_int() +G_ixmin[Rsph]<<"\n";
  return ((R_com(c)-G_xmin[Rsph])/CI.phys_per_int() +G_ixmin[Rsph]);
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in physical units, and outputs in 
/// physical units.
///
double uniform_grid_sph::distance(const double *p1, ///< position 1 (physical)
          const double *p2  ///< position 2 (physical)
          )
{
  return fabs(p1[Rsph]-p2[Rsph]);
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between two points, where the two position
/// are interpreted in the appropriate geometry.
/// This function takes input in code integer units, and outputs in
/// integer units (but obviously the answer is not an integer).
///
double uniform_grid_sph::idistance(const int *p1, ///< position 1 (integer)
           const int *p2  ///< position 2 (integer)
           )
{
  return fabs(static_cast<double>(p1[Rsph]-p2[Rsph]));
}
   


// ##################################################################
// ##################################################################



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in physical units (e.g. centimetres).
///
double uniform_grid_sph::distance_cell2cell(const cell *c1, ///< cell 1
              const cell *c2  ///< cell 2
              )
{
  return fabs(R_com(c1)-R_com(c2));
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between two cell--centres (will be between
/// centre-of-volume of cells if non-cartesian geometry).
/// Result returned in grid--integer units (one cell has a diameter
/// two units).
///
double uniform_grid_sph::idistance_cell2cell(const cell *c1, ///< cell 1
               const cell *c2  ///< cell 2
               )
{
  return fabs(R_com(c1)-R_com(c2))/CI.phys_per_int();
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are physical units.
///
double uniform_grid_sph::distance_vertex2cell(const double *v, ///< vertex (physical)
                const cell *c    ///< cell
                )
{
  return fabs(v[Rsph] - R_com(c));
}




// ##################################################################
// ##################################################################



///
/// As distance_vertex2cell(double[],cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
double uniform_grid_sph::difference_vertex2cell(
      const double *v,  ///< vertex (double)
      const cell *c, ///< cell
      const axes a   ///< Axis to calculate.
      )
{
  if      (a==Rsph) {
    return R_com(c) - v[Rsph];
  }
  else {
    cerr <<" Requested spherical distance in not-radial dir.\n";
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################



///
/// Calculate distance between a cell-vertex and a cell--centres
/// (will be between centre-of-volume of cells if non-cartesian
/// geometry).  Here both input and output are code-integer units.
///
double uniform_grid_sph::idistance_vertex2cell(const int *v, ///< vertex (integer)
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



///
/// As idistance_vertex2cell(int,cell) but for a single component
/// of the position vector, and not the absolute value.  It returns
/// the *cell* coordinate minus the *vertex* coordinate.
///
double uniform_grid_sph::idifference_vertex2cell(const int *v,  ///< vertex (integer)
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



///
/// As idifference_vertex2cell(int,cell,axis) but for the coordinate
/// difference between two cell positions along a given axis.
/// It returns *cell2* coordinate minus *cell1* coordinate.
///
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



int uniform_grid_sph::BC_assign_STWIND_add_cells2src(
        const int id, ///< source id
        struct boundary_data *b
        )
{
  //
  // this is for spherical geometry in 1D, with not necessarily equal
  // radius cells, so things are not as simple as before.  We run
  // through all cells, and if any are within the source's radius of
  // influence, then we add them to the lists.
  //
  if (G_ndim!=1) rep.error("Code winds for 2D-Sph!",G_ndim);
  //
  // We use physical locations and not integer units.
  //
  int err=0;
  int ncell=0;
  double srcpos[MAX_DIM];
  Wind->get_src_posn(id,srcpos);
  double srcrad;
  Wind->get_src_drad(id,&srcrad);
#ifdef TESTING
  cout <<"***SPH-1D srcrad="<<srcrad<<"\n";
  rep.printVec("***src", srcpos, G_ndim);
#endif

  //
  // For spherical geometry we need to add boundary cells, so we go
  // back to the boundary data before starting our 1D traverse across
  // the grid.
  //
  cell *c = FirstPt();
  while (NextPt(c,RNsph)!=0) c=NextPt(c,RNsph);

  do {
    //
    // GEOMETRY: This is to centre--of--volume of cell!
    //
    if (distance_vertex2cell(srcpos,c) <= srcrad) {
      ncell++;
#ifdef TESTING
      cout <<"grid-Sph1D: adding cell "<<c->id<<" to list.\n";
      rep.printVec("cell-pos", c->pos, G_ndim);
#endif
      //b->data.push_back(c); // don't need b to have a list too.
      err += Wind->add_cell(this, id,c);
      //cout <<"SPH: adding cell "<<c->id<<" to list. d=";
      //rep.printVec("src", srcpos, G_ndim);
      //rep.printVec("pos", c->pos, G_ndim);
    }
  } while ( (c=NextPt(c,RPsph))!=0);
  
  err += Wind->set_num_cells(id,ncell);

  return err;
}



// ##################################################################
// ##################################################################



//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------

