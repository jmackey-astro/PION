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

#include "constants.h"
#include "grid/uniform_grid.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include <fstream>
using namespace std;

//#define TEST_BC89FLUX
//#define BC_DEBUG

// ##################################################################
// ##################################################################

UniformGrid::UniformGrid(
    int nd,
    int nv,
    int eqt,
    int Nbc,  ///< Number of boundary cells to use.
    const std::array<double, MAX_DIM> &g_xn,  // this grid xmin
    const std::array<double, MAX_DIM> &g_xp,  // this grid xmax
    const std::array<int, MAX_DIM> &g_nc,
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM> &sim_xn,  // sim xmin
    std::array<double, MAX_DIM> &sim_xp   // sim xmax
    ) :
    VectorOps_Cart(nd),
    G_ndim(nd), G_nvar(nv), G_eqntype(eqt)  //,
// BC_nbc(Nbc)
{

  spdlog::debug(
      "Setting up UniformGrid with G_ndim={} and G_nvar={}", G_ndim, G_nvar);

#ifndef NDEBUG
  spdlog::debug("\tXmin : {}", g_xn);
  spdlog::debug("\tXmax : {}", g_xp);
  spdlog::debug("\tNpt  : {}", g_nc);
#endif
  G_coordsys = COORD_CRT;  // Cartesian Coordinate system

  //
  // Allocate arrays for dimensions of grid.
  //
  G_ng.fill(1);
  G_ng_all.fill(1);

  //
  // Assign values to grid dimensions and set the cell-size
  //
  G_ncell     = 1;
  G_ncell_all = 1;

  for (int i = 0; i < G_ndim; i++) {
    G_ng[i] = g_nc[i];
    G_ncell *= G_ng[i];  // Get total number of cells.
    G_xmin[i]  = g_xn[i];
    G_xmax[i]  = g_xp[i];
    G_range[i] = G_xmax[i] - G_xmin[i];
    if (G_range[i] < 0.)
      spdlog::error("{}: {}", "Negative range in direction i", i);

    //
    // For single static grid, Sim_xmin/xmax/range are the same as
    // for the grid (becuase there is only one grid).  But we set
    // them here anyway so that it is correctly placed in a multigrid
    // setup.
    //
    Sim_xmin[i]  = sim_xn[i];
    Sim_xmax[i]  = sim_xp[i];
    Sim_range[i] = sim_xp[i] - sim_xn[i];

    L_xmin[i]  = lev_xn[i];
    L_xmax[i]  = lev_xp[i];
    L_range[i] = lev_xp[i] - lev_xn[i];

    // if grid boundary is a level boundary, then add BC_nbc cells,
    // otherwise (i.e. MPI domain boundary) add 4.
    if (!pconst.equalD(G_xmin[i], L_xmin[i]))
      G_nbc[2 * i] = min(4, Nbc);
    else
      G_nbc[2 * i] = Nbc;
    if (!pconst.equalD(G_xmax[i], L_xmax[i]))
      G_nbc[2 * i + 1] = min(4, Nbc);
    else
      G_nbc[2 * i + 1] = Nbc;
    // all cells including boundary.
    G_ng_all[i] = g_nc[i] + G_nbc[2 * i] + G_nbc[2 * i + 1];
    G_ncell_all *= G_ng[i] + G_nbc[2 * i] + G_nbc[2 * i + 1];
  }
  for (int i = G_ndim; i < MAX_DIM; i++) {
    G_nbc[2 * i]     = 0;
    G_nbc[2 * i + 1] = 0;
  }

  spdlog::debug("MIN.MAX for x = {}\t{}", G_xmin[XX], G_xmax[XX]);
  spdlog::debug("Setting cell size..");

  //
  // Checks grid dimensions and discretisation is reasonable,
  // and that it gives cells with equal length in each dimension.
  //
  set_cell_size();

  // grid dimensions including boundary data
  for (int v = 0; v < G_ndim; v++) {
    G_xmin_all[v] = G_xmin[v] - G_nbc[2 * v] * G_dx;
  }
  for (int v = 0; v < G_ndim; v++) {
    G_xmax_all[v] = G_xmax[v] + G_nbc[2 * v + 1] * G_dx;
  }
  for (int v = 0; v < G_ndim; v++) {
    G_range_all[v] = G_xmax_all[v] - G_xmin_all[v];
  }

  // HERE WE SHOULD CALL A FUNCTION TO ALLOCATE A GRID OF CELLS, AND
  // ALLOCATE A BIG BLOCK OF MEMORY FOR ALL THE ARRAY POINTERS OF
  // EACH CELL.  PROBABLY I CAN RE-WRITE allocate_grid_data().
#ifndef NEWGRIDDATA
  //
  // Now create the first cell, and then allocate data from there.
  // Safe to assume we have at least one cell...
  //
  spdlog::debug("... done. Initialising first cell...");

  G_fpt_all     = CI.new_cell();
  G_fpt_all->id = 0;

  spdlog::debug(" done");

  if (G_fpt_all == 0) {
    spdlog::error(
        "{}: {}", "Couldn't assign memory to first cell in grid.", G_fpt_all);
  }
#endif  // NEWGRIDDATA

  // allocate memory for all cells, including boundary cells.

  spdlog::debug("allocating memory for grid");

  int err = allocate_grid_data();
  if (err != 0)
    spdlog::error("{}: {}", "Error setting up grid, allocate_grid_data", err);

  // assign grid structure on cells, setting positions and ngb
  // pointers.
  spdlog::debug("assigning pointers to neighbours");

  err += assign_grid_structure();
  if (err != 0)
    spdlog::error(
        "{}: {}", "Error setting up grid, assign_grid_structure", err);

#ifndef NDEBUG
  spdlog::debug(
      "\tFirst Pt. integer position : {}",
      std::vector<int>(FirstPt()->pos, FirstPt()->pos + nd));
  spdlog::debug(
      "\tLast  Pt. integer position : {}",
      std::vector<int>(LastPt()->pos, LastPt()->pos + nd));
#endif  // NDEBUG

  //
  // Leave boundaries uninitialised.
  //
  BC_nbd = -1;

  //
  // Set integer dimensions/location of grid.
  //
  CI.get_ipos_vec(G_xmin, G_ixmin);
  CI.get_ipos_vec(G_xmax, G_ixmax);
  for (int v = 0; v < G_ndim; v++)
    G_irange[v] = G_ixmax[v] - G_ixmin[v];
  G_idx = G_irange[XX] / G_ng[XX];
  for (int v = 0; v < G_ndim; v++)
    G_ixmin_all[v] = G_ixmin[v] - G_nbc[2 * v] * G_idx;
  for (int v = 0; v < G_ndim; v++)
    G_ixmax_all[v] = G_ixmax[v] + G_nbc[2 * v + 1] * G_idx;
  for (int v = 0; v < G_ndim; v++)
    G_irange_all[v] = G_ixmax_all[v] - G_ixmin_all[v];

  CI.get_ipos_vec(Sim_xmin, Sim_ixmin);
  CI.get_ipos_vec(Sim_xmax, Sim_ixmax);
  for (int v = 0; v < G_ndim; v++) {
    Sim_irange[v] = Sim_ixmax[v] - Sim_ixmin[v];
  }
  CI.get_ipos_vec(L_xmin, L_ixmin);
  CI.get_ipos_vec(L_xmax, L_ixmax);
  for (int v = 0; v < G_ndim; v++) {
    L_irange[v] = L_ixmax[v] - L_ixmin[v];
  }

  // set external boundary cells to be not part of the simulation
  // domain if they are outside of Sim_Xmin and Sim_Xmax.
  class cell *c = FirstPt_All();
  bool dom      = true;
  do {
    c->isdomain = true;
    dom         = true;
    for (int v = 0; v < G_ndim; v++) {
      if (c->pos[v] < Sim_ixmin[v] || c->pos[v] > Sim_ixmax[v]) dom = false;
    }
    if (!dom) c->isdomain = false;
    c->isleaf = true;  // assume all cells are leaves unless changed.
    if (dom) {
      c->timestep = true;
    }
  } while ((c = NextPt_All(c)) != 0);

#ifndef NDEBUG
  spdlog::debug("grid ixmin  : {}", G_ixmin);
  spdlog::debug("grid ixmax  : {}", G_ixmax);
  spdlog::debug("grid irange : {}", G_irange);
  spdlog::debug("grid Ncell : {}", G_ng);
  spdlog::debug("grid ixmin_all  : {}", G_ixmin_all);
  spdlog::debug("grid ixmax_all  : {}", G_ixmax_all);
  spdlog::debug("grid irange_all : {}", G_irange_all);
  spdlog::debug("grid Ncell_all : {}", G_ng_all);
  spdlog::debug("grid xmin  : {}", G_xmin);
  spdlog::debug("grid xmax  : {}", G_xmax);
  spdlog::debug("grid range : {}", G_range);
  spdlog::debug("grid xmin_all  : {}", G_xmin_all);
  spdlog::debug("grid xmax_all  : {}", G_xmax_all);
  spdlog::debug("grid range_all : {}", G_range_all);
  spdlog::debug("Sim xmin  : {}", Sim_xmin);
  spdlog::debug("Sim xmax  : {}", Sim_xmax);
  spdlog::debug("Sim range : {}", Sim_range);
  spdlog::debug("boundary depth : {}", G_nbc);
#endif

  spdlog::debug("Cartesian grid: dr={}", G_dx);
  RT = 0;

  spdlog::debug("UniformGrid Constructor done");
}  // UniformGrid Constructor



// ##################################################################
// ##################################################################



UniformGrid::~UniformGrid()
{
  //
  // Delete the grid data.
  //
#ifdef NEWGRIDDATA

  // deallocate grid data
  *griddata = mem.myfree(*griddata);
  griddata  = mem.myfree(griddata);

  for (size_t i = 0; i < G_ncell_all; i++) {
    (*gridcells)[i].ngb      = mem.myfree((*gridcells)[i].ngb);
    (*gridcells)[i].pos      = mem.myfree((*gridcells)[i].pos);
    (*gridcells)[i].isbd_ref = mem.myfree((*gridcells)[i].isbd_ref);
    if (!(*gridcells)[i].F.empty()) {
      for (int id = 0; id < G_ndim; id++) {
        if ((*gridcells)[i].F[id]) {
          (*gridcells)[i].F[id] = mem.myfree((*gridcells)[i].F[id]);
        }
      }
    }
  }
  gridcells = mem.myfree(gridcells);

#else

  cell *cpt = FirstPt_All();
  cell *npt = NextPt_All(cpt);
  do {
    CI.delete_cell(cpt);
    cpt = npt;
  } while ((npt = NextPt_All(cpt)) != 0);
  CI.delete_cell(cpt);

#endif  // NEWGRIDDATA

  BC_deleteBoundaryData();

  spdlog::debug("UniformGrid Destructor:\tdone");
}  // Destructor



// ##################################################################
// ##################################################################



int UniformGrid::allocate_grid_data()
{
  spdlog::debug("Allocating grid data... G_ncell={}", G_ncell);
#ifdef NEWGRIDDATA
  // allocate grid-data
  // * find out how many doubles are needed for each cell
  size_t nel = static_cast<size_t>(CI.get_Nel());
  // * multiply by number of cells to get total array size
  nel *= G_ncell_all;

  // set stride (number of bytes per cell in big array)
  UniformGrid::gdata_stride = nel;

  // * allocate data for big array
  griddata  = mem.myalloc(griddata, 1);
  *griddata = mem.myalloc(*griddata, nel);
  // * allocate array of cells
  gridcells    = mem.myalloc(gridcells, 1);
  *gridcells   = mem.myalloc(*gridcells, G_ncell_all);
  cell *carr   = *gridcells;
  double *data = *griddata;

  // * initialise first and last pointers
  G_fpt_all = &(carr[0]);
  G_lpt_all = &(carr[G_ncell_all - 1]);

  // * loop over cells, add pointers to elements in big data array
  // for each cell, and add an npt_all pointer to the next cell so
  // that the cells are a linked list.
  size_t ix = 0;
  for (size_t i = 0; i < G_ncell_all; i++) {
#ifndef NDEBUG
    // spdlog::debug("i={} of {}: ix={}, nel={}", i, G_ncell_all, ix, nel);
#endif
    carr[i].id = static_cast<long int>(i);
    ix         = CI.set_cell_pointers(&(carr[i]), data, ix);
    if (i < G_ncell_all - 1)
      carr[i].npt_all = &(carr[i + 1]);
    else
      carr[i].npt_all = 0;
  }
#else
  // add a npt_all pointer so cells are in a singly linked list.
  cell *c      = G_fpt_all;
  size_t count = 0;
  while (c->id < static_cast<long int>(G_ncell_all - 1)) {
    c->npt_all = CI.new_cell();
    c          = c->npt_all;
    count++;
    c->id = count;
  }
  c->npt_all = 0;
  G_lpt_all  = c;
#endif  // NEWGRIDDATA
  spdlog::debug("Finished Allocating Data");
  return 0;
}  // allocate_grid_data

// ##################################################################
// ##################################################################



// ##################################################################
// ##################################################################



int UniformGrid::assign_grid_structure()
{
  spdlog::debug("AssignGridStructure");

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
  // int ipos[MAX_DIM], offset[MAX_DIM];
  // for (int i=0; i<MAX_DIM; i++) ipos[i] = offset[i] = 0;
  // for (int i=0; i<G_ndim; i++) {
  //  offset[i] = 1+
  //  2*static_cast<int>(ONE_PLUS_EPS*(G_xmin[i]-Sim_xmin[i])/G_dx);
  //}

  // ----------------------------------------------------------------
  // ---------------------- SET CELL POSITIONS ----------------------
  // Go through full grid, setting positions, where the first cell is
  // BC_nbc cells negative of G_xmin in all directions.
  //
  int ix[MAX_DIM];  // ix[] is a counter for cell number.
  for (int i = 0; i < MAX_DIM; i++)
    ix[i] = -G_nbc[2 * i];

  std::array<double, MAX_DIM> dpos;
  for (int i = 0; i < MAX_DIM; i++)
    dpos[i] = 0.0;

  class cell *c = FirstPt_All();
  do {
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
    for (int i = 0; i < G_ndim; i++)
      dpos[i] = G_xmin[i] + G_dx * (ix[i] + 0.5);
    CI.set_pos(c, dpos);
#ifndef NDEBUG
    // rep.printVec("    pos", c->pos, G_ndim);
#endif  // NDEBUG

    //
    // Initialise the cell data to zero.
    //
    for (int v = 0; v < G_nvar; v++)
      c->P[v] = 0.0;
    if (!CI.query_minimal_cells()) {
      for (int v = 0; v < G_nvar; v++)
        c->Ph[v] = c->dU[v] = 0.;
    }

    //
    // See if cell is on grid or not, and set isgd/isbd accordingly.
    //
    CI.get_dpos(c, dpos);
    bool on_grid = true;
    for (int v = 0; v < G_ndim; v++)
      if (dpos[v] < G_xmin[v] || dpos[v] > G_xmax[v]) on_grid = false;
#ifndef NDEBUG
        // rep.printVec("    dpos",dpos,G_ndim);
        // rep.printVec("    xmax",G_xmax,G_ndim);
#endif  // NDEBUG

    if (on_grid) {
      c->isgd = true;
      c->isbd = false;
    }
    else {
      c->isgd = false;
      c->isbd = true;
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
    c->isedge = 0;
    if (ix[XX] == 0)
      c->isedge += 1;
    else if (ix[XX] == G_ng[XX] - 1)
      c->isedge += 2;
    if (G_ndim > 1) {
      if (ix[YY] == 0)
        c->isedge += 1 * 3;
      else if (ix[YY] == G_ng[YY] - 1)
        c->isedge += 2 * 3;
    }
    if (G_ndim > 2) {
      if (ix[ZZ] == 0)
        c->isedge += 1 * 3 * 3;
      else if (ix[ZZ] == G_ng[ZZ] - 1)
        c->isedge += 2 * 3 * 3;
    }

    //
    // Boundary values for isedge:  if the cell is off-grid, then set
    // isedge to equal the number of cells it is from the grid.
    //
    if (ix[XX] < 0) {
      c->isedge = ix[XX];
    }
    else if (ix[XX] >= G_ng[XX]) {
      c->isedge = G_ng[XX] - 1 - ix[XX];
    }
    if (G_ndim > 1) {
      if (ix[YY] < 0) {
        c->isedge = ix[YY];
      }
      else if (ix[YY] >= G_ng[YY]) {
        c->isedge = G_ng[YY] - 1 - ix[YY];
      }
    }
    if (G_ndim > 2) {
      if (ix[ZZ] < 0)
        c->isedge = ix[ZZ];
      else if (ix[ZZ] >= G_ng[ZZ])
        c->isedge = G_ng[ZZ] - 1 - ix[ZZ];
    }

    //
    // Increment counters
    //
    ix[XX]++;
    //
    // if we got to the end of a row, then loop back and increment iy
    //
    if (ix[XX] == G_ng[XX] + G_nbc[XP]) {
      ix[XX] = -G_nbc[XN];
      if (G_ndim > 1) {
        ix[YY]++;
        //
        // if we got to the end of a y-column, then loop back and
        // increment z.
        //
        if (ix[YY] == G_ng[YY] + G_nbc[YP]) {
          ix[YY] = -G_nbc[YN];
          if (G_ndim > 2) {
            ix[ZZ]++;
            if (ix[ZZ] > G_ng[ZZ] + G_nbc[ZP]) {
              spdlog::error("\tWe should be done by now");
              return (ix[ZZ]);
            }
          }  // If 3D.
        }    // If at end of YY row.
      }      // If 2D.
    }        // If at end of XX row.
  } while ((c = NextPt_All(c)) != 0);
  // ---------------------- SET CELL POSITIONS ----------------------
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // ----------------------  SET npt POINTERS  ----------------------
  // Now run through grid and set npt pointer for on-grid cells.
  // Also set fpt for the first on-grid point, and lpt for the last.
  //
  bool set_fpt = false, set_lpt = false;
  cell *ctemp = 0;
  c           = FirstPt_All();
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
        G_fpt   = c;
        set_fpt = true;
      }
      //
      // if next cell in list is also on-grid, this is npt.
      //
      if (NextPt_All(c) != 0 && NextPt_All(c)->isgd) {
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
        } while (ctemp != 0 && !ctemp->isgd);

        if (ctemp) {
          // there is another on-grid point.
          if (!ctemp->isgd)
            spdlog::error("{}: {}", "Setting npt error", ctemp->id);
          c->npt = ctemp;
        }
        else {
          // must be last point on grid.
          c->npt  = 0;
          G_lpt   = c;
          set_lpt = true;
        }
      }
    }  // if c->isgd
  } while ((c = NextPt_All(c)) != 0);
  if (!set_lpt) spdlog::error("{}: {}", "failed to find last on-grid point", 1);
  if (!set_fpt)
    spdlog::error("{}: {}", "failed to find first on-grid point", 1);
  // ----------------------  SET npt POINTERS  ----------------------
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // ----------------------  SET ngb POINTERS  ----------------------
  //
  // Set up pointers to neighbours in x-direction
  //
  c            = FirstPt_All();
  cell *c_prev = 0, *c_next = c;
  do {
    c_next     = NextPt_All(c_next);
    c->ngb[XN] = c_prev;

    if (!c_next) {
      c->ngb[XP] = 0;
    }
    else if (c_next->pos[XX] > c->pos[XX]) {
      c->ngb[XP] = c_next;
      c_prev     = c;
    }
    else {
      c->ngb[XP] = 0;
      c_prev     = 0;
    }

    c = c_next;
  } while (c != 0);

  //
  // Pointers to neighbours in the Y-Direction, if it exists.
  // This is not a very smart algorithm and is probably quite slow,
  // but it does the job.
  //
  if (G_ndim > 1) {
    c = FirstPt_All();
    do {
      c_next = c;
      do {
        c_next = NextPt_All(c_next);
      } while ((c_next != 0) && (c_next->pos[YY] >= c->pos[YY])
               && (c_next->pos[XX] != c->pos[XX]));
      //
      // So now maybe c_next is zero, in which case ngb[YP]=0
      //
      if (!c_next) {
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
        c->ngb[YP]      = c_next;
        c_next->ngb[YN] = c;
      }
    } while ((c = NextPt_All(c)) != 0);
  }  // if G_ndim>1

  //
  // Pointers to neighbours in the Z-Direction, if it exists.
  //
  if (G_ndim > 2) {
    //
    // First get two cells, one above the other:
    //
    c      = FirstPt_All();
    c_next = NextPt_All(c);
    // CI.print_cell(c);
    // CI.print_cell(c_next);
    while (c_next->pos[ZZ] == c->pos[ZZ]) {
      c_next = NextPt_All(c_next);
      // Let c_next loop through x-y plane until i get to above the
      // first point.
    }
    //
    // Now c_next should be c's ZP neighbour.
    //
    do {
      // rep.printVec("C-",c->pos,G_ndim);
      // rep.printVec("C+",c_next->pos,G_ndim);
      c->ngb[ZP]      = c_next;
      c_next->ngb[ZN] = c;
      c               = NextPt_All(c);
      c_next          = NextPt_All(c_next);
    } while (c_next != 0);
  }
  // ----------------------  SET ngb POINTERS  ----------------------
  // ----------------------------------------------------------------

  return 0;
}  // assign_grid_structure()



// ##################################################################
// ##################################################################



int UniformGrid::set_cell_size()
{
  //
  // Uniform Cartesian grid, with cells that have the same length in
  // each direction, so this is very easy...
  //
  spdlog::debug("Setting G_dx=constant for all cells");

  G_dx = G_range[0] / (G_ng[0]);

  if (G_ndim > 1) {
    if (!pconst.equalD(G_range[1] / G_dx, static_cast<double>(G_ng[1])))
      spdlog::error(
          "{}: {}",
          "Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          G_range[1] / G_dx / G_ng[1]);
  }

  if (G_ndim > 2) {
    if (!pconst.equalD(G_range[2] / G_dx, static_cast<double>(G_ng[2])))
      spdlog::error(
          "{}: {}",
          "Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          G_range[2] / G_dx / G_ng[2]);
  }

  //
  // Surface area of interface: It is assumed extra dimensions are
  // per unit length.
  //
  if (G_ndim == 1)
    G_dA = 1.;
  else if (G_ndim == 2)
    G_dA = G_dx;
  else
    G_dA = G_dx * G_dx;

  //
  // Volume of cell.
  //
  if (G_ndim == 1)
    G_dV = G_dx;
  else if (G_ndim == 2)
    G_dV = G_dx * G_dx;
  else
    G_dV = G_dx * G_dx * G_dx;

  return 0;
}  // set_cell_size



// ##################################################################
// ##################################################################



enum direction UniformGrid::OppDir(enum direction dir)
{
  // This should be in VectorOps, or somewhere like that...
  if (dir == XP)
    return (XN);
  else if (dir == XN)
    return (XP);
  else if (dir == YP)
    return (YN);
  else if (dir == YN)
    return (YP);
  else if (dir == ZP)
    return (ZN);
  else if (dir == ZN)
    return (ZP);
  else {
    spdlog::error("{}: {}", "Bad direction given to OppDir", dir);
    return (NO);
  }
}



// ##################################################################
// ##################################################################



class cell *UniformGrid::FirstPt()
{
  ///
  /// First point is always the at the negative corner of the grid in
  /// all directions.
  ///
  return (G_fpt);
}  // FirstPt



// ##################################################################
// ##################################################################



class cell *UniformGrid::FirstPt_All()
{
  ///
  /// First point is always the at the negative corner of the grid in
  /// all directions.
  ///
  return (G_fpt_all);
}  // FirstPt



// ##################################################################
// ##################################################################



class cell *UniformGrid::LastPt()
{
  return (G_lpt);
}  // LastPt



// ##################################################################
// ##################################################################



class cell *UniformGrid::LastPt_All()
{
  return (G_lpt_all);
}  // LastPt



// ##################################################################
// ##################################################################



class cell *UniformGrid::PrevPt(const class cell *p, enum direction dir)
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
  // This is going to be very inefficient...
  spdlog::debug(
      "This function is very inefficient and probably shouldn't be used");
  return (p->ngb[opp]);
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
      c = cy;
      // add boundary cells beside the grid.
      for (int v = 0; v < G_nbc[XN]; v++)
        c = NextPt(c, XN);
      if (!c) spdlog::error("{}: {}", "Got lost on grid! XN", cy->id);
      for (int v = 0; v < G_nbc[XN]; v++) {
        BC_bd[XN]->data.push_back(c);
        // spdlog::debug(" Adding cell {} to XN boundary", c->id);
        c = NextPt(c, XP);
      }
      if (G_ndim > 1) cy = NextPt(cy, YP);
    } while (G_ndim > 1 && cy != 0 && cy->isgd);
    if (G_ndim > 2) cz = NextPt(cz, ZP);
  } while (G_ndim > 2 && cz != 0 && cz->isgd);

  spdlog::debug("Setup XN boundary, got {} grid cells", BC_bd[XN]->data.size());

  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // Now do the same for XP, beginning at the most negative point
  // at the XP boundary and moving to the most positive.
  //
  cz = FirstPt();
  while (NextPt(cz, XP)->isgd)
    cz = NextPt(cz, XP);
  // loop in ZP direction
  do {
    cy = cz;
    // loop in YP direction
    do {
      c = cy;
      // add boundary cells beside the grid.
      for (int v = 0; v < G_nbc[XP]; v++) {
        c = NextPt(c, XP);
        if (!c) {
          CI.print_cell(cy);
          spdlog::error("{}: {}", "Got lost on grid! XP", cy->id);
        }
        BC_bd[XP]->data.push_back(c);
        // spdlog::debug(" Adding cell {} to XP boundary", c->id);
      }
      if (G_ndim > 1) cy = NextPt(cy, YP);
    } while (G_ndim > 1 && cy != 0 && cy->isgd);
    if (G_ndim > 2) cz = NextPt(cz, ZP);
  } while (G_ndim > 2 && cz != 0 && cz->isgd);
  spdlog::debug("Setup XP boundary, got {} grid cells", BC_bd[XP]->data.size());
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // The YN/YP boundaries are different because they have the
  // corner data too.
  //
  if (G_ndim > 1) {
    // --------------------------------------------------------------
    //
    // First YN boundary.
    // Get to most negative cell, and add cells, progressing to
    // the next one.
    //
    cz = FirstPt();
    while (NextPt(cz, XN) != 0)
      cz = NextPt(cz, XN);
    while (NextPt(cz, YN) != 0)
      cz = NextPt(cz, YN);
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
        cy = NextPt_All(cy);
      } while (cy->pos[YY] < G_ixmin[YY]);

      if (G_ndim > 2) cz = NextPt(cz, ZP);
    } while (G_ndim > 2 && cz != 0 && cz->pos[ZZ] < G_ixmax[ZZ]);
    spdlog::debug(
        "Setup YN boundary, got {} grid cells", BC_bd[YN]->data.size());
    // --------------------------------------------------------------

    // --------------------------------------------------------------
    //
    // Then YP boundary.
    // Get to the first row of YP boundary cells, move to the
    // first cell in the XY plane, and go from there.
    //
    cz = FirstPt();
    while (cz->isgd)
      cz = NextPt(cz, YP);
    while (NextPt(cz, XN) != 0)
      cz = NextPt(cz, XN);

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
        cy = NextPt_All(cy);
      } while ((cy != 0) && (cy->pos[YY] > G_ixmax[YY]));

      if (G_ndim > 2) cz = NextPt(cz, ZP);
    } while (G_ndim > 2 && cz != 0 && cz->pos[ZZ] < G_ixmax[ZZ]);

    spdlog::debug(
        "Setup YP boundary, got {} grid cells", BC_bd[YP]->data.size());
    // --------------------------------------------------------------
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // Now the Z-boundary
  //
  if (G_ndim > 2) {
    //
    // ZN is easy... all points until pos[Z] > xmin[Z]
    //
    cz = FirstPt_All();
    do {
      BC_bd[ZN]->data.push_back(cz);
      cz = NextPt_All(cz);
    } while (cz->pos[ZZ] < G_ixmin[ZZ]);
    spdlog::debug(
        "Setup ZN boundary, got {} grid cells", BC_bd[ZN]->data.size());
    //
    // ZP is also easy... all points with pos[Z] > xmax[Z]
    //
    cz = LastPt();
    while (cz->pos[ZZ] < G_ixmax[ZZ])
      cz = NextPt_All(cz);
    do {
      BC_bd[ZP]->data.push_back(cz);
      cz = NextPt_All(cz);
    } while (cz != 0);
    spdlog::debug(
        "Setup ZP boundary, got {} grid cells", BC_bd[ZP]->data.size());
  }
  // ----------------------------------------------------------------
  return 0;
}



// ##################################################################
// ##################################################################



int UniformGrid::BC_printBCdata(boundary_data *b)
{
  list<cell *>::iterator c = b->data.begin();
  for (c = b->data.begin(); c != b->data.end(); ++c) {
    CI.print_cell(*c);
  }
  return 0;
}



// ##################################################################
// ##################################################################



void UniformGrid::BC_deleteBoundaryData(boundary_data *b)
{
  if (b->refval != 0) {
    b->refval = mem.myfree(b->refval);
  }

  list<cell *>::iterator i = b->data.begin();
  if (b->data.empty()) {
    spdlog::debug("BC destructor: No boundary cells to delete");
  }
  else {
    do {
      b->data.erase(i);
      i = b->data.begin();
    } while (i != b->data.end());
  }

  if (b->send_data.size() > 0) {
    i = b->send_data.begin();
    do {
      b->send_data.erase(i);
      i = b->send_data.begin();
    } while (i != b->send_data.end());
  }

  for (unsigned int j = 0; j < b->NGrecvF2C.size(); j++) {
    i = b->NGrecvF2C[j].begin();
    do {
      b->NGrecvF2C[j].erase(i);
      i = b->NGrecvF2C[j].begin();
    } while (i != b->NGrecvF2C[j].end());
  }
  b->NGrecvF2C.clear();

  for (unsigned int j = 0; j < b->NGrecvC2F.size(); j++) {
    i = b->NGrecvC2F[j].begin();
    do {
      b->NGrecvC2F[j].erase(i);
      i = b->NGrecvC2F[j].begin();
    } while (i != b->NGrecvC2F[j].end());
  }
  b->NGrecvC2F.clear();

  for (unsigned int j = 0; j < b->NGsendC2F.size(); j++) {
    b->NGsendC2F[j]->c.clear();
    struct c2f *t = b->NGsendC2F[j];
    delete t;
    // b->NGsendC2F[j] = mem.myfree( b->NGsendC2F[j]);
  }
  b->NGsendC2F.clear();

  for (unsigned int j = 0; j < b->avg.size(); j++) {
    b->avg[j].avg_state = mem.myfree(b->avg[j].avg_state);
    b->avg[j].c.clear();
  }
  b->avg.clear();

  return;
}



// ##################################################################
// ##################################################################



void UniformGrid::BC_deleteBoundaryData()
{
  spdlog::debug("BC destructor: deleting Boundary data..");
  struct boundary_data *b;
  for (unsigned int ibd = 0; ibd < BC_bd.size(); ibd++) {
    b = BC_bd[ibd];
    BC_deleteBoundaryData(b);
  }  // loop over all boundaries.
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
    const std::array<double, MAX_DIM> &p1,  ///< position 1 (physical)
    const std::array<double, MAX_DIM> &p2   ///< position 2 (physical)
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += (p1[i] - p2[i]) * (p1[i] - p2[i]);
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
    const std::array<int, MAX_DIM> &p1,  ///< position 1 (integer)
    const std::array<int, MAX_DIM> &p2   ///< position 2 (integer)
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += static_cast<double>((p1[i] - p2[i]) * (p1[i] - p2[i]));
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
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += pow(CI.get_dpos(c1, i) - CI.get_dpos(c2, i), 2.0);
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
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp +=
        pow(static_cast<double>(CI.get_ipos(c1, i) - CI.get_ipos(c2, i)), 2.0);
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
    const std::array<double, MAX_DIM> &v,  ///< vertex (physical)
    const cell *c                          ///< cell
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += pow(v[i] - CI.get_dpos(c, i), 2.0);
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
    const cell *c,    ///< cell
    const axes a      ///< Axis to calculate.
)
{
  return (CI.get_dpos(c, a) - v[a]);
}



// ##################################################################
// ##################################################################



//
// Calculate distance between a cell-vertex and a cell--centres
// (will be between centre-of-volume of cells if non-cartesian
// geometry).  Here both input and output are code-integer units.
//
double UniformGrid::idistance_vertex2cell(
    const std::array<int, MAX_DIM> &v,  ///< vertex (integer)
    const cell *c                       ///< cell
)
{
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += pow(static_cast<double>(v[i] - CI.get_ipos(c, i)), 2.0);
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
    const int *v,   ///< vertex (integer)
    const cell *c,  ///< cell
    const axes a    ///< Axis to calculate.
)
{
  return (CI.get_ipos(c, a) - v[a]);
}



// ##################################################################
// ##################################################################



//
// As idifference_vertex2cell(int,cell,axis) but for the coordinate
// difference between two cell positions along a given axis.
// It returns *cell2* coordinate minus *cell1* coordinate.
//
double UniformGrid::idifference_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2,  ///< cell 2
    const axes a     ///< Axis.
)
{
  return (CI.get_ipos(c2, a) - CI.get_ipos(c1, a));
}



// ##################################################################
// ##################################################################



bool UniformGrid::point_on_grid(const double *pos  ///< position
)
{
  bool on = true;
  for (int v = 0; v < G_ndim; v++) {
    if (pos[v] < G_xmin[v] || pos[v] > G_xmax[v]) on = false;
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



uniform_grid_cyl::uniform_grid_cyl(
    int nd,   ///< ndim, length of position vector.
    int nv,   ///< nvar, length of state vectors.
    int eqt,  ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,  ///< Number of boundary cells to use.
    const std::array<double, MAX_DIM>
        &g_xn,  ///< array of minimum values of x,y,z for this grid.
    const std::array<double, MAX_DIM>
        &g_xp,  ///< array of maximum values of x,y,z for this grid.
    const std::array<int, MAX_DIM>
        &g_nc,  ///< array of number of cells in x,y,z directions.
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM>
        &sim_xn,  ///< array of min. x/y/z for full simulation.
    std::array<double, MAX_DIM>
        &sim_xp  ///< array of max. x/y/z for full simulation.
    ) :
    VectorOps_Cart(nd),
    UniformGrid(
        nd, nv, eqt, Nbc, g_xn, g_xp, g_nc, lev_xn, lev_xp, sim_xn, sim_xp),
    VectorOps_Cyl(nd)
{
  spdlog::debug(
      "Setting up cylindrical uniform grid with G_ndim={} and G_nvar={}",
      G_ndim, G_nvar);

  if (G_ndim != 2)
    spdlog::error("{}: {}", "Need to write code for !=2 dimensions", G_ndim);
  G_coordsys = COORD_CYL;  // Cylindrical Coordinate system

  spdlog::debug("cylindrical grid: dr={}", G_dx);
  return;
}



// ##################################################################
// ##################################################################



uniform_grid_cyl::~uniform_grid_cyl()
{
  spdlog::debug("uniform_grid_cyl destructor. Present and correct");
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
  spdlog::error("{}: {}", "redefine iR_cov() for parallel grid", "please");
#endif
  return (R_com(c, G_dx) - G_xmin[Rcyl]) / CI.phys_per_int() + G_ixmin[Rcyl];
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::distance(
    const std::array<double, MAX_DIM> &p1,  ///< position 1 (physical)
    const std::array<double, MAX_DIM> &p2   ///< position 2 (physical)
)
{
  //
  // This is the same as for cartesian as long as there is no 3rd
  // dimension.
  //
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::idistance(
    const std::array<int, MAX_DIM> &p1,  ///< position 1 (integer)
    const std::array<int, MAX_DIM> &p2   ///< position 2 (integer)
)
{
  //
  // This is the same as for cartesian as long as there is no 3rd
  // dimension.
  //
  double temp = 0.0;
  for (int i = 0; i < G_ndim; i++)
    temp += static_cast<double>((p1[i] - p2[i]) * (p1[i] - p2[i]));
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::distance_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d = 0.0, temp = 0.0;
  temp = CI.get_dpos(c1, Zcyl) - CI.get_dpos(c2, Zcyl);
  d += temp * temp;
  temp = R_com(c1, G_dx) - R_com(c2, G_dx);
  d += temp * temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::idistance_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d = 0.0, temp;
  // r1,r2;
  temp = static_cast<double>(CI.get_ipos(c1, Zcyl) - CI.get_ipos(c2, Zcyl));
  d += temp * temp;
  temp = iR_cov(c1) - iR_cov(c2);
  d += temp * temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::distance_vertex2cell(
    const std::array<double, MAX_DIM> &v,  ///< vertex (physical)
    const cell *c                          ///< cell
)
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  double d = 0.0, temp = 0.0;
  temp = v[Zcyl] - CI.get_dpos(c, Zcyl);
  d += temp * temp;
  temp = v[Rcyl] - R_com(c, G_dx);
  d += temp * temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::difference_vertex2cell(
    const double *v,  ///< vertex (double)
    const cell *c,    ///< cell
    const axes a      ///< Axis to calculate.
)
{
  if (a == Zcyl) {
    return (CI.get_dpos(c, a) - v[a]);
  }
  else if (a == Rcyl) {
    return R_com(c, G_dx) - v[Rcyl];
  }
  else {
    spdlog::error("Requested cylindrical distance in theta dir");
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::idistance_vertex2cell(
    const std::array<int, MAX_DIM> &v,  ///< vertex (integer)
    const cell *c                       ///< cell
)
{
  //
  // The z-direction is a simple cartesian calculation, but the radial
  // coordinate needs to be at the centre--of--volume of each cell.
  //
  // rep.printVec("idist vertex",v,G_ndim);
  // rep.printVec("idist cell  ",c->pos,G_ndim);
  double d = 0.0, temp = 0.0;
  temp = static_cast<double>(v[Zcyl] - CI.get_ipos(c, Zcyl));
  d += temp * temp;
  // iR_cov() gives coordinates of centre-of-volume radius in integer
  // units.
  temp = static_cast<double>(v[Rcyl]) - iR_cov(c);
  d += temp * temp;
  return sqrt(d);
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::idifference_vertex2cell(
    const int *v,   ///< vertex (integer)
    const cell *c,  ///< cell
    const axes a    ///< Axis to calculate.
)
{
  if (a == Zcyl)
    return (CI.get_ipos(c, a) - v[a]);
  else if (a == Rcyl) {
    return (iR_cov(c) - v[a]);
  }
  else
    spdlog::error(
        "{}: {}", "Bad axis for uniform_grid_cyl::idifference_vertex2cell()",
        a);

  return -1.0;
}



// ##################################################################
// ##################################################################



double uniform_grid_cyl::idifference_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2,  ///< cell 2
    const axes a     ///< Axis.
)
{
  if (a == Zcyl)
    return (CI.get_ipos(c2, a) - CI.get_ipos(c1, a));
  else if (a == Rcyl)
    return (iR_cov(c2) - iR_cov(c1));
  else
    spdlog::error(
        "{}: {}", "Bad axis for uniform_grid_cyl::idifference_cell2cell()", a);
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
    int nd,   ///< ndim, length of position vector.
    int nv,   ///< nvar, length of state vectors.
    int eqt,  ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,  ///< Number of boundary cells to use.
    const std::array<double, MAX_DIM>
        &g_xn,  ///< array of minimum values of x,y,z for this grid.
    const std::array<double, MAX_DIM>
        &g_xp,  ///< array of maximum values of x,y,z for this grid.
    const std::array<int, MAX_DIM>
        &g_nc,  ///< array of number of cells in x,y,z directions.
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM>
        &sim_xn,  ///< array of min. x/y/z for full simulation.
    std::array<double, MAX_DIM>
        &sim_xp  ///< array of max. x/y/z for full simulation.
    ) :
    VectorOps_Cart(nd),
    UniformGrid(
        nd, nv, eqt, Nbc, g_xn, g_xp, g_nc, lev_xn, lev_xp, sim_xn, sim_xp),
    VectorOps_Cyl(nd), VectorOps_Sph(nd)
{
  spdlog::debug(
      "Setting up spherical uniform grid with G_ndim={} and G_nvar={}", G_ndim,
      G_nvar);
  if (G_ndim != 1)
    spdlog::error("{}: {}", "Need to write code for >1 dimension", G_ndim);
  G_coordsys = COORD_SPH;  // Spherical Coordinate system

  spdlog::debug("spherical grid: dr=", G_dx);
  return;
}



// ##################################################################
// ##################################################################



uniform_grid_sph::~uniform_grid_sph()
{
  spdlog::debug("uniform_grid_sph destructor. Present and correct");
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
  spdlog::error("{}: {}", "redefine iR_cov() for parallel grid", "please");
#endif
  return ((R_com(c, G_dx) - G_xmin[Rsph]) / CI.phys_per_int() + G_ixmin[Rsph]);
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::distance(
    const std::array<double, MAX_DIM> &p1,  ///< position 1 (physical)
    const std::array<double, MAX_DIM> &p2   ///< position 2 (physical)
)
{
  return fabs(p1[Rsph] - p2[Rsph]);
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::idistance(
    const std::array<int, MAX_DIM> &p1,  ///< position 1 (integer)
    const std::array<int, MAX_DIM> &p2   ///< position 2 (integer)
)
{
  return fabs(static_cast<double>(p1[Rsph] - p2[Rsph]));
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::distance_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  return fabs(R_com(c1, G_dx) - R_com(c2, G_dx));
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::idistance_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2   ///< cell 2
)
{
  return fabs(R_com(c1, G_dx) - R_com(c2, G_dx)) / CI.phys_per_int();
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::distance_vertex2cell(
    const std::array<double, MAX_DIM> &v,  ///< vertex (physical)
    const cell *c                          ///< cell
)
{
  return fabs(v[Rsph] - R_com(c, G_dx));
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::difference_vertex2cell(
    const double *v,  ///< vertex (double)
    const cell *c,    ///< cell
    const axes a      ///< Axis to calculate.
)
{
  if (a == Rsph) {
    return R_com(c, G_dx) - v[Rsph];
  }
  else {
    spdlog::error(" Requested spherical distance in not-radial dir");
    return -1.0e99;
  }
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::idistance_vertex2cell(
    const std::array<int, MAX_DIM> &v,  ///< vertex (integer)
    const cell *c                       ///< cell
)
{
  return fabs(static_cast<double>(v[Rsph]) - iR_cov(c));
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::idifference_vertex2cell(
    const int *v,   ///< vertex (integer)
    const cell *c,  ///< cell
    const axes a    ///< Axis to calculate.
)
{
  if (a == Rsph)
    return (iR_cov(c) - static_cast<double>(v[a]));
  else
    spdlog::error(
        "{}: {}", "Bad axis for uniform_grid_sph::idifference_vertex2cell()",
        a);

  return -1.0;
}



// ##################################################################
// ##################################################################



double uniform_grid_sph::idifference_cell2cell(
    const cell *c1,  ///< cell 1
    const cell *c2,  ///< cell 2
    const axes a     ///< Axis.
)
{
  if (a == Rsph)
    return (iR_cov(c2) - iR_cov(c1));
  else
    spdlog::error(
        "{}: {}", "Bad axis for uniform_grid_sph::idifference_cell2cell()", a);
  return -1.0;
}



// ##################################################################
// ##################################################################

//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------
