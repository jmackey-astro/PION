///
/// \file NG_BC89flux.cpp
///
/// \author Jonathan Mackey
///
/// Function definitions for fluxes leaving/entering
/// levels on a nested grid, to make sure that fluxes between
/// adjacent levels are consistent with each other.
///
/// \author Jonathan Mackey
///
/// Modifications:n
/// - 2019.12.04 JM: moved functions from uniform_grid and
///   sim_control_NG classes to gather everything in one file.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_BC89flux.h"
#include "constants.h"
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>


using namespace std;

// ##################################################################
// ##################################################################

NG_BC89flux::NG_BC89flux() {}

// ##################################################################
// ##################################################################

NG_BC89flux::~NG_BC89flux()
{
  struct flux_interface *fi = 0;
  for (unsigned int l = 0; l < flux_update_recv.size(); l++) {
    for (unsigned int d = 0; d < flux_update_recv[l].size(); d++) {
      for (unsigned int f = 0; f < flux_update_recv[l][d].fi.size(); f++) {
        fi = flux_update_recv[l][d].fi[f];
        if (fi) {
          fi->c.clear();
          fi->area.clear();
          fi->flux = mem.myfree(fi->flux);
          flux_update_recv[l][d].fi[f] =
              mem.myfree(flux_update_recv[l][d].fi[f]);
        }
      }
    }
  }
  flux_update_recv.clear();
  for (unsigned int l = 0; l < flux_update_send.size(); l++) {
    for (unsigned int d = 0; d < flux_update_send[l].size(); d++) {
      for (unsigned int f = 0; f < flux_update_send[l][d].fi.size(); f++) {
        fi = flux_update_send[l][d].fi[f];
        if (fi) {
          fi->c.clear();
          fi->area.clear();
          fi->flux = mem.myfree(fi->flux);
          flux_update_send[l][d].fi[f] =
              mem.myfree(flux_update_send[l][d].fi[f]);
        }
      }
      flux_update_send[l].clear();
    }
  }
  flux_update_send.clear();
  return;
}

// ##################################################################
// ##################################################################

void NG_BC89flux::setup_flux_vectors(
    const int nl  ///< number of levels in nested grid.
)
{
  flux_update_recv.resize(nl);
  flux_update_send.resize(nl);
  return;
}

// ##################################################################
// ##################################################################

int NG_BC89flux::setup_flux_recv(
    class SimParams &par,       // simulation params (including BCs)
    class GridBaseClass *grid,  // pointer to coarse grid.
    const int lp1               // fine level to receive from
)
{
#ifdef TEST_BC89FLUX
  spdlog::info(" NG_BC89flux::setup_flux_recv() recv from level={}", lp1);
#endif
  //
  // Get size of interface region and number of cells.
  //
  size_t nc = 1;  // number of cells in each interface
  int l     = lp1 - 1;
  std::array<int, MAX_DIM> ixmin, ixmax, ncell;  // interface
  std::array<int, MAX_DIM> f_lxmin, f_lxmax;     // finer level
  std::array<int, MAX_DIM> cg_ixmin, cg_ixmax;   // coarser grid

  bool recv[2 * par.ndim];   // whether to get data in this direction
  size_t nel[2 * par.ndim];  // number of interfaces in each direction
  struct flux_interface *fi = 0;
  int idx                   = grid->idx();

  if (grid != par.levels[lp1].parent)
    spdlog::error("{}: {}", "level lp1 is not my child!", lp1);

  CI.get_ipos_vec(par.levels[lp1].Xmin, f_lxmin);
  CI.get_ipos_vec(par.levels[lp1].Xmax, f_lxmax);

  for (int v = 0; v < par.ndim; v++)
    cg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    cg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

  // define interface region of fine and coarse grids, and whether
  // each direction is to be included or not.  Note that this allows
  // for a fine grid that is not completely encompassed by the coarse
  // grid.
  for (int v = 0; v < par.ndim; v++) {
    ixmin[v] = std::max(cg_ixmin[v], f_lxmin[v]);
    if (cg_ixmin[v] < f_lxmin[v])
      recv[2 * v] = true;
    else
      recv[2 * v] = false;

    ixmax[v] = std::min(cg_ixmax[v], f_lxmax[v]);
    if (cg_ixmax[v] > f_lxmax[v])
      recv[2 * v + 1] = true;
    else
      recv[2 * v + 1] = false;

    ncell[v] = (ixmax[v] - ixmin[v]) / idx;
    if ((ixmax[v] - ixmin[v]) % idx != 0) {
      spdlog::error(
          "{}: {}", "interface region not divisible!", ixmax[v] - ixmin[v]);
    }
  }
  for (int v = 0; v < 2 * par.ndim; v++)
    nel[v] = 0;

  // different number of interfaces depending on dimensionality.
  switch (par.ndim) {
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
      if (recv[XN]) nel[XN] = ncell[YY] * ncell[ZZ];
      if (recv[XP]) nel[XP] = ncell[YY] * ncell[ZZ];
      if (recv[YN]) nel[YN] = ncell[XX] * ncell[ZZ];
      if (recv[YP]) nel[YP] = ncell[XX] * ncell[ZZ];
      if (recv[ZN]) nel[ZN] = ncell[XX] * ncell[YY];
      if (recv[ZP]) nel[ZP] = ncell[XX] * ncell[YY];
      break;
    default:
      spdlog::error("{}: {}", "bad ndim in setup_flux_recv", par.ndim);
      break;
  }

  // initialize arrays for class member data
  flux_update_recv[l].resize(2 * par.ndim);
  for (int d = 0; d < 2 * par.ndim; d++) {
    if (recv[d] == true) {
      flux_update_recv[l][d].fi.resize(nel[d]);
      flux_update_recv[l][d].dir    = d;
      flux_update_recv[l][d].ax     = d / 2;
      flux_update_recv[l][d].Ncells = nc;
      for (size_t i = 0; i < nel[d]; i++) {
        flux_update_recv[l][d].fi[i] =
            mem.myalloc(flux_update_recv[l][d].fi[i], 1);
        fi = flux_update_recv[l][d].fi[i];
        fi->c.resize(nc);
        fi->area.resize(nc);
        fi->flux = mem.myalloc(fi->flux, par.nvar);
        for (int v = 0; v < par.nvar; v++)
          fi->flux[v] = 0.0;
      }
    }
    else {
      flux_update_recv[l][d].fi.resize(1);
      flux_update_recv[l][d].fi[0]  = 0;
      flux_update_recv[l][d].dir    = d;
      flux_update_recv[l][d].ax     = d / 2;
      flux_update_recv[l][d].Ncells = nc;
    }
  }

  // For each interface, find the cell that is outside the fine grid
  // and that includes the interface.
  for (int d = 0; d < 2 * par.ndim; d++) {
    if (recv[d]) {
#ifdef TEST_BC89FLUX
      spdlog::debug(
          "UG:FRECV: d={}, nel={}... starting; l+1={}, l={}", d, nel[d], lp1,
          lp1 - 1);
#endif
      add_cells_to_face(
          par, grid, static_cast<enum direction>(d), ixmin, ixmax, ncell, 1,
          flux_update_recv[l][d]);
#ifdef TEST_BC89FLUX
      spdlog::debug("UG:FRECV: d={}, nel={}", d, nel[d]);
#endif
    }
  }
  return 0;
}

// ##################################################################
// ##################################################################

int NG_BC89flux::setup_flux_send(
    class SimParams &par,       ///< simulation params (including BCs)
    class GridBaseClass *grid,  // pointer to finer grid.
    const int lm1               ///< level to send to
)
{
#ifdef TEST_BC89FLUX
  spdlog::info(" NG_BC89flux::setup_flux_send() send to level={}", lm1);
#endif
  //
  // Get size of interface region and number of cells.
  //
  size_t nc = 1;  // number of cells in each interface
  int l     = lm1 + 1;
  std::array<int, MAX_DIM> ixmin, ixmax, ncell,
      nface;                                    // interface
  std::array<int, MAX_DIM> c_lxmin, c_lxmax;    // coarser grid
  std::array<int, MAX_DIM> fg_ixmin, fg_ixmax;  // finer grid
  bool send[2 * par.ndim];   // whether to send data in this direction
  size_t nel[2 * par.ndim];  // number of interfaces in each direction
  struct flux_interface *fi = 0;
  int idx                   = grid->idx();

#ifdef TEST_BC89FLUX
  if (grid != par.levels[lm1].child)
    spdlog::error("{}: {}", "level l is not my parent!", lm1);
#endif

  CI.get_ipos_vec(par.levels[lm1].Xmin, c_lxmin);
  CI.get_ipos_vec(par.levels[lm1].Xmax, c_lxmax);

  for (int v = 0; v < par.ndim; v++)
    fg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    fg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

  // define interface region of fine and coarse grids, and whether
  // each direction is to be included or not.  Note that this allows
  // for a fine grid that is not completely encompassed by the coarse
  // grid.
  for (int ax = 0; ax < par.ndim; ax++) {
    ixmin[ax] = std::max(fg_ixmin[ax], c_lxmin[ax]);
    if (fg_ixmin[ax] > c_lxmin[ax])
      send[2 * ax] = true;
    else
      send[2 * ax] = false;

    ixmax[ax] = std::min(fg_ixmax[ax], c_lxmax[ax]);
    if (fg_ixmax[ax] < c_lxmax[ax])
      send[2 * ax + 1] = true;
    else
      send[2 * ax + 1] = false;

    ncell[ax] = (ixmax[ax] - ixmin[ax]) / idx;
    nface[ax] = ncell[ax] / 2;  // # face elements on coarse grid
#ifdef TEST_BC89FLUX
    spdlog::debug(
        "axis {}, ncell={}, nface={}\n\tixmin={}, ixmax={}", ax, ncell[ax],
        nface[ax], ixmin[ax], ixmax);
    if ((ixmax[ax] - ixmin[ax]) % 2 * idx != 0) {
      spdlog::error(
          "{}: {}", "interface region not divisible (send)!",
          ixmax[ax] - ixmin[ax]);
    }
#endif
  }
  for (int d = 0; d < 2 * par.ndim; d++)
    nel[d] = 0;

  // This is only for parallel execution: if the boundary
  // of my grid is not at the boundary of my level, then set send
  // to false
  CI.get_ipos_vec(par.levels[lm1 + 1].Xmin, c_lxmin);
  CI.get_ipos_vec(par.levels[lm1 + 1].Xmax, c_lxmax);
  for (int ax = 0; ax < par.ndim; ax++) {
    if (fg_ixmin[ax] > c_lxmin[ax]) send[2 * ax] = false;
    if (fg_ixmax[ax] < c_lxmax[ax]) send[2 * ax + 1] = false;
  }

  // different number of interfaces depending on dimensionality.
  switch (par.ndim) {
    case 1:
      if (send[XN]) nel[XN] = 1;
      if (send[XP]) nel[XP] = 1;
      nc = 1;
      break;
    case 2:
      if (send[XN]) nel[XN] = nface[YY];
      if (send[XP]) nel[XP] = nface[YY];
      if (send[YN]) nel[YN] = nface[XX];
      if (send[YP]) nel[YP] = nface[XX];
      nc = 2;
      break;
    case 3:
      if (send[XN]) nel[XN] = nface[YY] * nface[ZZ];
      if (send[XP]) nel[XP] = nface[YY] * nface[ZZ];
      if (send[YN]) nel[YN] = nface[XX] * nface[ZZ];
      if (send[YP]) nel[YP] = nface[XX] * nface[ZZ];
      if (send[ZN]) nel[ZN] = nface[XX] * nface[YY];
      if (send[ZP]) nel[ZP] = nface[XX] * nface[YY];
      nc = 4;
      break;
    default:
      spdlog::error("{}: {}", "bad ndim in setup_flux_send", par.ndim);
      break;
  }

  // initialize arrays
  flux_update_send[l].resize(2 * par.ndim);
  for (int d = 0; d < 2 * par.ndim; d++) {
    flux_update_send[l][d].dir = d;
    flux_update_send[l][d].ax  = d / 2;

    if (send[d] == true) {
#ifdef TEST_BC89FLUX
      spdlog::debug("d={}, nel={}", d, nel[d]);
#endif
      flux_update_send[l][d].fi.resize(nel[d]);
#ifdef TEST_BC89FLUX
      spdlog::debug(
          "d={}, nel={}... allocating memory for {} elements", d, nel[d],
          nel[d]);
#endif
      flux_update_send[l][d].Ncells = nc;
      for (size_t i = 0; i < nel[d]; i++) {
        flux_update_send[l][d].fi[i] =
            mem.myalloc(flux_update_send[l][d].fi[i], 1);
        fi = flux_update_send[l][d].fi[i];
        fi->c.resize(nc);
        fi->area.resize(nc);
        fi->flux = mem.myalloc(fi->flux, par.nvar);
        for (int v = 0; v < par.nvar; v++)
          fi->flux[v] = 0.0;
      }
    }
    else {
      flux_update_send[l][d].fi.resize(1);
      flux_update_send[l][d].fi[0]  = 0;
      flux_update_send[l][d].Ncells = nc;
    }
  }

  // For each interface, find the cell that is outside the fine grid
  // and that includes the interface.
  for (int d = 0; d < 2 * par.ndim; d++) {
    if (send[d]) {
#ifdef TEST_BC89FLUX
      spdlog::debug(
          "FLUX_SEND: d={}, nel={}... adding cells to face", d, nel[d]);
#endif
      add_cells_to_face(
          par, grid, static_cast<enum direction>(d), ixmin, ixmax, nface, 2,
          flux_update_send[l][d]);
    }
  }

  return 0;
}

// ##################################################################
// ##################################################################

int NG_BC89flux::add_cells_to_face(
    class SimParams &par,       ///< simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    enum direction d,           ///< which direction we're facing
    std::array<int, MAX_DIM>
        &ixmin,  ///< xmin of interface region (integer units)
    std::array<int, MAX_DIM>
        &ixmax,  ///< xmax of interface region (integer units)
    std::array<int, MAX_DIM>
        &nface,               ///< number of elements in interface region
    const int ncell,          ///< number of cells per face.
    struct flux_update &flux  ///< list to populate
)
{
  int a   = static_cast<int>(d) / 2;
  cell *c = grid->FirstPt_All();
  int idx = grid->idx();

#ifdef TEST_BC89FLUX
  spdlog::debug(
      "add_cells_to_face({}, {}, {}, {}, {}), list-size={}", d, ixmin[a],
      ixmax[a], nface[a], ncell, flux.fi.size());
  spdlog::debug("ixmin : {}", ixmin);
  spdlog::debug("ixmax : {}", ixmax);
  spdlog::debug("nface : {}", nface);
#endif

  //
  // Split code for different dimensionality of grid:
  // For all dimensions, every cell on the face needs to have an
  // extra Flux vector allocated, c->F[a], to store it.
  //
  if (par.ndim == 1) {
    if (d == XN) {
      while (c->pos[XX] < ixmin[XX] - idx)
        c = grid->NextPt(c, XP);
    }
    else if (d == XP) {
      while (c->pos[XX] < ixmax[XX])
        c = grid->NextPt(c, XP);
    }
    flux.fi[0]->c[0]    = c;
    c->F[a]             = mem.myalloc(c->F[a], par.nvar);
    c->isbd_ref[d]      = true;
    flux.fi[0]->area[0] = grid->CellInterface(c, grid->OppDir(d), 0);
  }  // 1D

  else if (par.ndim == 2) {
    enum direction perpdir = XP;
    enum axes perpaxis     = XX;
    //
    // get to first cell at interface, most negative position in all
    // directions.
    //
    switch (d) {
      case XN:
        while (c->pos[XX] < ixmin[XX] - idx)
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        perpdir  = YP;
        perpaxis = YY;
        break;
      case XP:
        while (c->pos[XX] < ixmax[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        perpdir  = YP;
        perpaxis = YY;
        break;
      case YN:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY] - idx)
          c = grid->NextPt(c, YP);
        perpdir  = XP;
        perpaxis = XX;
        break;
      case YP:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmax[YY])
          c = grid->NextPt(c, YP);
        perpdir  = XP;
        perpaxis = XX;
        break;
      default:
        spdlog::error("{}: {}", "bad direction in add_cells_to_face 2D", d);
    }

      // loop over cells in interface:
#ifdef TEST_BC89FLUX
    spdlog::debug("{}, {}", nface[perpaxis], flux.fi.size());
    if (!c) {
      spdlog::error("got lost on grid!");
      spdlog::debug("ixmin : {}", ixmin);
      spdlog::debug("ixmax : {}", ixmax);
      // rep.printVec("G_xmin",G_xmin,par.ndim);
      // rep.printVec("G_xmax",G_xmax,par.ndim);
      spdlog::error("{}: {}", "lost on grid", c);
    }
#endif
    if (nface[perpaxis] != static_cast<int>(flux.fi.size())) {
      spdlog::error(
          "add_cells_to_face({}, {}, {}, {}, {}), list-size={}", d, ixmin[a],
          ixmax[a], nface[a], ncell, flux.fi.size());
      spdlog::debug("ixmin : {}", ixmin);
      spdlog::debug("ixmax : {}", ixmax);
      spdlog::debug("nface : {}", nface);
      spdlog::error(
          "ERROR: nface[perpaxis]={}, and flux struct size is {}",
          nface[perpaxis], flux.fi.size());
      spdlog::error(
          "{}: {}", "wrong number of cells 2D interface", flux.fi.size());
    }

    for (int i = 0; i < nface[perpaxis]; i++) {
      for (int ic = 0; ic < ncell; ic++) {
        if (!c)
          spdlog::error(
              "{}: {}", "Cell is null in BC89 add_cells", fmt::ptr(c));
        flux.fi[i]->c[ic]    = c;
        c->F[a]              = mem.myalloc(c->F[a], par.nvar);
        c->isbd_ref[d]       = true;
        flux.fi[i]->area[ic] = grid->CellInterface(c, grid->OppDir(d), 0);
#ifdef TEST_BC89FLUX
        // cout <<"area["<<ic<<"] = "<<flux.fi[i]->area[ic]<<": adding
        // cell: "; rep.printVec("pos",c->pos,par.ndim);
        // CI.print_cell(c);
#endif
        c = grid->NextPt(c, perpdir);
      }
    }
  }  // 2D

  else {
    cell *marker = c, *m2 = c;
    enum direction perpdir1 = XP, perpdir2 = XP;
    enum axes perpaxis1 = XX, perpaxis2 = XX;
    // get to first cell at interface, most negative position in all
    // directions.
    switch (d) {
      case XN:
        while (c->pos[XX] < ixmin[XX] - idx)
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmin[ZZ])
          c = grid->NextPt(c, ZP);
        perpdir1  = YP;
        perpdir2  = ZP;
        perpaxis1 = YY;
        perpaxis2 = ZZ;
        break;
      case XP:
        while (c->pos[XX] < ixmax[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmin[ZZ])
          c = grid->NextPt(c, ZP);
        perpdir1  = YP;
        perpdir2  = ZP;
        perpaxis1 = YY;
        perpaxis2 = ZZ;
        break;
      case YN:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY] - idx)
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmin[ZZ])
          c = grid->NextPt(c, ZP);
        perpdir1  = XP;
        perpdir2  = ZP;
        perpaxis1 = XX;
        perpaxis2 = ZZ;
        break;
      case YP:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmax[YY])
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmin[ZZ])
          c = grid->NextPt(c, ZP);
        perpdir1  = XP;
        perpdir2  = ZP;
        perpaxis1 = XX;
        perpaxis2 = ZZ;
        break;
      case ZN:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmin[ZZ] - idx)
          c = grid->NextPt(c, ZP);
        perpdir1  = XP;
        perpdir2  = YP;
        perpaxis1 = XX;
        perpaxis2 = YY;
        break;
      case ZP:
        while (c->pos[XX] < ixmin[XX])
          c = grid->NextPt(c, XP);
        while (c->pos[YY] < ixmin[YY])
          c = grid->NextPt(c, YP);
        while (c->pos[ZZ] < ixmax[ZZ])
          c = grid->NextPt(c, ZP);
        perpdir1  = XP;
        perpdir2  = YP;
        perpaxis1 = XX;
        perpaxis2 = YY;
        break;
      default:
        spdlog::error("{}: {}", "bad direction in add_cells_to_face 3D", d);
    }

#ifdef TEST_BC89FLUX
    spdlog::debug(
        "d={}, perpdir1={}, perpdir2={}, perpaxis1={}, perpaxis2={}", d,
        perpdir1, perpdir2, perpaxis1, perpaxis2);
#endif

    // loop over cells in interface:
    if (nface[perpaxis1] * nface[perpaxis2]
        != static_cast<int>(flux.fi.size())) {
      spdlog::error(
          "add_cells_to_face({}, {}, {}, {}, {}), list-size={}", d, ixmin[a],
          ixmax[a], nface[a], ncell, flux.fi.size());
      spdlog::debug("ixmin : {}", ixmin);
      spdlog::debug("ixmax : {}", ixmax);
      spdlog::debug("nface : {}", nface);
      spdlog::error(
          "ERROR: nfaces={}, and flux struct size is {}",
          nface[perpaxis1] * nface[perpaxis2], flux.fi.size());
      spdlog::error(
          "{}: {}", "wrong number of cells 3D interface", flux.fi.size());
    }

    marker = c;
    for (int i = 0; i < nface[perpaxis2]; i++) {
      for (int j = 0; j < nface[perpaxis1]; j++) {
        if (ncell == 1) {
#ifdef TEST_BC89FLUX
          // cout <<"i="<<i<<", j="<<j<<", id="<<i*nface[perpaxis1]+j;
          // cout <<", fisize="<<flux.fi.size()<<"\n";
          // CI.print_cell(c); cout.flush();
#endif
          flux.fi[i * nface[perpaxis1] + j]->c[0] = c;
          c->F[a] = mem.myalloc(c->F[a], par.nvar);
          for (int v = 0; v < par.nvar; v++)
            c->F[a][v] = 0.0;
          c->isbd_ref[d] = true;
          flux.fi[i * nface[perpaxis1] + j]->area[0] =
              grid->CellInterface(c, grid->OppDir(d), 0);
        }
        else {
          // need to get 4 cells onto this face.

          flux.fi[i * nface[perpaxis1] + j]->c[0] = c;
          c->F[a] = mem.myalloc(c->F[a], par.nvar);
          for (int v = 0; v < par.nvar; v++)
            c->F[a][v] = 0.0;
          c->isbd_ref[d] = true;
          flux.fi[i * nface[perpaxis1] + j]->area[0] =
              grid->CellInterface(c, grid->OppDir(d), 0);
#ifdef TEST_BC89FLUX
          // cout <<"c1="<<c->id<<" ";
          // rep.printVec("pos",c->pos,par.ndim);
#endif

          m2                                      = grid->NextPt(c, perpdir1);
          flux.fi[i * nface[perpaxis1] + j]->c[1] = m2;
          m2->F[a] = mem.myalloc(m2->F[a], par.nvar);
          for (int v = 0; v < par.nvar; v++)
            m2->F[a][v] = 0.0;
          m2->isbd_ref[d] = true;
          flux.fi[i * nface[perpaxis1] + j]->area[1] =
              grid->CellInterface(m2, grid->OppDir(d), 0);
#ifdef TEST_BC89FLUX
          // cout <<"c2="<<m2->id<<" ";
          // rep.printVec("pos",m2->pos,par.ndim);
#endif

          m2                                      = grid->NextPt(c, perpdir2);
          flux.fi[i * nface[perpaxis1] + j]->c[2] = m2;
          m2->F[a] = mem.myalloc(m2->F[a], par.nvar);
          for (int v = 0; v < par.nvar; v++)
            m2->F[a][v] = 0.0;
          m2->isbd_ref[d] = true;
          flux.fi[i * nface[perpaxis1] + j]->area[2] =
              grid->CellInterface(m2, grid->OppDir(d), 0);
#ifdef TEST_BC89FLUX
          // cout <<"c3="<<m2->id<<" ";
          // rep.printVec("pos",m2->pos,par.ndim);
#endif

          m2                                      = grid->NextPt(m2, perpdir1);
          flux.fi[i * nface[perpaxis1] + j]->c[3] = m2;
          m2->F[a] = mem.myalloc(m2->F[a], par.nvar);
          for (int v = 0; v < par.nvar; v++)
            m2->F[a][v] = 0.0;
          m2->isbd_ref[d] = true;
          flux.fi[i * nface[perpaxis1] + j]->area[3] =
              grid->CellInterface(m2, grid->OppDir(d), 0);
#ifdef TEST_BC89FLUX
          // cout <<"c4="<<m2->id<<" ";
          // rep.printVec("pos",m2->pos,par.ndim);
#endif
        }
        for (int ic = 0; ic < ncell; ic++)
          c = grid->NextPt(c, perpdir1);
      }
      for (int ic = 0; ic < ncell; ic++)
        marker = grid->NextPt(marker, perpdir2);
      c = marker;
    }
  }  // 3D

  return 0;
}

// ##################################################################
// ##################################################################

void NG_BC89flux::save_fine_fluxes(
    class SimParams &par,  ///< simulation parameters
    const int l            ///< level
)
{
#ifdef SKIP_BC89_FLUX
  return;
#endif
#ifdef TEST_BC89FLUX
  spdlog::info(
      "save_fine_fluxes() \nsize of flux_update_send = {}",
      flux_update_send[l].size());
#endif
  int level_step            = par.levels[l].step;
  double dt                 = par.levels[l].dt;
  struct flux_interface *fi = 0;

  for (unsigned int d = 0; d < flux_update_send[l].size(); d++) {

#ifdef TEST_BC89FLUX
    spdlog::debug(
        "size of flux_update_send[{}] = {}", d,
        flux_update_send[l][d].fi.size());
#endif
    int a = flux_update_send[l][d].ax;  // normal axis to face

    for (unsigned int f = 0; f < flux_update_send[l][d].fi.size(); f++) {
      // these faces have 2^(ndim-1) cells.
      fi = flux_update_send[l][d].fi[f];

#ifdef TEST_BC89FLUX
      // cout <<"save_fine_fluxes["<<d<<"]["<<f<<"] = "<<fi<<"\n";
#endif

      if (fi == 0 && f == 0)
        continue;
      else if (fi == 0)
        spdlog::error("{}: {}", "save_fine_fluxes fi=0", d);

      // zero the fine fluxes only every 2nd step, because the 2
      // steps sum to a full step on the coarser grid.
      if (level_step % 2 == 0) {
        for (int v = 0; v < par.nvar; v++)
          fi->flux[v] = 0.0;
      }
      for (int i = 0; i < flux_update_send[l][d].Ncells; i++) {
        if (!fi->c[i]->F[a])
          spdlog::error("{}: {}", "fine flux is not allocated!", f);
          // Add cell flux to the full flux for this face over 2
          // steps.
#ifdef TEST_BC89FLUX
          // CI.print_cell(fi->c[i]);
          // CI.print_cell(grid->NextPt(fi->c[i],YP));
          // cout <<"PRE: flux="; rep.printVec("",fi->flux,par.nvar);
          // cout <<"save_fine_fluxes["<<d<<"]["<<f<<"]: i="<<i;
          // cout <<", f0="<<fi->c[i]->F[a][0]<<",
          // area="<<fi->area[i]<<", dt="<<dt<<"\n";
#endif
        for (int v = 0; v < par.nvar; v++) {
          fi->flux[v] += fi->c[i]->F[a][v] * fi->area[i] * dt;
        }
#ifdef TEST_BC89FLUX
        // rep.printVec("cell fine flux",fi->c[i]->F[a],par.nvar);
        // rep.printVec("cell fine Ph",fi->c[i]->Ph,par.nvar);
        // cout <<"save_fine_fluxes["<<d<<"]["<<f<<"]: i="<<i;
        // cout <<", flux="; rep.printVec("",fi->flux,par.nvar);
#endif
      }
#ifdef TEST_BC89FLUX
      spdlog::debug("save_fine_fluxes[{}][{}]: , flux=", d, f);
      spdlog::debug(" : {}", fi->flux);
#endif
    }
  }
  return;
}

// ##################################################################
// ##################################################################

void NG_BC89flux::save_coarse_fluxes(
    class SimParams &par,  ///< simulation parameters
    const int l            ///< level
)
{
#ifdef SKIP_BC89_FLUX
  return;
#endif
#ifdef TEST_BC89FLUX
  spdlog::info(
      "save_coarse_fluxes() \nsize of flux_update_recv = {}",
      flux_update_recv[l].size());
#endif
  double dt                 = par.levels[l].dt;
  struct flux_interface *fi = 0;

  for (unsigned int d = 0; d < flux_update_recv[l].size(); d++) {
    int a = flux_update_recv[l][d].ax;  // normal axis to face
    for (unsigned int f = 0; f < flux_update_recv[l][d].fi.size(); f++) {
      // these elements have only one cell per face.
      fi = flux_update_recv[l][d].fi[f];
#ifdef TEST_BC89FLUX
      // cout <<"save_coarse_fluxes["<<d<<"]["<<f<<"] = "<<fi<<"\n";
#endif

      // if first element is null then list is empty.
      if (fi == 0 && f == 0)
        continue;
      else if (fi == 0)
        spdlog::error("{}: {}", "save_fine_fluxes fi=0", d);
      if (!fi->c[0]->F[a])
        spdlog::error("{}: {}", "coarse flux is not allocated!", f);

        // set face flux to be the negative of the intercell flux
#ifdef TEST_BC89FLUX
        // cout <<"save_coarse_fluxes["<<d<<"]["<<f<<"]: i="<<0;
        // cout <<", f0="<<fi->c[0]->F[a][0]<<", area="<<fi->area[0];
        // cout <<", dt="<<dt<<"\n";
#endif
      for (int v = 0; v < par.nvar; v++) {
        fi->flux[v] = -fi->c[0]->F[a][v] * fi->area[0] * dt;
      }
    }
  }
  return;
}

// ##################################################################
// ##################################################################

int NG_BC89flux::recv_BC89_fluxes_F2C(
    class FV_solver_base *spatial_solver,  ///< solver, for gradients
    class SimParams &par,                  ///< simulation parameters
    const int l,                           ///< My level in grid hierarchy.
    const int step,                        ///< OA2 or OA1
    const int ooa  ///< Full order of accuracy of simulation
)
{
#ifdef SKIP_BC89_FLUX
  return 0;
#endif
  if (step != ooa) {
    spdlog::error("don't receive fluxes on half step");
    return 1;
  }
  if (l == par.grid_nlevels - 1) {
    spdlog::error("{}: {}", "finest level trying to receive data from l+1", l);
  }
#ifdef TEST_BC89FLUX
  spdlog::debug("level {}: correcting fluxes from finer grid", l);
#endif

  int err                   = 0;
  class GridBaseClass *grid = par.levels[l].grid;
  enum axes ax;
  double dt = par.levels[l].dt;

  // loop over boundaries, and, for each direction that has a
  // non-zero boundary, correct the coarse fluxes.
  for (unsigned int d = 0; d < flux_update_recv[l].size(); d++) {
#ifdef TEST_BC89FLUX
    spdlog::debug("flux update: {}, looping through faces", d);
#endif
    ax = static_cast<enum axes>(d / 2);

    if (flux_update_recv[l][d].fi[0] == 0) {
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      spdlog::debug("flux update: {}, receiving.", d);
#endif
      err += recv_BC89_flux_boundary(
          spatial_solver, par, grid, dt, flux_update_send[l + 1][d],
          flux_update_recv[l][d], d, ax);
    }

#ifdef TEST_BC89FLUX
    spdlog::debug("Direction: {}, finished.", d);
#endif
  }  // loop over directions
  return err;
}

// ##################################################################
// ##################################################################

int NG_BC89flux::recv_BC89_flux_boundary(
    class FV_solver_base *spatial_solver,  ///< solver, for gradients
    class SimParams &par,                  ///< pointer to simulation parameters
    class GridBaseClass *grid,             ///< pointer to coarse grid
    const double dt,                       ///< timestep
    struct flux_update &send,              ///< data for fine grid
    struct flux_update &recv,              ///< data for coarse grid
    const unsigned int d,                  ///< direction of outward normal
    const axes ax                          ///< axis of normal direction.
)
{
#ifdef SKIP_BC89_FLUX
  return 0;
#endif
  struct flux_interface *fc = 0;
  struct flux_interface *ff = 0;
  double ftmp[par.nvar], utmp[par.nvar];
  for (int v = 0; v < par.nvar; v++)
    ftmp[v] = 0.0;
  for (int v = 0; v < par.nvar; v++)
    utmp[v] = 0.0;

  if (send.fi.size() != recv.fi.size()) {
    spdlog::debug("send={}, recv={}", send.fi.size(), recv.fi.size());
    spdlog::error("{}: {}", "fine and parent face arrays r different size", 2);
  }

  for (unsigned int f = 0; f < recv.fi.size(); f++) {
    fc = recv.fi[f];
    ff = send.fi[f];

    for (int v = 0; v < par.nvar; v++)
      fc->flux[v] += ff->flux[v];
    for (int v = 0; v < par.nvar; v++)
      fc->flux[v] /= fc->area[0];
    for (int v = 0; v < par.nvar; v++)
      ftmp[v] = 0.0;

#ifdef TEST_BC89FLUX
      // for (int v=0;v<par.nvar;v++) fc->flux[v] *= fc->area[0];
#endif
    //
    // fc->flux is now the error in dU made for both coarse cells.
    // Correct dU in outer cell only, because inner cell is on
    // top of fine grid and gets overwritten anyway.
    //
    // If we are at a negative boundary, then it is the positive
    // face of the cell that we correct, and vice versa for a
    // positive boundary.  So we call the DivStateVectorComponent
    // function with the flux in different order, depending on
    // what direction the boundary is in.
    // This function is aware of geometry, so it calculates the
    // divergence correctly for curvilinear grids (important!!!).
    // The other face of the cell is set to zero flux.
    if (d % 2 == 0) {
      spatial_solver->DivStateVectorComponent(
          fc->c[0], grid, ax, par.nvar, ftmp, fc->flux, utmp);
    }
    else {
      spatial_solver->DivStateVectorComponent(
          fc->c[0], grid, ax, par.nvar, fc->flux, ftmp, utmp);
    }
    for (int v = 0; v < par.nvar; v++)
      fc->c[0]->dU[v] += utmp[v];

  }  // loop over interfaces.
  return 0;
}

// ##################################################################
// ##################################################################
