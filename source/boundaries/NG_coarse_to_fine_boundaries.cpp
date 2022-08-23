/// \file NG_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_coarse_to_fine boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "boundaries/NG_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

//#define TEST_C2F
//#define C2F_TAU

// ##################################################################
// ##################################################################

int NG_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
    class SimParams &par,        ///< simulation parameters
    class GridBaseClass *grid,   ///< pointer to grid.
    boundary_data *b,            ///< boundary data
    class GridBaseClass *parent  ///< pointer to parent grid.
)
{
  //
  // Make a list of pointers to cells in the coarser grid that map
  // onto this (finer) grid external boundary, and then write an
  // alogrithm to interpolate the coarse data onto the finer grid.
  //
  if (b->data.empty())
    spdlog::error(
        "{}: {}", "BC_assign_COARSE_TO_FINE: empty boundary data", b->itype);

    //
    // If we are doing raytracing, then also send the column densities
    // from coarse to fine grids.  Here we set up the number of them.
    //
#ifdef C2F_TAU
  struct rad_src_info *s;
  C2F_Nxd = 0;
  C2F_tauoff.resize(par.RS.Nsources);
  for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
    s                = &(par.RS.sources[isrc]);
    C2F_tauoff[isrc] = C2F_Nxd;
    C2F_Nxd += 2 * s->NTau;  // Need col2cell and cell_col for each Tau.
#ifdef TEST_MPI_NG
    spdlog::debug(
        "C2F_BC: RT Source {}: adding {} cols, running total = {}", isrc,
        2 * s->NTau, C2F_Nxd);
#endif
  }
#endif  // C2F_TAU

#ifdef TEST_C2F
  spdlog::debug("assigning C2F boundary for dir {}", b->dir);
#endif
  list<cell *>::iterator bpt = b->data.begin();
  int gidx                   = grid->idx();
  cell *pc                   = parent->FirstPt_All();  // parent cell.
  double distance            = 0.0;

  // Match parent cells with child grid boundary cells.
  do {
    pc = parent->FirstPt_All();

    for (int d = 0; d < par.ndim; d++) {
      // Find parent cell that covers this boundary cell.  It should be
      // G_idx/2 away from the boundary cell in each direction.
      enum axes ax       = static_cast<axes>(d);
      enum direction pos = static_cast<direction>(2 * d + 1);
      while (fabs(distance = grid->idifference_cell2cell(**bpt, *pc, ax))
             > 0.75 * gidx)
        pc = parent->NextPt(*pc, pos);
      if (!pc) spdlog::error("{}: {}", "C2F boundaries setup", distance);
    }
#ifdef TEST_C2F
    // rep.printVec("bpt pos",(*bpt)->pos,par.ndim);
    // cout <<"found parent: ";
    // rep.printVec("pc->pos",pc->pos,par.ndim);
#endif
    // set boundary cell's 'npt' pointer to point to the parent cell.
    (*bpt)->npt = pc;
    ++bpt;
  } while (bpt != b->data.end());

  return 0;
}

// ##################################################################
// ##################################################################

int NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(
    class SimParams &par,          ///< simulation parameters
    class FV_solver_base *solver,  ///< pointer to equations
    const int level,               ///< level in the NG grid structure
    struct boundary_data *b,
    const int step)
{
  if (level == 0) return 0;
#ifdef C2F_FULLSTEP
  if ((step) % 2 != 0) return 0;
#endif
#ifdef TEST_C2F
  spdlog::debug(
      "C2F: updating boundary data from coarse grid to level {}", level);
#endif

  //
  // Use linear interpolation (or
  // higher order), must conserve mass, momentum and
  // energy between the two levels.  It should also be monotonic and
  // produce cell-averaged values. Flash and pluto seem
  // to use very old Fortran code to accomplish this, and the code is
  // almost impenetrable.  AMRVAC uses linear/bilinear/trilinear
  // interpolation with renormalisation to conserve mass/P/E.
  // PION is doing what AMRVAC does, as far as I can tell.
  //

  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].parent;
  class GridBaseClass *fine   = par.levels[level].grid;
#ifdef TEST_C2F
  spdlog::debug("C2F: fine={}, parent={}, level={}", fine, coarse, level);
#endif

  list<cell *>::iterator f_iter = b->data.begin();
  cell *f, *c;

  // ----------------------------------------------------------------
  //
  // if on an odd-numbered step, then need to update the data on the
  // coarse grid to half way through a coarse step.  Assume dU has
  // been calculated for the coarse grid, but the state vector not
  // updated, so we can convert U+0.5*dU into a new primitive
  // state for the half step.
  // NB: We overwrite Ph[] in the coarse cell, assuming it is not
  // needed anymore on the coarse grid because dU[] is already
  // calculated.
  //
  if (step % 2 != 0) {
#ifdef TEST_C2F
    spdlog::info("C2F: odd step, interpolating coarse data in time.");
    int ct = 0;
#endif
    std::vector<double> U(par.nvar);
    for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
      f = (*f_iter);
      c = f->npt;
      solver->PtoU(c->P.data(), U.data(), par.gamma);
      for (int v = 0; v < par.nvar; v++)
        U[v] += 0.5 * c->dU[v];
      solver->UtoP(U.data(), c->Ph, par.EP.MinTemperature, par.gamma);
#ifdef TEST_C2F
      ct++;
      // cout <<"updated cell "<<ct<<" of "<<b->data.size()<<"\n";
      for (int v = 0; v < par.nvar; v++) {
        if (!isfinite(c->Ph[v])) {
          spdlog::debug("NAN c->P  : {}", c->P);
          spdlog::debug("NAN c->Ph : {}", c->Ph);
        }
      }
#endif
    }
  }
  // ----------------------------------------------------------------

#ifdef C2F_TAU
  // ----------------------------------------------------------------
  // save coarse cell optical depths onto fine cells (Experimental!)
  for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
    cell *f1 = 0, *f2 = 0;
    struct rad_src_info *s;
    int off = 0, ncells = 0;
    std::array<pion_flt, MAX_DIM> Tau, T, cpos;
    std::vector<cell *> fcl;
    c = (*f_iter)->npt;
    CI.get_dpos(c, cpos);

    // save Tau and dTau for each source in array T.
    for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
      s   = &(par.RS.sources[isrc]);
      off = C2F_tauoff[isrc];
      CI.get_col(c, s->id, Tau);
      for (int iT = 0; iT < s->NTau; iT++)
        T[off + iT] = Tau[iT];
      off += s->NTau;
      CI.get_cell_col(c, s->id, Tau);
      for (int iT = 0; iT < s->NTau; iT++)
        T[off + iT] = Tau[iT];
    }
    // rep.printVec("T coarse",T,C2F_Nxd);

    // add fine cells to list.
    if (par.ndim == 1) {
      fcl.push_back(*f_iter);
      f_iter++;
      fcl.push_back(*f_iter);
      ncells = 2;
    }
    else if (par.ndim == 2) {
      if (!fine->NextPt((*f_iter), YP)
          || fine->NextPt((*f_iter), YP)->npt != c) {
        // cout <<"skipping column\n";
        ncells = 0;
      }
      else {
        // get list of four fine cells.
        fcl.push_back(*f_iter);
        f1 = (*f_iter);
        f_iter++;
        fcl.push_back(*f_iter);
        f2 = (*f_iter);
        fcl.push_back(fine->NextPt(f1, YP));
        fcl.push_back(fine->NextPt(f2, YP));
        ncells = 4;
      }
    }
    else {
      spdlog::error("{}: {}", "C2F RT 3D not implemented", par.ndim);
    }
    if (ncells == 0) {
      // cout <<"ncells=0: skipping column\n";
      continue;  // must be on odd-numbered column
    }

    // rep.printVec("f0 pos",fcl[0]->pos,2);
    // rep.printVec("f1 pos",fcl[1]->pos,2);
    // rep.printVec("f2 pos",fcl[2]->pos,2);
    // rep.printVec("f3 pos",fcl[3]->pos,2);
    // rep.printVec("c  pos",cpos,2);

    // This function assigns interpolated Tau,dTau to fine cells.
    get_C2F_Tau(par, fcl, cpos, T);
  }     // loop over cells and assign optical depths
        // ----------------------------------------------------------------
#endif  // C2F_TAU

  // ----------------------------------------------------------------
  // if spatial order of accuracy is 1, then we have piecewise
  // constant data, so there is no interpolation to be done.
  if (par.spOOA == OA1) {
    for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
      f = (*f_iter);
      c = f->npt;
      // ----- constant data -----
      for (int v = 0; v < par.nvar; v++)
        f->Ph[v] = c->Ph[v];
      for (int v = 0; v < par.nvar; v++)
        f->P[v] = f->Ph[v];
      for (int v = 0; v < par.nvar; v++)
        f->dU[v] = 0.0;
      // ----- constant data -----
    }
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  else if (
      par.spOOA == OA2
      // || par.spOOA == OA1
  ) {
    std::vector<pion_flt> sx(par.nvar), sy(par.nvar), sz(par.nvar);
    for (int v = 0; v < par.nvar; v++)
      sx[v] = 0.0;
    for (int v = 0; v < par.nvar; v++)
      sy[v] = 0.0;
    for (int v = 0; v < par.nvar; v++)
      sz[v] = 0.0;
    //
    // Dimensions is sufficiently different that we have an if/else
    // loop for each dimension, and then do linear/bilinear/trilinear
    // interpolation as needed.
    //
    if (par.ndim == 1) {
      pion_flt c_vol = 0.0;
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        cell *f1, *f2, *c;
        f1 = (*f_iter);
        f_iter++;
        f2 = (*f_iter);
        // coarse cell properties:
        c     = f1->npt;
        c_vol = coarse->CellVolume(*c, 0);
        solver->SetSlope(*c, XX, par.nvar, sx.data(), OA2, coarse);
        interpolate_coarse2fine1D(
            par, fine, solver, c->Ph, c_vol, sx.data(), *f1, *f2);
      }
    }  // 1D

    else if (par.ndim == 2) {
#ifdef TEST_C2F
      spdlog::debug(
          "interpolating coarse to fine 2d: ncells={}", b->data.size());
#endif
      pion_flt c_vol = 0.0;

      //
      // Need to do bilinear interpolation, 4 cells at a time.
      //
      // First get the values at the 4 corners of coarse-cell c.
      // Then interpolate with bilinear interp. to the fine-cell
      // positions as 0.25*dx in each direction from the corners.
      //
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        cell *f1, *f2, *f3, *f4, *c;
        c     = (*f_iter)->npt;
        c_vol = coarse->CellVolume(*c, 0);
        solver->SetSlope(*c, XX, par.nvar, sx.data(), OA2, coarse);
        solver->SetSlope(*c, YY, par.nvar, sy.data(), OA2, coarse);
        // for (int v=0;v<par.nvar;v++) sx[v] = 0.0;
        // for (int v=0;v<par.nvar;v++) sy[v] = 0.0;
        // only do this on every second row because we update 4
        // cells at a time.
        if (!fine->NextPt(**f_iter, YP)
            || fine->NextPt(**f_iter, YP)->npt != c) {
#ifdef TEST_C2F
          if (level == 3) {
            spdlog::debug(
                "skipping for f_iter id={} c={}", (*f_iter)->id, c->id);
          }
#endif
          continue;
        }
        else {
#ifdef TEST_C2F
          if (level == 3) {
            spdlog::debug(
                "NOT skipping for f_iter id={} c={}", (*f_iter)->id, c->id);
          }
#endif
        }

        // get list of four fine cells.
        f1 = (*f_iter);
        f_iter++;
        f2 = (*f_iter);
        f3 = fine->NextPt(*f1, YP);
        f4 = fine->NextPt(*f2, YP);

        interpolate_coarse2fine2D(
            par, fine, solver, c->Ph, c->pos.data(), c_vol, sx.data(),
            sy.data(), *f1, *f2, *f3, *f4);

      }  // loop over fine cells
    }    // 2D

    else if (par.ndim == 3) {
      //
      // Need to do trilinear interpolation on 8 fine cells.
      //
#ifdef TEST_C2F
      spdlog::debug(
          "interpolating coarse to fine 2d: ncells={}", b->data.size());
#endif
      pion_flt c_vol = 0.0;
      cell *fch[8];  // fine-level cells (children)
      cell *c;       // coarse-level cell.
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        c = (*f_iter)->npt;

        // only do this on every second plane/row because we update 8
        // cells at a time.
        if (!fine->NextPt(**f_iter, YP)
            || fine->NextPt(**f_iter, YP)->npt != c) {
          continue;
        }
        if (!fine->NextPt(**f_iter, ZP)
            || fine->NextPt(**f_iter, ZP)->npt != c) {
          continue;
        }

        c_vol = coarse->CellVolume(*c, 0);
        solver->SetSlope(*c, XX, par.nvar, sx.data(), OA2, coarse);
        solver->SetSlope(*c, YY, par.nvar, sy.data(), OA2, coarse);
        solver->SetSlope(*c, ZZ, par.nvar, sz.data(), OA2, coarse);

        // get list of 8 fine cells.
        fch[0] = (*f_iter);
        f_iter++;
        fch[1] = (*f_iter);
        fch[2] = fine->NextPt(*fch[0], YP);
        fch[3] = fine->NextPt(*fch[1], YP);
        fch[4] = fine->NextPt(*fch[0], ZP);
        fch[5] = fine->NextPt(*fch[1], ZP);
        fch[6] = fine->NextPt(*fch[2], ZP);
        fch[7] = fine->NextPt(*fch[3], ZP);

#ifdef TEST_C2F
        for (int v = 0; v < 8; v++) {
          if (fch[v] == 0) spdlog::error("{}: {}", "fine-cell list", v);
        }
#endif

        interpolate_coarse2fine3D(
            par, fine, solver, c->Ph, c->pos.data(), c_vol, sx.data(),
            sy.data(), sz.data(), fch);

      }  // loop over fine cells
    }    // 3D
  }      // 2nd-order accuracy
  // ----------------------------------------------------------------

  return 0;
}

// ##################################################################
// ##################################################################

void NG_coarse_to_fine_bc::interpolate_coarse2fine1D(
    class SimParams &par,          ///< simulation parameters
    class GridBaseClass *fine,     ///< pointer to fine grid
    class FV_solver_base *solver,  ///< pointer to equations
    const pion_flt *P,             ///< state vector of coarse cell.
    const pion_flt c_vol,          ///< volume of coarse cell.
    pion_flt *sx,                  ///< dP/dx in coarse cell.
    cell &f1,                      ///< pointer to first fine cell  (XN)
    cell &f2                       ///< pointer to second fine cell (XP)
)
{
  std::vector<double> fU(par.nvar), f1U(par.nvar), f2U(par.nvar), cU(par.nvar);
  double f_vol[2];
  double dx = fine->DX();  // dx
  // Set fine-cell values based on slope, linear interpolation to cell centre
  for (int v = 0; v < par.nvar; v++)
    f1.P[v] = P[v] - 0.5 * dx * sx[v];
  for (int v = 0; v < par.nvar; v++)
    f2.P[v] = P[v] + 0.5 * dx * sx[v];

  // Check mass/momentum/energy conservation between coarse and fine levels
  solver->PtoU(f1.P.data(), f1U.data(), par.gamma);
  solver->PtoU(f2.P.data(), f2U.data(), par.gamma);
  f_vol[0] = fine->CellVolume(f1, 0);
  f_vol[1] = fine->CellVolume(f2, 0);

  for (int v = 0; v < par.nvar; v++)
    fU[v] = f1U[v] * f_vol[0] + f2U[v] * f_vol[1];
  // compare with coarse cell.
  solver->PtoU(P, cU.data(), par.gamma);
  for (int v = 0; v < par.nvar; v++)
    cU[v] *= c_vol;

#ifdef TEST_C2F
  spdlog::debug("1D coarse : {}", cU);
  spdlog::debug("1D fine   : {}", fU);
#endif

  // scale f1U, f2U by ratio of coarse to fine energy.
  // scale fine conserved vec by adding the difference between
  // conserved quantities on the fine and coarse grids.
  for (int v = 0; v < par.nvar; v++)
    cU[v] = 0.5 * (cU[v] - fU[v]) / c_vol;
  for (int v = 0; v < par.nvar; v++)
    f1U[v] += cU[v];
  for (int v = 0; v < par.nvar; v++)
    f2U[v] += cU[v];
#ifdef TEST_C2F
  for (int v = 0; v < par.nvar; v++)
    fU[v] = f1U[v] + f2U[v];
  spdlog::debug("1D fine 2 : {}", fU);
#endif
  solver->UtoP(f1U.data(), f1.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f1.P[v] = f1.Ph[v];
  solver->UtoP(f2U.data(), f2.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f2.P[v] = f2.Ph[v];

  for (int v = 0; v < par.nvar; v++)
    f1.dU[v] = 0.0;
  for (int v = 0; v < par.nvar; v++)
    f2.dU[v] = 0.0;

#ifdef DEBUG_NG
  for (int v = 0; v < par.nvar; v++) {
    if (!isfinite(f1.P[v]) || !isfinite(f2.P[v])) {
      spdlog::error(
          "{}: {} {}", "1D C2F: fine 1,2 not finite", f1.P[v], f2.P[v]);
      exit(1);
    }
  }
#endif
}

// ##################################################################
// ##################################################################

void NG_coarse_to_fine_bc::interpolate_coarse2fine3D(
    class SimParams &par,          ///< pointer to simulation parameters
    class GridBaseClass *fine,     ///< pointer to fine grid
    class FV_solver_base *solver,  ///< pointer to equations
    const pion_flt *P,             ///< state vector of coarse cell.
    const int *,                   ///< position of coarse cell.
    const pion_flt c_vol,          ///< volume of coarse cell.
    pion_flt *sx,                  ///< dP/dx in coarse cell.
    pion_flt *sy,                  ///< dP/dy in coarse cell.
    pion_flt *sz,                  ///< dP/dz in coarse cell.
    cell **fch                     ///< pointer to array of 8 fine cells
)
{
  double **fU = 0;
  std::vector<double> Utot(par.nvar), cU(par.nvar);
  fU = mem.myalloc(fU, 8);
  for (int i = 0; i < 8; i++)
    fU[i] = mem.myalloc(fU[i], par.nvar);
  double dxo2 = 0.5 * fine->DX();  // dx
  double f_vol[8];
  // int idx = 2*fine->idx(); // idx of coarse cell
  // double f_psi[4];
  //
  // Need to do trilinear interpolation, 4 cells at a time.
  // use slopes in each direction to get cell-centred values in each
  // of the fine cells.
  // These positions are half-way to the coarse cell-corners.
  //
  for (int v = 0; v < par.nvar; v++)
    sx[v] *= dxo2;  // centres of fine cell
  for (int v = 0; v < par.nvar; v++)
    sy[v] *= dxo2;
  for (int v = 0; v < par.nvar; v++)
    sz[v] *= dxo2;

  for (int v = 0; v < par.nvar; v++)
    fU[0][v] = P[v] - sx[v] - sy[v] - sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[1][v] = P[v] + sx[v] - sy[v] - sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[2][v] = P[v] - sx[v] + sy[v] - sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[3][v] = P[v] + sx[v] + sy[v] - sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[4][v] = P[v] - sx[v] - sy[v] + sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[5][v] = P[v] + sx[v] - sy[v] + sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[6][v] = P[v] - sx[v] + sy[v] + sz[v];
  for (int v = 0; v < par.nvar; v++)
    fU[7][v] = P[v] + sx[v] + sy[v] + sz[v];

  for (int i = 0; i < 8; i++)
    for (int v = 0; v < par.nvar; v++)
      fch[i]->P[v] = fU[i][v];

  // Need to check mass/momentum/energy conservation between
  // coarse and fine levels.  Re-use fU[][] array for this
  //
  for (int i = 0; i < 8; i++)
    solver->PtoU(fch[i]->P.data(), fU[i], par.gamma);
  for (int i = 0; i < 8; i++)
    f_vol[i] = fine->CellVolume(*fch[i], 0);

  for (int v = 0; v < par.nvar; v++)
    Utot[v] = 0.0;
  for (int i = 0; i < 8; i++) {
    for (int v = 0; v < par.nvar; v++)
      Utot[v] += fU[i][v] * f_vol[i];
  }
  // compare with coarse cell.
  solver->PtoU(P, cU.data(), par.gamma);
  for (int v = 0; v < par.nvar; v++)
    cU[v] *= c_vol;

#ifdef DEBUG_NG
  for (int i = 0; i < 8; i++) {
    for (int v = 0; v < par.nvar; v++) {
      if (!isfinite(fU[i][v])) {
        spdlog::debug("error in 3D C2F interpolation: i={}, v={}", i, v);
        spdlog::debug("Unscaled fine : {}", fU[i][v]);
      }
    }
  }
#endif

  // scale fine conserved vec by adding the difference between
  // conserved quantities on the fine and coarse grids.
  for (int v = 0; v < par.nvar; v++)
    cU[v] = 0.125 * (cU[v] - Utot[v]) / c_vol;
  for (int i = 0; i < 8; i++) {
    for (int v = 0; v < par.nvar; v++)
      fU[i][v] += cU[v];
  }

  // put scaled conserved variable vectors back into fine cells
  for (int i = 0; i < 8; i++) {
    solver->UtoP(fU[i], fch[i]->Ph, par.EP.MinTemperature, par.gamma);
    for (int v = 0; v < par.nvar; v++)
      fch[i]->P[v] = fch[i]->Ph[v];
  }

  for (int i = 0; i < 8; i++)
    for (int v = 0; v < par.nvar; v++)
      fch[i]->dU[v] = 0.0;

#ifdef DEBUG_NG
  for (int i = 0; i < 8; i++) {
    for (int v = 0; v < par.nvar; v++) {
      if (!isfinite(fU[i][v])) {
        spdlog::debug("error in 3D C2F interpolation: i={}, v={}", i, v);
        spdlog::error("{}: {}", "C2F fine cell not finite", fU[i][v], par.nvar);
      }
    }
  }
#endif

  for (int i = 0; i < 8; i++)
    fU[i] = mem.myfree(fU[i]);
  mem.myfree(fU);
}

// ##################################################################
// ##################################################################

void NG_coarse_to_fine_bc::interpolate_coarse2fine2D(
    class SimParams &par,          ///< pointer to simulation parameters
    class GridBaseClass *fine,     ///< pointer to fine grid
    class FV_solver_base *solver,  ///< pointer to equations
    const pion_flt *P,             ///< state vector of coarse cell.
    const int *,                   ///< position of coarse cell.
    const pion_flt c_vol,          ///< volume of coarse cell.
    pion_flt *sx,                  ///< dP/dx in coarse cell.
    pion_flt *sy,                  ///< dP/dy in coarse cell.
    cell &f1,                      ///< pointer to first fine cell  (XN,YN)
    cell &f2,                      ///< pointer to second fine cell (XP,YN)
    cell &f3,                      ///< pointer to third fine cell  (XN,YP)
    cell &f4                       ///< pointer to fourth fine cell (XP,YP)
)
{
  std::vector<double> fU(par.nvar), f1U(par.nvar), f2U(par.nvar);
  std::vector<double> f3U(par.nvar), f4U(par.nvar), cU(par.nvar);
  double dxo2 = 0.5 * fine->DX();  // dx
  double f_vol[4];

  //
  // Need to do bilinear interpolation, 4 cells at a time.
  // use slopes in each direction to get corner values for the
  // coarse cell.
  //
  for (int v = 0; v < par.nvar; v++)
    sx[v] *= dxo2;  // fine (dx/2)
  for (int v = 0; v < par.nvar; v++)
    sy[v] *= dxo2;  // fine (dx/2)
  for (int v = 0; v < par.nvar; v++)
    f1.P[v] = P[v] - sx[v] - sy[v];
  for (int v = 0; v < par.nvar; v++)
    f2.P[v] = P[v] + sx[v] - sy[v];
  for (int v = 0; v < par.nvar; v++)
    f3.P[v] = P[v] - sx[v] + sy[v];
  for (int v = 0; v < par.nvar; v++)
    f4.P[v] = P[v] + sx[v] + sy[v];

  // Need to check mass/momentum/energy conservation between
  // coarse and fine levels
  //
  solver->PtoU(f1.P.data(), f1U.data(), par.gamma);
  solver->PtoU(f2.P.data(), f2U.data(), par.gamma);
  solver->PtoU(f3.P.data(), f3U.data(), par.gamma);
  solver->PtoU(f4.P.data(), f4U.data(), par.gamma);
  f_vol[0] = fine->CellVolume(f1, 0);
  f_vol[1] = fine->CellVolume(f2, 0);
  f_vol[2] = fine->CellVolume(f3, 0);
  f_vol[3] = fine->CellVolume(f4, 0);

  for (int v = 0; v < par.nvar; v++)
    fU[v] = f1U[v] * f_vol[0] + f2U[v] * f_vol[1] + f3U[v] * f_vol[2]
            + f4U[v] * f_vol[3];
  // compare with coarse cell.
  solver->PtoU(P, cU.data(), par.gamma);
  for (int v = 0; v < par.nvar; v++)
    cU[v] *= c_vol;

#ifdef DEBUG_NG
  for (int v = 0; v < par.nvar; v++) {
    if (!isfinite(f1U[v]) || !isfinite(f1U[v]) || !isfinite(f3U[v])
        || !isfinite(f4U[v])) {
      spdlog::debug("Unscaled fine00 : {}", f1U);
      spdlog::debug("Unscaled fine10 : {}", f2U);
      spdlog::debug("Unscaled fine01 : {}", f3U);
      spdlog::debug("Unscaled fine11 : {}", f4U);
    }
  }
#endif

  // scale fine conserved vec by adding the difference between
  // conserved quantities on the fine and coarse grids.
  for (int v = 0; v < par.nvar; v++)
    cU[v] = 0.25 * (cU[v] - fU[v]) / c_vol;
  for (int v = 0; v < par.nvar; v++)
    f1U[v] += cU[v];
  for (int v = 0; v < par.nvar; v++)
    f2U[v] += cU[v];
  for (int v = 0; v < par.nvar; v++)
    f3U[v] += cU[v];
  for (int v = 0; v < par.nvar; v++)
    f4U[v] += cU[v];

  // put scaled conserved variable vectors back into fine cells
  solver->UtoP(f1U.data(), f1.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f1.P[v] = f1.Ph[v];
  solver->UtoP(f2U.data(), f2.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f2.P[v] = f2.Ph[v];
  solver->UtoP(f3U.data(), f3.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f3.P[v] = f3.Ph[v];
  solver->UtoP(f4U.data(), f4.Ph, par.EP.MinTemperature, par.gamma);
  for (int v = 0; v < par.nvar; v++)
    f4.P[v] = f4.Ph[v];

  for (int v = 0; v < par.nvar; v++)
    f1.dU[v] = 0.0;
  for (int v = 0; v < par.nvar; v++)
    f2.dU[v] = 0.0;
  for (int v = 0; v < par.nvar; v++)
    f3.dU[v] = 0.0;
  for (int v = 0; v < par.nvar; v++)
    f4.dU[v] = 0.0;

#ifdef DEBUG_NG
  for (int v = 0; v < par.nvar; v++) {
    if (!isfinite(f3.P[v]) || !isfinite(f4.P[v])) {
      spdlog::error(
          "{}: {} {}", "2D C2F fine 3,4 not finite", f3.P[v], f4.P[v]);
      exit(1);
    }
  }
#endif
}

// ##################################################################
// ##################################################################

void NG_coarse_to_fine_bc::get_C2F_Tau(
    class SimParams &par,    ///< pointer to simulation parameters
    std::vector<cell *> &c,  ///< list of cells
    const pion_flt *cpos,    ///< centre of coarse cell.
    double *T                ///< coarse-cell optical depths
)
{
  // data needed for getting coarse cell Taus onto coarse grid.
  double Tau[MAX_TAU], dTau[MAX_TAU];
  class cell *f0, *f1, *f2, *f3;
  struct rad_src_info *s;
  double diffx, diffy;

  if (par.ndim == 1) {
    f0 = c[0];
    f1 = c[1];

    for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      for (int iT = 0; iT < s->NTau; iT++)
        Tau[iT] = T[C2F_tauoff[isrc] + iT];
      for (int iT = 0; iT < s->NTau; iT++)
        dTau[iT] = T[C2F_tauoff[isrc] + s->NTau + iT];
      // column from source through cell depends on
      // which direction is to the source
      if (s->pos[XX] > cpos[XX]) {
        // f0 has same column as coarse, f1 has 0.5*dTau subtracted.
        // Both cells have dTau divided by 2.
        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
      }
      else {
        // Other way around.  Now f0 is closer to source, so we
        // subtract 0.5*dTau from it.
        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
      }
    }  // loop over sources
  }    // if 1D

  else if (par.ndim == 2) {
#ifdef RT_TESTING
    spdlog::info("2D C2F RT Routine");
#endif
    f0 = c[0];
    f1 = c[1];
    f2 = c[2];
    f3 = c[3];
#ifdef RT_TESTING
    spdlog::debug(
        "f: [{},{}], c: [{},{}] : ", f0->pos[XX], f0->pos[YY], cpos[XX],
        cpos[YY]);
#endif
    for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      for (int iT = 0; iT < s->NTau; iT++)
        Tau[iT] = T[C2F_tauoff[isrc] + iT];
      for (int iT = 0; iT < s->NTau; iT++)
        dTau[iT] = T[C2F_tauoff[isrc] + s->NTau + iT];
      diffx = s->pos[XX] - cpos[XX];
      diffy = s->pos[YY] - cpos[YY];
#ifdef RT_TESTING
      spdlog::debug(
          "diffxy = {}  {}\nt1={}, t2={}, t3={}, t4={}", diffx, diffy, *Tau1,
          *Tau2, *Tau3, *Tau4);
#endif
      // column from source through cell depends on
      // which direction is to the source.
      //
      /// We assume we are in the far-field limit, so that the
      /// rays through the coarse and 4 fine cells are essentially
      /// parallel to each other, and the ray segments throug the
      /// fine cells have the same length.
      //
      if (diffx > 0 && fabs(diffx) >= fabs(diffy)) {
        // Source in Q1 coming from dir XP (-45 < theta < 45 deg)
        // f2,f0 have Tau same as coarse cell, and f1,f3 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
        CI.set_col(*f2, s->id, Tau);
        CI.set_cell_col(*f2, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
        CI.set_col(*f3, s->id, Tau);
        CI.set_cell_col(*f3, s->id, dTau);
      }
      else if (diffy > 0 && fabs(diffx) < fabs(diffy)) {
        // source in Q2, coming from dir YP (45 < theta < 135 deg)
        // f0,f1 have Tau same as coarse cell, and f2,f3 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f2, s->id, Tau);
        CI.set_cell_col(*f2, s->id, dTau);
        CI.set_col(*f3, s->id, Tau);
        CI.set_cell_col(*f3, s->id, dTau);
      }
      else if (diffx < 0 && fabs(diffx) >= fabs(diffy)) {
        // source in Q3, coming from XN (135 < theta < 225 deg)
        // f1,f3 have Tau same as coarse cell, and f0,f2 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        // cout <<"in Q3: pos = ";
        // rep.printVec("cpos",cpos,par.ndim);
        // rep.printVec("Tau",Tau,s->NTau);
        // rep.printVec("dTau",dTau,s->NTau);

        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
        CI.set_col(*f3, s->id, Tau);
        CI.set_cell_col(*f3, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
        CI.set_col(*f2, s->id, Tau);
        CI.set_cell_col(*f2, s->id, dTau);

        // double fpos[MAX_DIM]; CI.get_dpos(f0,fpos);
        // if (fpos[YY]<-1.9e18) {
        //  rep.printVec("fpos",fpos,par.ndim);
        //  CI.print_cell(f0);
        //  CI.print_cell(f1);
        //  CI.print_cell(f2);
        //  CI.print_cell(f3);
        //}
      }
      else {
        // source in Q4, coming from YN (225 < theta < 315 deg)
        // f2,f3 have Tau same as coarse cell, and f0,f1 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT = 0; iT < s->NTau; iT++)
          dTau[iT] *= 0.5;
        CI.set_col(*f2, s->id, Tau);
        CI.set_cell_col(*f2, s->id, dTau);
        CI.set_col(*f3, s->id, Tau);
        CI.set_cell_col(*f3, s->id, dTau);
        for (int iT = 0; iT < s->NTau; iT++)
          Tau[iT] -= dTau[iT];
        CI.set_col(*f0, s->id, Tau);
        CI.set_cell_col(*f0, s->id, dTau);
        CI.set_col(*f1, s->id, Tau);
        CI.set_cell_col(*f1, s->id, dTau);
      }
#ifdef RT_TESTING
      spdlog::debug("  tc={}", *tmp);
#endif
    }  // loop over sources.
  }    // if 2D

  else {
    spdlog::error("{}: {}", "3D C2F RT not implemented in NG grid", par.ndim);
  }  // if 3D
}  // if there is a finer level

// ##################################################################
// ##################################################################
