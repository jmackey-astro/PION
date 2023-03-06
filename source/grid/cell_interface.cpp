///
/// \file   cell_interface.cpp
/// \author Jonathan Mackey
/// \date   12.11.2010
///
/// Purpose: Define a set of routines for accessing cell data; these
/// would be too memory intensive if each cell had to have a copy of
/// the functions.
///
/// History: Used to be in global.cc (up to SVN rev. 236).
///
/// Modifications:
/// - 2010.11.12 JM: Added support for H-correction speeds.  Changed
///   'col' to be an access function 'monochromatic_tau(cell *c)'.
/// - 2010.11.15/19 JM: Debugged.
/// - 2010.12.30 JM: Added DivV() set/get functions.
/// - 2011.02.17 JM: Enhanced optical depth storage for multiple sources.
/// - 2011.02.24 JM: Simplified optical depth storage a little. (tidied up
/// 02.25).
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2011.03.21 JM: Updated optical-depth info, multiple variables per source.
/// - 2011.04.18 JM: Added storage for dS, path length through a cell for
/// raytracing.
/// - 2011.10.17 JM: Updated RT storage.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.08.20 JM: Moved raytracing set/get functions to header and
///    made them inline.
/// - 2013.09.20 JM: Changed initialisation of unsigned ints to zero.
/// - 2015.01.10 JM: New include statements for new file structure.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).  Changed int_converter depending on pion_flt.
/// - 2017.12.09 JM: Added ndim, nvar to get rid of SimPM references.

#include "cell_interface.h"
#include "constants.h"
#include "raytracing/rad_src_data.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

using namespace std;

#ifdef COUNT_ENERGETICS
struct energetics *GLOBAL_CE = 0;
///< for tracking rates in microphysics/raytracing.
#endif

/// global class for accessing cell data, positions, neigbours.
class cell_interface CI;

/************************* CELL INTERFACE ***********************/

// ##################################################################
// ##################################################################



cell_interface::cell_interface()
{
  minimal_cell = false;
  dxo2         = -HUGEVALUE;
  ndim         = -1;
  nvar         = -1;
  //
  // NG grid parameters
  //
  nlevels = 1;

  if (sizeof(pion_flt) == sizeof(double)) {
    int_converter = ONE_PLUS_EPS;
  }
  else {
    // this means I can have grids with up to 5e5 zones before it
    // fails with floating point variables...
    int_converter = 1.0 + 2.0e-6;
  }
  cell_size_int_units = 2;

  /// This must be set to true to create a cell.
  have_setup_extra_data = false;
  /// Flag: 0=not doing RT,  N>=1 = N tau values.
  using_RT = 0;
  /// Flag: 0=no Hcorr, N=need N variables (Lapidus=1,Hcorr=Ndim).
  using_Hcorr        = 0;
  using_DivV         = 0;
  using_GradP        = 0;
  using_compton_cool = 0;
  using_wind_acc     = 0;
  star_data.clear();
  /// Size of extra_data array (can be zero).
  N_extra_data = 0;
  //
  // index arrays initialise to zero.
  //
  iDivV  = 0;
  iGradP = 0;
  /// indices of Hcorrection values in extra_data [XX,YY,ZZ].
  fill(iHcorr.begin(), iHcorr.end(), 0);
}



// ##################################################################
// ##################################################################



void cell_interface::set_minimal_cell_data()
{
  minimal_cell = true;
}



// ##################################################################
// ##################################################################



void cell_interface::unset_minimal_cell_data()
{
  minimal_cell = false;
}



// ##################################################################
// ##################################################################



void cell_interface::set_ndim(const int nd)
{
  ndim = nd;
}



// ##################################################################
// ##################################################################



void cell_interface::set_nvar(const int nv)
{
  nvar = nv;
}



// ##################################################################
// ##################################################################



void cell_interface::set_xmin(const std::array<double, MAX_DIM> &xm)
{
  copy(xm.begin(), xm.end(), xmin.begin());
}



// ##################################################################
// ##################################################################



//
// Set variables for extra_data based on what we need per cell.
// Currently monochromatic radiation needs one double, the
// H-correction needs Ndim doubles, and Lapidus viscosity one double.
//
void cell_interface::setup_extra_data(
    const struct rad_sources &rsi,      ///< Flag for ray-tracing
    const int hc_flag,                  ///< Flag for H-correction
    const int dv_flag,                  ///< Flag for Div(V).
    const int gp_flag,                  ///< Flag for |grad(P)|.
    const struct which_physics &ep,     ///< some extra parameters
    const struct stellarwind_list &swp  ///< list of wind sources.
)
{
  spdlog::debug("cell_interface::setup_extra_data()");
  //
  // Set up a 1D array for all of the extra data that a cell needs,
  // and set indices to access the required elements.
  //
  // Start with no extra data:
  //
  N_extra_data = 0;

  //
  // Now add ray-tracing optical depth variable(s).
  //
  if (rsi.Nsources > 0) {
    //
    // A number of sources, each with a number of variables determined by
    // the source type.
    //
    using_RT = 1;
    //
    // Each source, s, numbered from 0 to Nsources-1, needs a number
    // of variables:
    // - NTau[s]:  the number of Tau and DTau variables for source s.
    // - iTau[s]:  index of first Tau variable in s (others follow).
    // - iDTau[s]: index of first DTau variable in s (others follow).
    // - iVsh[s]:  index of Vshell variable in extra_data[] for s.
    // - idS[s]:   index of dS variable in extra_data[] for s.
    //
    NTau.resize(rsi.Nsources);
    iTau.resize(rsi.Nsources);
    iDTau.resize(rsi.Nsources);
    iVsh.resize(rsi.Nsources);
    idS.resize(rsi.Nsources);

    for (int s = 0; s < rsi.Nsources; ++s) {
      //
      // Number of quantities traced from source:
      //
      NTau[s] = rsi.sources[s].NTau;
      //
      // New update with rates: need NTau[s] vars for Tau and DTau,
      // one var for Vshell and dS.
      //
      iTau[s] = N_extra_data;
      N_extra_data += NTau[s];
      iDTau[s] = N_extra_data;
      N_extra_data += NTau[s];
      iVsh[s] = N_extra_data++;
      idS[s]  = N_extra_data++;
    }  // loop over radiation sources.
    spdlog::debug("Adding RT: N={}", N_extra_data);
  }

  //
  // Now add H-correction/Lapidus viscosity variables
  //
  if (hc_flag) {
    using_Hcorr = hc_flag;
    if (hc_flag > MAX_DIM)
      spdlog::error(
          "{}: {}", "Can't ask for more than MAX_DIM H-corr variables",
          hc_flag);

    for (int v = 0; v < hc_flag; v++) {
      iHcorr[v] = N_extra_data++;
    }
    spdlog::debug("Adding HCORR: N={}", N_extra_data);
  }

  //
  // DivV just needs one extra variable.
  //
  if (dv_flag) {
    using_DivV = dv_flag;
    iDivV      = N_extra_data++;
    spdlog::debug("Adding DIVV: N={}", N_extra_data);
  }

  //
  // GradP just needs one extra variable.
  //
  if (gp_flag) {
    using_GradP = dv_flag;
    iGradP      = N_extra_data++;
    spdlog::debug("Adding GRADP: N={}", N_extra_data);
  }

  //
  // Compton cooling needs a distance and a radiation flux
  //
  if (ep.compton_cool != 0) {
    spdlog::info("CI-xd: compton cooling, Nsrc: {}", swp.Nsources);
    using_compton_cool = 1;
    star_data.resize(swp.Nsources);
    for (int v = 0; v < swp.Nsources; v++) {
      star_data[v].star_dist = N_extra_data++;
      star_data[v].star_urad = N_extra_data++;
    }
  }
  else {
    using_compton_cool = 0;
  }

  //
  // wind acceleration also needs distance, and Ndim accelerations
  //
  if (ep.wind_acceleration != 0) {
    using_wind_acc = 1;
    if (!using_compton_cool) {
      // need to setup the distance if not already done.
      spdlog::info("CI-xd: wind acc: nsrc: {}, ndim {}", swp.Nsources, ndim);
      star_data.resize(swp.Nsources);
      for (int v = 0; v < swp.Nsources; v++) {
        star_data[v].star_dist = N_extra_data++;
      }
    }
    // now add accelerations
    for (int v = 0; v < swp.Nsources; v++) {
      star_data[v].wind_acc.resize(ndim);
      for (int d = 0; d < ndim; d++)
        star_data[v].wind_acc[d] = N_extra_data++;
    }
  }
  else {
    using_wind_acc = 0;
  }

  // finished
  have_setup_extra_data = true;
  spdlog::debug("cell_interface::setup_extra_data() finished");
}


// ##################################################################
// ##################################################################


double cell_interface::get_wind_acceleration_el(
    const cell &c,     ///< cell pointer
    const int src,     ///< star id.
    const int element  ///< which axis to return acceleration for
)
{
  // spdlog::info("wa-el: {} {} {}, sd[] {} : DATA
  // {:12.6e}",c.id,src,element,star_data.size(),c.extra_data[star_data[src].wind_acc[element]]);
  return c.extra_data[star_data[src].wind_acc[element]];
}

// ##################################################################
// ##################################################################



// returns true if using minimal cells.
bool cell_interface::query_minimal_cells()
{
  return minimal_cell;
}



// ##################################################################
// ##################################################################



void cell_interface::set_cell_pointers(
    cell &c,                    ///< cell to add pointers to
    std::vector<double> &d_Ph,  ///< data array
    size_t ix_Ph,               ///< index of first free element in array
    std::vector<double> &d_xd,  ///< data array
    size_t ix_xd                ///< index of first free element in array
)
{
  if (!have_setup_extra_data)
    spdlog::error(
        "{}: {}", "Setup extra data before calling new_cell", using_RT);
  if (dxo2 < 0.0) spdlog::error("{}: {}", "Cell Interface: set dx", dxo2);
  if (ndim < 0) spdlog::error("{}: {}", "Cell Interface: set ndim", ndim);
  if (nvar < 0) spdlog::error("{}: {}", "Cell Interface: set nvar", nvar);

  c.pos.resize(ndim, 0);

  c.P.resize(nvar, 0.0);
  fill(c.ngb.begin(), c.ngb.end(), nullptr);
  c.npt      = 0;
  c.npt_all  = 0;
  c.isedge   = -999;
  c.isbd     = false;
  c.isgd     = true;
  c.isdomain = false;
  c.rt       = false;
  c.timestep = false;
  c.hll      = false;  // Flag to temporarily use HLL solver on this cell.
  fill(c.isbd_ref.begin(), c.isbd_ref.end(), false);

  //
  // If we need all the [dU,Ph] arrays, initialise them, but if we have
  // set "minimal_cell" to true, then skip them to save memory in
  // analysis code.
  //
  if (!minimal_cell) {
    c.dU.resize(nvar, 0.0);
    c.Ph = &(d_Ph[ix_Ph]);
    fill(c.Ph, c.Ph + nvar, 0.0);
  }

  if (N_extra_data > 0) {
    c.extra_data = &(d_xd[ix_xd]);
    fill(c.extra_data, c.extra_data + N_extra_data, 0.0);
  }
}



// ##################################################################
// ##################################################################



void cell_interface::delete_cell(cell &c)
{
  c.Ph = mem.myfree(c.Ph);

  if (N_extra_data >= 1) c.extra_data = mem.myfree(c.extra_data);
}



// ##################################################################
// ##################################################################



void cell_interface::set_pos(
    cell &c,  ///< pointer to cell
    const std::array<double, MAX_DIM>
        &p_in  ///< double array of size ndim, containing cell position.
)
{
  //
  // Set position integer according to Xmin+i*DX/2=x
  //
  for (int v = 0; v < ndim; v++) {
    c.pos[v] = static_cast<int>(int_converter * ((p_in[v] - xmin[v]) / dxo2));
  }
#ifdef DEBUG
  spdlog::debug("int-pos from double : {}", c->pos);
#endif
}



// ##################################################################
// ##################################################################



void cell_interface::set_pos(
    cell &c,  ///< pointer to cell
    const std::array<int, MAX_DIM>
        &p_in  ///< integer array of size ndim, containing cell position.
)
{
  //
  // Set position integer according to Xmin+i*DX/2=x
  // This function assumes a clever person has set p_in to have these values!
  // If not, the code may fail catastrophically.
  //
  copy(p_in.begin(), p_in.end(), c.pos.begin());
}



// ##################################################################
// ##################################################################



void cell_interface::get_dpos(
    const cell &c,                      ///< pointer to cell
    std::array<double, MAX_DIM> &p_out  ///< array to write position into.
)
{
  for (int v = 0; v < ndim; v++)
    p_out[v] = xmin[v] + c.pos[v] * dxo2;
}



// ##################################################################
// ##################################################################



double cell_interface::get_dpos(
    const cell &c,  ///< pointer to cell
    const int v     ///< element of position vector we want
)
{
  return xmin[v] + c.pos[v] * dxo2;
}



// ##################################################################
// ##################################################################



void cell_interface::get_ipos(
    const cell &c,  ///< pointer to cell
    int *ipos_out   ///< array to write integer position into.
)
{
  for (int v = 0; v < ndim; v++)
    ipos_out[v] = c.pos[v];
}



// ##################################################################
// ##################################################################



int cell_interface::get_ipos(
    const cell &c,  ///< pointer to cell
    const int v     ///< element of position we want.
)
{
  return c.pos[v];
}



// ##################################################################
// ##################################################################



void cell_interface::get_ipos_vec(
    const std::array<double, MAX_DIM> &p_in,  ///< physical position (input)
    std::array<int, MAX_DIM> &p_out           ///< integer position (output)
)
{
  if (dxo2 < 0.0)
    spdlog::error(
        "{}: {}", "set up grid before trying to get integer positions!!!",
        dxo2);
  for (int v = 0; v < ndim; v++) {
    if (fabs(p_in[v]) > VERY_LARGE_VALUE)
      p_out[v] = -1234567;
    else
      p_out[v] = static_cast<int>(int_converter * ((p_in[v] - xmin[v]) / dxo2));
  }
}



// ##################################################################
// ##################################################################



void cell_interface::get_ipos_as_double(
    const std::array<double, MAX_DIM> &p_in,  ///< physical position (input)
    std::array<double, MAX_DIM> &p_out        ///< integer position (output)
)
{
  if (dxo2 < 0.0)
    spdlog::error(
        "{}: {}", "set up grid before trying to get integer positions!!!",
        dxo2);
  for (int v = 0; v < ndim; v++) {
    if (fabs(p_in[v]) > VERY_LARGE_VALUE)
      p_out[v] = -VERY_LARGE_VALUE;
    else
      p_out[v] = (p_in[v] - xmin[v]) / dxo2;
  }
}



// ##################################################################
// ##################################################################



void cell_interface::get_dpos_vec(
    const std::array<int, MAX_DIM> &p_in,  ///< integer position (input)
    std::array<double, MAX_DIM> &p_out     ///< physical position (output)
)
{
  for (int v = 0; v < ndim; v++)
    p_out[v] = xmin[v] + (p_in[v]) * dxo2;
}



// ##################################################################
// ##################################################################



void cell_interface::copy_cell(const cell &c1, cell &c2)
{
  c2.pos = c1.pos;
  c2.P   = c1.P;
  c2.dU  = c1.dU;
  copy(c2.Ph, c2.Ph + nvar, c1.Ph);
  c2.F   = c1.F;
  c2.ngb = c1.ngb;
  copy(c1.extra_data, c1.extra_data + N_extra_data, c2.extra_data);
  c2.npt      = c1.npt;
  c2.npt_all  = c1.npt_all;
  c2.id       = c1.id;
  c2.isedge   = c1.isedge;
  c2.isbd     = c1.isbd;
  c2.isgd     = c1.isgd;
  c2.isdomain = c1.isdomain;
  c2.isbd_ref = c1.isbd_ref;
}



// ##################################################################
// ##################################################################



void cell_interface::print_cell(const cell &c)
{
  spdlog::info(
      "cell:\t id = {}\tisedge: {}\tisbd: {}\tisgd: {}\tisdomain: {}\tisleaf: {}\ttimestep: {}\trt: {}",
      c.id, c.isedge, c.isbd, c.isgd, c.isdomain, c.isleaf, c.timestep, c.rt);
  if (c.npt != 0)
    spdlog::info("\tnpt[id]: {}", c.npt->id);
  else
    spdlog::info("\tnpt is not addressed (last point?)");
  if (c.npt_all != 0)
    spdlog::info("\tnpt_all[id]: {}", c.npt_all->id);
  else
    spdlog::warn("\tnpt_all is not addressed (last point?)");
  if (N_extra_data > 0) {
    spdlog::info(
        "extra_data[] : {}",
        std::vector<double>(c.extra_data, c.extra_data + N_extra_data));
  }
  spdlog::info("pos[] : {}", c.pos);
  std::array<double, MAX_DIM> p = {0.0, 0.0, 0.0};
  get_dpos(c, p);
  spdlog::info("dpos[] : {}", p);
  spdlog::info("P[]   : {}", c.P);
  if (!minimal_cell) {
    spdlog::info("Ph[]  : {}", std::vector<double>(c.Ph, c.Ph + nvar));
    spdlog::info("dU[]  : {}", c.dU);
  }
  for (int i = 0; i < ndim; i++) {
    if (!c.F[i].empty()) {
      spdlog::info("F[i]  : {}", c.F[i]);
    }
  }
  // TODO spdlog::info("ngb[] : {}", std::vector<cell *>(c.ngb, c.ngb + 2 *
  // ndim));
  spdlog::info("isbd_ref[] : {}\n", c.isbd_ref);
}



// ##################################################################
// ##################################################################



// ----------------------------------------------------------------
// *** Methods for a NG grid ***
// ----------------------------------------------------------------



// ##################################################################
// ##################################################################



void cell_interface::set_nlevels(
    const double dx,  ///< dx on coarsest grid.
    const int n       ///< number of levels in NG grid.
)
{
  nlevels = n;
  n_idx.resize(n);
  n_dxo2.resize(n);
  n_dx.resize(n);

  n_idx[n - 1] = 2;  // cell diameter is 2 units on finest level.
  // each coarser level has 2x larger cells.
  for (int l = n - 2; l >= 0; l--)
    n_idx[l] = 2 * n_idx[l + 1];

  n_dx[0] = dx;
  for (int l = 1; l < n; l++)
    n_dx[l] = n_dx[l - 1] * 0.5;

  n_dxo2[0] = 0.5 * dx;
  for (int l = 1; l < n; l++)
    n_dxo2[l] = n_dxo2[l - 1] * 0.5;

  cell_diameter = n_dx[n - 1];
  dxo2          = n_dxo2[n - 1];  // refers to the finest grid now.
}



// ##################################################################
// ##################################################################

/************************* CELL INTERFACE ***********************/
