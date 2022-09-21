///
/// \file stellar_wind_BC.h
/// \author Jonathan Mackey
/// \date 2010.10.05
///
/// Header for the stellar wind boundary condition class.
/// Moved from global.h b/c it got too big.
///
/// Modifications:
///
/// - 2010.10.05 JM: Added Rstar parameter to calculate wind pressure.
///
/// - 2011.02.14 JM: Added stellar_wind_evolution class for winds with
///   evolving properties, determined by a stellar evolution model, and
///   fitted with spline functions. 02.15 JM: Debugged.
/// - 2011.11.22 JM: Added t_scalefactor parameter for stellar winds.
/// - 2013.09.06 JM: Removed integer position functions, because the
///    rounding errors created potential errors in 1D.
/// - 2015.01.10 JM: New include statements for new file structure.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.10.19 JM: Fixed wind-tracer to always use pion_flt.
/// - 2017.07.21 RK: Removed Vinf as argument for wind_source, added
///    v_rot and v_esc. Created stellar_wind_angle class for
///    angle-dependent BSG winds.

#ifndef STELLAR_WIND_BC_H
#define STELLAR_WIND_BC_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_constants.h"
#include "sim_params.h"
#include <array>


// Defines for type of wind:
#define WINDTYPE_CONSTANT 0  // spherically symmetric, constant in time
#define WINDTYPE_EVOLVING 1  // spherically symmetric, with evolution
#define WINDTYPE_ANGLE 2     // Langer et al. (1999) model implemented
#define WINDTYPE_LATDEP 3    // Modification of LGM99, original work

// ##################################################################
// ##################################################################

///
/// Struct to hold a grid cell, the wind properties for the point, and the
/// distance from the grid cell-centre to the source position.
///
struct wind_cell {
  class cell *c;  ///< pointer to cell we are interested in.
  pion_flt *p;    ///< primitive vector with wind properties for this cell.
  double dist;    ///< distance of cell centre to wind src.
  double theta;   ///< polar angle in radians of cell from source.
  double cfac;    ///< wind correction factor for density (curvilinear coords)
};

// ##################################################################
// ##################################################################

///
/// Struct for a wind source, containing the source position and the wind
/// properties, and also an STL vector list of cells whose properties are
/// set automatically to be in the freely expanding wind.
///
struct wind_source {
  struct stellarwind_params *pars;
  int ncell,   ///< number of cells in the artificially fixed region.
      Hplus,   ///< index of H+ tracer variable, if present.
      iHplus;  ///< index of H+ tracer variable, offset from 1st tracer.
  double Md0,  ///< Mdot equiv. non-rotating star (lat-dep.wind) (g/s)
      v_esc,   ///< wind escape velocity (cm/s)
      dpos_init[MAX_DIM];  /// Initial position of the source (calculated)
  bool cells_added;        ///< false until we add all cells to the source list.
  std::vector<struct wind_cell *> wcells;
};

// ##################################################################
// ##################################################################

///
/// Stellar wind class, to keep a list of stellar wind sources and to
/// help update the internal boundary conditions by maintaining a list
/// of cells and their values.
///
class stellar_wind {
public:
  stellar_wind(
      const int,            ///< ndim
      const int,            ///< nvar
      const int,            ///< ntracer
      const int,            ///< ftr
      const std::string *,  ///< List of tracer variable names.
      const int,            ///< coord_sys
      const int,            ///< eqn_type
      const double          ///< minimum temperature allowed
  );

  virtual ~stellar_wind();

  ///
  /// Add a wind source, returns source id (count from zero).
  /// Note the temperature is in Kelvin if we have a pure neutral atomic
  /// hydrogen gas, otherwise it will be modified accordingly.
  /// This is for a spherically symmetric wind.
  ///
  int add_source(
      struct stellarwind_params *  ///< pointer to wind parameters struct
  );

  ///
  /// This function is only used in a derived class.
  ///
  virtual int add_evolving_source(
      const double,                ///< current time.
      struct stellarwind_params *  ///< pointer to wind parameters struct
  )
  {
    spdlog::error("{}: {}", "Don't call add_evolving_source from here.", 99);
    return 99;
  }

  ///
  /// Add a wind source, returns source id (count from zero).
  /// Note the temperature is in Kelvin if we have a pure neutral atomic
  /// hydrogen gas, otherwise it will be modified accordingly.
  ///
  virtual int add_rotating_source(
      struct stellarwind_params *  ///< pointer to wind parameters struct
  )
  {
    spdlog::error("{}: {}", "Don't call add_rotating_source from here.", 99);
    return 99;
  }

  ///
  /// Return number of sources
  ///
  int Nsources();

  ///
  /// Add a cell to the list of boundary cells and assign it a
  /// (possibly fixed) boundary value.  Returns non-zero on error.
  ///
  int add_cell(
      class GridBaseClass *,
      const int,  ///< src id
      cell &      ///< cell to add to list.
  );

  ///
  /// Remove all cells from the list of boundary cells.
  /// Returns non-zero on error.
  ///
  int remove_cells(
      class GridBaseClass *,
      const int  ///< src id
  );

  ///
  /// Set the total number of cells for this source; function
  /// checks with internal counting to make sure we got all the
  /// cells in the internal list.
  ///
  int set_num_cells(
      const int,  ///< src id
      const int   ///< number of cells.
  );

  ///
  /// Return total number of cells for this source.
  ///
  int get_num_cells(const int  ///< src id
  );

  ///
  /// For the given cell, reset P[],Ph[] to the boundary values,
  /// or return a positive value if the cell is not in the list.
  ///
  virtual int set_cell_values(
      class GridBaseClass *,
      const int,    ///< src id
      const double  ///< simulation time
  );

  // --------------------------------------------------------------
  // Now functions for data output -- return all the values we need
  // --------------------------------------------------------------

  void get_src_posn(
      const int,                     ///< src id
      std::array<double, MAX_DIM> &  ///< position vector (output)
  );

#ifdef ANALYTIC_ORBITS
  void get_src_orbit(
      const int,  ///< src id
      int *,      ///< is star moving? 1=yes, 0=no (output)
      double *,   ///< eccentricity
      double *,
      double *,
      double *,
      double *);
#endif  // ANALYTIC_ORBITS

  void set_src_posn(
      const int,  ///< src id
      double *    ///< position vector (output)
  );

  /// allows you to set a pointer to a microphysics class
  void SetMicrophysics(
      class microphysics_base *ptr  ///< pointer to microphysics class
  )
  {
    MP = ptr;
    return;
  }

  // --------------------------------------------------------------

protected:
  const int ndim;     ///< number of dimensions on grid
  const int nvar;     ///< number of variables in state vec.
  const int ntracer;  ///< number of tracer variables in state vec.
  std::vector<std::string> tracers;
  const int ftr;                ///< first tracer index in state vec.
  const int coordsys;           ///< identifier of coordinate system
  const int eqntype;            ///< Type of equations to solve
  const double Tmin;            ///< Minimum Temperature allowed on grid.
  class microphysics_base *MP;  ///< pointer to microphysics class

  ///
  /// Set values of wind_cell reference state based on Wind-Source properties
  /// and the cell-to-source distance.
  ///
  void set_wind_cell_reference_state(
      class GridBaseClass *,
      struct wind_cell *,
      const struct wind_source *,
      const double  ///< EOS gamma
  );

  //
  // Eldridge et al. (2006, MN, 367, 186).
  // v_inf = sqrt(beta)*v_esc
  //
  double beta(const double);  ///< Teff

  std::vector<struct wind_source *> wlist;  ///< list of sources.
  int nsrc;  ///< number of sources (size of wlist vector)
};

// ##################################################################
// ##################################################################

///
/// Struct for holding evolving stellar wind data (spline arrays, Npt).
///
struct evolving_wind_data {
  int Npt;
  vector<double> time_evo;
  vector<double> M_evo;
  vector<double> L_evo;
  vector<double> R_evo;
  vector<double> Teff_evo;
  vector<double> Mdot_evo;
  vector<double> vinf_evo;
  vector<double> vrot_evo;
  vector<double> vcrt_evo;
  vector<double> vesc_evo;
  vector<double> X_H_evo;
  vector<double> X_He_evo;
  vector<double> X_C_evo;
  vector<double> X_N_evo;
  vector<double> X_O_evo;
  vector<double> X_Z_evo;
  vector<double> X_D_evo;
  // double *t;
  // double *mdot;
  // double *vinf;
  // double *Teff;
  struct wind_source *ws;
  int i_XH, i_XHe, i_XC, i_XN, i_XO, i_XZ, i_XD;
  bool is_active;  ///< Set to false initially, then to true once its time has
                   ///< come.
  double offset,   ///< sim-time minus spline-time.
      tstart,      ///< When source switches on, in sim-time units.
      tfinish,  ///< When source switches off, or sim ends, in sim-time units.
      t_next_update,  ///< Next time we need to update the source properties.
      update_freq;    ///< Time interval between calls to update wind
                      ///< properties.
};

// ##################################################################
// ##################################################################

///
/// \brief stellar_wind_evolution: Derived class for time-varying
/// stellar winds, with values set by interpolation of a table of
/// values read in from a text file.
///
/// \date 2011.02.14
/// \author Jonathan Mackey
///
class stellar_wind_evolution : virtual public stellar_wind {
public:
  ///
  /// Constructor:
  ///
  stellar_wind_evolution(
      const int,            ///< ndim
      const int,            ///< nvar
      const int,            ///< ntracer
      const int,            ///< ftr
      const std::string *,  ///< List of tracer variable names.
      const int,            ///< coord_sys
      const int,            ///< eqn_type
      const double,         ///< minimum temperature allowed
      const double,         ///< Simulation start time.
      const double          ///< Simulation finish time.
  );

  ///
  /// Destructor (Delete spline arrays, etc.)
  ///
  ~stellar_wind_evolution();

  ///
  /// Add a wind source, returns source id (count from zero).
  /// This just wraps the stellar_wind version.
  ///
  int add_source(
      struct stellarwind_params *  ///< pointer to wind parameters struct
  );

  ///
  /// New add source function: takes in the filename of the source data, and
  /// a time offset between the start of the simulation and the time in the
  /// stellar model (may need to be <0 so that wind feedback starts
  /// immediately).
  ///
  virtual int add_evolving_source(
      const double t_now,          ///< current time.
      struct stellarwind_params *  ///< pointer to wind parameters struct
  );

  ///
  /// For the given cell, reset P[],Ph[] to the boundary values,
  /// or return a positive value if the cell is not in the list.
  ///
  /// If the source is active, it checks if it needs to update the source
  /// properties. If it does, then it updates the properties, and the
  /// reference state vectors for all of the wind cells (by calling
  /// update_source()). If it is not active, and it needs to be activated, it
  /// will update properties and set it to be active.
  ///
  /// The state vectors P[] and Ph[] are set to wind values, if the source is
  /// active.
  ///
  int set_cell_values(
      class GridBaseClass *,
      const int,    ///< src id
      const double  ///< simulation time
  );

protected:
  const double sim_start;   ///< start time of simulation.
  const double sim_finish;  ///< finish time of simulation.

  ///
  /// Read data from evolution file for evolving stellar wind
  /// All columns must be in CGS units already.
  /// All columns must be present, but may not be used depending
  /// on the type of wind source.
  int read_evolution_file(
      const string,                ///< file name to read data from.
      struct evolving_wind_data *  ///< where to put the data
  );

  /// Look through list of tracers and find elements and map them to
  /// tracer indices.
  void set_element_indices(
      struct evolving_wind_data *  ///< pointer to evolving wind data
  );

  ///
  /// If it is time to update the wind properties then this function does it,
  /// updating both the wind properties and the state vectors of all of the
  /// wind cells.
  ///
  virtual int update_source(
      class GridBaseClass *,
      struct evolving_wind_data *,  ///< source to update.
      const double,                 ///< current simulation time.
      const double                  ///< EOS Gamma
  );

  /// Vink et al. (2000) mass-loss recipe for the hot side of the
  /// bistability jump.
  double Mdot_Vink_hot(
      const double,  ///< luminosity (Lsun)
      const double,  ///< mass (Msun)
      const double,  ///< T_eff (K)
      const double,  ///< Metallicity wrt solar
      const double   ///< beta of wind on hot side of BSJ
  );

  /// Vink et al. (2000) mass-loss recipe for the cool side of the
  /// bistability jump.
  double Mdot_Vink_cool(
      const double,  ///< luminosity (Lsun)
      const double,  ///< mass (Msun)
      const double,  ///< T_eff (K)
      const double,  ///< Metallicity wrt solar
      const double   ///< beta of wind on cool side of BSJ
  );

  ///  Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
  double Mdot_Nieuwenhuijzen(
      const double,  ///< luminosity (Lsun)
      const double,  ///< mass (Msun)
      const double,  ///< Radius (Rsun)
      const double   ///< Metallicity wrt solar
  );

  /// Implementation of the Brott et al. (2011) mass-loss recipe.
  /// This uses beta=2.6 for hot side of bistability jump, and 1.3
  /// for the cool side, and Nieuwenhuijzen & de Jager for RSG.
  /// Returns mass-loss rate in g/s
  double Mdot_Brott(
      const double,  ///< luminosity (Lsun)
      const double,  ///< mass (Msun)
      const double,  ///< T_eff (K)
      const double,  ///< Radius (Rsun)
      const double   ///< Metallicity wrt solar
  );

  ///
  /// List of wind sources, which may or may not be active... There is
  /// a one-to-one correspondence between sources in this list and
  /// the constant wind list wlist[].
  ///
  std::vector<struct evolving_wind_data *> wdata_evol;
};

#endif  // STELLAR_WIND_BC_H
