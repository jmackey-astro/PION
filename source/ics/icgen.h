/// \file icgen.h
///
/// Headers for the various initial conditions classes.
/// Jonathan Mackey
///
/// Modifications.\n
/// - 2009-12-18 JM: Added in the laser ablation problem.
/// - 2010.10.05 JM: Modified basic_tests::setup_uniformgrid()
///
/// - 2010.12.15 JM: Added extra ambient medium functions to
///    photoevap_multi_clumps() class.
/// - 2011.03.24 JM: Added spherically symmetric clump generator.
/// - 2012.02.07 JM: Added class for Harpreet's 1D to 2D mapping.
/// - 2012.02.25 JM: Added optional velocity vector for photoevaporated clumps.
/// - 2012.09.16 JM: Added new photoevaporating clump function+var.
/// - 2013.01.10 JM: Added new StarBench workshop test probs class.
/// - 2013.03.23 JM: Added another StarBench test.
/// - 2013.03.24 JM: Added another StarBench test.
/// - 2015.02.03 JM: changed to use IC_base class sub_domain pointer.
/// - 2016.05.02 JM: A planar ionisation-front Test

#ifndef ICGEN_H
#define ICGEN_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "ics/icgen_base.h"
#include "ics/inside_sphere.h"

#include "dataIO/readparams.h"
#ifdef PARALLEL
#include "sub_domain/sub_domain.h"
#endif  // PARALLEL

// ##################################################################
// ##################################################################

class IC_basic_tests : public ICsetup_base {
public:
  IC_basic_tests();
  ~IC_basic_tests();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int eqns;  ///< =1 for Euler equations, =2 for MHD equations.
  /** \brief Sets up a completely uniform grid, as a sanity check that the
   * integration algorithm is stable.  Noise at some level can be added if
   * needed (and probably is needed to asess stability).
   * */
  int setup_uniformgrid(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );
  /** \brief Sets up a constant state with an overdense clump, and a sine
   * wave variation in the Y-velocity across the grid. */
  int setup_sinewave_velocity();
  /** \brief Sets up constant state with an overdense clump, with everything
   * moving linearly across the grid, at some angle to the coordinate axes.
   * */
  int setup_advection();
  /** \brief Function that sets up a constant state with divB non-zero near
   * the centre of the domain, used to test divB cleaning algorithms.
   *
   * From Dedner et al. 2002, JCP, 175, 645.
   * */
  int setup_divBpeak();
  /** \brief Set up Advection of a Field Loop test. */
  int setup_FieldLoop(double  ///< Z-velocity of fluid
  );
  /** \brief Set up an Orszag-Tang Vortex initial condition file. */
  int setup_OrszagTang();
  /** \brief Set up Woodward \& Colella (1984) Double Mach Reflection. */
  int setup_DoubleMachRef();
  /** \brief Set up Kelvin Helmholtz Instability, params from Stone's
   * Code Test page. */
  int setup_KelvinHelmholtz_Stone();

  /// Setup Kelvin Helmholtz Instability, params from Frank et al. 1996,
  /// ApJ, 460, 777.
  int setup_KelvinHelmholtz();

  /// Setup Liska & Wendroff (2003) implosion test problem.
  int setup_LWImplosion();
};

// ##################################################################
// ##################################################################

///
/// Set up one of a number of blast wave test problems.
///
class IC_blastwave : public ICsetup_base {
public:
  IC_blastwave();
  ~IC_blastwave();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int eqns;         ///< set to 1 for euler eqns, 2 for mhd equations
                    ///< (including all divB cleaning methods)
  double gam,       ///< eos gamma
      bw_PG,        ///< ambient pressure.
      bw_RO,        ///< ambient density.
      bw_energy,    ///< Blast Wave energy (CGS).
      bw_nzones,    ///< Number of zones to put energy in.
      bw_blastRO,   ///< density in blast wave region.
      bw_BX,        ///< ambient B_x field.
      bw_BY,        ///< ambient B_y field.
      bw_BZ;        ///< ambient B_z field.
  double *ambient;  ///< 2nd ambient medium vector (optional).
  double *BW_tr;    ///< ambient medium tracers.
  double
      interface;  ///< optional position of interface between 2 ambient media.
  void get_amb2_params();
  /** \brief Set up a blast wave in cylindrical symmetry, centred on the
   * origin.
   */
  int setup_cyl_bw();
  /** \brief Set up a blast wave in cartesian geometry, at the centre of the
   * grid. */
  int setup_cart_bw();
  /// Set up a blast wave in spherical symmetry, centred on the origin.
  int setup_sph_bw();
  ///
  /// Set up a 1D spherically-symmetric blast-wave into a pre-calculated
  /// density field, to be read in from an initial conditions file.
  ///
  int setup_sph_bw_File(
      const string,  ///< filename to read.
      const int      ///< index of tracer variable for H+ fraction.
  );
};

// ##################################################################
// ##################################################################

class IC_shock_cloud : public ICsetup_base {
public:
  IC_shock_cloud();
  ~IC_shock_cloud();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int ndim;    ///< grid dimensionality.
  int coords;  ///< coord sys used
  int eqns;    ///< set to 1 for euler eqns, 2 for mhd equations (including all
               ///< divB cleaning methods)
  double shockpos,  ///< position of shock.
      clrad,        ///< cloud radius.
      cltr0,        ///< value of optional colour tracer for cloud (e.g.-1.0)
      cltr1,        ///< value of optional tracer for cloud (e.g.-1.0)
      gam,          ///< eos gamma.
      dratio,       ///< density ratio rho(cloud)/rho(preshock)
      pratio,       ///< pressure ratio p(cloud)/p(preshock)
      Bratio,       ///< B-field ratio B(cloud)/B(preshock)
      *preshock,    ///< preshock state vector.
      *postshock;   ///< postshock state vector.
  std::array<double, MAX_DIM> cloudcentre;  ///< centre of cloud.
  int setup_shockcloud();  ///< setup data now that we have all infor from file.
};

// ##################################################################
// ##################################################################

class IC_shocktube : public ICsetup_base {
public:
  IC_shocktube();
  ~IC_shocktube(){};
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int number,  ///< shocktube test problem number
      ndim;    ///< grid dimensionality.
  int coords;  ///< coord sys used
  int eqns;    ///< set to 1 for euler eqns, 2 for mhd equations (including all
               ///< divB cleaning methods)
  double gam,  ///< eos gamma;
      shockpos,  ///< position of shock.
      angleXY,   ///< angle shock normal makes with x-axis, in XY plane.
      angleXZ;   ///< angle shock normal makes with x-axis, in XZ plane.

  std::vector<double> preshock,  ///< preshock state vector.
      postshock;                 ///< postshock state vector.

  int assign_data(
      double *,  ///< left state.
      double *,  ///< right state.
      double     ///< shock position (in units of X-range)
  );             ///< given angles,shockpos,l,r, set data.
  int get_riemann_ics(
      int,       ///< test number to run
      double *,  ///< left state.
      double *,  ///< right state.
      double *   ///< shock position
  );             ///< get appropriate test from a list.
};

// ##################################################################
// ##################################################################

class IC_jet : public ICsetup_base {
public:
  IC_jet();
  ~IC_jet();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int ndim,    ///< grid dimensionality
      coords,  ///< coord sys used
      eqns,    ///< set to 1 for euler eqns, 2 for mhd equations (including all
               ///< divB cleaning methods)
      jetrad;  ///< jet radius in units of cellsize
  double *ambient,  ///< ambient state vector
      jdens,        ///< jet density.
      jpres,        ///< jet pressure.
      jvel,         ///< jet velocity.
      j_bax,        ///< jet axial B-field
      j_btor,       ///< jet toriodal B-field
      jtr0;         ///< value of colour tracer in jet
};

class IC_radiative_shock : public ICsetup_base {
public:
  IC_radiative_shock();
  ~IC_radiative_shock();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int eqns;     ///< 1=euler; 2=mhd
  double vsh;   ///< shock velocity (cgs)
  double rho0;  ///< initial gas density (cgs)
  double T0;    ///< initial gas temperature (K)
  double B0;    ///< optional transverse field (BY)
  double gam;   ///< eos gamma.
  /** \brief Set up a radiative shock problem with Wall B.C. (cgs units). */
  int setup_RadiativeShock();

  /** \brief Set up a radiative shock problem with outflow BC (cgs units). */
  int setup_OutflowRadiativeShock();
};

// ##################################################################
// ##################################################################

///
/// Laser ablation problem (for Turlough).
///
class IC_laser_ablation : public ICsetup_base {
public:
  IC_laser_ablation();
  ~IC_laser_ablation();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int eqns;          ///< 1=euler; 2=mhd
  double vel0;       ///< shock velocity (cgs)
  double rho0;       ///< initial gas density (cgs)
  double Pressure0;  ///< initial gas pressure (cgs)
  double Dratio;     ///< inside/outside density ratio.
  double Pratio;     ///< inside/outside pressure ratio.
  double BT0;        ///< optional transverse field (BY)
  double BX0;        ///< optional longitudinal field.
  double gam;        ///< eos gamma.
  /// Set up a laser ablation problem with Wall B.C. (cgs units) in
  /// axisymmetry.
  int setup_LaserAblationAxi();

  /// Set up a  laser ablation problem with Wall B.C. (cgs units) in 3D.
  int setup_LaserAblation3D();
};

/** \brief Photoevaporating clump problem in 2D or 3D cartesian. */
class IC_photoevaporatingclump : public ICsetup_base {
public:
  IC_photoevaporatingclump();
  ~IC_photoevaporatingclump();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int ndim;    ///< grid dimensionality.
  int coords;  ///< coord sys used
  int eqns;    ///< set to 1 for euler eqns, 2 for mhd equations (including all
               ///< divB cleaning methods)
  double clrad,      ///< cloud radius.
      *cltr,         ///< value of optional tracers for cloud.
      gam,           ///< eos gamma.
      dratio,        ///< density ratio rho(cloud)/rho(amb)
      pratio,        ///< pressure ratio p(cloud)/p(amb)
      Bratio,        ///< B-field ratio B(cloud)/B(amb)
      radial_slope,  ///< power law density profile, if needed.
      core_radius,   ///< core-radius of ISM (so rho(0)!=infty)
      *ambient;      ///< ambient state vector.
  std::array<double, MAX_DIM> cloudcentre;  ///< centre of cloud.
  int setup_pec();   ///< setup data now that we have all info from file.
  int setup_pec2();  ///< setup two PECs now that we have all info from file.
  int setup_radialprofile();     ///< setup data for a density profile out from
                                 ///< source.
  int setup_powerlaw_density();  ///< data for a slab/axi-symmetric power law
                                 ///< density profile.
  int setup_paralleltest();  ///< setup data for parallel rays with differing
                             ///< density in each ray.
  ///
  /// Setup a problem with an ISM that has a core radius and then
  /// drops off with a power-law density profile outside that, and
  /// also has a top-hat dense clump at some (other) position.
  ///
  int setup_cloud_clump();
};

// ##################################################################
// ##################################################################

// RANDOM CLUMPS BEING PHOTO-EVAPORATED BY A SOURCE

/** \brief struct to hold the data for the clump properties. */
struct clump {
  double mass;         ///< total clump mass (not used in variable mass setup).
  double overdensity;  ///< overdensity.
  double
      centre[MAX_DIM];   ///< Position of centre of clump (same units as grid).
  double size[MAX_DIM];  ///< Scale radius of clump in each direction in its
                         ///< principal axes frame.
  double ang[MAX_DIM];   ///< Rotation vector [alpha,beta,gamma] from grid axes
                         ///< to principal axes.
  double rm[MAX_DIM]
           [MAX_DIM];  ///< Rotation matrix associated with orientation[].
  double tracer_vals[MAX_NVAR];  ///< values of tracer for this clump.
  double Vel[MAX_DIM];           ///< Optional velocity of clumps w.r.t. grid.
};

// ##################################################################
// ##################################################################

/** \brief Setup a 2D or 3D cartesian domain with a number of random clumps
 * (random in position, size, orientation, overdensity), with a radiation source
 * either on or off grid, or at infinity.
 */
class IC_photevap_random_clumps : public ICsetup_base {
public:
  IC_photevap_random_clumps();
  ~IC_photevap_random_clumps();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int ndim;     ///< grid dimensionality.
  int coords;   ///< coord sys used
  int eqns;     ///< set to 1 for euler eqns, 2 for mhd equations (including all
                ///< divB cleaning methods)
  int Nclumps;  ///< Number of clumps in grid.
  int profile;  ///< flag to say what kind of radial density profile to use
                ///< (0=top-hat, 1=Gaussian)
  double min_overdensity,  ///< minimum overdensity of clumps (multiple of
                           ///< background)
      max_overdensity,     ///< maximum overdensity of clumps
      clump_mass,    ///< total mass in clumps (used in setup_perc_fixedmass()).
      min_size,      ///< minimum size of clumps (fraction of grid size)
      max_size,      ///< maximum size of clumps (fraction of grid size)
      *cltr,         ///< value of optional tracers for cloud.
      radial_slope,  ///< optional slope in background density.
      gam,           ///< eos gamma.
      *ambient,      ///< ambient state vector.
      ambdens;  ///< ambient gas density (after taking clump mass out of total
                ///< mass).
  struct clump *cl;
  int setup_perc();  ///< setup data now that we have all info from file.
  int setup_perc_fixedmass();  ///< setup clumps for fixed total mass on grid.
  int clumps_random_setup();   ///< setup clumps (original version, random
                               ///< total mass).
  int clumps_random_setup_fixedmass();  ///< setup clumps (fixed total mass).
#ifdef PARALLEL
  int clumps_random_setup_pllel();  ///< Setup clumps in parallel (orig.
                                    ///< variable total mass).
  int clumps_random_setup_pllel_fixedmass();  ///< setup clumps in parallel
                                              ///< (fixed total mass).
#endif                                        // PARALLEL
  int clumps_random_set_dens(cell &           ///< cell to operate on.
  );
  void print_clump(struct clump *);
  double
  random_frac();  ///< returns random value on [0,1], assuming seed is set.
};

struct random_clump_data {
  bool used;    ///< set to true so destructor knows to delete dynamic memory.
  int Nclumps;  ///< Number of random clumps in grid.
  int profile;  ///< flag to say what kind of radial density profile to use
                ///< (0=top-hat, 1=Gaussian)
  int random_seed;   ///< Seed to use for generating clumps.
  double *border;    ///< 2*ndim doubles saying what fraction of the domain in
                     ///< each direction has no clumps.
  double density,    ///< mean number density of gas to go into clumps (within
                     ///< clumpy region).
      total_mass,    ///< total mass to draw clumps from.
      min_mass,      ///< min. fraction of total mass that a clump can have (for
                     ///< fixed mass clumps).
      max_mass,      ///< max. fraction of total mass that a clump can have (for
                     ///< fixed mass clumps).
      min_size,      ///< min fraction of Y-Range for clump radius.
      max_size;      ///< max fraction of Y-Range for clump radius.
  struct clump *cl;  ///< pointer to info for each of the N random clumps.
};

// ##################################################################
// ##################################################################

struct ambient_data {
  bool used;  ///< set to true so destructor knows to delete dynamic memory.
  double *ambient;      ///< Ambient gas state vector
  int radial_profile;   ///< Optional radial profile for ambient gas (radial
                        ///< out from star).
  double cloudradius;   ///< Radius of core of cloud with given radial profile.
  bool xscale;          ///< if scale density in x-direction via power law?
  double xscale_x0;     ///< rho(x) = rho_0 * [(x-x_0)/l_0)]^alpha
  double xscale_l0;     ///< rho(x) = rho_0 * [(x-x_0)/l_0)]^alpha
  double xscale_alpha;  ///< rho(x) = rho_0 * [(x-x_0)/l_0)]^alpha
};

// ##################################################################
// ##################################################################

struct strategic_clump_data {
  bool used;    ///< set to true so destructor knows to delete dynamic memory.
  int Nclumps;  ///< Number of strategically placed clumps.
  struct clump *cl;  ///< Pointer to list of clumps.
  int profile;       ///< What kind of radial density profile to use (0=top-hat,
                     ///< 1=Gaussian)
};

// ##################################################################
// ##################################################################

/** \brief Setup a 2D or 3D cartesian domain with a number of random clumps
 * (random in position, size, orientation, overdensity), with a radiation source
 * either on or off grid, or at infinity, and a number of strategically
 * positioned clumps with specified properties.
 */
class IC_photevap_multi_clumps : public ICsetup_base {
public:
  IC_photevap_multi_clumps();
  ~IC_photevap_multi_clumps();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int ndim;    ///< grid dimensionality.
  int coords;  ///< coord sys used
  int eqns;    ///< set to 1 for euler eqns, 2 for mhd equations (including all
               ///< divB cleaning methods)
  double ambdens;     ///< ambient gas density (after taking clump mass out of
                      ///< total mass).
  double ambdivider;  ///< dividing x-value between the two ambient medium
                      ///< states.
  double gamma;       ///< eos gamma;
  struct random_clump_data rc_data;
  struct strategic_clump_data sc_data;
  struct ambient_data amb_data;

  //
  // Ambient Medium Routines
  //
  int get_ambient_params(class ReadParams *, struct ambient_data *);
  int add_ambient_data_to_grid(class GridBaseClass *, struct ambient_data *);

  //
  // Ambient Medium Routines (the alternate values are only written
  // to the first 10 per cent of the x-domain).
  //
  int get_alternate_ambient_params(class ReadParams *, struct ambient_data *);
  int add_alternate_ambient_data_to_grid(
      class GridBaseClass *, struct ambient_data *);
  //
  // Generic Clump Routines
  //
  double random_frac();
  void print_clump(struct clump *);
  int clumps_set_dens(
      class cell &,    ///< cell to operate on.
      const int,       ///< number of clumps (length of clump array following)
      struct clump *,  ///< pointer to list of clump structs.
      const int        ///< clump profile.
  );

  //
  // Random Clumps Routines:
  //
  int get_random_clump_params(class ReadParams *, struct random_clump_data *);
  int rc_fixed_mass_range_params(
      class ReadParams *, struct random_clump_data *);
  int rc_fixed_number_params(class ReadParams *, struct random_clump_data *);
  /** \brief This function is hardcoded for Gaussian density distributions
   * (for total mass). It uses an analytic result to go from mass and size to
   * peak overdensity.
   */
  int rc_set_clump_properties(struct random_clump_data *);
  int add_random_clumps_to_grid(
      class GridBaseClass *, struct random_clump_data *);
  //
  // Strategically Placed Clumps Routines
  //
  int get_strategic_clump_params(
      class ReadParams *, struct strategic_clump_data *);
  int add_strategic_clumps_to_grid(
      class GridBaseClass *, struct strategic_clump_data *);
};

// ##################################################################
// ##################################################################

class IC_spherical_clump : public ICsetup_base {
public:
  IC_spherical_clump();
  ~IC_spherical_clump();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  int setup_clump();   ///< allocate data to grid cells.
  int eqns;            ///< set to 1 for euler eqns, 2 for mhd equations.
  double gam,          ///< eos gamma
      AMB_density,     ///< ambient density.
      AMB_pressure,    ///< ambient pressure
      SC_rad,          ///< cloud radius
      SC_overdensity,  ///< cloud overdensity = rho_max/rho_ambient.
      SC_BX, SC_BY,
      SC_BZ;                ///< magnetic field (optional).
  int SC_pressure_profile,  ///< cloud pressure scaling 1=isothermal,
                            ///< 2=constant pressure.
      SC_density_profile;   ///< cloud density profile: 0=top-hat,
                            ///< 1=isothermal, 2=Gaussian.
};

// ##################################################################
// ##################################################################

///
/// Read 1D data onto a 2D grid.  The 1D data should be in a text
/// file with format for each non-comment line as follows:
/// radius density pressure v_r v_theta v_phi <tracer variables>
/// Also works for 3D data, despite the name.
///
class IC_read_1Dto2D : public ICsetup_base {
public:
  IC_read_1Dto2D();
  ~IC_read_1Dto2D();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  void get_data_vals(
      std::array<double, MAX_DIM> &,  ///< Cell centre
      vector<double> &,               ///< radius vector
      vector<vector<double> > &,      ///< arrays of variable data.
      const int,                      ///< number of variables.
      double *                        ///< array for output data values at pos.
  );
  void get_3D_data_vals(
      std::array<double, MAX_DIM> &,  ///< Cell centre
      vector<double> &,               ///< radius vector
      vector<vector<double> > &,      ///< arrays of variable data.
      const int,                      ///< number of variables.
      double *                        ///< array for output data values at pos.
  );
};

// ##################################################################
// ##################################################################

#ifdef HARPREETS_CODE_EXT
#ifndef EXCLUDE_HD_MODULE
//
// Class for Harpreet to map her 1D Supernova/cloud sims onto a 2D axisymmetric
// grid.
//
class IC_HD_2D_ShockCloud : public ICsetup_base {
public:
  IC_HD_2D_ShockCloud();
  ~IC_HD_2D_ShockCloud();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  void get_data_vals(
      double *,                   ///< Cell centre
      const double,               ///< Shock position.
      vector<double> &,           ///< radius vector
      vector<vector<double> > &,  ///< arrays of variable data.
      const int,                  ///< number of variables.
      double *                    ///< array for output data values at pos.
  );
};
#endif  // don't EXCLUDE_HD_MODULE
#endif  // HARPREETS_CODE_EXT

// ##################################################################
// ##################################################################

#ifdef BBTURBULENCE_CODE_EXT
class IC_read_BBurkhart_data : public ICsetup_base {
public:
  IC_read_BBurkhart_data();
  ~IC_read_BBurkhart_data();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  void read_file(
      const string,          ///< input fits file to read
      const int,             ///< variable in state vector
      const double,          ///< optional scaling of value.
      const double,          ///< optional offset of value, after scaling.
      class GridBaseClass *  ///< pointer to grid
  );

  std::string File_RO, File_PG, File_VX, File_VY, File_VZ, File_BX, File_BY,
      File_BZ;
  double VXoffset, VYoffset, VZoffset;
};
#endif  // BBTURBULENCE_CODE_EXT

// ##################################################################
// ##################################################################

//
// Class for StarBench Workshop test problems
//
class IC_StarBench_Tests : public ICsetup_base {
public:
  IC_StarBench_Tests();
  ~IC_StarBench_Tests();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );

protected:
  ///
  /// The contact discontinuity advection tests from Dale and
  /// Wuensch.
  ///
  int setup_ContactDiscontinuity(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );

  ///
  /// The planar ionisation-front instability tests A-C from Tom
  /// Haworth.
  ///
  int setup_StarBench_IFI(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );

  ///
  /// A planar ionisation-front Test.
  ///
  int setup_StarBench_planarIF(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );

  ///
  /// The Irradiated Cloud test from Tom Haworth.
  ///
  int setup_StarBench_IrrCl(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );

  ///
  /// The shadowed region cooling/recomb test from Pascal Tremblin.
  ///
  int setup_StarBench_TremblinCooling(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );

  ///
  /// conical HII region expansion.
  ///
  int setup_StarBench_Cone(
      class ReadParams *,     ///< pointer to parameter list.
      class GridBaseClass *,  ///< pointer to grid
      string &test            ///< String with which test to run.
  );
};


///
/// This reads in a snapshot and adds a supernova to the data.
///
class IC_overwrite_snapshot : public ICsetup_base {
public:
  IC_overwrite_snapshot();
  ~IC_overwrite_snapshot();
  int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
  );
};



#endif  //  ICGEN_H
