///
/// \file sim_projection.h
///
/// - 2010-03-22 JM: Moved to MHD_ET_2010, modified to also get
///   projected magnetic field Q and U components
/// - 2010.12.13 JM: Added  NEW_STOKES_CALC ifdef to Makefile; the new
///    code in the ifdef does a different Stokes Q,U calculation and
///    replaces the projected Bx,By with values calculated from Q,U.
/// - 2013.10.15 JM: Updated to use microphysics classes to get the
///    gas temperature and n(H), and added [N II] forbidden line
///    emission.
/// - 2015.07.03 JM: updated for pion_dev: uses MCMD, SimSetup,
///    constants.h


#ifndef SIM_PROJECTION_H
#define SIM_PROJECTION_H

#define I_DENSITY    0
#define I_NEUTRAL_NH 1
#define I_VEL_LOS    2
#define I_VX         3
#define I_EMISSION   4
#define I_B_STOKESQ  5
#define I_B_STOKESU  6
#define I_ALL_SCALARS 7
#define I_BXabs 8
#define I_BYabs 9
#define I_NII6584    10

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
//#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"

#include "grid/cell_interface.h"
#include "grid/grid_base_class.h"
#include "microphysics/microphysics_base.h"

using namespace std;

/** \brief Basic class for converting between positive/negative directions and
 * associated axes, and back again.  Also cross product between axes.
 * */
class axes_directions {
public:
  virtual ~axes_directions() {}
  enum axes get_axis_from_dir(const enum direction ///< direction to convert.
			      );
  enum direction cross_product(const enum direction, ///< first direction
			       const enum direction  ///< second direction
			       );
  /** \brief Returns positive direction along an axis. */
  enum direction get_posdir(const enum axes);
  /** \brief Returns negative direction along an axis. */
  enum direction get_negdir(const enum axes);
};




// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------




/** \brief integration point along line of sight of a pixel.  This point 
 * is to be placed along a plane of cells, so the state vector at the 
 * point can be obtained by linear interpolation between two cells to the 
 * left and right of the point (in that plane).
 * */
struct point_4cellavg {
  double pos[3]; ///< position of point, in image coords.
  cell *ngb[4];  ///< pointers to four surrounding cells.
  double wt[4];  ///< weight of right cell value (wt of left=1-wt).
};


/** \brief Struct containing relevant information for a line of sight
 * integration through a pixel.
 * N.B. The 2-point average is not very good, and potentially misses half
 * a cell at each boundary, so it would be much better to do a 4pt average, 
 * I think.
 * */
struct integration_points {
  int npt; ///< number of integration points.
  double dx;      ///< interval between points (in image units)
  double dx_phys; ///< interval in physical units.
  //  double theta, tantheta, costheta; ///< angle w.r.t. normal in x-z plane.
  struct point_4cellavg *p; ///< array of integration points.
};

/** \brief Pixel struct -- contains a list of integration points, ipix, ix[2]. */
struct pixel {
  cell *inpixel;
  int ncells;
  int ipix, ix[2];

  struct integration_points int_pts; ///< list of points to integrate.
};




// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------



class coordinate_conversion : public axes_directions {
public:
  coordinate_conversion(const enum direction, ///< Line of sight direction
			const int,            ///< Angle of LOS w.r.t. los direction.
			const enum direction, ///< vertical direction, which stays in image plane.
			class GridBaseClass * ///< pointer to grid of data.
			);
  virtual ~coordinate_conversion();
  /** \brief Given that we know the grid size and the viewing angle, return
   * the number of pixels in each direction, based on the assumption that
   * pixels have the same surface area (square) as a simulation cell.
   */
  void get_npix(int * ///< 2D array to put number of pixels in each direction.
		);
  void get_image_Ipos(const int *, ///< integer position in sim coords.
		      double *     ///< converted position in image coords.
		      );
  void get_image_Dpos(const double *, ///< integer position in sim coords.
		      double *        ///< converted position in image coords.
		      );
  void get_sim_Dpos(const double *, ///< position in image coordinates.
		    double *        ///< position in sim coords (dx=2).
		    );
  /** \brief Given a pixel centre, calculates the position of the first point
   * to use for the line of sight integration, and then the interval between 
   * points and the number of points.
   */
  void set_integration_points(const double *, ///< pixel centre (x,y,zmin=0)
			      double *,       ///< dx between points (image units)
			      double *,       ///< position of nearest point.
			      int *           ///< number of points.
			      );
  bool point_in_Isim_domain(const double * /// Point in sim coords (integer)
			    );
  enum axes get_normal_axis()     {return sa[ZZ];}
  enum axes get_vertical_axis()   {return sa[YY];}
  enum axes get_horizontal_axis() {return sa[XX];}
  void get_image_axes(enum axes out[3]) {for (int i=0;i<3;i++) out[i]=sa[i];return;}
  void get_image_dir_signs(int out[3]) {for (int i=0;i<3;i++) out[i]=ss[i];return;}
  
protected:
  class GridBaseClass *gptr;
  //
  // Simulation coordinates:
  //
  /** \brief Simulation directions in image coords.  e.g. if image YP points
   * along ZN in simulation coords, then sd[YY]=ZN. */
  enum direction sd[3]; ///< direction of 3 simulation axes in image coordinates.
  enum axes sa[3]; ///< Axes associated with sim closest to each image axis.
  /** \brief Sign function.  e.g. sd[0]=XN/YN/ZN, then ss[0]=-1; else ss[0]=+1. */
  int ss[3];  ///< sign of sim direction in image direction.
  double sim_xminP[3], ///< xmin of sim:  sim coords, physical
    sim_xmaxP[3],      ///< xmax of sim:  sim coords, physical
    sim_rangeP[3],     ///< range of sim in each direction: sim coords, physical.
    sim_dxP;           ///< cell size, sim coords, physical.
  int sim_xminI[3], ///< xmin of sim: sim coords, integer.
    sim_xmaxI[3],   ///< xmax of sim: sim coords, integer.
    sim_rangeI[3],  ///< range of sim: sim coords, integer.
    sim_ncell[3],   ///< number of cells along each sim axis.
    sim_dxI;        ///< Cell size, sim coords, integer.
  //
  // Image coordinates:
  //
  int im_dx,     ///< image pixel size: image coords, image units (should =1).
    im_npix[2],  ///< number of pixels in each direction in the image plane.
    im_npixels;  ///< total number of pixels.
  int theta_deg; ///< Angle of LOS wrt z-axis, in degrees obviously!
  double theta,  ///< angle of LOS wrt z-axis, in radians.
    costheta,
    tantheta,
    sintheta;
  bool zero_angle; ///< true if we are at zero angle to a grid axis.
  double s_xmin_img[3], ///< most negative corner of sim, in image coords, so [0,0,0].
    s_xmax_img[3],      ///< most positive corner of sim, in image coords.
    s_origin_img[3];  ///< [XN,0,ZN] corner of grid, in image coords and units.
  //
  // Setup Functions.
  //
  void set_npix(); ///< set the number of pixels in each direction.
  /** \brief Called by constructor to get max/min bounds of simulation in 
   * image coordinates. */
  void set_sim_extents_in_image_coords();
};


// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------

/** \brief Set of Functions and Data for calculating the velocity at a point, and 
 * smoothing/broadening this into a binned profile array.
 *
 * We can either do Doppler Broadening at each point, based on the gas temperature,
 * or else calculate a profile along a line of sight and then smooth that by some
 * fixed velocity dispersion.  The latter is much faster, but also produces some 
 * numerical FFT artefacts (ringing, and excess signal at the ends of the profile).
 *
 * */
class point_velocity {
protected:
  //
  // Geometry variables
  //
  int vx, ///< velocity component perp. to LOS direction (contributing)
    vz,   ///< velocity componenet along LOS
    sx, ///< +1 if looking along +ve x axis; -1 otherwise
    sz; ///< +1 if looking along +ve vz axis; -1 otherwise
  double ct, ///< Cosine of angle to LOS
    st;      ///< Sine of angle to LOS
  //
  // Velocity profile variables
  //
  double 
    v_min,   ///< minimum velocity to include in profile.
    v_max,   ///< maximum velocity to include in profile.
    v_binsize; ///< size of each velocity bin in profile.
  int v_Nbins; ///< Number of velocity bins in profile.
  int broaden; ///< [0=none], 1=constant Gaussian broadening.
  double sigma; ///< width of gaussian to smooth with.
  /** \brief This computes the forward and inverse Fast Fourier Transform on
   * a 1D array of data.  It is effectively the NR version, for zero offset arrays.
   *
   * Note for the reverse transform, you need to divide each output element by nn to get
   * the actual inverse transform.
   * */
  void four1(double *,              ///< data array (Cx. data, where el. 2i=i-th real, 2i+1=i-th imag.)
	     unsigned long int nn,  ///< half the length of the array (i.e. number of Cx. values).
	     int isign              ///< =1 for forward transform, =-1 for inverse transform.
	     );

  int FT_Ng; ///< number of elements in data array for FFTing.
public:

  /** \brief constructor sets geometry info. */
  point_velocity(const int,    ///< velocity component perp. to LOS direction (contributing)
		 const int,    ///< velocity componenet along LOS
		 const int,    ///< +1 if looking along +ve vx axis; -1 otherwise
		 const int,    ///< +1 if looking along +ve vz axis; -1 otherwise
		 const double, ///< Angle to LOS (radians)
		 const double, ///< minimum velocity in range
		 const double, ///< maximum velocity in range
		 const int     ///< Number of bins.
		 );

  ~point_velocity() {}

  /** \brief Set the line broadening technique. */
  void set_broadening(
          const int,    ///< Type: 1=constant Gaussian broadening.
          const double  ///< FWHM of Gaussian to smooth profile by.
          );

  /** \brief Returns the LOS velocity profile at the point in question, subject to
   * the parameters imposed in the constructor.
   *
   * Adds the point's mass to the right velocity bin in the
   * temp_profile[] array, and if smooth==2 it will smooth this with a
   * Gaussian corresponding to the Doppler broadening from the point's
   * temperature.  This is the primary useful function call for this
   * class. */
  void get_point_v_los_profile(
          const struct point_4cellavg *, ///< point to add to profile.
          double *, ///< Array of velocity bins to put profile into.
          const int ///< index of ion fraction in prim.var.
          );

  /** \brief Returns the X-velocity componoent in a velocity profile at the
   * point in question, subject to the parameters imposed in the constructor.
   *
   * Adds the point's mass to the right velocity bin in the
   * temp_profile[] array, and if smooth==2 it will smooth this with a
   * Gaussian corresponding to the Doppler broadening from the point's
   * temperature.  This is for diagnostics rather than to mimic a real 
   * observation, since Vx is not a straightforward observation.
   */
  void get_point_VX_profile(
          const struct point_4cellavg *, ///< point to add to profile.
          double *, ///< Array of velocity bins to put profile into.
          const int ///< index of ion fraction in prim.var.
          );

  /** \brief Smooth the profile with whatever smoothing is required.  This only has any
   * effect if using fixed-width smoothing, when it does an FFT-based smoothing with a
   * Gaussian function.
   * */
  void smooth_profile_FFT(
          double * ///< Array of velocity bins to smooth.
          );

   ///
   /// Get the absorption and emission coefficients for H-alpha
   /// radiation, according to a fit to Osterbrock's data table for
   /// photoionised nebulae.
   ///
   void get_point_Halpha_params(
          const struct point_4cellavg *, ///< point in question.
          const int, ///< ifrac index in prim.vec.
          double *,  ///< absorption coefficient (/cm)
          double *   ///< emission coeff (erg/cm^3/s/sq.arcsec)
          );

   ///
   /// Get the absorption and emission coefficients for [N II] 6584AA
   /// radiation, according to a fit from Dopita (1973,A&A,29,387).
   ///
   void get_point_NII6584_params(
          const struct point_4cellavg *, ///< point in question.
          const int, ///< ifrac index in prim.vec.
          double *,  ///< absorption coefficient (/cm)
          double *   ///< emission coeff (erg/cm^3/s/sq.arcsec)
          );



 protected:

  /** \brief Get the density at the point, based on 4 cell bilinear interpolation. */
  double get_point_density(const struct point_4cellavg *);

  ///
  /// Get the H^0 number density, using SimPM.MP.H_MassFrac to
  /// convert from mass density to number density.
  ///
  double get_point_neutralH_numberdensity(
          const struct point_4cellavg *, // point
          const int // ifrac
          );

  ///
  /// Get the temperature at a point, based on 4-cell bilinear
  /// interpolation.  This requires a microphysics class.
  ///
  double get_point_temperature(
          const struct point_4cellavg *,
          const int // ifrac index in prim.vec.
          );

  /** \brief Returns the LOS velocity at the point, based on a 4 cell bilinear average. */
  double get_point_los_velocity(
          const struct point_4cellavg *
          );

  /** \brief Returns the X-velocity component at the point, based on a 4 cell bilinear average. */
  double get_point_VX(
          const struct point_4cellavg *
          );

  double get_point_perp_velocity(
          const struct point_4cellavg *
          );

  int get_velocity_bin_number(
          const double ///< point's velocity
          );
  
  /** \brief This does Doppler broadening of a single velocity point, into a 
   * temporary profile array.  Convolution is easy because the velocity is
   * a delta function.
   * */
  void broaden_profile(const struct point_4cellavg *, ///< point to add to profile.
		       double *,  ///< velocity bins.
		       const int, ///< index of i-fraction in P.V.
		       double,    ///< LOS velocity of point.
		       double     ///< Normalisation of profile.
		       );
};


// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------

//
// This is just a collection of bits of info the velocity profiling needs for each pixel.
//
struct vel_prof_stuff {
  int npix[3]; ///< number of pixels in image; 3rd element is Nbins, the number of vel. bins in profile.
  double v_min; ///< min velocity in bins.
  double v_max; ///< max velocity in bins.
  int smooth; ///< flag for what kind of smoothing to do.
  double broadening; ///< if constant smoothing, this tells us the FWHM.
};
  


/** \brief The Image Class.  This is the main class used to create
 * projected images through a simulation.
 * 
 * This sets up the coordinate system for the image, and all the lines
 * of sight for calculating projected quantities.  It also has the
 * driver function for calculating the image value for each pixel. */
class image : public coordinate_conversion {
public:
  image(const enum direction, ///< Line of sight direction
	const int,            ///< Angle of LOS w.r.t. los direction.
	const enum direction, ///< vertical direction, which stays in image plane.
	class GridBaseClass * ///< pointer to grid of data.
	);
  ~image();

  //
  // Setting cell positions in image space
  //
  void set_cell_positions_in_image();
  void delete_cell_positions(); ///< delete positions arrays; called by destructor.

  //
  // Image pixels
  //
  struct pixel *pix; ///< 1D array, pixel [i,j] = pix[Nx*j +i]

  //
  // Associating Cells and Pixels
  //
  /** \brief Tests if a cell's integer position is in a pixel. */
  bool cell_is_in_pixel(double *, ///< Cell position (in image coordinates).
			pixel  *  ///< pixel in question.
			);
  /** \brief Add cells to pixels, in a list. */
  void add_cells_to_pixels();

  //
  // Setting up points along the pixels' lines of sight
  //
  void add_integration_pts_to_pixels();

  //
  // Calculate a pixel
  //
  void calculate_pixel(struct pixel *, ///< pointer to pixel
		       const struct vel_prof_stuff *, ///< struct with info for velocity binning.
		       const int,      ///< flag for what to integrate.
		       double *,       ///< array of pixel data.
		       double *       ///< general purpose counter for stuff.
		       );

  /** \brief THESE FUNCTIONS ARE CARBON COPIES OF THOSE IN
   *  POINT_VELOCITY CLASS!!! MIGHT WANT TO CHANGE THIS SOMEDAY! */
  double get_pt_density(struct point_4cellavg *);
  /** \brief THESE FUNCTIONS ARE CARBON COPIES OF THOSE IN
   *  POINT_VELOCITY CLASS!!! MIGHT WANT TO CHANGE THIS SOMEDAY! */
  double get_pt_neutral_numberdensity(struct point_4cellavg *, ///< *pt,
				      int                      ///< ifrac
				      );
  ///
  /// Get the Stokes Q component for the perpendicular magnetic field
  ///
  double get_pt_StokesQ(struct point_4cellavg *, ///< pt
			       const int, ///< ifrac
			       const int, ///< bx index (image coords)
			       const int, ///< by index (image coords)
			       const int, ///< bz index (image coords)
			       const int, ///< sign(xx)
			       const int, ///< sign(yy)
			       const int, ///< sign(zz)
			       const double, ///< sin(theta)
			       const double  ///< cos(theta)
			       );
  ///
  /// Get the Stokes U component for the perpendicular magnetic field
  ///
  double get_pt_StokesU(struct point_4cellavg *, ///< pt
			       const int, ///< ifrac
			       const int, ///< bx index (image coords)
			       const int, ///< by index (image coords)
			       const int, ///< bz index (image coords)
			       const int, ///< sign(xx)
			       const int, ///< sign(yy)
			       const int, ///< sign(zz)
			       const double, ///< sin(theta)
			       const double  ///< cos(theta)
			       );
  ///
  /// Get the |BX| component for the perpendicular magnetic field.
  /// This returns sqrt(n_H)*Bx^2/|B|, so it's a density weighted
  /// value and also diluted by the proportion of B on the los
  /// direction
  ///
  double get_pt_BXabs(struct point_4cellavg *, ///< pt
		      const int, ///< ifrac
		      const int, ///< bx index (image coords)
		      const int, ///< by index (image coords)
		      const int, ///< bz index (image coords)
		      const int, ///< sign(xx)
		      const int, ///< sign(yy)
		      const int, ///< sign(zz)
		      const double, ///< sin(theta)
		      const double  ///< cos(theta)
		      );
  ///
  /// Get the |BY| component for the perpendicular magnetic field.
  /// This returns sqrt(n_H)*By^2/|B|, so it's a density weighted
  /// value and also diluted by the proportion of B on the los
  /// direction
  ///
  double get_pt_BYabs(struct point_4cellavg *, ///< pt
		      const int, ///< ifrac
		      const int, ///< bx index (image coords)
		      const int, ///< by index (image coords)
		      const int, ///< bz index (image coords)
		      const int, ///< sign(xx)
		      const int, ///< sign(yy)
		      const int, ///< sign(zz)
		      const double, ///< sin(theta)
		      const double  ///< cos(theta)
		      );
protected:
  bool cell_positions_set; ///< set to true if cell positions have been set.
  void initialise_pixels();        ///< allocate memory data in each pixel.
  void delete_pixel_data(pixel *); ///< Delete allocated memory in pixel.
  void find_surrounding_cells(const double *, ///< position of point, in simI coordinates.
			      cell *,         ///< cell in plane, move from here to point.
			      cell **,        ///< OUTPUT: list of 4 cells surrounding point
			      double *        ///< OUTPUT: list of weights for each cell.
			      );

};


// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------



#endif // SIM_PROJECTION_H
