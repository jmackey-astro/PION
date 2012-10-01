///
/// \file   cell_interface.h
/// \author Jonathan Mackey
/// \date   12.11.2010
///
/// Purpose: Declare a set of routines for accessing cell data; these
/// would be too memory intensive if each cell had to have a copy of
/// the functions.
/// 
/// History: Used to be in global.h (up to SVN rev. 236).
///
/// Modifications:\n
/// - 2010.11.12 JM: Added support for H-correction speeds.  Changed
///   'col' to be an access function 'monochromatic_tau(cell *c)'.
///
/// - 2010.11.15 JM: Added include global.h!
/// - 2010.12.30 JM: Added DivV() set/get functions.
/// - 2011.02.17 JM: Enhanced optical depth storage for multiple sources.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2011.03.21 JM: Updated optical-depth info, multiple variables per source.
/// - 2011.04.18 JM: Added storage for dS, path length through a cell for raytracing.

#ifndef CELL_INTERFACE_H
#define CELL_INTERFACE_H

#include "../defines/functionality_flags.h"
#include "../defines/testing_flags.h"

class cell_interface;

#ifdef COUNT_ENERGETICS
struct energetics {
  double 
  ci_cooling,   // in erg/cm3/s
    rr_cooling, // in erg/cm3/s
    fn_cooling, // in erg/cm3/s
    pi_heating, // in erg/cm3/s
    ci_rate,    // in 1/cm3/s
    rr_rate,    // in 1/cm3/s
    pi_rate,    // in 1/cm3/s
    tot_heating, // in erg/cm3/s
    tot_cooling, // in erg/cm3/s
    net_heating, // in erg/cm3/s
    cooling_time, // in years
    recomb_time;  // in years
};

#endif

///
/// Data Point for finite-volume grid. 
///
/// This gives the datapoint an id, a position (which is zone centred
/// for finite volume grid), a Primitive variable vector, P, a conserved update vector, dU,
/// and a half-timestep vector, Ph.
///
class cell {
 public:
  friend class cell_interface;
  double *P;   ///< Primitive Variables.
  double *Ph;  ///< Primitive State vector at half timestep.
  double *dU;  ///< Update vector to be added to U when we do the time update.
  cell **ngb;  ///< Pointers to cell's Neigbours.
  cell *npt;   ///< Pointer to next point in list, or some other cell.
  int id;      ///< Grid point's id.
  int isedge;  ///< Integer specifying if it is an edge cell, and if so, which edge.
  bool isbd;   ///< True if cell is boundary data, false if not.
  bool isgd;   ///< True if cell is grid data, false if not.
#ifdef COUNT_ENERGETICS
  struct energetics e; ///< to count up the energetics properties of the current cell.
#endif
  int *pos;
  private:
  double *extra_data; ///< General purpose data (Tau in ray-tracing, eta for H-correction)
};

  
#include "../global.h"
///
/// Cell Interface class, so each cell doesn't need to have a list of pointers
/// to all these functions.  This should contain non-geometry-specific and 
/// non-grid-specific functions.
///
class cell_interface {
 public:
  cell_interface();
  ~cell_interface();
  /** \brief Set code so Ph,dU are not initialised, to save memory for
   *  analysis code. */
  void set_minimal_cell_data();
  void unset_minimal_cell_data(); ///< Back to normal cells.
  bool query_minimal_cells(); ///< returns true if using minimal cells.
  class cell * new_cell();  ///< Create a new cell.
  void delete_cell(cell *); ///< Delete all the dynamically allocated memory in a cell.
  /** \brief Set cell position by passing in a double precision vector.  */
  void set_pos(cell *,        ///< pointer to cell
	       const double * ///< double array of size ndim, containing cell position.
	       );
  ///
  /// Set cell position using an integer vector, which is offset from
  /// SimPM.Xmin[] and with a cell width equal to 2 units.
  ///
  void set_pos(cell *,     ///< pointer to cell
	       const int * ///< int array of size ndim, containing cell integer position.
	       );
  /** \brief Returns double precision position of cell centre. */
  void get_dpos(const cell *, ///< pointer to cell
		double * ///< array to write position into.
		);
  /** \brief Returns double precision position of cell centre in a specified coordinate direction. */
  double get_dpos(const cell *, ///< pointer to cell
		  const int ///< element of position we want.
		  );
  /** \brief Returns integer position of cell centre. */
  void get_ipos(const cell *, ///< pointer to cell
		int * ///< array to write integer position into.
		);
  /** \brief Returns integer position of cell centre along a specified coordinate direction. */
  int get_ipos(const cell *, ///< pointer to cell
	       const int ///< element of position we want.
	       );
  /** \brief Converts a double precision physical position into an integer position. */
  void get_ipos_vec(const double *, ///< physical position (input)
		    int *           ///< integer position (output)
		    );
  ///
  /// Converts a double precision physical position into its
  /// equivalent value in internal (integer) units.
  ///
  void get_ipos_as_double(const double *, ///< physical position (input)
			  double *        ///< floating pt. position in integer position units (output)
			  );
  /** \brief Converts from an integer position to a double precision position. */
  void get_dpos_vec(const int *, ///< integer position  (input)
		    double *     ///< physical position (output)
		    );
  void copy_cell(const cell *, ///< input : copy from this cell.
		 cell *        ///< output: copy to this cell.
		 );
  void print_cell(const cell * ///< cell whose info we want to display.
		 );
  /** \brief returns the size of a cell in the internal integer units. */
  inline int get_integer_cell_size() {return cell_size_int_units;}
  ///
  /// Return physical size of one internal unit
  ///
  inline double phys_per_int() {return dxo2;}

  ///
  /// Set variables for extra_data based on what we need per cell.
  /// The H-correction needs Ndim doubles, and radiation sources need either 1, 3,
  /// or 5 doubles each.
  /// THIS DOES NOT ALLOW THE NUMBER OF SOURCES TO CHANGE OVER TIME.
  /// THIS *MUST* BE SET BEFORE CREATING ANY CELLS ON THE GRID.
  ///
  void setup_extra_data(const struct rad_sources &,  ///< Pointer to (SimPM.RS) ray-tracing struct.
			const int,  ///< Flag for H-correction requirements
			const int   ///< Flag for Div(V) variable required.
			);

  //
  // Return the optical depth along a ray between the source and the
  // cell centre, from the source to the point where the ray exits the
  // cell.  The indexing can be whatever the calling class wants it to
  // be, running from zero to Ntau-1.
  //
  inline double get_col(const cell *c,  ///< current cell.
                        const int v     ///< index of source.
                        )
  {return c->extra_data[iTau0[v]];}

  ///
  /// Set the optical depth along the ray (see get_col for details)
  ///
  inline void   set_col(cell *c,
                        const int v,  ///< index of source.
                        const double tau
                        )
  {c->extra_data[iTau0[v]] = tau; return;}


  ///
  /// Get cell optical depth
  ///
  double get_cell_col(const cell *c,  ///< current cell.
                      const int v     ///< index of source.
                      );

  ///
  /// Set cell optical depth
  ///
  void   set_cell_col(cell *c,
                      const int v,  ///< index of source.
                      const double tau
                      );
 
  ///
  /// Set cell Vshell value (for raytracing).
  ///
  void set_cell_Vshell(cell *,
                       const int ,  ///< index of source.
                       const double
                       );
  
  ///
  /// Get cell Vshell value (for raytracing).
  ///
  double get_cell_Vshell(const cell *,  ///< current cell.
                         const int      ///< index of source.
                         );
  ///
  /// Set raytracing path length through cell for source i.
  ///
  void set_cell_deltaS(cell *,
                       const int,   ///< index of source.
                       const double ///< path length dS
                       );
  
  ///
  /// Get raytracing path length through cell for source i.
  ///
  double get_cell_deltaS(const cell *,  ///< current cell.
                         const int      ///< index of source.
                         );



  
  //
  // Get the H-correction coefficient eta for the requested cell in
  // the requested direction (at the positive direction cell
  // interface).
  //
  inline double get_Hcorr(const cell *c, const axes a)
  {return c->extra_data[iHcorr[a]];}

  //
  // Set the H-correction coefficient eta for cell in direction a
  //
  inline void   set_Hcorr(cell *c, const axes a, double eta)
  {c->extra_data[iHcorr[a]] = eta; return;}

  //
  // Get div(v) for a cell.
  //
  inline double get_DivV(const cell *c)
  {return c->extra_data[iDivV];}

  //
  // Set div(v) for a cell.
  //
  inline void   set_DivV(cell *c, double divv)
  {c->extra_data[iDivV] = divv; return;}

 private:
  bool minimal_cell; ///< default is false, set to true if analysing sims.
  double dxo2; ///< Half a cell width.
  double *xmin; ///< The global Xmin of the domain, for counting integer positions from.
  double int_converter; ///< Number of integers per cell width.
  int cell_size_int_units; ///< size of a cell in integer units (==2)

  bool have_setup_extra_data; ///< Flag checked when creating a cell!
  int using_RT;     ///< Flag: 0=not doing RT, 1=doing RT.
  int using_Hcorr;  ///< Flag: 0=no Hcorr, N=need N variables (Hcorr=Ndim).
  int using_DivV;   ///< Flag: 0=don't need div(v), 1=do need div(v).
  int N_extra_data; ///< Size of extra_data array (can be zero).
  int *iTau0;  ///< index of First  Optical depth for source(s) in extra_data;
  //int *iTau1;  ///< index of Second Optical depth for source(s) in extra_data;
  int *iDTau0; ///< index of 1st cell Optical depth for source(s) in extra_data;
  int *iVsh;  ///< index in extra_data of Vshell parameter for photon-conserving-RT.
  int *idS;   ///< index in extra_data of path length through cell for source(s). 
  int iHcorr[MAX_DIM]; ///< indices of Hcorrection values in extra_data [XX,YY,ZZ].
  int iDivV;           ///< Index of Div(V) variable in extra data.

};

#endif // CELL_INTERFACE_H
