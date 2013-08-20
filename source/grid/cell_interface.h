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
/// - 2011.02.17 JM: Enhanced optical depth storage for multiple
///    sources.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid
///    now.
/// - 2011.03.21 JM: Updated optical-depth info, multiple variables
///    per source.
/// - 2011.04.18 JM: Added storage for dS, path length through a cell
///    for raytracing.
/// - 2013.01.14 JM: Updated comments and re-arranged things.
///    Added cell *npt_all for GRIDV2 (linked list including B.C.s).
/// - 2013.08.20 JM: Moved raytracing set/get functions to header and
///    made them inline.  Added option for radiation sources to be
///    associated with NTau>=1 optical depths, for flexibility.

#ifndef CELL_INTERFACE_H
#define CELL_INTERFACE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
// NEED AN INCLUDE TO GET struct rad_sources DECLARATION

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
#ifdef GRIDV2
  cell *npt_all;  ///< Pointer to next point in grid, including
                  ///< boundary data.
#endif
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


#ifndef GRIDV2
#include "global.h"
#endif
///
/// Cell Interface class, so each cell doesn't need to have a list of
/// pointers to all these functions.  This should contain
/// non-geometry-specific and non-grid-specific functions.
///
class cell_interface {
 public:
  cell_interface();
  ~cell_interface();
  ///
  /// Set flag so Ph,dU are not initialised, to save memory for
  /// analysis code.
  ///
  void set_minimal_cell_data();

  void unset_minimal_cell_data(); ///< Back to normal cells.

  bool query_minimal_cells(); ///< true if using minimal cells.

  class cell * new_cell();  ///< Create a new cell.

  /// Delete all the dynamically allocated memory in a cell.
  void delete_cell(cell *);

  ///
  /// Copy one cell's contents to another.
  ///
  void copy_cell(
    const cell *, ///< input : copy from this cell.
    cell *        ///< output: copy to this cell.
    );

  ///
  /// Print most of the information held by a cell.
  ///
  void print_cell(
    const cell * ///< cell whose info we want to display.
    );

  ///
  /// Returns the size of a cell in the internal integer units (2).
  ///
  inline int get_integer_cell_size() {return cell_size_int_units;}

  ///
  /// Return physical size of one internal unit
  ///
  inline double phys_per_int() {return dxo2;}

  ///
  /// Set cell position by passing in a double precision vector.
  ///
  void set_pos(
    cell *,        ///< pointer to cell
    const double * ///< double array of size ndim, cell position.
    );

  ///
  /// Set cell position using an integer vector, which is offset from
  /// SimPM.Xmin[] and with a cell width equal to 2 units.
  ///
  void set_pos(
    cell *,     ///< pointer to cell
    const int * ///< int array of size ndim, cell integer position.
    );

  ///
  /// Returns double precision position of cell centre (code units).
  ///
  void get_dpos(
    const cell *, ///< pointer to cell
    double *      ///< array to write position into.
    );

  ///
  /// Returns double precision position of cell centre in a specified
  /// coordinate direction (code units).
  ///
  double get_dpos(
    const cell *, ///< pointer to cell
    const int     ///< element of position we want.
    );

  ///
  /// Converts a double precision physical position into an integer
  /// position.
  ///
  void get_ipos_vec(
    const double *, ///< physical position (input)
    int *           ///< integer position (output)
    );

  ///
  /// Converts from an integer position to a double precision pos.
  ///
  void get_dpos_vec(
    const int *, ///< integer position  (input)
    double *     ///< physical position (output)
    );

  ///
  /// Returns integer position of cell centre.
  ///
  void get_ipos(
    const cell *, ///< pointer to cell
    int * ///< array to write integer position into.
    );

  ///
  /// Returns integer position of cell centre along a specified
  /// coordinate direction.
  ///
  int get_ipos(
    const cell *, ///< pointer to cell
    const int ///< element of position we want.
    );

  ///
  /// Converts a double precision physical position into its
  /// equivalent value in internal (integer) units.
  ///
  void get_ipos_as_double(
    const double *, ///< physical position (input)
    double *        ///< floating pt. position in integer position units (output)
    );



// ##################################################################
// ##################################################################


  // --------- EXTRA DATA FOR RAYTRACING AND H-CORRECTION -----------
  ///
  /// Set variables for extra_data based on what we need per cell.
  /// The H-correction needs Ndim doubles, and radiation sources need
  /// either 1, 3, or 5 doubles each.
  /// THIS DOES NOT ALLOW THE NUMBER OF SOURCES TO CHANGE OVER TIME.
  /// THIS *MUST* BE SET BEFORE CREATING ANY CELLS ON THE GRID.
  ///
  void setup_extra_data(
    const struct rad_sources &, ///< Pointer to (SimPM.RS)
                                ///< ray-tracing struct.
    const int,  ///< Flag for H-correction requirements
    const int   ///< Flag for Div(V) variable required.
    );


// ##################################################################
// ##################################################################

  ///
  /// Return the optical depth along a ray between the source and the
  /// cell centre, from the source to the point where the ray exits
  /// the cell.  Radiation source v should already be setup at the
  /// start of the simulation.
  ///
  inline void get_col(
    const cell *c,  ///< current cell.
    const int s,    ///< index of source.
    double *tau     ///< pointer to array for optical depths.
    )
  {
    for (short unsigned int v=0; v<NTau[s]; v++) {
      tau[v] = c->extra_data[iTau[s]+v];
    }
    return;
  }


// ##################################################################
// ##################################################################

  ///
  /// Set the optical depth along the ray (see get_col for details)
  ///
  inline void set_col(
    cell *c,
    const int s,       ///< index of source.
    const double tau[] ///< value to set
    )
  {
    for (short unsigned int v=0; v<NTau[s]; v++) {
      c->extra_data[iTau[s]+v] = tau[v];
    }
    return;
  }


// ##################################################################
// ##################################################################

  ///
  /// Get cell optical depth
  ///
  inline void get_cell_col(
    const cell *c,
    const int s,    ///< index of source.
    double *dtau    ///< pointer to array for optical depths.
    )
  {
    for (short unsigned int v=0; v<NTau[s]; v++) {
      dtau[v] = c->extra_data[iDTau[s]+v];
    }
    return;
  }


// ##################################################################
// ##################################################################

  ///
  /// Set cell optical depth
  ///
  inline void set_cell_col(
    cell *c,
    const int s,        ///< index of source.
    const double dtau[] ///< value to set
    )
  {
    for (short unsigned int v=0; v<NTau[s]; v++) {
      c->extra_data[iDTau[s]+v] = dtau[v];
    }
    return;
  }
 

// ##################################################################
// ##################################################################

  ///
  /// Set cell Vshell value (for raytracing).
  ///
  inline void set_cell_Vshell(
    cell *c,
    const int s,    ///< index of source.
    const double Vs ///< value to set.
    )
  {
#ifdef RT_TESTING
    if (iVsh[s] <0) { cout <<"source "<<s<<": ";
      rep.error("Source has no Vhsell variable", iVsh[s]); }
#endif // RT_TESTING
    c->extra_data[iVsh[s]] = Vs;
    return;
  }
  

// ##################################################################
// ##################################################################

  ///
  /// Get cell Vshell value (for raytracing).
  ///
  inline double get_cell_Vshell(
    const cell *c, ///< current cell.
    const int s    ///< index of source.
    )
  {
#ifdef RT_TESTING
    if (iVsh[s] <0) { cout <<"source "<<s<<": ";
      rep.error("Source has no Vhsell variable", iVsh[s]); }
#endif // RT_TESTING
    return c->extra_data[iVsh[s]];
  }


// ##################################################################
// ##################################################################

  ///
  /// Set raytracing path length through cell for source i.
  ///
  inline void set_cell_deltaS(
        cell *c,
        const int s,  ///< index of source.
        const double deltaS
        )
  {
#ifdef RT_TESTING
    if (idS[s] <0) { cout <<"source "<<s<<": ";
      rep.error("Source has no deltaS variable", idS[s]); }
#endif // RT_TESTING
    c->extra_data[idS[s]] = deltaS;
    return;
  }
  

// ##################################################################
// ##################################################################

  ///
  /// Get raytracing path length through cell for source i.
  ///
  inline double get_cell_deltaS(
    const cell *c,  ///< current cell.
    const int s     ///< index of source.
    )
  {
#ifdef RT_TESTING
    if (idS[s] <0) { cout <<"source "<<s<<": ";
      rep.error("Source has no Vhsell variable", idS[s]); }
#endif // RT_TESTING
    return c->extra_data[idS[s]];
  }
    

// ##################################################################
// ##################################################################

  ///
  /// Get the H-correction coefficient eta for the requested cell in
  /// the requested direction (at the positive direction cell
  /// interface).
  ///
  inline double get_Hcorr(
    const cell *c,
    const axes a
    )
  {return c->extra_data[iHcorr[a]];}


// ##################################################################
// ##################################################################

  ///
  /// Set the H-correction coefficient eta for cell in direction a
  ///
  inline void   set_Hcorr(
    cell *c,
    const axes a,
    double eta
    )
  {c->extra_data[iHcorr[a]] = eta; return;}


// ##################################################################
// ##################################################################

  ///
  /// Get div(v) for a cell.
  ///
  inline double get_DivV(
    const cell *c
    )
  {return c->extra_data[iDivV];}


// ##################################################################
// ##################################################################

  ///
  /// Set div(v) for a cell.
  ///
  inline void   set_DivV(
    cell *c,
    double divv
    )
  {c->extra_data[iDivV] = divv; return;}

// ##################################################################
// ##################################################################

  // --------- EXTRA DATA FOR RAYTRACING AND H-CORRECTION -----------

 private:
  bool minimal_cell; ///< default is false, set to true if analysing sims.
  double dxo2; ///< Half a cell width.
  double *xmin; ///< The global Xmin of the domain, for counting integer positions from.
  double int_converter; ///< Number of integers per cell width.
  int cell_size_int_units; ///< size of a cell in integer units (==2)

  bool have_setup_extra_data; ///< Flag checked when creating a cell!
  short unsigned int using_RT;     ///< Flag: 0=not doing RT, 1=doing RT.
  short unsigned int using_Hcorr;  ///< Flag: 0=no Hcorr, N=need N variables (Hcorr=Ndim).
  short unsigned int using_DivV;   ///< Flag: 0=don't need div(v), 1=do need div(v).
  unsigned int N_extra_data; ///< Size of extra_data array (can be zero).

  ///
  /// array where the value of the 'i'th element is the number of
  /// optical depths associated with radiation source i.
  ///
  short unsigned int *NTau;

  ///
  /// First optical depth for source v is contained in the iTau[v]'th
  /// element of extra_data, i.e. extra_data[iTau[v]].
  /// Any further optical depths are in the following elements e.g.
  /// extra_data[iTau[v]+1], extra_data[iTau[v]+2]
  ///
  short unsigned int *iTau;

  ///
  /// First cell optical depth for source v is contained in the
  /// iDTau[v]'th element of extra_data, i.e. extra_data[iDTau[v]].
  /// Any further optical depths are in the following elements e.g.
  /// extra_data[iDTau[v]+1], extra_data[iDTau[v]+2]
  ///
  short unsigned int *iDTau;

  ///
  /// Vshell parameter for photon-conserving-RT is given by value in
  /// extra_data[iVsh[v]] for radiation source v.
  ///
  short unsigned int *iVsh;

  ///
  /// Path length through cell for radiation is given by the value in
  /// extra_data[idS[v]] for radiation source v.
  ///
  short unsigned int *idS;

  ///
  /// indices of Hcorrection values in extra_data [XX,YY,ZZ] are
  /// given by extra_data[iHcorr[XX/YY/ZZ]].
  ///
  short unsigned int iHcorr[MAX_DIM];

  ///
  /// If an algorithm needs the divergence of the velocity field in
  /// the cell data, it is in extra_data[iDivV]].
  ///
  short unsigned int iDivV;

};

#endif // CELL_INTERFACE_H
