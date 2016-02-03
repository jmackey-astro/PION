/// \file raytracer_SC.h
/// \author: Jonathan Mackey
/// 
/// Base raytracing class, to define the interface functions.
///
/// Modifications:
/// - 2015.03.13 JM: Moved from global.h
///

#ifndef RAYTRACER_BASE_H
#define RAYTRACER_BASE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "microphysics/microphysics_base.h"

// ##################################################################
// ##################################################################



#ifndef CELL_CENTRED_SRC
#ifndef NON_CELL_CENTRED_SRC
#error "MUST define one of CELL_CENTRED_SRC and NON_CELL_CENTRED_SRC in global.h"
#endif
#endif

#ifdef CELL_CENTRED_SRC
#error "CELL_CENTRED_SRC has not been used for a LONG time and may be buggy."
#error "USE AT YOUR OWN RISK.  DELETE THIS ERROR IN source/raytracing/raytracer_SC.h"
#endif

#define NO_SOURCE_CELL_GEOMETRY ///< this stops the code changing tau in the
                      ///< source cell for different angles (makes no difference)

#define INTERPOLATE_METHOD 0 /// 0=C2Ray mintau=0.7, 1=Geometric,
        /// 2=Geometric in exp(-tau), 3=quadratic/exp(-tau) !!!ONLY C2-RAY IS ANY GOOD!!!



// ##################################################################
// ##################################################################


///
/// Related quadrants in x-y plane to indices for the order
/// they are to be traversed.  The source is in Q1, so it comes
/// first, then Q2, then Q4, and finally Q3.  Q2 and Q4 could be
/// interchanged, but this is the way I'm doing it.
///
enum Quadrants {Q1=0, Q2=1, Q4=2, Q3=3};
#define NQUADS 4


// ##################################################################
// ##################################################################



///
/// This arranges the octants in integers according to how I
/// am tracing them. The first four octants are in ZP dir, last
/// four in ZN, and each group of four traces an anticlockwise
/// circle in the x-y plane, starting with (XP,YP).
///
enum Octants {OCT1=0, OCT2=1, OCT4=2, OCT3=3, OCT5=4, OCT6=5, OCT8=6, OCT7=7};


// ##################################################################
// ##################################################################




///
/// Struct to hold info on radiation sources.  An extension of rad_src_info in global.h
///
struct rad_source {
  ///
  /// pointer to source (set to SimPM.RS.source[id]) with the basic
  /// info about the source.
  ///
  struct rad_src_info *s;

  class cell *sc; ///< nearest cell to source.
  bool src_on_grid; ///< true if source is at a grid cell.
  int ipos[MAX_DIM];    ///< source position in integer form (grid units, dx=2).

  ///
  /// This struct is used by the code to pass cell and source data
  /// to the microphysics integrator.  It contains the relevant source
  /// information, and also some cell-source geometry information which 
  /// must be set on a cell-by-cell basis as rays are traced.
  /// The struct is declared in source/microphysics/microphysics_base.h
  ///
  struct rt_source_data data;
};



// ##################################################################
// ##################################################################





///
/// Pure virtual ray-tracer base class. 
///
/// This provides the interface to the raytracer, regardless of the 
/// implementation.  So far these are the functions I'll need for 
/// a uniform grid with a single source, so I may need to add more
/// later.
///
class RayTracingBase {
 public:
  virtual ~RayTracingBase() {} ///< non-virtual destructor.

  ///
  /// \brief Adds a source to the list of sources to trace.
  /// 
  /// Returns the sources id, which starts at zero and increases 
  /// linearly.  So if we have 10 sources, and add another source,
  /// its id will be 10.  Note the ID is contained in the rad_src_info
  /// struct though, so each particular class may not have consecutively 
  /// numbered sources.
  ///
  virtual int Add_Source(struct rad_src_info * ///< ptr to source info.
       )=0;

  ///
  /// Processes a single source's effect on the grid over a timestep.
  ///
  virtual int RayTrace_SingleSource(const int, ///< Source id
                                    const double, ///< Timestep
                                    const double  ///< EOS Gamma.
                                    )=0;
  ///
  /// Just calculate the column densities required for RT.
  ///
  virtual int RayTrace_Column_Density(
                const int,    ///< Source id
                const double, ///< Timestep
                const double  ///< EOS Gamma.
                )=0;

  /// \brief Prints list of sources with id, location, strength.
  virtual void Print_SourceList()=0;
  
  /// \brief Returns the number of sources to track.
  virtual int NSources()=0;

  ///
  /// Returns whether we are doing an implicit (==0) or an explicit (==1)
  /// integration of the raytracing/microphysics.
  ///
  virtual int type_of_RT_integration()=0;

  ///
  /// This sets the number of ionising and UV heating sources of radiation,
  /// and makes sure the rt_source_data structs are populated correctly.
  /// It can be used to change the radiation sources if e.g. the luminosity
  /// changes over time.
  ///
  virtual void update_RT_source_properties(
                  const struct rad_src_info * ///< ptr to source info.
                  )=0;

  /// Returns the number of ionising sources
  virtual int N_ionising_sources()=0;

  /// Returns the number of UV heating sources
  virtual int N_heating_sources()=0;

  ///
  /// This function copies the ionising source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_ionising_src_list(
                std::vector<struct rt_source_data> & ///< list of data for ionising sources
                )=0;

  ///
  /// This function copies the UV-heating source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_UVheating_src_list(
                std::vector<struct rt_source_data> & ///< list of data for UV-heating sources
                )=0;

};


// ##################################################################
// ##################################################################



extern class RayTracingBase *RT; ///< Raytracer for all radiation sources
// ************************** RAY TRACER ***************************
// *****************************************************************



#endif // RAYTRACER_BASE_H
