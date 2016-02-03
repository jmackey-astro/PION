///
/// \file raytracer_shielding.h
/// \author Jonathan Mackey
///
/// \brief Declaration of ray-tracing class which only calculates optical depth
///        along grid axes to every cell.
///
/// Modifications:\n
/// - 2011.02.17 JM: File created.
/// - 2011.02.24 JM: removed trailing backslash on a comment line.
/// - 2011.03.02 JM: Added parallelised class.
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function. 
///

#ifndef RAYTRACING_SHIELDING_H
#define RAYTRACING_SHIELDING_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "raytracer_SC.h" // short-characteristics ray-tracer.

///
/// This class is designed to just calculate optical depths to a source based
/// on a diffuse interstellar radiation field (ISRF).  A source is placed at 
/// infinity in every coordinate direction (as long as the boundary is not 
/// reflecting), and this function will add a column density along each ray to
/// the current cell, or more precisely an A_v extinction value.
///
/// The algorithm is as described by:
/// - Nelson \& Langer, 1997, ApJ, 524, 923.
/// - Glover, Federrath, Mac-Low, \& Klessen, 2010, 404, 2.
///
class raytracer_shielding
	: virtual public raytracer_USC_infinity // uniform grid, source at infinity.
{
  public:
  ///
  /// Constructor (just sets pointers to local grid and to Microphysics class).
  ///
  raytracer_shielding(
                  class GridBaseClass *,   ///< Pointer to grid
                  class MicroPhysicsBase * ///< Pointer to MicroPhysics Class.
                  );
  ///
  /// Destructor really doesn't do anything.
  ///
  ~raytracer_shielding();

  ///
  /// Just calculate the column densities required for RT.
  ///
  virtual int RayTrace_Column_Density(
                const int,    ///< Source id
                const double, ///< Timestep
                const double  ///< EOS Gamma.
                );

  protected:
   ///
   /// Here we just calculate the optical depth to the cell in the current
   /// direction, and assign this using the cell_interface to the appropriate
   /// Av variable for the face the ray is entering the cell through.  The
   /// cell column density variable is used for the total column along the ray
   /// to the point where it exits the current cell (i.e. this is the column
   /// density to the point entering the next cell).
   ///
   virtual int ProcessCell(
                  class cell *,       ///< Current cell.
                  double,             ///< Column to cell.
                  double,             ///< Path Length through cell.
                  const rad_source *, ///< pointer to source struct.
                  const double        ///< Timestep
                  );

   private:
   std::vector<bool> have_set_Vshell; ///< Set to false initially, and then to true once it is set.
};

#ifdef PARALLEL
///
/// Parallel self-shielding calculator differs from the serial version only
/// in that it sets up send and receive boundaries for RT in uniform_gridMPI
/// when adding a new source.  Receive and send functions are then called 
/// before and after the serial ray-tracing algorithm.
///
class raytracer_shielding_pllel :
  virtual public raytracer_shielding
{
  public:
  ///
  /// constructor.  Doesn't do all that much.
  ///
  raytracer_shielding_pllel(class GridBaseClass *, ///< Pointer to grid
                            class MicroPhysicsBase * ///< Pointer to MicroPhysics Class.
                            );

  ///
  /// Destructor. Doesn't do all that much.
  ///
  ~raytracer_shielding_pllel();

  /// \brief Adds a source to the list of sources to trace.
  ///
  /// Differs from serial version in that it also calls 
  /// gridptr->setup_RT_boundaries() so that the grid can decide which (if any)
  /// extra corner boundaries it needs to set up for sending and receiving data.
  ///
  virtual int Add_Source(const struct rad_src_info * ///< ptr to source info.
                        );

  ///
  /// Parallel version of the ray-tracing algorithm.  This first receives any
  /// required boundary data from other processors, then performs the serial
  /// raytracing algorithm, and finally sends boundary data to any processes 
  /// which need it.
  ///
  virtual int RayTrace_SingleSource(const int,    ///< Source id
                                    const double, ///< Timestep
                                    const double  ///< eos gamma.
                                    );

};

#endif // PARALLEL

#endif // define RAYTRACING_SHIELDING_H
