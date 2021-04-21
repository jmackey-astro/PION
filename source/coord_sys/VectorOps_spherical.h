///
/// \file VectorOps_spherical.h
///
/// \author Jonathan Mackey
///
/// Class for spherical coordinates, equal spacing in r.
/// This is only 1D for now.
///
/// Created 2010.10.01
///
/// Modifications:
/// - 2010.11.03 JM: Updated R_com() and introduced R3() so that the
///  Finite Volume formulation is self-consistent between source terms
///  and divergence components.
///
/// - 2010.12.04 JM: Added constructor with only one argument.  Also
///   a set_dx() function.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.

#ifndef VECTOROPS_SPHERICAL_H
#define VECTOROPS_SPHERICAL_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "VectorOps.h"

class VectorOps_Sph : virtual public VectorOps_Cyl {
public:
  ///
  /// Constructor, sets number of spatial dimensions.
  ///
  VectorOps_Sph(const int  ///< number of spatial dimensions of grid
  );

  ///
  /// Destructor
  ///
  ~VectorOps_Sph();

  /// Returns Volume of cell.
  virtual double CellVolume(
      const cell *,  ///< cell pointer
      const double   ///< cell diameter
  );

  /// Returns Surface area of interface.
  virtual double CellInterface(
      const cell *,     ///< Cell
      const direction,  ///< outward normal to interface.
      const double      ///< cell diameter
  );

  ///
  /// Returns maximum of all gradients with neighbouring cells for Spherical
  /// Coordinates.
  ///
  virtual double maxGradAbs(
      const cell *,  ///< current point.
      const int,     ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,     ///< Index of Primitive Variable to get gradient of.
      class GridBaseClass *  ///< pointer to computational grid.
  );

  ///
  /// Vector Gradient Operator for Spherical Coordinates.
  /// Uses second order central difference formula (modified for
  /// spherical geometry)
  ///
  /// \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} =
  ///     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
  ///
  virtual void Gradient(
      const cell *,  ///< Cell to calculate gradient of
      const int,     ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,     ///< Index of Scalar (in state vector) to calculate
                     ///< gradient of.
      class GridBaseClass *,  ///< pointer to computational grid.
      pion_flt *              ///< Pointer to array to put gradient vector.
  );

  ///
  /// Calculate Vector Diff. Op. Divergence of a vector at a point for
  /// Spherical Coordinates.
  /// This is a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] for the x-component of the divergence.
  ///
  /// It assumes spherical coordinates in the order (r,theta,phi), and
  /// includes the relevant scale factors for each derivative.
  ///
  virtual double Divergence(
      cell *,       ///< point for which to calculate div(B)
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int *,  ///< Indices of vector quantity (in state vector) to
                    ///< calculate divergence of. Should contain 3
                    ///< elements, ordered as x,y,z components.
      class GridBaseClass *  ///< pointer to computational grid.
  );

  ///
  /// Calculate Curl of Vector Quantity at a point for Spherical Coordinates.
  /// Note this is identically zero in 1D.
  virtual void Curl(
      const cell *,  ///< point for which to calculate curl
      const int,     ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int *,   ///< Indices of vector quantity (in state vector) to
                     ///< calculate divergence of. Should contain 3 elements,
                     ///< ordered as x,y,z components.
      class GridBaseClass *,  ///< pointer to computational grid.
      pion_flt *              ///< Pointer to array to put curl vector.
  );

  ///
  /// Given a state, and a slope (dP/dx), construct an edge state.
  ///
  /// For first order accuracy, left and right edge states are just the
  /// cell centred value, P0.\n
  /// For second order accuracy, The vector of slopes is used to construct
  /// edge states \f$ P \pm (dP/dx)dx/2 \f$.
  ///
  /// This function assumes that the slope array contains an already
  /// calculated gradient for the grid point in question.  The pivot point is
  /// the centre of gravity of the cell.
  ///
  int SetEdgeState(
      const cell *,      ///< Current Cell.
      const direction,   ///< Add or subtract the slope depending on direction.
      const int,         ///< length of state vectors.
      const pion_flt *,  ///< Slope vector.
      pion_flt *,        ///< vector for edge state.
      const int,         ///< Order of spatial Accuracy.
      class GridBaseClass *  ///< pointer to computational grid.
  );

  ///
  /// Calculate the slope (dP/dx) in cell i.
  ///
  int SetSlope(
      const cell *,          ///< Current Cell.
      const axes,            ///< Which direction to calculate slope in.
      const int,             ///< length of state vectors.
      pion_flt *,            ///< Slope vector to be written to.
      const int,             ///< Order of spatial Accuracy.
      class GridBaseClass *  ///< pointer to computational grid.
  );

  ///
  /// Calculates the i-th term of the divergence for a vector of
  /// variables, given a cell and positive and negative flux vectors.
  ///
  int DivStateVectorComponent(
      const cell *,           ///< current cell.
      class GridBaseClass *,  ///< pointer to computational grid.
      const axes,             ///< current coordinate axis we are looking along.
      const int,              ///< length of state vectors.
      const pion_flt *,       ///< Negative direction flux.
      const pion_flt *,       ///< Positive direction flux.
      pion_flt *              ///< Vector to assign divergence component to.
  );

protected:
  ///
  /// Returns \f[ \int r^2 dr \f] from \f[ r_i-\delta r/2\f] to
  /// \f[r_i+\delta r/2\f], normalised by \f[r_i\delta r \f].
  /// This is used for calculating the centre-of-volume of a cell, and
  /// for evaluating volume averages of source terms.
  ///
  virtual inline double R3(
      const cell *c,   ///< Cell to operate on.
      const double dR  ///< cell diameter
  )
  {
    return (CI.get_dpos(c, Rsph) + dR * dR / 12.0 / CI.get_dpos(c, Rsph));
  }

  ///
  /// Returns the centre-of-volume of a cell in spherically symmetry:
  /// This is \f[ (\int r dV) / (\int dV ) = \int r^3dr/\int r^2dr \f]
  /// This point is not the midpoint of a cell for non-cartesian
  /// geometries, and is used as the pivot point for the slope
  /// calculation in the 2nd order algorithm.  Obviously for
  /// \f[ r>>\delta r \f] it approaches the midpoint.
  ///
  virtual inline double R_com(
      const cell *c,   ///< cell to operate on
      const double dR  ///< cell diameter
  )
  {
    double delta2 = dR / CI.get_dpos(c, Rsph);
    delta2 *= delta2;
    return CI.get_dpos(c, Rsph) * (1.0 + 0.25 * delta2) / (1.0 + delta2 / 12.0);
  }
};

#endif  // VECTOROPS_SPHERICAL_H
