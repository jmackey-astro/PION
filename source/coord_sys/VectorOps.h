/// \file VectorOps.h
///
/// \brief Declares VectorOps utility classes
///
/// \author Jonathan Mackey
///
///
/// This declares classes of functions which act on a grid of cells, but which
/// depend on the coordinate system being used (for example the non-local vector
/// differential operators grad, div, and curl).  The base class is
/// incomplete, and the derived classes have names VectorOps_Cart,
/// VectorOps_Cyl, etc. for the different coordinate systems.
///
/// Modified:\n
/// - 2007-08-01 File Created
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.12.04 JM: Added constructor with only one argument.  Also
///   a set_dx() function.
/// - 2013.07.19 JM: Added TRACER_SLOPES_CONSERVED_VARS option, but
///    it is more diffusive than primitive variables, so it is not
///    used.
/// - 2015.01.13 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#ifndef VECTOROPS_H
#define VECTOROPS_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "grid/cell_interface.h"
#include "grid/grid_base_class.h"

///
/// Class to hold vector differential operations which depend on
/// the coordinate system.
///
class BaseVectorOps {
public:
  virtual ~BaseVectorOps();

  virtual double CellVolume(
      const cell*,  ///< cell pointer
      const double  ///< cell diameter
      ) = 0;        ///< Returns Volume of cell.

  virtual double CellInterface(
      const cell*,      ///< Cell
      const direction,  ///< outward normal to interface.
      const double      ///< cell diameter
      ) = 0;            ///< Returns Surface area of interface.

  ///
  /// Returns maximum of all gradients with neighbouring cells
  ///
  /// Pass in a point, and a primitive variable to calculate on, and it
  /// will return the maximum of all (absolute value) gradients between
  /// the cell and its neighbours.  These are one-sided first order
  /// differences, calculated w.r.t. all 2n neighbours.
  ///
  virtual double max_grad_abs(
      const cell*,  ///< current point.
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Primitive Variable to get gradient of.
      class GridBaseClass*  ///< pointer to computational grid.
      ) = 0;

  ///
  /// Vector Gradient Operator for Cartesian Coordinates.
  ///
  /// Uses second order central differenc formula:
  /// \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} =
  ///     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
  ///
  virtual void Gradient(
      const cell*,  ///< Cell to calculate gradient of
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Scalar (in state vector) to calculate
                    ///< gradient of.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put gradient vector.
      ) = 0;

  ///
  /// Calculate Vector Diff. Op. Divergence of a vector at a point.
  ///
  /// This is a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] for the x-component of the divergence.
  ///
  virtual double Divergence(
      cell*,       ///< point for which to calculate div(B)
      const int,   ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int*,  ///< Indices of vector quantity (in state vector) to
                   ///< calculate divergence of. Should contain 3
                   ///< elements, ordered as x,y,z components.
      class GridBaseClass*  ///< pointer to computational grid.
      ) = 0;

  ///
  /// Calculate Curl of Vector Quantity at a point.
  ///
  /// This uses a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] to calculate derivatives.
  ///
  virtual void Curl(
      const cell*,  ///< point for which to calculate curl
      const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int*,   ///< Indices of vector quantity (in state vector) to
                    ///< calculate curl of. Should contain 3 elements, ordered
                    ///< as x,y,z components.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put curl vector.
      ) = 0;

  /// central difference operator for a given variable.
  virtual double CentralDiff(
      class GridBaseClass*,  ///< pointer to computational grid.
      class cell*,           ///< point for which to calculate curl
      const int,             ///< axis along which to take difference
      const int,  ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int   ///< index in state vector of variable
      ) = 0;

  virtual double GradZone(
      class GridBaseClass*,  ///< pointer to computational grid.
      class cell*,           ///< pointer to cell
      const int,             ///< axis along which to take difference
      const int,  ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int   ///< index in state vector of variable
      ) = 0;
  ///
  /// Given a state, and a slope (dP/dx), construct an edge state.
  ///
  /// For first order accuracy, left and right edge states are just the
  /// cell centred value, P0.\n
  /// For second order accuracy, The vector of slopes is used to construct
  /// edge states \f$ P \pm (dP/dx)dx/2 \f$.
  ///
  /// \section Slope
  /// This function assumes that the slope array contains an already
  /// calculated gradient for the grid point in question.
  ///
  virtual int SetEdgeState(
      const cell*,      ///< Current Cell.
      const direction,  ///< Add or subtract the slope depending on direction.
      const int,        ///< length of state vectors.
      const pion_flt*,  ///< Slope vector.
      pion_flt*,        ///< vector for edge state.
      const int,        ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
      ) = 0;

  /// Calculate the slope (dP/dx) in cell i.
  virtual int SetSlope(
      const cell*,          ///< Current Cell.
      const axes,           ///< Which direction to calculate slope in.
      const int,            ///< length of state vectors.
      pion_flt*,            ///< Slope vector to be written to.
      const int,            ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
      ) = 0;

  ///
  /// Calculates the i-th term of the divergence for a vector of
  /// variables, given a cell and positive and negative flux vectors.
  ///
  virtual int DivStateVectorComponent(
      const cell*,           ///< current cell.
      class GridBaseClass*,  ///< pointer to computational grid.
      const axes,            ///< current coordinate axis we are looking along.
      const int,             ///< length of state vectors.
      const pion_flt*,       ///< Negative direction flux.
      const pion_flt*,       ///< Positive direction flux.
      pion_flt*              ///< Vector to assign divergence component to.
      ) = 0;

  ///
  /// Calculates the average of a left and right slope, using limiters.
  ///
  /// This uses Sam Falle's averaging procedure, which is generally
  /// known as the symmetric Van Albada limiter.
  ///
  double AvgFalle(
      const double,  ///< Left Slope.
      const double   ///< Right Slope.
  );

  ///
  /// Calculates Dot product of two vectors.
  ///
  /// This function takes in two vectors and their length, and returns
  /// their scalar or dot product.  Note this is a local quantity, so
  /// is independent of the coordinates used.
  ///
  double DotProduct(
      const pion_flt*,  ///< Vector 1.
      const pion_flt*,  ///< Vector 2.
      const int         ///< length of vectors.
  );

  ///
  /// Calculates Cross (aka vector) Product of two vectors.
  ///
  /// This function takes in two vectors and an int with their length,
  /// and calculates their vector or cross product, putting the result into
  /// the result vector array.  Note that the cross product is a local
  /// quantity, independent of the coordinate system used (as long as
  /// it is an orthonormal system).
  ///
  int CrossProduct(
      const pion_flt*,  ///< Vector 1.
      const pion_flt*,  ///< Vector 2.
      const int,        ///< length of vectors.
      pion_flt*         ///< Result vector
  );
};

///
/// Class to hold vector differential operations for Cartesian geometry
///
/// It is assumed that vector quantities are three-dimensional
/// (i.e. \f$ \vec{v} = (v_1,v_2,v_3) \f$), and the number of spatial dimensions
/// on the grid is passed into the constructor.
///
class VectorOps_Cart : virtual public BaseVectorOps {
protected:
  const int VOnd;  ///< Number of spatial dimensions in grid.
                   // double VOdx;      ///< Length of cell side.
                   // double VOdA;      ///< Cell surface area.
                   // double VOdV;      ///< Cell volume.

public:
  ///
  /// Constructor, sets number of spatial dimensions.
  ///
  VectorOps_Cart(const int  ///< number of spatial dimensions of grid
  );

  ///
  /// Destructor
  ///
  ~VectorOps_Cart();

  virtual double CellVolume(
      const cell*,  ///< cell pointer
      const double  ///< cell diameter
  );                ///< Returns Volume of cell.

  virtual double CellInterface(
      const cell*,      ///< Cell
      const direction,  ///< outward normal to interface.
      const double      ///< cell diameter
  );                    ///< Returns Surface area of interface.

  ///
  /// Returns maximum of all gradients with neighbouring cells
  ///
  /// Pass in a point, and a primitive variable to calculate on, and it
  /// will return the maximum of all (absolute value) gradients between
  /// the cell and its neighbours.  These are one-sided first order
  /// differences, calculated w.r.t. all 2n neighbours.
  ///
  virtual double max_grad_abs(
      const cell*,  ///< current point.
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Primitive Variable to get gradient of.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Vector Gradient Operator for Cartesian Coordinates.
  ///
  /// Uses second order central differenc formula:
  /// \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} =
  ///     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
  ///
  virtual void Gradient(
      const cell*,  ///< Cell to calculate gradient of
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Scalar (in state vector) to calculate
                    ///< gradient of.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put gradient vector.
  );

  ///
  /// Calculate Vector Diff. Op. Divergence of a vector at a point.
  ///
  /// This is a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] for the x-component of the divergence. \note It assumes cartesian
  /// geometry.
  ///
  virtual double Divergence(
      cell*,       ///< point for which to calculate div(B)
      const int,   ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int*,  ///< Indices of vector quantity (in state vector) to
                   ///< calculate divergence of. Should contain 3
                   ///< elements, ordered as x,y,z components.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Calculate Curl of Vector Quantity at a point.
  ///
  /// This uses a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] to calculate derivatives. \note It assumes cartesian geometry, and
  /// square/cubic cells.
  ///
  virtual void Curl(
      const cell*,  ///< point for which to calculate curl
      const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int*,   ///< Indices of vector quantity (in state vector) to
                    ///< calculate divergence of. Should contain 3 elements,
                    ///< ordered as x,y,z components.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put curl vector.
  );

  /// central difference operator for a given variable.
  virtual double CentralDiff(
      class GridBaseClass*,  ///< pointer to computational grid.
      class cell*,           ///< pointer to cell
      const int,             ///< axis along which to take difference
      const int,  ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int   ///< index in state vector of variable
  );

  virtual double GradZone(
      class GridBaseClass*,  ///< pointer to computational grid.
      class cell*,           ///< pointer to cell
      const int,             ///< axis along which to take difference
      const int,  ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int   ///< index in state vector of variable
  );
  ///
  /// Given a state, and a slope (dP/dx), construct an edge state.
  ///
  /// For first order accuracy, left and right edge states are just the
  /// cell centred value, P0.\n
  /// For second order accuracy, The vector of slopes is used to construct
  /// edge states \f$ P \pm (dP/dx)dx/2 \f$.
  ///
  /// \section Slope
  /// This function assumes that the slope array contains an already
  /// calculated gradient for the grid point in question.
  ///
  virtual int SetEdgeState(
      const cell*,      ///< Current Cell.
      const direction,  ///< Add or subtract the slope depending on direction.
      const int,        ///< length of state vectors.
      const pion_flt*,  ///< Slope vector.
      pion_flt*,        ///< vector for edge state.
      const int,        ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  /// Calculate the slope (dP/dx) in cell i.
  virtual int SetSlope(
      const cell*,          ///< Current Cell.
      const axes,           ///< Which direction to calculate slope in.
      const int,            ///< length of state vectors.
      pion_flt*,            ///< Slope vector to be written to.
      const int,            ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Calculates the i-th term of the divergence for a vector of
  /// variables, given a cell and positive and negative flux vectors.
  ///
  virtual int DivStateVectorComponent(
      const cell*,           ///< current cell.
      class GridBaseClass*,  ///< pointer to computational grid.
      const axes,            ///< current coordinate axis we are looking along.
      const int,             ///< length of state vectors.
      const pion_flt*,       ///< Negative direction flux.
      const pion_flt*,       ///< Positive direction flux.
      pion_flt*              ///< Vector to assign divergence component to.
  );
};

///
/// Class to hold vector differential operations for Cylindrical Coordinates.
///
/// It is assumed that vector quantities are three-dimensional
/// (i.e. \f$ \vec{v} = (v_1,v_2,v_3) \f$), and the number of spatial dimensions
/// on the grid is passed into the constructor.
///
/// \section com Centre Of Mass vs. Centre of Cell
/// The coordinates of each cell are at the midpoint of the cell in each
/// direction, but the state vectors are assumed to be at the centre of mass of
/// the cell, so for calculating gradients, the pivot points in the radial
/// direction should be the calculated centre of mass.  Centre of mass in z and
/// theta directions is of course the geometric centre of the cell.
///
class VectorOps_Cyl : virtual public VectorOps_Cart {
protected:
  ///
  /// Returns centre of mass of cell in R-direction.
  ///
  /// Centre of mass is \f[\int Rdm / \int dm = R_i+\frac{h^2}{12R_i} \;,\f]
  /// where h is the size of the cell in the R-direction, R is the radial
  /// coordinate, and \f$R_i = 0.5(R_{\mbox{min}} + R_{\mbox{max}})\f$ is the
  /// midpoint of cell i.
  ///
  virtual inline double R_com(
      const cell* c,  ///< Cell to operate on.
      const double dR)
  {
    return (CI.get_dpos(c, Rcyl) + dR * dR / 12. / CI.get_dpos(c, Rcyl));
  }

public:
  ///
  /// Constructor, sets number of spatial dimensions.
  ///
  VectorOps_Cyl(const int  ///< number of spatial dimensions of grid
  );

  ///
  /// Destructor
  ///
  ~VectorOps_Cyl();

  virtual double CellVolume(
      const cell*,  ///< cell pointer
      const double  ///< cell diameter
  );                ///< Returns Volume of cell.

  virtual double CellInterface(
      const cell*,      ///< Cell
      const direction,  ///< outward normal to interface.
      const double      ///< cell diameter
  );                    ///< Returns Surface area of interface.

  ///
  /// Returns maximum of all gradients with neighbouring cells for
  /// Cylindrical Coordinates.
  ///
  /// Pass in a point, and a primitive variable to calculate on, and it
  /// will return the maximum of all (absolute value) gradients between
  /// the cell and its neighbours.  These are one-sided first order
  /// differences, calculated w.r.t. all 2n neighbours.
  ///
  virtual double max_grad_abs(
      const cell*,  ///< current point.
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Primitive Variable to get gradient of.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Vector Gradient Operator for Cylindrical Coordinates.
  ///
  /// Uses second order central differenc formula:
  /// \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} =
  ///     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
  ///
  virtual void Gradient(
      const cell*,  ///< Cell to calculate gradient of
      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int,    ///< Index of Scalar (in state vector) to calculate
                    ///< gradient of.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put gradient vector.
  );

  ///
  /// Calculate Vector Diff. Op. Divergence of a vector at a point for
  /// Cylindrical Coordinates.
  ///
  /// This is a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] for the x-component of the divergence.
  ///
  /// It assumes cylindrical coordinates in the order (z,R,theta), and
  /// includes the relevant scale factors for each derivative.
  /// */
  virtual double Divergence(
      cell*,       ///< point for which to calculate div(B)
      const int,   ///< Which State Vector to take scalar from (P=0,Ph=1)
      const int*,  ///< Indices of vector quantity (in state vector) to
                   ///< calculate divergence of. Should contain 3
                   ///< elements, ordered as x,y,z components.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Calculate Curl of Vector Quantity at a point for Cylindrical
  /// Coordinates.
  ///
  /// This uses a second order cell-centred finite difference
  /// \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x}
  /// \f] to calculate derivatives.
  ///
  /// It assumes cylindrical coordinates in the order (z,R,theta), and
  /// includes the relevant scale factors for each derivative. THIS FUNCTION
  /// DOESN'T WORK!!!
  ///
  virtual void Curl(
      const cell*,  ///< point for which to calculate curl
      const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int*,   ///< Indices of vector quantity (in state vector) to
                    ///< calculate divergence of. Should contain 3 elements,
                    ///< ordered as x,y,z components.
      class GridBaseClass*,  ///< pointer to computational grid.
      pion_flt*              ///< Pointer to array to put curl vector.
  );

  ///
  /// Given a cell-centred state vector, and a slope (dP/dx),
  /// construct an edge state vector.
  ///
  /// For first order accuracy, left and right edge states are just the
  /// cell centred value, P0.\n
  /// For second order accuracy, The vector of slopes is used to construct
  /// edge states \f$ P \pm (dP/dx)dx/2 \f$.
  ///
  /// \section Slope
  /// This function assumes that the slope array contains an already
  /// calculated gradient for the grid point in question.
  ///
  int SetEdgeState(
      const cell*,      ///< Current Cell.
      const direction,  ///< Add or subtract the slope depending on direction.
      const int,        ///< length of state vectors.
      const pion_flt*,  ///< Slope vector.
      pion_flt*,        ///< vector for edge state.
      const int,        ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  /// Calculate the slope (dP/dx) in cell i.
  int SetSlope(
      const cell*,          ///< Current Cell.
      const axes,           ///< Which direction to calculate slope in.
      const int,            ///< length of state vectors.
      pion_flt*,            ///< Slope vector to be written to.
      const int,            ///< Order of spatial Accuracy.
      class GridBaseClass*  ///< pointer to computational grid.
  );

  ///
  /// Calculates the i-th term of the divergence for a vector of
  /// variables, given a cell and positive and negative flux vectors.
  ///
  int DivStateVectorComponent(
      const cell*,           ///< current cell.
      class GridBaseClass*,  ///< pointer to computational grid.
      const axes,            ///< current coordinate axis we are looking along.
      const int,             ///< length of state vectors.
      const pion_flt*,       ///< Negative direction flux.
      const pion_flt*,       ///< Positive direction flux.
      pion_flt*              ///< Vector to assign divergence component to.
  );
};

#endif  // VECTOROPS_H
