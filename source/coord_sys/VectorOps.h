/** \file VectorOps.h
 * 
 * \brief Declares VectorOps utility classes
 * 
 * \author Jonathan Mackey
 * 
 * 
 * This declares classes of functions which act on a grid of cells, but which depend
 * on the coordinate system being used (for example the non-local vector
 * differential operators grad, div, and curl).  The base class is pure
 * virtual, and the derived classes have names VectorOps_Cart, VectorOps_Cyl,
 * etc. for the different coordinate systems.
 * 
 * Modified:\n
 *  - 2007-08-01 File Created
 * */
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///
/// - 2010.12.04 JM: Added constructor with only one argument.  Also
///   a set_dx() function.
///


#ifndef VECTOROPS_H
#define VECTOROPS_H

//
// These tells code what to compile and what to leave out.
//
#include "../defines/functionality_flags.h"
#include "../defines/testing_flags.h"


#include "../global.h"

//#define FIX_TRACER_SLOPES_TO_DENSITY ///< This sets tracers to have zero slope across cell.

#define VOdz VOdx

/** \brief Class to hold vector differential operations which depend on 
 * the coordinate system.  This is a virtual base class. */
class BaseVectorOps {
  public:
   virtual ~BaseVectorOps();
   virtual double CellVolume(const cell *)=0; ///< Returns Volume of cell.
   virtual double CellInterface(const cell *, ///< Cell
				const direction ///< outward normal to interface.
				)=0; ///< Returns Surface area of interface.
   /** \brief Returns maximum of all gradients with neighbouring cells
    * 
    * Pass in a point, and a primitive variable to calculate on, and it
    * will return the maximum of all (absolute value) gradients between 
    * the cell and its neighbours.  These are one-sided first order 
    * differences, calculated w.r.t. all 2n neighbours.
    */
   virtual double maxGradAbs(const cell *,      ///< current point.
			     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
			     const int     ///< Index of Primitive Variable to get gradient of.
			     )=0;
   /** \brief Vector Gradient Operator for Cartesian Coordinates.
    * 
    * Uses second order central differenc formula:
    * \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} = 
    *     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
    */
   virtual void Grad(const cell *, ///< Cell to calculate gradient of
		     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		     const int,    ///< Index of Scalar (in state vector) to calculate gradient of.
		     double *      ///< Pointer to array to put gradient vector.
		     )=0;
	     
   /** \brief Calculate Vector Diff. Op. Divergence of a vector at a point.
    * 
    * This is a second order cell-centred finite difference
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * for the x-component of the divergence.  
    * */
   virtual double Div(const cell *, ///< point for which to calculate div(B)
		      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		      const int *   ///< Indices of vector quantity (in state vector) to calculate divergence of.
		                    ///< Should contain 3 elements, ordered as x,y,z components.
		      )=0;
   /** \brief Calculate Curl of Vector Quantity at a point.
    * 
    * This uses a second order cell-centred finite difference 
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * to calculate derivatives.
    * */
   virtual void Curl(const cell *, ///< point for which to calculate curl
		     const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
		     const int *,  ///< Indices of vector quantity (in state vector) to calculate curl of.
		                   ///< Should contain 3 elements, ordered as x,y,z components.
		     double *      ///< Pointer to array to put curl vector.
		     )=0;
   /**  \brief Given a state, and a slope (dP/dx), construct an edge state.
    * 
    * For first order accuracy, left and right edge states are just the
    * cell centred value, P0.\n
    * For second order accuracy, The vector of slopes is used to construct
    * edge states \f$ P \pm (dP/dx)dx/2 \f$.
    * 
    * \section Slope
    * This function assumes that the slope array contains an already calculated
    * gradient for the grid point in question.
    */
   virtual int SetEdgeState(const cell *,    ///< Current Cell.
			    const direction, ///< Add or subtract the slope depending on direction.
			    const int,  ///< length of state vectors.
			    const double *,  ///< Slope vector.
			    double *,        ///< vector for edge state. 
			    const int        ///< Order of spatial Accuracy.
			    )=0;
   /** \brief Calculate the slope (dP/dx) in cell i. */
   virtual int SetSlope(const cell *, ///< Current Cell.
			const axes,   ///< Which direction to calculate slope in.
			const int,  ///< length of state vectors.
			double *,     ///< Slope vector to be written to.
			const int     ///< Order of spatial Accuracy.
			)=0;
   /** \brief Calculates the i-th term of the divergence for a vector of 
    * variables, given a cell and positive and negative flux vectors.
    * */
   virtual int DivStateVectorComponent(const cell *, ///< current cell.
				       const axes, ///< current coordinate axis we are looking along.
				       const int,  ///< length of state vectors.
				       const double *, ///< Negative direction flux.
				       const double *, ///< Positive direction flux.
				       double * ///< Vector to assign divergence component to.
				       )=0;
   /** \brief Calculates the average of a left and right slope, using limiters.
    * 
    * This uses Sam Falle's averaging procedure.
    * */
   double AvgFalle(const double, ///< Left Slope.
		   const double  ///< Right Slope.
		   );
   /** \brief Calculates Dot product of two vectors.
    * 
    * This function takes in two vectors and their length, and returns
    * their scalar or dot product.  Note this is a local quantity, so 
    * is independent of the coordinates used.
    */
   double DotProduct(const double *, ///< Vector 1.
		     const double *, ///< Vector 2.
		     const int       ///< length of vectors.
		     );
   /** \brief Calculates Cross (aka vector) Product of two vectors. 
    * 
    * This function takes in two vectors and an int with their length,
    * and calculates their vector or cross product, putting the result into the
    * result vector array.  Note that the cross product is a local 
    * quantity, independent of the coordinate system used (as long as
    * it is an orthonormal system).
    */
   int CrossProduct(const double *, ///< Vector 1.
		    const double *, ///< Vector 2.
		    const int,      ///< length of vectors.
		    double *        ///< Result vector
		    );
};


/** \brief Class to hold vector differential operations for Cartesian geometry 
 * 
 * It is assumed that vector quantities are three-dimensional 
 * (i.e. \f$ \vec{v} = (v_1,v_2,v_3) \f$), and the number of spatial dimensions
 * on the grid is passed into the constructor.
 */
class VectorOps_Cart : virtual public BaseVectorOps
{
  protected:
   const int VOnd;   ///< Number of spatial dimensions in grid.
   double VOdx;      ///< Length of cell side.
   bool have_set_dx; ///< set to true once VOdx has been set.
   double VOdA;      ///< Cell surface area.
   double VOdV;      ///< Cell volume.
  public:
   ///
   /// constructor, sets number of spatial dimensions of grid.
   ///
   VectorOps_Cart(int, ///< Number of spatial dimensions of grid
		  double ///< Length of cell in all directions (cells must be cubic).
		  );

   ///
   /// Alternate constructor which doesn't take delta-x as a parameter.
   ///
   VectorOps_Cart(const int ///< number of spatial dimensions of grid
		  );

   ///
   /// Specify delta-x for the grid (if a uniform grid!)
   ///
   virtual void set_dx(const double ///< size of grid cell.
		       );
   
   ///
   /// Destructor
   ///
   ~VectorOps_Cart();

   virtual double CellVolume(const cell *); ///< Returns Volume of cell.
   virtual double CellInterface(const cell *,   ///< Cell
				const direction ///< outward normal to interface.
				); ///< Returns Surface area of interface.
   /** \brief Returns maximum of all gradients with neighbouring cells
    * 
    * Pass in a point, and a primitive variable to calculate on, and it
    * will return the maximum of all (absolute value) gradients between 
    * the cell and its neighbours.  These are one-sided first order 
    * differences, calculated w.r.t. all 2n neighbours.
    */
   virtual double maxGradAbs(const cell *,      ///< current point.
			     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
			     const int     ///< Index of Primitive Variable to get gradient of.
			     );
   /** \brief Vector Gradient Operator for Cartesian Coordinates.
    * 
    * Uses second order central differenc formula:
    * \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} = 
    *     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
    */
   virtual void Grad(const cell *, ///< Cell to calculate gradient of
		     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		     const int,    ///< Index of Scalar (in state vector) to calculate gradient of.
		     double *      ///< Pointer to array to put gradient vector.
		     );
	     
   /** \brief Calculate Vector Diff. Op. Divergence of a vector at a point.
    * 
    * This is a second order cell-centred finite difference
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * for the x-component of the divergence.  
    * \note It assumes cartesian geometry.
    * */
   virtual double Div(const cell *, ///< point for which to calculate div(B)
		      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		      const int *   ///< Indices of vector quantity (in state vector) to calculate divergence of.
		                    ///< Should contain 3 elements, ordered as x,y,z components.
		      );
   /** \brief Calculate Curl of Vector Quantity at a point.
    * 
    * This uses a second order cell-centred finite difference 
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * to calculate derivatives.
    * \note It assumes cartesian geometry, and square/cubic cells.
    * */
   virtual void Curl(const cell *, ///< point for which to calculate curl
		     const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
		     const int *,  ///< Indices of vector quantity (in state vector) to calculate divergence of.
		                   ///< Should contain 3 elements, ordered as x,y,z components.
		     double *      ///< Pointer to array to put curl vector.
		     );
   /**  \brief Given a state, and a slope (dP/dx), construct an edge state.
    * 
    * For first order accuracy, left and right edge states are just the
    * cell centred value, P0.\n
    * For second order accuracy, The vector of slopes is used to construct
    * edge states \f$ P \pm (dP/dx)dx/2 \f$.
    * 
    * \section Slope
    * This function assumes that the slope array contains an already calculated
    * gradient for the grid point in question.
    */
   virtual int SetEdgeState(const cell *,    ///< Current Cell.
		    const direction, ///< Add or subtract the slope depending on direction.
		    const int,  ///< length of state vectors.
		    const double *,  ///< Slope vector.
		    double *,        ///< vector for edge state. 
		    const int        ///< Order of spatial Accuracy.
		    );
   /** \brief Calculate the slope (dP/dx) in cell i. */
   virtual int SetSlope(const cell *, ///< Current Cell.
		const axes,   ///< Which direction to calculate slope in.
		const int,  ///< length of state vectors.
		double *,     ///< Slope vector to be written to.
		const int     ///< Order of spatial Accuracy.
		);
   /** \brief Calculates the i-th term of the divergence for a vector of 
    * variables, given a cell and positive and negative flux vectors.
    * */
   virtual int DivStateVectorComponent(const cell *, ///< current cell.
			       const axes, ///< current coordinate axis we are looking along.
			       const int,  ///< length of state vectors.
			       const double *, ///< Negative direction flux.
			       const double *, ///< Positive direction flux.
			       double * ///< Vector to assign divergence component to.
			       );
};

/** \brief Class to hold vector differential operations for Cylindrical Coordinates.
 * 
 * It is assumed that vector quantities are three-dimensional 
 * (i.e. \f$ \vec{v} = (v_1,v_2,v_3) \f$), and the number of spatial dimensions
 * on the grid is passed into the constructor.
 * 
 * \section com Centre Of Mass vs. Centre of Cell
 * The coordinates of each cell are at the midpoint of the cell in each direction,
 * but the state vectors are assumed to be at the centre of mass of the cell, so
 * for calculating gradients, the pivot points in the radial direction should be
 * the calculated centre of mass.  Centre of mass in z and theta directions is of
 * course the geometric centre of the cell.
 *
 * \todo Since switching cell positions to integer values, I think I could gain 
 * a good bit in terms of efficiency by not calling double precision positions 
 * all the time in the div/curl calculations.  So I would like to do that eventually,
 * but it's a bit low on priorities at the moment.
 */
class VectorOps_Cyl : virtual public VectorOps_Cart
{
  protected:
   double VOdR; ///< length of cell in R-dir.
   /** \brief Returns centre of mass of cell in R-direction.
    * 
    * Centre of mass is \f[\int Rdm / \int dm = R_i+\frac{h^2}{12R_i} \;,\f]
    * where h is the size of the cell in the R-direction, R is the radial
    * coordinate, and \f$R_i = 0.5(R_{\mbox{min}} + R_{\mbox{max}})\f$ is the
    * midpoint of cell i.
    */
   virtual inline double R_com(const cell *c ///< Cell to operate on.
			       )
   {return(CI.get_dpos(c,Rcyl) + VOdR*VOdR/12./CI.get_dpos(c,Rcyl));}

  public:
   ///
   /// constructor, sets number of spatial dimensions of grid.
   ///
   VectorOps_Cyl(int, ///< Number of spatial dimensions of grid
		 double ///< Length of cell in all directions (cells must be cubic).
		 );

   ///
   /// Alternate constructor which doesn't take delta-x as a parameter.
   ///
   VectorOps_Cyl(const int ///< number of spatial dimensions of grid
		 );

   ///
   /// Specify delta-x for the grid (if a uniform grid!)
   ///
   virtual void set_dx(const double ///< size of grid cell.
		       );

   ///
   /// Destructor
   ///
   ~VectorOps_Cyl();

   virtual double CellVolume(const cell *); ///< Returns Volume of cell.
   virtual double CellInterface(const cell *,   ///< Cell
				const direction ///< outward normal to interface.
				); ///< Returns Surface area of interface.
   /** \brief Returns maximum of all gradients with neighbouring cells for Cylindrical Coordinates.
    * 
    * Pass in a point, and a primitive variable to calculate on, and it
    * will return the maximum of all (absolute value) gradients between 
    * the cell and its neighbours.  These are one-sided first order 
    * differences, calculated w.r.t. all 2n neighbours.
    */
   virtual double maxGradAbs(const cell *, ///< current point.
			     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
			     const int     ///< Index of Primitive Variable to get gradient of.
			     );
   /** \brief Vector Gradient Operator for Cylindrical Coordinates.
    * 
    * Uses second order central differenc formula:
    * \f[ \left( \frac{\partial \phi}{\partial x} \right)_{x_i} = 
    *     \frac{\phi(x_{i+1})-\phi(x_{i-1})}{2\delta x} \f]
    */
   virtual void Grad(const cell *, ///< Cell to calculate gradient of
		     const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		     const int,    ///< Index of Scalar (in state vector) to calculate gradient of.
		     double *      ///< Pointer to array to put gradient vector.
		     );
	     
   /** \brief Calculate Vector Diff. Op. Divergence of a vector at a point for Cylindrical Coordinates.
    * 
    * This is a second order cell-centred finite difference
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * for the x-component of the divergence.
    * 
    * It assumes cylindrical coordinates in the order (z,R,theta), and includes
    * the relevant scale factors for each derivative.
    * */
   virtual double Div(const cell *, ///< point for which to calculate div(B)
		      const int,    ///< Which State Vector to take scalar from (P=0,Ph=1)
		      const int *   ///< Indices of vector quantity (in state vector) to calculate divergence of.
		                    ///< Should contain 3 elements, ordered as x,y,z components.
		      );
   /** \brief Calculate Curl of Vector Quantity at a point for Cylindrical Coordinates.
    * 
    * This uses a second order cell-centred finite difference 
    * \f[ \frac{\partial B_i}{\partial x} = \frac{B_{i+1}-B_{i-1}}{2\delta x} \f]
    * to calculate derivatives.
    * 
    * It assumes cylindrical coordinates in the order (z,R,theta), and includes
    * the relevant scale factors for each derivative. THIS FUNCTION DOESN'T WORK!!!
    * */
   virtual void Curl(const cell *, ///< point for which to calculate curl
		     const int,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
		     const int *,  ///< Indices of vector quantity (in state vector) to calculate divergence of.
		                   ///< Should contain 3 elements, ordered as x,y,z components.
		     double *      ///< Pointer to array to put curl vector.
		     );
   /**  \brief Given a state, and a slope (dP/dx), construct an edge state.
    * 
    * For first order accuracy, left and right edge states are just the
    * cell centred value, P0.\n
    * For second order accuracy, The vector of slopes is used to construct
    * edge states \f$ P \pm (dP/dx)dx/2 \f$.
    * 
    * \section Slope
    * This function assumes that the slope array contains an already calculated
    * gradient for the grid point in question.
    */
   int SetEdgeState(const cell *,    ///< Current Cell.
		    const direction, ///< Add or subtract the slope depending on direction.
		    const int,       ///< length of state vectors.
		    const double *,  ///< Slope vector.
		    double *,        ///< vector for edge state. 
		    const int        ///< Order of spatial Accuracy.
		    );
   /** \brief Calculate the slope (dP/dx) in cell i. */
   int SetSlope(const cell *, ///< Current Cell.
		const axes,   ///< Which direction to calculate slope in.
		const int ,   ///< length of state vectors.
		double *,     ///< Slope vector to be written to.
		const int     ///< Order of spatial Accuracy.
		);
   /** \brief Calculates the i-th term of the divergence for a vector of 
    * variables, given a cell and positive and negative flux vectors.
    * */
   int DivStateVectorComponent(const cell *, ///< current cell.
			       const axes, ///< current coordinate axis we are looking along.
			       const int,  ///< length of state vectors.
			       const double *, ///< Negative direction flux.
			       const double *, ///< Positive direction flux.
			       double * ///< Vector to assign divergence component to.
			       );
};

#endif // VECTOROPS_H
