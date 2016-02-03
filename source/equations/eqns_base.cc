///
/// \file eqns_base.cc
/// \author Jonathan Mackey
/// Some functions in the base class for equations which are common to all solvers.
///
/// Modifications:
///
///  - 2009-12-22 JM: Added comments to rotation to make clear that we
///    are moving the vector quantity through +theta, NOT rotating the
///    coord.sys. by +theta.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to
///    flux functions).
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.12.23 JM: Added eq_refvec[] to eqns_base.  Deleted
///   riemann_base altogether.  Moved flux_solver_base class
///   definitions to flux_base.cc
/// - 2011.01.03 JM: Added eq_posdir and eq_negdir direction variables.
/// - 2015.01.14 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"

#include "eqns_base.h"
using namespace std;


// ##################################################################
// ##################################################################



eqns_base::eqns_base(const int n ///< Number of Variables in State Vector
		     )
  : eq_nvar(n)
{
  //
  // We assume at least 5 variables in the state vector -- density, pressure, and
  // a 3d velocity
  //
  if (eq_nvar<5) rep.error("Bad eq_nvar; must be >=5!",eq_nvar);

  //
  // These are constants regardless of what way we are looking
  //
  eqRO = RO;  eqRHO=RHO;
  eqPG = PG;  eqERG=ERG;
  //
  // Initial direction is set to X
  //
  eq_dir = XX;
  eq_posdir = XP; eq_negdir = XN;

  eqVX = VX; eqVY = VY; eqVZ = VZ;
  eqMMX=MMX; eqMMY=MMY; eqMMZ=MMZ;
  eqBX = BX; eqBY = BY; eqBZ = BZ;
  eqBBX=BBX; eqBBY=BBY; eqBBZ=BBZ;

  eq_gamma=-1.0e99;

  //
  // Allocate memory for reference vector:
  //
  eq_refvec = mem.myalloc(eq_refvec,eq_nvar);
  for (int v=0;v<eq_nvar;v++)
    eq_refvec[v] = 0.0;

  return;
}


// ##################################################################
// ##################################################################



eqns_base::~eqns_base()
{
  eq_refvec = mem.myfree(eq_refvec);
  return;
}


// ##################################################################
// ##################################################################




void eqns_base::SetDirection(const enum axes d)
{
  //
  // If we don't need to change direction, just return:
  //
  if (d==eq_dir) return;

  //
  // Else reset all the direction--dependent quantities.
  //
  eq_dir = d;
  switch (eq_dir) {
  case XX:
    eq_posdir = XP; eq_negdir = XN;
    eqVX = VX; eqVY = VY; eqVZ = VZ;
    eqMMX=MMX; eqMMY=MMY; eqMMZ=MMZ;
    eqBX = BX; eqBY = BY; eqBZ = BZ;
    eqBBX=BBX; eqBBY=BBY; eqBBZ=BBZ;
    break;
  case YY:
    eq_posdir = YP; eq_negdir = YN;
    eqVX = VY; eqVY = VZ; eqVZ = VX;
    eqMMX=MMY; eqMMY=MMZ; eqMMZ=MMX;
    eqBX = BY; eqBY = BZ; eqBZ = BX;
    eqBBX=BBY; eqBBY=BBZ; eqBBZ=BBX;
    break;
  case ZZ:
    eq_posdir = ZP; eq_negdir = ZN;
    eqVX = VZ; eqVY = VX; eqVZ = VY;
    eqMMX=MMZ; eqMMY=MMX; eqMMZ=MMY;
    eqBX = BZ; eqBY = BX; eqBZ = BY;
    eqBBX=BBZ; eqBBY=BBX; eqBBZ=BBY;
    break;
   default:
     rep.error("bad direction in eqns_base::SetDirection",d);
     break;
  }
  return;
}


// ##################################################################
// ##################################################################



enum axes eqns_base::GetDirection()
{
  return eq_dir;
}


// ##################################################################
// ##################################################################



void eqns_base::rotate(
      pion_flt *vec, ///< State vector
      enum axes initdir, ///< Initial orientation.
      enum axes finaldir ///< Final Orientation.
      )
{
  ///
  /// \section Directions
  /// This only rotates between the three positive directions, XP,YP,ZP, 
  /// so it never introduces sign changes to the elements, just reordering
  /// of the vector components.
  ///
  /// \section Equations
  /// This only rotates the velocity components, since all of my equations
  /// classes I use have at least a 3D velocity vector.
  ///
  if (initdir==finaldir) return;
  double temp[eq_nvar];
  for (int i=0;i<eq_nvar;i++) temp[i] = vec[i];
  int offset = (static_cast<int>(finaldir-initdir+3))%3;
  //  if(offset!=1 && offset!=2) rep.error("rotate function broken.",offset);

  //
  // Only two options, a positive permutation of indices, or a
  // negative permutation.  We keep a right-handed system, so 
  // X->Y or X->Z define the only possible permutations.
  //
  if (offset ==1) {
    temp[eqVX] = vec[eqVY];
    temp[eqVY] = vec[eqVZ];
    temp[eqVZ] = vec[eqVX];
  }
  else if (offset==2) {
    temp[eqVX] = vec[eqVZ];
    temp[eqVY] = vec[eqVX];
    temp[eqVZ] = vec[eqVY];
  }

  //
  // Now copy the new ordering back into vec and return
  //
  vec[eqVX] = temp[eqVX];
  vec[eqVY] = temp[eqVY];
  vec[eqVZ] = temp[eqVZ];
  return;
}


// ##################################################################
// ##################################################################



void eqns_base::rotateXY(
      pion_flt *v, ///< State vector
      double theta ///< rotation angle
      )
{
  //
  // This rotates the 3D vector fields by an angle +theta around the z-axis
  // This is equivalent to rotating the coordinate system by an angle -theta
  //
  double ct=cos(theta); double st=sin(theta);
  double vx = v[eqVX]*ct - v[eqVY]*st;
  double vy = v[eqVX]*st + v[eqVY]*ct;
  v[eqVX] = vx; v[eqVY]=vy;
  return;
}


// ##################################################################
// ##################################################################



void eqns_base::PtoFlux(
      const pion_flt *p,
      pion_flt *f,
      const double gamma)
{
   pion_flt u[eq_nvar];
   PtoU(p, u, gamma);
   PUtoFlux(p,u,f);
   return;
}

