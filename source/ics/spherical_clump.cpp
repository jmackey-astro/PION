///
/// \file spherical_clump.cc
/// \author Jonathan Mackey
/// \date 2011.03.24
///
/// Sets up a spherically symmetric clump at the origin, in a 1D spherical grid,
/// or 2D axisymmetric grid.
///
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "coord_sys/VectorOps.h"
#include "ics/icgen.h"
#include "ics/icgen_base.h"

IC_spherical_clump::IC_spherical_clump() {}
IC_spherical_clump::~IC_spherical_clump() {}

int IC_spherical_clump::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!", ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams", rp);

  string seek, str;

  //
  // Ambient medium density (cgs)
  //
  seek = "AMB_density";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::AMB_density = atof(str.c_str());

  //
  // Ambient medium pressure (cgs)
  //
  seek = "AMB_pressure";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::AMB_pressure = atof(str.c_str());

  //
  // Radius of clump, in cgs units.  For Top-hat this is the radius,
  // for 1/r^2 this is the scale radius of the constant density core
  // (rho/(1+(r/rc)^2)) for gaussian this is sigma in the exp(-r^2/(2sigma^2))
  //
  seek = "SC_radius";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::SC_rad = atof(str.c_str());

  //
  // type of pressure profile :1=isothermal, 2=constant pressure.
  //
  seek = "SC_pressure";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::SC_pressure_profile = atoi(str.c_str());

  //
  // overdensity at centre of clump.
  //
  seek = "SC_overdensity";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::SC_overdensity = atof(str.c_str());

  //
  // Density profile, 0=top-hat, 1=1/r^2, 2=gaussian.
  //
  seek = "SC_density_profile";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  IC_spherical_clump::SC_density_profile = atoi(str.c_str());

  seek = "SC_magfieldX";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_spherical_clump::SC_BX = 0.0;  // could have no field.
  else
    IC_spherical_clump::SC_BX = atof(str.c_str());
  seek = "SC_magfieldY";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_spherical_clump::SC_BY = 0.0;  // could have no field.
  else
    IC_spherical_clump::SC_BY = atof(str.c_str());
  seek = "SC_magfieldZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_spherical_clump::SC_BZ = 0.0;  // could have no field.
  else
    IC_spherical_clump::SC_BZ = atof(str.c_str());

#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  SC_BX /= sqrt(4.0 * M_PI);
  SC_BY /= sqrt(4.0 * M_PI);
  SC_BZ /= sqrt(4.0 * M_PI);
#endif

  IC_spherical_clump::gam = SimPM->gamma;

  IC_spherical_clump::eqns = SimPM->eqntype;
  if (eqns == EQEUL)
    eqns = 1;
  else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
    eqns = 2;
  else
    rep.error("Bad equations", eqns);

  // now make sure we are to do a blast wave sim.
  string ics = rp->find_parameter("ics");

  if (ics == "")
    rep.error("didn't get any ics to set up.", ics);
  else if (ics == "Clump_Spherical" && SimPM->coord_sys == COORD_SPH) {
    cout << "Setting up Spherically symmetric 1D cloud.\n";
    err += setup_clump();
  }
  else if (ics == "Clump_Axisymmetric" && SimPM->coord_sys == COORD_CYL) {
    cout << "Setting up Axisymmetric 2D cloud.\n";
    err += setup_clump();
  }
  else
    rep.error("Don't know what Initial Condition is!", ics);

  // Add noise to data?  Smooth data?
  double noise = 0.0;
  int smooth   = 0;
  ics          = rp->find_parameter("noise");
  if (ics != "")
    noise = atof(ics.c_str());
  else
    noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number", noise);
  if (noise > 0) err += AddNoise2Data(gg, *SimPM, 2, noise);

  ics = rp->find_parameter("smooth");
  if (ics != "")
    smooth = atoi(ics.c_str());
  else
    smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number", smooth);
  if (smooth > 0) err += SmoothData(smooth);

  return err;
}

int IC_spherical_clump::setup_clump()
{
  int ndim = gg->Ndim();
  if (ndim > 2) rep.error("Bad ndim in setup spherical Clump 1D", ndim);
  cout << "Setting up spherically symmetric clump with radius " << SC_rad;
  cout << " with an overdensity of " << SC_overdensity << " in an ambient";
  cout << " medium with rho=" << AMB_density << ", p=" << AMB_pressure << "\n";

  //
  // Centred at [0]
  //
  double centre[ndim];
  for (int i = 0; i < ndim; i++)
    centre[0] = 0.0;

  //
  // Data.
  //
  cout << "Assigning primitive vectors.\n";
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    cpt->P[RO] = AMB_density;
    cpt->P[PG] = AMB_pressure;
    cpt->P[VX] = cpt->P[VY] = cpt->P[VZ] =
        0.0;  // Stationary gas to start with.
    if (eqns == 2) {
      cpt->P[BX] = SC_BX;
      cpt->P[BY] = SC_BY;
      cpt->P[BZ] = SC_BZ;
    }
    for (int i = 0; i < SimPM->ntracer; i++) {
      cpt->P[SimPM->ftr + i] = 0.0;
    }
    //
    // This is where I set the state inside the clump.  We overlay the clump
    // on top of the ambient medium.
    //
    // CI.get_dpos(cpt,dpos);

    switch (SC_density_profile) {
      case 0:  // top-hat
        if (gg->distance_vertex2cell(centre, cpt) <= SC_rad)
          cpt->P[RO] *= SC_overdensity;
        break;

      case 1:  // 1/r^2 with a core.
        cpt->P[RO] *=
            SC_overdensity
            / (1.0 + pow(gg->distance_vertex2cell(centre, cpt) / SC_rad, 2.0));
        break;

      case 2:  // Gaussian
        cpt->P[RO] *=
            SC_overdensity
            * exp(-0.5
                  * pow(gg->distance_vertex2cell(centre, cpt) / SC_rad, 2.0));
        break;

      default:
        rep.error("Bad density profile in spherical_clump", SC_density_profile);
        break;
    }

    switch (SC_pressure_profile) {
      case 1:
        // isothermal, so p/rho = constant.
        // so if the density is 100x higher, so must the pressure be.
        cpt->P[PG] *= cpt->P[RO] / AMB_density;
        break;

      case 2:
        // constant pressure equal to the ambient pressure.  So do
        // nothing.
        break;

      default:
        rep.error("Bad pressure profile", SC_pressure_profile);
        break;
    }
    //       cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
  } while ((cpt = gg->NextPt(cpt)) != 0);
  cout << "Got through data successfully.\n";
  // Data done.
  return (0);
}
