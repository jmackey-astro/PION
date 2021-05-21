/// \file photoevaporating_multiclumps.cc
/// \author Jonathan Mackey
///
/// File for setting up photo-evaporation of many random clumps and some
/// strategically positioned clumps.
///
/// Modifications:\n
///  - 2010-01-20 JM: line 714, corrected mass to be rho*(1+delta) [from
///  rho*delta]
///  - 2010-01-27 JM: comments
///  - 2010-02-08 JM: Changed tracer calculation in clumps_set_dens() so that
///      it gives the right answer (if all clumps have the same tracer value).
///  - 2010-02-10 JM: Added isothermal euler equations.
/// - 2010.12.15 JM: Added get_alternate_ambient_params(rp, &amb_data) and
///    add_alternate_ambient_data_to_grid(ggg, &amb_data) to allow for the
///    nearest 10th of the grid to the negative x-boundary to have a
///    different ambient medium property.  Note this only works well if all
///    clumps are well away from this region (otherwise the tracer values
///    get overwritten with the other ambient medium values which can have
///    unintended consequences!)
/// - 2011.01.06 JM: fixed memory leak.
/// - 2012.02.25 JM: Added optional velocity vector for clump.
/// - 2012.04.23 JM: Added optional density gradient in ambient medium.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "ics/icgen.h"
#include "ics/icgen_base.h"
#include <algorithm>
#include <sstream>

// ##################################################################
// ##################################################################

IC_photevap_multi_clumps::IC_photevap_multi_clumps()
{
  ndim = coords = eqns = -1;
  ambdens = gamma = ambdivider = -1.0e99;
  rc_data.used                 = false;
  sc_data.used                 = false;
  amb_data.used                = false;
}

// ##################################################################
// ##################################################################

IC_photevap_multi_clumps::~IC_photevap_multi_clumps()
{
  if (rc_data.used) {
    rc_data.border = mem.myfree(rc_data.border);
    rc_data.cl     = mem.myfree(rc_data.cl);
  }

  if (sc_data.used) {
    sc_data.cl = mem.myfree(sc_data.cl);
  }

  if (amb_data.used) {
    amb_data.ambient = mem.myfree(amb_data.ambient);
  }
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  //
  // Set pointers to grid and read-params classes.
  //
  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!", ggg);
  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams", rp);

  //
  // Get general sim parameters:
  //
  ndim = SimPM->ndim;
  if (ndim != 1 && ndim != 2 && ndim != 3)
    rep.error("Photoevaporation problem must be 1-3D", ndim);
  coords = SimPM->coord_sys;
  if (coords != COORD_CRT) {
    cout << "WARNING! using cylindrical coords, so MAKE SURE SOURCE AND ALL "
            "CLUMPS ARE ON-AXIS!!\n";
    // rep.error("Bad coord sys",coords);
  }
  eqns = SimPM->eqntype;
  if (eqns == EQEUL || eqns == EQEUL_ISO)
    eqns = 1;
  else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
    eqns = 2;
  else
    rep.error("Bad equations", eqns);

  //
  // Get ambient data from p-file, and write it to grid cells.
  //
  err += get_ambient_params(rp, &amb_data);
  err += add_ambient_data_to_grid(ggg, &amb_data);

  //
  // Get info for a portion of the domain which may have a different
  // ambient medium value (e.g. inflowing boundary from source).
  // Then reset the ambient data to the original value (in case it is
  // used later for anything).
  //
  err += get_alternate_ambient_params(rp, &amb_data);
  err += add_alternate_ambient_data_to_grid(ggg, &amb_data);
  err += get_ambient_params(rp, &amb_data);

  //
  // Get Random data from p-file, and add random clumps to grid.
  //
  err += get_random_clump_params(rp, &rc_data);
  err += add_random_clumps_to_grid(ggg, &rc_data);

  //
  // Get Strategic clumps from p-file, and add them to the grid.
  //
  err += get_strategic_clump_params(rp, &sc_data);
  err += add_strategic_clumps_to_grid(ggg, &sc_data);

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::get_ambient_params(
    class ReadParams *rparams, struct ambient_data *amb)
{
  int err = 0;
  if (!amb || !rparams)
    rep.error("Null pointer  passed to get_ambient_params!", amb);

  cout << "\t**** Getting Ambient Params...";
  if (!amb->used) {
    amb->used    = true;
    amb->ambient = mem.myalloc(amb->ambient, SimPM->nvar);
    for (int v = 0; v < SimPM->nvar; v++)
      amb->ambient[v] = 0.0;
  }

  string seek, str;
  ostringstream temp;
  string v;

  v = "RO";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[RO] = atof(str.c_str());
  else
    amb->ambient[RO] = -1.0e99;
  ambdens = amb->ambient[RO];

  v = "PG";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[PG] = atof(str.c_str());
  else
    amb->ambient[PG] = -1.0e99;

  v = "VX";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VX] = atof(str.c_str());
  else
    amb->ambient[VX] = -1.0e99;

  v = "VY";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VY] = atof(str.c_str());
  else
    amb->ambient[VY] = -1.0e99;

  v = "VZ";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VZ] = atof(str.c_str());
  else
    amb->ambient[VZ] = -1.0e99;

  if (eqns == 2) {  // mhd sim
    v = "BX";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BX] = atof(str.c_str());
    else
      amb->ambient[BX] = -1.0e99;

    v = "BY";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BY] = atof(str.c_str());
    else
      amb->ambient[BY] = -1.0e99;

    v = "BZ";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BZ] = atof(str.c_str());
    else
      amb->ambient[BZ] = -1.0e99;

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    amb->ambient[BX] /= sqrt(4.0 * M_PI);
    amb->ambient[BY] /= sqrt(4.0 * M_PI);
    amb->ambient[BZ] /= sqrt(4.0 * M_PI);
#endif

    if (SimPM->eqntype == EQGLM) amb->ambient[SI] = 0.;
  }  // if mhd vars

  // tracer variables
  for (int t = 0; t < SimPM->ntracer; t++) {
    temp.str("");
    temp << "PERC_ambTR" << t;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[t + SimPM->ftr] = atof(str.c_str());
    else
      amb->ambient[t + SimPM->ftr] = -1.0e99;
  }

  //
  // Radial slope is outwards from the origin, with a core radius of
  // "cloudradius" where density is contant, and a power-law slope
  // of "radialslope" beyond this distance.
  //
  temp.str("");
  temp << "PERCradialslope";
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->radial_profile = atoi(str.c_str());
  else
    amb->radial_profile = 0;
  cout << "\t\tradial_profile=" << amb->radial_profile << "\n";

  temp.str("");
  temp << "PERCcloudradius";
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "") {
    amb->cloudradius = atof(str.c_str());
    // radius is given in units of y-dir range.
    if (ndim > 1)
      amb->cloudradius *= SimPM->Range[YY];
    else
      amb->cloudradius *= SimPM->Range[XX];
    cout << "Cloud Radius in cm is " << amb->cloudradius << endl;
  }
  else
    amb->cloudradius = 0.0;
  cout << "\t\tcloudradius=" << amb->cloudradius << "\n";

  //
  // power-law slope in x-direction:
  //  rho(x) = rho_0 * [(x-x_0)/l_0)]^alpha
  //
  temp.str("");
  temp << "PERC_amb_xscale";
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str == "true") {
    amb->xscale = true;
    // cout <<"Getting X0, L0, A0 for x-scale: ";

    // get parameters rho_0, x_0, l_0, alpha (rho_0 = PERC_ambRO
    // all in cgs units.
    temp.str("");
    temp << "PERC_xscale_x0";
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->xscale_x0 = atof(str.c_str());
    else
      amb->xscale_x0 = 1.0;

    temp.str("");
    temp << "PERC_xscale_l0";
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->xscale_l0 = atof(str.c_str());
    else
      amb->xscale_l0 = 1.0;

    temp.str("");
    temp << "PERC_xscale_alpha";
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->xscale_alpha = atof(str.c_str());
    else
      amb->xscale_alpha = 1.0;

    // cout <<amb->xscale_x0<<", ";
    // cout <<amb->xscale_l0<<", ";
    // cout <<amb->xscale_alpha<<"\n";
  }
  else {
    amb->xscale = false;
  }

  cout << "\tGot Ambient Params. ****\n";
  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::get_alternate_ambient_params(
    class ReadParams *rparams, struct ambient_data *amb)
{
  int err = 0;
  if (!amb || !rparams)
    rep.error("Null pointer  passed to get_ambient_params!", amb);

  string seek, str;
  ostringstream temp;
  string v;

  v = "RO";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str == "") {
    cout << "\t**** No Alternate Ambient Params...\n";
    return 0;
  }

  cout << "\t**** Getting Alternate Ambient Params...";
  if (!amb->ambient) {
    amb->ambient = mem.myalloc(amb->ambient, SimPM->nvar);
    amb->used    = true;
    for (int v = 0; v < SimPM->nvar; v++)
      amb->ambient[v] = 0.0;
  }

  v = "RO";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[RO] = atof(str.c_str());
  else
    amb->ambient[RO] = -1.0e99;
  ambdens = amb->ambient[RO];

  v = "PG";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[PG] = atof(str.c_str());
  else
    amb->ambient[PG] = -1.0e99;

  v = "VX";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VX] = atof(str.c_str());
  else
    amb->ambient[VX] = -1.0e99;

  v = "VY";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VY] = atof(str.c_str());
  else
    amb->ambient[VY] = -1.0e99;

  v = "VZ";
  temp.str("");
  temp << "PERC_ALTamb" << v;
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    amb->ambient[VZ] = atof(str.c_str());
  else
    amb->ambient[VZ] = -1.0e99;

  if (eqns == 2) {  // mhd sim
    v = "BX";
    temp.str("");
    temp << "PERC_ALTamb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BX] = atof(str.c_str());
    else
      amb->ambient[BX] = -1.0e99;

    v = "BY";
    temp.str("");
    temp << "PERC_ALTamb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BY] = atof(str.c_str());
    else
      amb->ambient[BY] = -1.0e99;

    v = "BZ";
    temp.str("");
    temp << "PERC_ALTamb" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[BZ] = atof(str.c_str());
    else
      amb->ambient[BZ] = -1.0e99;

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    amb->ambient[BX] /= sqrt(4.0 * M_PI);
    amb->ambient[BY] /= sqrt(4.0 * M_PI);
    amb->ambient[BZ] /= sqrt(4.0 * M_PI);
#endif

    if (SimPM->eqntype == EQGLM) amb->ambient[SI] = 0.;
  }  // if mhd vars

  // tracer variables
  for (int t = 0; t < SimPM->ntracer; t++) {
    temp.str("");
    temp << "PERC_ALTambTR" << t;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    if (str != "")
      amb->ambient[t + SimPM->ftr] = atof(str.c_str());
    else
      amb->ambient[t + SimPM->ftr] = -1.0e99;
  }

  //
  // Divider between the first and (optional) second ambient medium.
  //
  temp.str("");
  temp << "PERC_ALTdivider";
  seek = temp.str();
  str  = rparams->find_parameter(seek);
  if (str != "")
    ambdivider = atof(str.c_str());
  else
    ambdivider = SimPM->Xmin[XX] + 0.1 * SimPM->Range[XX];
  cout << "\t\tdivision between two media at " << ambdivider << "\n";

  cout << "\tGot Alternate Ambient Params. ****\n";
  rep.printVec("\tAlt.Amb.", amb->ambient, SimPM->nvar);
  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::add_ambient_data_to_grid(
    class GridBaseClass *ggg, struct ambient_data *amb)
{
  int err = 0;
  if (!ggg || !amb)
    rep.error("Null pointer passed to add_ambient_data_to_grid()", ggg);

  class cell *c = ggg->FirstPt();
  double cloudcentre[MAX_DIM];
  for (int v = 0; v < MAX_DIM; v++)
    cloudcentre[v] = 0.0;
  double dist = 0.0;
  double dpos[ggg->Ndim()];

  do {
    // Set values of primitive variables.
    for (int v = 0; v < SimPM->nvar; v++)
      c->P[v] = amb->ambient[v];
    //
    // Now if we have a radial profile in the slope, we need to adjust rho
    // Note the radial profile is measured from cloudcentre, not
    // neccessarily centred on the source (srcpos).  cloudcentre is set by
    // (PEC_xpos,PEC_ypos,PEC_zpos).
    //
    if (amb->radial_profile != 0) {
      dist = gg->distance_vertex2cell(cloudcentre, c);
      // cout <<"dist="<<dist<<", rad="<<amb->cloudradius<<",
      // rho="<<c->P[RO]; cout <<",
      // multiplier="<<exp(amb->radial_profile*log(amb->cloudradius/dist))<<"\n";
      //
      // Following the Iliev et al 2009 test 6, we use rho=rho0(r0/r)^n if
      // r>r0 We also change the pressure so there is a constant
      // temperature state.
      //
      if (dist > amb->cloudradius) {
        c->P[RO] *= exp(amb->radial_profile * log(amb->cloudradius / dist));
        c->P[PG] *= exp(amb->radial_profile * log(amb->cloudradius / dist));
      }
    }

    if (amb->xscale == true) {
      // cout <<"*#*#*#*#*#*# adding scaled X-data: ";
      // cout <<amb->xscale_x0<<", "<<amb->xscale_l0<<",
      // "<<amb->xscale_alpha<<": ";
      CI.get_dpos(c, dpos);
      // scale factor:
      dist = exp(-pow(
          fabs(dpos[XX] - amb->xscale_x0) / amb->xscale_l0, amb->xscale_alpha));
      // cout <<"     x="<<dpos[XX]<<", and scale factor="<<dist<<"\n";
      c->P[RO] = amb->ambient[RO] * dist;
      c->P[PG] = amb->ambient[PG] * dist;
    }

  } while ((c = ggg->NextPt(c)) != 0);

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::add_alternate_ambient_data_to_grid(
    class GridBaseClass *ggg, struct ambient_data *amb)
{
  int err = 0;
  if (!ggg || !amb)
    rep.error("Null pointer passed to add_ambient_data_to_grid()", ggg);

  //
  // The default is to set only the most negative tenth of the
  // x-domain to the alternate ambient values, but it can be set by the
  // parameter PERC_ALTdivider
  //
  double cutoff = ambdivider;

  class cell *c = ggg->FirstPt();
  //  int xmin=;
  do {
    // Set values of primitive variables.
    if (CI.get_dpos(c, XX) < cutoff) {
      // cout <<"cutoff="<<cutoff<<"\tix="<<c->pos[XX]<<"\n";
      // rep.printVec("P old",c->P,SimPM->nvar);
      for (int v = 0; v < SimPM->nvar; v++)
        c->P[v] = amb->ambient[v];
      // rep.printVec("P new",c->P,SimPM->nvar);
      // if (amb->radial_profile!=0) c->P[RO] *=
      // exp(amb->radial_profile*log((c->pos[XX]-xmin)
    }
  } while ((c = ggg->NextPt(c)) != 0);

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::get_random_clump_params(
    class ReadParams *rparams, struct random_clump_data *rcd)
{
  int err = 0;
  string seek, str;
  ostringstream temp;

  //
  // First see if we are doing random clumps:
  //
  seek = "PERC_addrandomclumps";
  str  = rparams->find_parameter(seek);
  if (str == "" || str == "NO" || str == "no" || str == "No") {
    cout << seek << "=" << str << ", so not doing random clumps...\n";
    return 0;
  }
  else
    cout << seek << "=" << str << ", so continuing to add random clumps.\n";
  rcd->used = true;

  //
  // Initialise border pointer, and get values for the non-clump regions.
  // Border values are as a fraction of the total domain in that dimension.
  //
  rcd->border          = mem.myalloc(rcd->border, 2 * MAX_DIM);
  double default_empty = 0.05;
  seek                 = "PERCempty_xn";
  str                  = rparams->find_parameter(seek);
  if (str == "")
    rcd->border[XN] = default_empty;
  else
    rcd->border[XN] = atof(str.c_str()) * SimPM->Range[XX];
  seek = "PERCempty_xp";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rcd->border[XP] = default_empty;
  else
    rcd->border[XP] = atof(str.c_str()) * SimPM->Range[XX];
  if (ndim > 1) {
    seek = "PERCempty_yn";
    str  = rparams->find_parameter(seek);
    if (str == "")
      rcd->border[YN] = default_empty;
    else
      rcd->border[YN] = atof(str.c_str()) * SimPM->Range[YY];
    seek = "PERCempty_yp";
    str  = rparams->find_parameter(seek);
    if (str == "")
      rcd->border[YP] = default_empty;
    else
      rcd->border[YP] = atof(str.c_str()) * SimPM->Range[YY];
  }
  if (ndim > 2) {
    seek = "PERCempty_zn";
    str  = rparams->find_parameter(seek);
    if (str == "")
      rcd->border[ZN] = default_empty;
    else
      rcd->border[ZN] = atof(str.c_str()) * SimPM->Range[ZZ];
    seek = "PERCempty_zp";
    str  = rparams->find_parameter(seek);
    if (str == "")
      rcd->border[ZP] = default_empty;
    else
      rcd->border[ZP] = atof(str.c_str()) * SimPM->Range[ZZ];
  }

  //
  // Get mean density of matter to go into clumps, and convert this
  // to a total mass with the box volume.
  //
  seek = "PERC_mean_clump_density";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("Failed to get mass to put into clumps", seek);
  else
    rcd->density = atof(str.c_str());

  //
  // Can set the volume to be the total sim.volume, or just the fraction
  // of it we are putting clumps into.  I think the 'R.C. mean density'
  // means more if it is the total volume.
  //
  double volume = 1.0;
  volume *= SimPM->Range[XX];  //-rcd->border[XN]-rcd->border[XP];
  // cout <<"\tvolume="<<volume;
  if (ndim > 1) volume *= SimPM->Range[YY];  //-rcd->border[YN]-rcd->border[YP];
  // cout <<"\tvolume="<<volume;
  if (ndim > 2) volume *= SimPM->Range[ZZ];  //-rcd->border[ZN]-rcd->border[ZP];
  // cout <<"\tvolume="<<volume;
  rcd->total_mass = rcd->density * volume;
  cout << "Mean number density for Random Clumps:" << rcd->density
       << ", giving total mass=" << rcd->total_mass << " grams.\n";

  //
  // Get Random Seed from File
  //
  seek = "PERCrandomseed";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("didn't find parameter", seek);
  else
    rcd->random_seed = atoi(str.c_str());
  srand(rcd->random_seed);

  //
  // profile of clumps, and min/max size
  //
  seek = "PERCclump_profile";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("didn't find parameter", seek);
  else
    rcd->profile = atoi(str.c_str());

  seek = "PERCmin_size";
  str  = rparams->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  rcd->min_size = atof(str.c_str());

  seek = "PERCmax_size";
  str  = rparams->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  rcd->max_size = atof(str.c_str());

  //
  // Find out if "Fixed Mass Range Clumps" or "Fixed Number of Clumps".
  // Get Nclumps, and clump parameters, based on different criteria, so call
  // different functions.
  //
  // N.B. These functions will bug out if less than one random clump is
  // generated, so safe for following code to assume at least one random clump
  // exists.
  //
  seek = "PERC_selection_criterion";
  str  = rparams->find_parameter(seek);
  if (str == "") {
    cout << "WARNING! Parameter PERC_selection_criterion not found, assuming "
            "fixed mass range clumps.\n";
    err += rc_fixed_mass_range_params(rparams, rcd);
  }
  else if (str == "Fixed_Mass_Range") {
    err += rc_fixed_mass_range_params(rparams, rcd);
  }
  else if (str == "Fixed_Number") {
    err += rc_fixed_number_params(rparams, rcd);
  }
  else
    rep.error("Bad parameter for PERC_selection_criterion", str);

  //
  // Tracer Values for Random Clumps:
  //
  for (int v = 0; v < SimPM->ntracer; v++) {
    temp.str("");
    temp << "PERCcloudTR" << v;
    seek = temp.str();
    str  = rparams->find_parameter(seek);
    for (int c = 0; c < rcd->Nclumps; c++) {
      if (str == "")
        rcd->cl[c].tracer_vals[v] = 0.0;
      else
        rcd->cl[c].tracer_vals[v] = atof(str.c_str());
    }
    cout << "***********tracer[" << v << "] = " << rcd->cl[0].tracer_vals[v]
         << endl;
  }

  //
  // Now we have Nclumps, and the clump properties, so select their
  // properties:
  //
  err += rc_set_clump_properties(rcd);

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::rc_fixed_number_params(
    class ReadParams *rparams, struct random_clump_data *rcd)
{
  int err = 0;
  string seek, str;
  ostringstream temp;

  //
  // Get Number of clumps, and initialise *cl with that number.
  //
  seek = "PERCnumclumps";
  str  = rparams->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  rcd->Nclumps = atoi(str.c_str());

  if (rcd->Nclumps < 1) {
    cout << "\tFound " << rcd->Nclumps
         << " Random Clumps, so not adding any clumps.\n";
    rep.error(
        "Need at least one clump for Fixed number of clumps.", rcd->Nclumps);
  }
  else {
    rcd->cl = mem.myalloc(rcd->cl, rcd->Nclumps);
  }

  //
  // choose N+1 random values on the interval [0,M_cl], and use them to get
  // clump masses. Note for large N this will give roughly equal mass clumps
  // because random values have uniform probability (this is why I moved away
  // from this method...).
  //
  std::vector<double> masses;
  masses.push_back(0.0);
  for (int j = 1; j < rcd->Nclumps; j++) {
    masses.push_back(random_frac() * rcd->total_mass);
  }
  masses.push_back(rcd->total_mass);
  std::sort(masses.begin(), masses.end());
  for (int j = 0; j < rcd->Nclumps + 1; j++)
    cout << "masses: " << masses[j] << endl;
  for (int j = 0; j < rcd->Nclumps; j++) {
    rcd->cl[j].mass = masses[j + 1] - masses[j];
  }

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::rc_fixed_mass_range_params(
    class ReadParams *rparams, struct random_clump_data *rcd)
{
  int err = 0;

  //
  // Get min and max clump masses, as fraction of total mass.
  //
  string seek, str;
  seek = "PERCmin_mass";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("Failed to get param", seek);
  else
    rcd->min_mass = atof(str.c_str()) * rcd->total_mass;
  seek = "PERCmax_mass";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("Failed to get param", seek);
  else
    rcd->max_mass = atof(str.c_str()) * rcd->total_mass;
  cout << "min/max masses: " << rcd->min_mass << ", " << rcd->max_mass
       << ",  and total mass: " << rcd->total_mass << endl;

  //
  // Take a random chunk (between min and max) out of total mass, as store as
  // a clump mass
  //
  std::vector<double> masses;
  double mass_remaining = rcd->total_mass, this_mass = 0.0;
  int nnn = 0;
  do {
    this_mass = rcd->min_mass + random_frac() * (rcd->max_mass - rcd->min_mass);
    this_mass = std::min(this_mass, mass_remaining);
    masses.push_back(this_mass);
    mass_remaining -= this_mass;
    nnn++;
  } while (mass_remaining > (rcd->total_mass * SMALLVALUE));
  cout << "Used up all the mass in " << nnn << " clumps (Nclumps was "
       << rcd->Nclumps << ")  Mass left over=" << mass_remaining << "\n";
  rcd->Nclumps = nnn;

  if (rcd->Nclumps < 1)
    rep.error(
        "Got less than one clump in fixed total mass routine!", rcd->Nclumps);

  //
  // set up clumps struct and add masses for each clump.
  //
  rcd->cl = mem.myalloc(rcd->cl, rcd->Nclumps);
  for (int j = 0; j < rcd->Nclumps; j++) {
    rcd->cl[j].mass = masses[j];
  }

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::rc_set_clump_properties(
    struct random_clump_data *rcd)
{
  int err     = 0;
  double xmax = SimPM->Range[XX], ymax = SimPM->Range[YY],
         zmax = SimPM->Range[ZZ];

  for (int j = 0; j < rcd->Nclumps; j++) {
    //
    // Set clump centres, allowing for empty regions.
    // If periodic BCs, we only have empty regions in the X-dir, since Y,Z
    // are periodic
    //
    if (SimPM->BC_YN == "periodic") {
      // cout <<"find = "<<SimPM->typeofbc.find("YNper")<<endl;
      rcd->cl[j].centre[XX] =
          SimPM->Xmin[XX] + rcd->border[XN]
          + (xmax - rcd->border[XP] - rcd->border[XN]) * random_frac();
      rcd->cl[j].centre[YY] = SimPM->Xmin[YY] + ymax * random_frac();
      rcd->cl[j].centre[ZZ] = SimPM->Xmin[ZZ] + zmax * random_frac();
    }
    else {
      // cout <<"not periodic bcs!\n";
      rcd->cl[j].centre[XX] =
          SimPM->Xmin[XX] + rcd->border[XN]
          + (xmax - rcd->border[XP] - rcd->border[XN]) * random_frac();
      rcd->cl[j].centre[YY] =
          SimPM->Xmin[YY] + rcd->border[YN]
          + (ymax - rcd->border[YP] - rcd->border[YN]) * random_frac();
      rcd->cl[j].centre[ZZ] =
          SimPM->Xmin[ZZ] + rcd->border[ZN]
          + (zmax - rcd->border[ZP] - rcd->border[ZN]) * random_frac();
    }

    //
    // radius of clumps
    //
    rcd->cl[j].size[XX] =
        rcd->min_size * ymax
        + (rcd->max_size - rcd->min_size) * ymax * random_frac();
    rcd->cl[j].size[YY] =
        rcd->min_size * ymax
        + (rcd->max_size - rcd->min_size) * ymax * random_frac();
    rcd->cl[j].size[ZZ] =
        rcd->min_size * ymax
        + (rcd->max_size - rcd->min_size) * ymax * random_frac();
    // orientation from grid axes to principal axes.
    rcd->cl[j].ang[XX] = random_frac() * M_PI;
    rcd->cl[j].ang[YY] = random_frac() * M_PI;
    rcd->cl[j].ang[ZZ] = random_frac() * M_PI;
    //
    // zero the z-dir if 2d
    if (SimPM->ndim == 2) {
      rcd->cl[j].ang[ZZ] = rcd->cl[j].ang[YY] = rcd->cl[j].size[ZZ] =
          rcd->cl[j].centre[ZZ]               = 0.0;
    }

    //
    // set overdensity (FIX THIS SO THAT IT CAN DO OTHER THAN GAUSSIAN
    // PROFILES!!!!)
    //
    rcd->cl[j].overdensity = rcd->cl[j].mass / ambdens;
    for (int v = 0; v < SimPM->ndim; v++)
      rcd->cl[j].overdensity /= sqrt(2.0 * M_PI) * rcd->cl[j].size[v];

    //
    // print clump...
    print_clump(&rcd->cl[j]);
    //
    // now set up the rotation matrix.

    rcd->cl[j].rm[XX][XX] = cos(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[YY])
                                * cos(rcd->cl[j].ang[ZZ])
                            - sin(rcd->cl[j].ang[XX]) * sin(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[XX][YY] = sin(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[YY])
                                * cos(rcd->cl[j].ang[ZZ])
                            + cos(rcd->cl[j].ang[XX]) * sin(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[XX][ZZ] = -sin(rcd->cl[j].ang[YY]) * cos(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[YY][XX] = -cos(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[YY])
                                * sin(rcd->cl[j].ang[ZZ])
                            - sin(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[YY][YY] = -sin(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[YY])
                                * sin(rcd->cl[j].ang[ZZ])
                            + cos(rcd->cl[j].ang[XX]) * cos(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[YY][ZZ] = sin(rcd->cl[j].ang[YY]) * sin(rcd->cl[j].ang[ZZ]);
    rcd->cl[j].rm[ZZ][XX] = cos(rcd->cl[j].ang[XX]) * sin(rcd->cl[j].ang[YY]);
    rcd->cl[j].rm[ZZ][YY] = sin(rcd->cl[j].ang[XX]) * sin(rcd->cl[j].ang[YY]);
    rcd->cl[j].rm[ZZ][ZZ] = cos(rcd->cl[j].ang[YY]);

    // set velocity
    rcd->cl[j].Vel[0] = amb_data.ambient[VX];
    rcd->cl[j].Vel[1] = amb_data.ambient[VY];
    rcd->cl[j].Vel[2] = amb_data.ambient[VZ];

    // all done, move on to next clump.
  }
  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::add_random_clumps_to_grid(
    class GridBaseClass *ggg, struct random_clump_data *rcd)
{
  int err = 0;
  if (!rcd->used) {
    cout << "\tnot using random clumps, so not adding any to grid.\n";
    return 0;
  }
  cout << "Adding random clumps to grid...";

  cell *c = ggg->FirstPt();
  do {
    err += clumps_set_dens(c, rcd->Nclumps, rcd->cl, rcd->profile);
  } while ((c = ggg->NextPt(c)) != 0);

  //
  // Set tracer values, assuming all random clumps have the same tracer value,
  // and that it is based on density compared to ambient density.  Need to
  // change it in future if I want to change this.
  //
  // cout <<"ntracer = "<<SimPM->ntracer<<"\t
  // amb[tr1]="<<ambient[SimPM->ftr+1]<<"\tcl[tr1]="<<cltr[1];
  c = ggg->FirstPt();
  cout << "\tamb-density=" << ambdens << " clump_mass=" << rcd->total_mass
       << endl;
  do {
    for (int v = SimPM->ftr; v < SimPM->nvar; v++) {
      // cout <<"trval before="<<c->P[v];
      c->P[v] =
          (ambdens / c->P[RO]) * amb_data.ambient[v]
          + (1.0 - ambdens / c->P[RO]) * rcd->cl[0].tracer_vals[v - SimPM->ftr];
      // cout <<"\tafter="<<c->P[v]<<endl;
    }
  } while ((c = ggg->NextPt(c)) != 0);

  cout << "\tFinished adding random clumps to grid.\n";

  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::get_strategic_clump_params(
    class ReadParams *rparams, struct strategic_clump_data *scd)
{
  int err = 0;
  string seek, str;
  ostringstream temp;

  //
  // First see if we have strategic clumps to place:
  //
  seek = "PE_SC_addstrategicclumps";
  str  = rparams->find_parameter(seek);
  if (str == "" || str == "NO" || str == "no" || str == "No") {
    cout << seek << "=" << str << ", so not doing any strategic clumps...\n";
    return 0;
  }
  else
    cout << seek << "=" << str << ", so continuing to add strategic clumps.\n";
  scd->used = true;

  //
  // See how many clumps we need:
  //
  seek = "PE_SC_numclumps";
  str  = rparams->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  scd->Nclumps = atoi(str.c_str());

  if (scd->Nclumps < 1) {
    cout << "\tFound " << scd->Nclumps
         << " Strategic Clumps, so not adding any clumps.\n";
    rep.error("Need at least one clump for strategic clumps.", scd->Nclumps);
  }
  else {
    scd->cl = mem.myalloc(scd->cl, scd->Nclumps);
  }

  //
  // profile of clumps, and min/max size
  //
  seek = "PE_SC_clump_profile";
  str  = rparams->find_parameter(seek);
  if (str == "")
    rep.error("didn't find parameter", seek);
  else
    scd->profile = atoi(str.c_str());

  //
  // Now loop over all clumps, reading in data for each.
  //
  for (int c = 0; c < scd->Nclumps; c++) {
    //
    // Clump Overdensity.
    //
    temp.str("");
    temp << "PE_SC" << c << "_overdensity";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "")
      rep.error("no parameter found for clump", seek);
    else
      scd->cl[c].overdensity = atof(str.c_str());

    //
    // Clump Position (read in as fraction of Y-Range).
    //
    temp.str("");
    temp << "PE_SC" << c << "_xpos";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("didn't find parameter", seek);
    scd->cl[c].centre[XX] =
        SimPM->Xmin[XX] + atof(str.c_str()) * SimPM->Range[YY];
    if (ndim > 1) {
      temp.str("");
      temp << "PE_SC" << c << "_ypos";
      seek = temp.str();
      str  = rp->find_parameter(seek);
      if (str == "") rep.error("didn't find parameter", seek);
      scd->cl[c].centre[YY] =
          SimPM->Xmin[YY] + atof(str.c_str()) * SimPM->Range[YY];
    }
    if (ndim > 2) {
      temp.str("");
      temp << "PE_SC" << c << "_zpos";
      seek = temp.str();
      str  = rp->find_parameter(seek);
      if (str == "") rep.error("didn't find parameter", seek);
      scd->cl[c].centre[ZZ] =
          SimPM->Xmin[ZZ] + atof(str.c_str()) * SimPM->Range[YY];
    }

    //
    // Clump size: (note this is hard-coded to be circular, this param is as
    // fraction of the Y-Range).
    //
    temp.str("");
    temp << "PE_SC" << c << "_radius";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("didn't find parameter", seek);
    for (int v = 0; v < MAX_DIM; v++)
      scd->cl[c].size[v] = atof(str.c_str()) * SimPM->Range[YY];

    //
    // Clump Velocity (cm/s).
    //
    temp.str("");
    temp << "PE_SC" << c << "_vx";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "")
      scd->cl[c].Vel[0] = 0.0;
    else {
      scd->cl[c].Vel[0] = atof(str.c_str());
      cout << "VX=" << scd->cl[c].Vel[0] << "\n";
    }

    temp.str("");
    temp << "PE_SC" << c << "_vy";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "")
      scd->cl[c].Vel[1] = 0.0;
    else
      scd->cl[c].Vel[1] = atof(str.c_str());

    temp.str("");
    temp << "PE_SC" << c << "_vz";
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str == "")
      scd->cl[c].Vel[2] = 0.0;
    else
      scd->cl[c].Vel[2] = atof(str.c_str());

    //
    // Tracer Values for each cloud (can be different for each!)
    //
    for (int v = 0; v < SimPM->ntracer; v++) {
      temp.str("");
      temp << "PE_SC" << c << "_cloudTR" << v;
      seek = temp.str();
      str  = rp->find_parameter(seek);
      if (str == "")
        scd->cl[c].tracer_vals[v] = 0.0;
      else
        scd->cl[c].tracer_vals[v] = atof(str.c_str());
    }

    // orientation from grid axes to principal axes.
    scd->cl[c].ang[XX] = 0.0;
    scd->cl[c].ang[YY] = 0.0;
    scd->cl[c].ang[ZZ] = 0.0;
    //
    // zero the z-dir if 2d
    if (SimPM->ndim == 2) {
      scd->cl[c].ang[ZZ] = scd->cl[c].ang[YY] = scd->cl[c].size[ZZ] =
          scd->cl[c].centre[ZZ]               = 0.0;
    }

    //
    // Set Mass, based on profile.  Note this is total mass, including
    // ambient gas.
    //
    scd->cl[c].mass = ambdens * (1.0 + scd->cl[c].overdensity);
    if (scd->profile == 0) {
      // top hat profile: M=Volume*density
      scd->cl[c].mass *= 4.0 * M_PI * scd->cl[c].size[XX] * scd->cl[c].size[XX]
                         * scd->cl[c].size[XX] / 3.0;
    }
    else if (scd->profile == 1) {
      // Gaussian profile
      for (int v = 0; v < SimPM->ndim; v++)
        scd->cl[c].mass *= sqrt(2.0 * M_PI) * scd->cl[c].size[v];
      //
      // If axisymmetric, then 3rd dimension is same as radial dimension,
      // so multiply again.
      //
      if (SimPM->coord_sys == COORD_CYL)
        scd->cl[c].mass *= sqrt(2.0 * M_PI) * scd->cl[c].size[Rcyl];
    }
    else
      rep.error("Bad profile in get_strategic_clump_params()", scd->profile);

    //
    // print clump...
    print_clump(&scd->cl[c]);
    //
    // now set up the rotation matrix.
    //
    scd->cl[c].rm[XX][XX] = cos(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[YY])
                                * cos(scd->cl[c].ang[ZZ])
                            - sin(scd->cl[c].ang[XX]) * sin(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[XX][YY] = sin(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[YY])
                                * cos(scd->cl[c].ang[ZZ])
                            + cos(scd->cl[c].ang[XX]) * sin(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[XX][ZZ] = -sin(scd->cl[c].ang[YY]) * cos(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[YY][XX] = -cos(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[YY])
                                * sin(scd->cl[c].ang[ZZ])
                            - sin(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[YY][YY] = -sin(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[YY])
                                * sin(scd->cl[c].ang[ZZ])
                            + cos(scd->cl[c].ang[XX]) * cos(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[YY][ZZ] = sin(scd->cl[c].ang[YY]) * sin(scd->cl[c].ang[ZZ]);
    scd->cl[c].rm[ZZ][XX] = cos(scd->cl[c].ang[XX]) * sin(scd->cl[c].ang[YY]);
    scd->cl[c].rm[ZZ][YY] = sin(scd->cl[c].ang[XX]) * sin(scd->cl[c].ang[YY]);
    scd->cl[c].rm[ZZ][ZZ] = cos(scd->cl[c].ang[YY]);
    //
    // all done, move on to next clump.
    //
  }
  cout << "\tSet parameters for " << scd->Nclumps << " clumps. returning.\n";
  return err;
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::add_strategic_clumps_to_grid(
    class GridBaseClass *ggg, struct strategic_clump_data *scd)
{
  int err = 0;
  if (!scd->used) {
    cout << "\tnot using strategic clumps, so not adding any to grid.\n";
    return 0;
  }
  cout << "Adding strategic clumps to grid...";

  cell *c = ggg->FirstPt();
  do {
    err += clumps_set_dens(c, scd->Nclumps, scd->cl, scd->profile);

  } while ((c = ggg->NextPt(c)) != 0);

  cout << "\tFinished adding strategic clumps to grid.\n";

  return err;
}

// ##################################################################
// ##################################################################

double IC_photevap_multi_clumps::random_frac()
{
  //
  // Returns a random value between 0 and 1.
  //
  return ((double)rand()) / ((double)RAND_MAX);
}

// ##################################################################
// ##################################################################

int IC_photevap_multi_clumps::clumps_set_dens(
    class cell *c, const int Nclumps, struct clump *cl, const int profile)
{
  int err = 0;
  double x0[ndim], x1[ndim], dpos[ndim];
  CI.get_dpos(c, dpos);
  for (int j = 0; j < Nclumps; j++) {
    //
    // get distance from centre to clump.
    // If we have periodic BCs, need to do a bit more work...
    //
    // if (SimPM->typeofbc.find("YNper") != std::string::npos) {
    if (SimPM->BC_YN == "periodic") {
      x0[XX] = dpos[XX] - cl[j].centre[XX];
      double temp;
      for (int i = 1; i < ndim; i++) {
        temp = dpos[i] - cl[j].centre[i];
        if (temp > 0.0) {
          if (temp > SimPM->Range[i] / 2.)
            x0[i] = temp - SimPM->Range[i];
          else
            x0[i] = temp;
        }
        else {
          if (temp < -SimPM->Range[i] / 2.)
            x0[i] = temp + SimPM->Range[i];
          else
            x0[i] = temp;
        }
      }
    }
    else {
      // not periodic bcs.
      for (int i = 0; i < ndim; i++)
        x0[i] = dpos[i] - cl[j].centre[i];
    }

    //
    // rotate clump according to rotation matrix generated from random
    // angles.
    //
    for (int u = 0; u < ndim; u++) {
      x1[u] = 0.0;
      for (int v = 0; v < ndim; v++)
        x1[u] += x0[v] * cl[j].rm[u][v];
    }

    //
    // add density from clump onto existing density in cell.
    //
    double add_rho = 0.0;
    if (profile == 0) {
      // top-hat profile
      bool inside = true;
      double dist = 0.0;
      for (int k = 0; k < ndim; k++)
        dist += x1[k] * x1[k] / cl[j].size[k] / cl[j].size[k];
      if (dist > 1.0) inside = false;
      //      for (int k=0;k<ndim;k++) if (fabs(x1[k])> cl[j].size[k])
      //      inside=false;
      if (inside) {
        add_rho = ambdens * cl[j].overdensity;
        c->P[RO] += add_rho;
      }
    }
    else if (profile == 1) {
      //
      // Gaussian Profile: Set exponential factor based on how far from
      // centre we are. Then set density based on that. Note profile has
      // factor of 2, so that size is easily related to normal
      // distribution. Set tracers according to fraction of mass that is
      // ambient and clump:
      //
      double ef = 0.0;
      for (int v = 0; v < ndim; v++)
        ef += x1[v] * x1[v] / cl[j].size[v] / cl[j].size[v] / 2.0;
      add_rho = ambdens * cl[j].overdensity * exp(-ef);
      c->P[RO] += add_rho;
    }
    else
      rep.error("Bad profile id in parameter-file", profile);

    //
    // Set tracer values:
    // (random clumps over-writes this later, because it is sketchy if
    // clumps overlap) Criterion: If we have added density>0.5*ambient to
    // cell, label it with
    //            clump's tracer values.
    // if (add_rho >= 0.1*ambdens*cl[j].overdensity) {
    //  for (int v=SimPM->ftr; v<SimPM->nvar;v++) {
    //	c->P[v] = cl[j].tracer_vals[v-SimPM->ftr];
    // }
    //}

    //
    // Update 2010-02-08 JM: compare density to ambient and apply
    // tracer value appropriately.  Note this only works if all the
    // clumps have the same tracer value (and also for ambient data!).
    //
    if (add_rho / c->P[RO] > 0.001) {
      for (int v = SimPM->ftr; v < SimPM->nvar; v++) {
        c->P[v] =
            cl[j].tracer_vals[v - SimPM->ftr]
            + ambdens / c->P[RO]
                  * (amb_data.ambient[v] - cl[j].tracer_vals[v - SimPM->ftr]);
      }
      //
      // Update 2012.02.25 JM: Add clump velocity w.r.t. grid, smoothly
      // varying from the max value at the centre to the ambient value
      // near the edge.
      //
      // cout
      // <<"V=["<<cl[j].Vel[0]<<","<<cl[j].Vel[1]<<","<<cl[j].Vel[2]<<"\n";
      c->P[VX] = cl[j].Vel[0]
                 + ambdens / c->P[RO] * (amb_data.ambient[VX] - cl[j].Vel[0]);
      c->P[VY] = cl[j].Vel[1]
                 + ambdens / c->P[RO] * (amb_data.ambient[VY] - cl[j].Vel[1]);
      c->P[VZ] = cl[j].Vel[2]
                 + ambdens / c->P[RO] * (amb_data.ambient[VZ] - cl[j].Vel[2]);
    }

    //
    // Done with this clump, loop onto next...
    //
  }

  //
  // If we are using isothermal equations, then ambient[PG] is set
  // as the sound speed everywhere (assuming everywhere is neutral
  // initially), so we don't need to adjust it at all.
  //

  return err;
}


// ##################################################################
// ##################################################################



void IC_photevap_multi_clumps::print_clump(struct clump *rc)
{
  cout << "--clump overdensity:" << rc->overdensity
       << "  mass/Msun:" << rc->mass / pconst.Msun() << endl;
  rep.printVec("Centre", rc->centre, SimPM->ndim);
  rep.printVec("Radius", rc->size, SimPM->ndim);
  rep.printVec("Angles", rc->ang, SimPM->ndim);
  cout << "-------------------\n";
  return;
}
