/// \file raytracer_SC.cc
/// \brief Contains Function Definitions of Short Characteristic Tracer Class.
/// 
/// Author: Jonathan Mackey
/// 
/// Modifications:
///  - 2008-03-19 Wrote file.
///  - 2010-01-19 JM: put weighting schemes in new function to avoid code duplication in parallel version.
///  - 2010-01-20 JM: changed ds/2 to 0.5*ds on line 1334, and similar in process_source_cell(). Tidied up code.
///  - 2010-01-22 JM: worked on sources at cell-corners and got it working in 2D.
///  - 2010-01-23 JM: changed for corner-centred sources: find_src_on_grid(), centre_source_on_cell()
///  - 2010-01-24 JM: changed for corner-centred sources: cell_cols_3D()
///  - 2010-01-26 JM: Debugging sources at cell-vertices in 3D.
///  - 2010-02-10 JM: Removed obselete treatment of RN boundary for axisymmetry.
///  - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.02.17 JM: Moved to raytracing/ subdir, so change include for uniform_grid.
///     Allowed parallel rays class to have a source at infinity in the +R direction.
/// - 2011.02.25 JM: Tidied up parallel rays class a little.  Allowed col2cell_1d()
///    to read data from boundary cells if appropriate; so I now *MUST* make sure 
///    their column data are initialised to zero!
///    Removed NEW_RT_MP_INTERFACE flags. Should handle multiple sources now. Of course I 
///    need a "rates-only" update function to actually use multiple emitting sources.
///    Removed HCORR ifdef around new code.
/// - 2011.02.28 JM: Undid an error I introduced last week.
///    Put more comments in RT_TESTING.
/// - 2011.03.02 JM: parallel rays, col1d now only treats src-cell differently if
///    the source is not at infinity.
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function.  It is 
///    just a dummy function for now here.
/// - 2011.04.15 JM: Updated the Raytracing code so that it can do just Column-density 
///    calculations, depending on the source->update and source->opacity_src settings.
///    This is decided in the ProcessCell() function, so RayTrace_Column_Density() is
///    redundant, but it does make the code clearer in gridMethods, so I might leave it.
/// - 2011.04.17 JM: Debugged the Vshell setting.  Seems to be working now.
/// - 2011.04.18 JM: Added interpolate_2D_RHO() function.  Fixed other bugs.  Working now!
/// - 2011.04.22 JM: Put the setting-Vshell subroutine in an ifdef for serial only.
///    Parallel version calls it in its own Add_Source() function.
/// - 2011.04.23 JM: Added TauMin[] array of values of TauMin for each source, ordered
///    by source id.  This removes the need for interpolate_2D_RHO() so I got rid of it.
///    Now the interpolate functions read TauMin[src_id] and apply that to the interpolation.
/// - 2011.05.02 JM: Added multifreq ionisation option to TauMin[] chooser (same as mono).
/// - 2011.05.04 JM: Changed TauMin for multifreq to be larger (still trying to get best value).
/// - 2011.06.21 JM: Changed it so that all cell data is read from Ph[] not P[].  This means
///    I can get time-centred optical depths from the half-step state.  The C2Ray update is
///    still for the full timestep only (fully operator-split).
///
/// - 2011.10.14 JM: Decreased the quantity and improved the formatting of
///    information printed out when adding radiation sources.  Added new interface
///    functions for the number of ionising and heating sources, and to get the 
///    vector of data structures.
/// - 2011.10.17 JM: Debugging.  Getting new implicit integrator working.
/// - 2011.10.22 JM: Debugging. New implicit integrator now seems to work as it should.
///    Deleted the NEW_MP_INTERFACE flags, so that the old code is now gone and
///    replaced with the NEW_MP_INTERFACE stuff.
/// - 2011.10.24 JM: Put more comments in RT_TESTING.
/// - 2012.01.16 JM: Gave update_RT_source_properties() a real function.
/// - 2012.03.12 JM: Added a hack TauMin for Harpreet's H2 opacity (unused for
///    1D simulations, so it is not important).
/// - 2012.03.31 JM: Moved setting Vshell for sources into its own function (it
///    was causing trouble with the parallel code, which needs to call it 
///    somewhere else).  Also split the Add_Source() function into 2 functions
///    so that I can get rid of the serial ifdefs.
/// - 2012.04.23 JM: Added option changing TauMin[] for primordial gas tests.
/// - 2012.05.16 JM: Added OPACITY_HALPHA to calculate H-alpha emission along a
///    ray, including absorption by dust (which is assumed to be everywhere).
/// - 2013.01.17 JM: Put some error checking in cell_cols_3D() inside
///    an RT_TESTING flag.  Enforced column density >= 0 with max().
/// - 2013.08.12 JM: Added raytracing of the recombination rate for a
///    photoionisation equilibrium calculation.
/// - 2013.08.20 JM: Modifications so a given source can have an
///    array of optical depths rather than a single value per source.
///    This is tough going, and is only half-way done so far.
/// - 2013.08.23 JM: Debugging new HHe module code.
/// - 2014.06.10 JM: Fixed bug in 1D parallel grid in spherical
///    coordinates, where when run in parallel, the column densities
///    from internal boundary cells were not correctly picked up.
/// - 2014.08.01 JM: Set ipos[] to have the correct sign in the fn
///    raytracer_USC_infinity::add_source_to_list()
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2017.08.24 JM: Set code to ignore internal boundary regions
///    when calculating absorption.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "raytracing/raytracer_SC.h"
#include "constants.h"
#include "grid/uniform_grid.h"
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

//#define RT_TESTING

#define LARGE_TAU_FINISHED -10

/*********************************************************************/
/************************* PARALLEL RAYS *****************************/



// ##################################################################
// ##################################################################



raytracer_USC_infinity::raytracer_USC_infinity(
    class GridBaseClass *ggg,   ///< Pointer to grid
    class microphysics_base *mmm, ///< Pointer to MicroPhysics.
    int nd,         ///< number of dimensions of grid
    int csys,  ///< coordinate system
    int nv,         ///< number of variables in state vector
    int ftr         ///< index of first tracer variable in state vector
    )
: ndim(nd), coord_sys(csys), nvar(nv), first_tr(ftr)
{
#ifdef RT_TESTING
  cout <<"SC src-at-infinity (parallel-rays) raytracer class constructor!\n";
#endif
  if (XN!=0 || XP!=1 || YN!=2 || YP!=3 || ZN!=4 || ZP!=5 || NO!=-1)
    rep.error("direction enums not what is expected, bad stuff will happen.",XN);
  if (XX!=0  || YY!=1 || ZZ!=2)
    rep.error("axes enums not what is expected, bad stuff will happen.",XX);

  raytracer_USC_infinity::gridptr=0; gridptr = ggg;
  raytracer_USC_infinity::mpptr=0;   mpptr   = mmm;
  if (!gridptr || !mpptr) rep.error("grid or microphysics null pointers",gridptr);

#ifdef RT_TESTING
  cout <<"  ndim = "<<ndim<<" and grid pointer = "<<gridptr<<"\n";
#endif
  if (ndim <1 || ndim>3) rep.error("bad ndim for raytracing!",ndim);
  if (SourceList.size() != 0)
    rep.error("SourceList badly initialised",SourceList.size());

  N_ion_srcs = 0;
  N_uvh_srcs = 0;
  type_of_RT_int = RT_UPDATE_EXPLICIT; // default is explicit integration.

  return;
}




// ##################################################################
// ##################################################################



raytracer_USC_infinity::~raytracer_USC_infinity()
{
#ifdef RT_TESTING
  cout <<"SC src at infinity raytracer class destructor!\n";
#endif
  for (vector<rad_source>::iterator i=SourceList.begin(); i!=SourceList.end(); ++i) {
    //    cout <<(*i).pos[0]<<"\n";
    //if ((*i).pos !=0) (*i).pos  = mem.myfree((*i).pos);
    //if ((*i).ipos!=0) (*i).ipos = mem.myfree((*i).ipos);
  }
  SourceList.clear();
  gridptr=0;
  mpptr=0;
  return;
}



// ##################################################################
// ##################################################################




int raytracer_USC_infinity::Add_Source(
      struct rad_src_info *src ///< source info.
      )
{
  //
  // This happens in two stages: first we check all of the source properties
  // and add it to the list of sources using the add_source_to_list() function.
  // Then we set Vshell in every cell for the source.  It has to be in two
  // steps because the parallelised raytracer needs to connect the domains in
  // between these two steps.
  //

#ifdef RT_TESTING
  cout <<"\n--BEGIN-----raytracer_USC_infinity::AddSource()------------\n";
#endif
  add_source_to_list(src);

  //
  // Now we need to set Vshell in every grid point for this source. 
  // 
  set_Vshell_for_source(&SourceList.back());

#ifdef RT_TESTING
  cout <<"Add_Source() all done.\n";
  cout <<"\n--END-----raytracer_USC_infinity::AddSource()------------\n";
#endif
  return (SourceList.back()).s->id;
}



// ##################################################################
// ##################################################################


void raytracer_USC_infinity::add_source_to_list(
      struct rad_src_info *src ///< source info.
      )
{
  rad_source rs;
  
  //
  // Set pointer in rad_source to *src.  This sets the position, id,
  // and all the parameters necessary to set up the source.
  //
  rs.s = src;

  //rs.pos=0; rs.ipos=0;
  //rs.pos  = mem.myalloc(rs.pos, ndim);
  //rs.ipos = mem.myalloc(rs.ipos,ndim);
  //for (int i=0; i<ndim; i++) rs.pos[i] = src->position[i];
  for (int i=0; i<ndim; i++) rs.ipos[i] = -100; // ipos not needed here.

  //rs.id          = src->id;
  //rs.strength    = src->strength;
  //rs.type        = src->type;
  //rs.at_infinity = (src->at_infinity>0) ? true : false;
  //if (!rs.at_infinity) {
  //  rep.error("Source is not at infinity",src->id);
  //}
  //rs.effect      = src->effect;
  //rs.opacity_src = src->opacity_src;
  //rs.opacity_var = src->opacity_var;
  //rs.update      = src->update;

  //
  // Set the source-dependent parts of the rt_source_data struct.
  //
  rs.data.id       = rs.s->id;
  rs.data.type     = rs.s->type;
  rs.data.strength = rs.s->strength;
  rs.data.Vshell   = 0.0;
  rs.data.dS       = 0.0;
  rs.data.NTau     = rs.s->NTau;
  for (unsigned short int iT=0; iT<src->NTau; iT++) {
    rs.data.DelCol[iT]   = 0.0;
    rs.data.Column[iT]   = 0.0;
  }
  
  //
  // Check that opacity_var is not running off the end of the state vector.
  //
  if (rs.s->opacity_var+first_tr >nvar-1) {
    cout <<"opacity_var="<<rs.s->opacity_var<<" and ftr="<<first_tr;
    cout <<", but state-vec has only "<<nvar;
    cout <<" elements.  The opacity";
    cout <<" var will run off then end of the state vector array.\n";
    cout <<"OPACITY_VAR IS OFFSET - 1ST TRACER HAS OPACITY_VAR=0.\n";
    rep.error("Bad opacity var for source",rs.s->id);
  }

  //
  // now find the cell closest to source so we can set up a pointer to it.
  //
#ifdef RT_TESTING
  cout <<"Add_Source() finding source.\n";
#endif
  cell *c=gridptr->FirstPt();
  //
  // first find direction to source at infinity, and error if not at infinity.
  //
  enum direction dir=NO;
  for (int i=0;i<ndim;i++) {
    if (rs.s->pos[i] <= gridptr->Xmin(static_cast<axes>(i)) ||
	rs.s->pos[i] >= gridptr->Xmax(static_cast<axes>(i))) {
      if      (rs.s->pos[i] < -1.e99) {
        dir=static_cast<direction>(2*i);
        rs.ipos[i] = -10000;
      }
      else if (rs.s->pos[i] >  1.e99) {
        dir=static_cast<direction>(2*i+1);
        rs.ipos[i] = 10000;
      }
      else {
#ifdef RT_TESTING
	cout <<"Source off grid, but not at infinity!";
        cout <<" pos="<<rs.s->pos[i];
        cout <<", resetting to centre\n";
#endif
	rs.s->pos[i] = 0.5*(gridptr->Xmin(static_cast<axes>(i))+
                            gridptr->Xmax(static_cast<axes>(i)));
      }
    }
    else {
      double centre = 0.5*(gridptr->Xmin(static_cast<axes>(i))+
                           gridptr->Xmax(static_cast<axes>(i)));
      if (!pconst.equalD(centre,rs.s->pos[i])) {
#ifdef RT_TESTING
	cout <<"source not at infinity, or at centre of grid in dir";
        cout <<" "<<i<<", resetting to centre.\n";
	cout <<"old pos: "<<rs.s->pos[i]<<"  centre: "<<centre<<"\n";
#endif
	rs.s->pos[i] = centre;
      }
    }
  }
  if (dir==NO)
    rep.error("source not at infinity in any direction",dir);

  //
  // Check that coord sys allows source at infinity in the direction
  // found.
  // 2011.02.17 JM: I'm adding the possibility for a source at
  // infinity in the positive radial direction, since I need it for
  // diffuse field calculations.
  //
  if (coord_sys==COORD_CYL) {
    if (dir!=ZNcyl && dir!=ZPcyl && dir!=RPcyl) {
      rep.error("Cylindrical coords, src at R<0!",dir);
    }
  }

  //
  // now go from first point on grid in direction dir until we get 
  // to the edge.  Then go back one cell.  This is rs.sc.
  //
  while (gridptr->NextPt(c,dir)) c=gridptr->NextPt(c,dir);
  c = gridptr->NextPt(c,gridptr->OppDir(dir));

  rs.sc = c;
#ifdef RT_TESTING
  cout <<"Add_Source() source->sc = "<<rs.sc<<", id="<<rs.sc->id<<"\n";
#endif
  rs.src_on_grid = false;

  SourceList.push_back(rs);
  update_local_variables_for_new_source(rs);
  return;
}



// ##################################################################
// ##################################################################


void raytracer_USC_infinity::set_Vshell_for_source(
      struct rad_source *this_src
      )
{
  //
  // Now we need to set Vshell in every grid point for this source. 
  // Change opacity flag to 'Vshell' so ProcessCell() will just set Vshell in each 
  // cell.  Then call RayTrace_Column_Density() where Vshell is set, and then 
  // change the opacity flag back to its original value.
  //
//#ifdef RT_TESTING
  cout <<"\t\tSetting Vshell for source.\n";
//#endif

  int temp = this_src->s->opacity_src;
  int upd  = this_src->s->update;
  this_src->s->opacity_src = RT_OPACITY_VSHELL;
  this_src->s->update      = RT_UPDATE_EXPLICIT;

  int err = RayTrace_Column_Density(this_src->s->id,1.0,1.0);
  if (err)
    rep.error("raytracer_USC_infinity::RayTrace_Column_Density() Vshell",err);

  this_src->s->opacity_src = temp; // revert opacity type
  this_src->s->update      = upd;  // revert update type.
  return;
}


// ##################################################################
// ##################################################################



void raytracer_USC_infinity::update_local_variables_for_new_source(
      const struct rad_source rs_new
      )
{
  //
  // Set local class variables:
  // - Increment either the number of ionising or uv-heating sources.
  // - Check if we are doing an implicit raytracing/integration.
  //
  if (rs_new.s->effect==RT_EFFECT_UV_HEATING) {
    N_uvh_srcs++;
    UVH_data.push_back(rs_new.data);
  }
  else {
    N_ion_srcs++;
    ION_data.push_back(rs_new.data);
  }

  if (rs_new.s->update==RT_UPDATE_IMPLICIT)
    type_of_RT_int=RT_UPDATE_IMPLICIT;

  return;
}



// ##################################################################
// ##################################################################


void raytracer_USC_infinity::update_RT_source_properties(
      const struct rad_src_info *rs ///< ptr to source info.
      )
{
  //
  // rs contains updated strength, Rstar, and Tstar data, so update strength
  // here.
  // First find the source in the list and make sure it is a multi-frequency
  // photoionising source.
  //
  struct rad_source *src=0;
  for (vector<rad_source>::iterator i=SourceList.begin();
                                    i!=SourceList.end(); ++i)
    if ( (*i).s->id==rs->id ) src=&(*i);
  if (!src) {
    rep.error("update_RT_source_properties() Couldn't find source in source list.",rs->id);
  }
  //if (src->s->effect != RT_EFFECT_PION_MULTI) {
  //  rep.error("update_RT_source_properties() Don't know how to update source type for this id",src->s->id);
  //}
  //
  // copy parameters.
  //
  for (int v=0; v<MAX_DIM;v++)
    src->s->pos[v] = rs->pos[v];
  src->s->strength = rs->strength;
  src->s->Rstar = rs->Rstar;
  src->s->Tstar = rs->Tstar;
  //src->s->id = rs->id;
  //src->s->type = rs->type;
  //src->s->update = rs->update;
  //src->s->at_infinity = rs->at_infinity;
  //src->s->effect = rs->effect;
  //src->s->NTau = rs->NTau;
  //src->s->opacity_src = rs->opacity_src;
  //src->s->opacity_var = rs->opacity_var;
  //src->s->EvoFile = rs->EvoFile;

    
  return;
}


// ##################################################################
// ##################################################################


int raytracer_USC_infinity::populate_ionising_src_list(
      std::vector<struct rt_source_data> &ion_list
      ///< list of data for ionising sources
      )
{
#ifdef RT_TESTING
  if (ion_list.size() != static_cast<unsigned int>(N_ion_srcs)) {
    cout <<"Wrong sized list passed to populate_ionising_src_list:";
    cout <<" list size: "<<ion_list.size();
    cout <<", N_ion_srcs: "<<N_ion_srcs<<"\n";
    rep.error("Wrong size rt_USC_infty::populate_ionising_src_list",
              N_ion_srcs);
  }
#endif // RT_TESTING

  for (int v=0; v<N_ion_srcs; v++) ion_list[v] = ION_data[v];

  return 0;
}

// ##################################################################
// ##################################################################


int raytracer_USC_infinity::populate_UVheating_src_list(
      std::vector<struct rt_source_data> &uvh_list
      ///< list of data for UV-heating sources
      )
{
#ifdef RT_TESTING
  if (uvh_list.size() != static_cast<unsigned int>(N_uvh_srcs)) {
    cout <<"Wrong sized list passed to populate_uvheating_src_list:";
    cout <<" list size: "<<uvh_list.size();
    cout <<", N_uvh_srcs: "<<N_uvh_srcs<<"\n";
    rep.error("Wrong size rt_USC_infty::populate_uvheat_src_list",
              N_uvh_srcs);
  }
#endif // RT_TESTING

  for (int v=0; v<N_uvh_srcs; v++) uvh_list[v] = UVH_data[v];

  return 0;
}




// ##################################################################
// ##################################################################



int raytracer_USC_infinity::NSources()
{
  return SourceList.size();
}



// ##################################################################
// ##################################################################


void raytracer_USC_infinity::Print_SourceList()
{
  vector<rad_source>::iterator i;
  for (i=SourceList.begin(); i!=SourceList.end(); ++i) {
    cout <<"Source "<<(*i).s->id<<" at position ["<<(*i).s->pos[XX];
    if (ndim>1) cout <<", "<<(*i).s->pos[YY];
    if (ndim>2) cout <<", "<<(*i).s->pos[ZZ];
    cout <<"] has strength "<<(*i).s->strength;
    if ( (*i).sc !=0 )
      cout <<" and cell id "<<(*i).sc->id;
    cout <<" on grid? (Y=1,N=0) "<<(*i).src_on_grid;
    cout <<"\n";
  }
  return;
}



// ##################################################################
// ##################################################################


int raytracer_USC_infinity::RayTrace_Column_Density(
      const int s_id,  ///< Source id
      const double dt, ///< Timestep
      const double g   ///< EOS Gamma.
      )
{
#ifdef RT_TESTING
  cout <<"raytracer_USC_infinity::RayTrace_Column_Density()";
  cout <<" calling RayTrace_SingleSource().\n";
#endif // RT_TESTING

  struct rad_source *source=0;
  for (vector<rad_source>::iterator i=SourceList.begin();
                                    i!=SourceList.end(); ++i) {
    if ( (*i).s->id==s_id ) source=&(*i);
  }
  if (!source) {
    cerr <<"Couldn't find source "<<s_id<<" in source list.\n";
    return 1;
  }

  // now trace the source.
  int err = RayTrace_SingleSource(s_id,dt,g);
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC_infinity::RayTrace_SingleSource(
      const int s_id, ///< Source id
      const double dt, ///< Timestep
      const double g /// eos gamma
      )
{
  //cout <<"*******RT single source. \n";
  raytracer_USC_infinity::gamma = g;
  raytracer_USC_infinity::delt = dt;

#ifdef RT_TESTING
  //
  // Testing, to make sure we assign column densities to all cells.
  //
  cell *c = gridptr->FirstPt();
  double Tau[MAX_TAU];
  do {
    cout <<"cell "<<c->id<<" setting Tau to -1\n";
    //for (short unsigned int iT=0;
    //     iT<SourceList[s_id].s->NTau; iT++)
    for (short unsigned int iT=0; iT<MAX_TAU; iT++) Tau[iT]=-1.0;
    CI.set_col(c, s_id, Tau);
  } while ( (c=gridptr->NextPt(c)) !=0);
#endif // RT_TESTING

  struct rad_source *source=0;
  for (vector<rad_source>::iterator i=SourceList.begin();
                                    i!=SourceList.end(); ++i) {
    if ( (*i).s->id==s_id ) source=&(*i);
  }
  if (!source)
  {
    cerr <<"Couldn't find source "<<s_id<<" in source list.\n";
    return 1;
  }
  //  cout <<"found source. moving to it.\n";

  enum direction dir=NO;
  for (int i=0;i<ndim;i++) {
    if      (source->s->pos[i] < -1.e99)
      dir=static_cast<direction>(2*i);
    else if (source->s->pos[i] >  1.e99)
      dir=static_cast<direction>(2*i+1);
  }
  if (dir==NO) rep.error("source not at infinity in any direction",dir);

  int err = trace_parallel_rays(source, dir);
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC_infinity::trace_parallel_rays(
      const rad_source *source, ///< source we are dealing with.
      const enum direction dir  ///< direction to source at infinity.
      )
{
  int err = 0;
  cell *c = source->sc; 
  // make sure source cell is at the edge of the grid, in direction dir.
  if (gridptr->NextPt(c,dir)!=0 && gridptr->NextPt(c,dir)->isgd)
    rep.error("source cell not set up right",source->s->id);
#ifdef RT_TESTING
  cout <<"raytracer_USC_infinity::trace_parallel_rays()";
  cout <<": start_cell: ";
  CI.print_cell(c);
#endif
  
  enum direction oppdir=gridptr->OppDir(dir);
  //
  // Now start at each cell on this boundary, and go in a column in oppdir.
  //
  if      (ndim==1) {
    err += trace_column_parallel(source,c,oppdir);
  }

  else if (ndim==2) {
    enum direction perpdir = NO;
    if      (dir==XP || dir==XN) perpdir=YP;
    else if (dir==YP || dir==YN) perpdir=XP;
    else rep.error("Source not in xy plane in 2d sim!",dir);
    do {
      err += trace_column_parallel(source,c,oppdir);
    } while ( (c=gridptr->NextPt(c,perpdir))!=0);
  } // 2D
  
  else if (ndim==3) {
    enum direction pdir1=NO, pdir2=NO;
    if      (dir==XP || dir==XN) {pdir1=YP; pdir2=ZP;}
    else if (dir==YP || dir==YN) {pdir1=XP; pdir2=ZP;}
    else if (dir==ZP || dir==ZN) {pdir1=XP; pdir2=YP;}
    else rep.error("bad dir in trace_parallel_rays()",dir);
    cell *cp;
    do {
      cp=c;
      do {err += trace_column_parallel(source,cp,oppdir);}
      while ( (cp=gridptr->NextPt(cp,pdir1))!=0);
    } while ( (c =gridptr->NextPt(c, pdir2))!=0);
  } // 3D
  
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC_infinity::trace_column_parallel(
      const rad_source *source, ///< source we are tracing from.
      cell *c,                  ///< cell to start from.
      const enum direction dir  ///< direction we are looking.
      )
{
  double ds=0.0, Nc[MAX_TAU];
  for (unsigned short int iT=0; iT<source->s->NTau; iT++)
    Nc[iT] = 0.0;

  int    err=0;
  enum direction oppdir=gridptr->OppDir(dir);

#ifdef RT_TESTING
  //
  // If we are at the source cell, we don't get the column to it, so
  // we treat it specially.
  //
  if (source->src_on_grid) {
    rep.error("Can't have source on grid for parallel rays!",
              source->s->id);
  }
#endif // RT_TESTING
  
  do {
#ifdef TESTING
    //      cout <<"setting cell pointer.\n";
    dp.c = c;
#endif
    err += cell_cols_1d(source, c, oppdir, Nc, &ds);
    err += ProcessCell(c,Nc,ds,source,delt);
    
    //cout <<"testing...\n";
    if (c->Ph[PG]<0.0 || !isfinite(c->Ph[PG]))
      rep.error("ProcessCell() gives negative energy!",c->Ph[PG]);

  } while ( (c=gridptr->NextPt(c,dir))!=0);

  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC_infinity::cell_cols_1d(
      const rad_source *src,
      cell *c,
      const enum direction sdir, // direction to source
      double Nc[],
      double *ds
      )
{
  if (!c) rep.error("cell_cols_1d() null cell",c);

  if (c==src->sc && !src->s->at_infinity && src->src_on_grid) {
    //
    // Cell is source cell (can happen for 1D spherical grids)
    //
    cout <<"source cell!source cell!source cell!source cell!\n";
    for (unsigned short int iT=0; iT<src->s->NTau; iT++) {
      Nc[iT] = 0.0;
    }
  }
  else{
    cell *ngb = gridptr->NextPt(c,sdir);
    // assume if neighbour doesn't exist, that the source is coming
    // in from off grid to cell c.
    if (ngb) {
      // boundary cells should have zero column (or some value if MPI
      // boundary)
      CI.get_col(ngb, src->s->id, Nc);
    }
    else {
      for (unsigned short int iT=0; iT<src->s->NTau; iT++)
        Nc[iT] = 0.0;
    }
  }

  //
  // column through cell is very easy in 1D with uniform grid.
  // (ds is back in physical units now!)
  //
  *ds = gridptr->DX();
  return 0;
}



// ##################################################################
// ##################################################################



int raytracer_USC_infinity::ProcessCell(
      cell *c,                  ///< Current cell.
      double col2cell[],        ///< Column to cell (optical depth)
      double ds,                ///< Path length through cell.
      const rad_source *source, ///< pointer to source struct.
      const double dt          ///< Timestep
      )
{
  //
  // TODO: Switch this to actually use optical depths, because
  //   otherwise TauMin is not dimensionless, and it is confusing
  //   to use it.
  //
  double cell_col[MAX_TAU];
  for (unsigned short int iT=0; iT<source->s->NTau; iT++)
    cell_col[iT] = ds;
  double temp=0.0;
  short unsigned int ix=0;
#ifdef RT_TESTING
  cout <<"\tsrc "<<source->s->id<<": Updating cell id="<<c->id;
  cout <<" with opacity-only step.  N[0]="<<col2cell[0]<<"\n";
#endif // RT_TESTING

  // set vshell in cell (setting up)
  if (source->s->opacity_src == RT_OPACITY_VSHELL) {
    //cout <<"set vshell\n";
    set_Vshell_in_cell(c, ds, source);
    cell_col[0] = ds;
    col2cell[0] += cell_col[0];
    //CI.set_cell_col(c, source->s->id, cell_col);
    //CI.set_col     (c, source->s->id, col2cell);
    return 0;
  }

  // 2 special cases: cell is not part of domain, so nothing is set,
  // or else cell is not a leaf, in which case col2cell/cell_col are
  // received from finer grid.
  if (!c->isdomain) {
    // if cell is not in the domain, set its column to be zero,
#ifdef RT_TESTING
    cout <<"off domain: "<<c->id<<", "; double dpos[ndim];
    CI.get_dpos(c,dpos); rep.printVec("pos",dpos,ndim);
#endif
    cell_col[0] = 0.0;
    CI.set_cell_col(c, source->s->id, cell_col);
    col2cell[0]  = 0.0;
    CI.set_col     (c, source->s->id, col2cell);
    return 0;
  }

  // if cell is not a leaf, then cell optical depths have already
  // been set and we can just return.
  // Just check that the interpolated d-Tau is less than Tau (which
  // is the sum of d-Tau and the optical depth from source to cell).
  //
  if (!c->isleaf) {
    CI.get_col(c,source->s->id, col2cell);
    CI.get_cell_col(c,source->s->id, cell_col);
#ifdef RT_TESTING
    cout <<"NG grid cell: "<<c->id <<" isdomain="<<c->isdomain;
    cout <<", col2cell="<<*col2cell<<", cell_col="<<*cell_col<<"  ";
    rep.printVec("pos",c->pos,ndim);
#endif
    if (col2cell[0]<cell_col[0]) {
      //cout <<"NG grid cell: "<<c->id <<" isleaf="<<c->isleaf;
      //cout <<" isdom="<<c->isdomain;
      //cout <<", col2cell="<<*col2cell<<", cell_col="<<*cell_col<<" ";
      //rep.printVec("pos",c->pos,ndim);
      //cout <<"\n";
      //CI.print_cell(c);
      //rep.error("negative col2cell",col2cell[0]-cell_col[0]);
      for (unsigned short int iT=0; iT<source->s->NTau; iT++)
        col2cell[iT] = cell_col[iT];
      CI.set_col(c, source->s->id, col2cell);
    }
    return 0;
  }

  // Now cell is a leaf (or boundary data), so we have to set the
  // optical depth to and through the cell.
  //
  // set cell_col based on opacity flag.  cell_col[]=ds from above,
  // so multiply it by the cell-values here.
  switch (source->s->opacity_src) {

  case RT_OPACITY_MINUS:
    // generic flag for integrating rho*(1-y_i)*ds
    cell_col[0] = ds*c->Ph[RO]*
                 (1.0-c->Ph[source->s->opacity_var+first_tr]);
    break;

  case RT_OPACITY_TOTAL:
    cell_col[0] = ds*c->Ph[RO];
    break;

  case RT_OPACITY_TRACER:
    // this is a generic flag for integrating rho*y_i*ds
    cell_col[0] = ds*c->Ph[RO]*
                   c->Ph[source->s->opacity_var+first_tr];
    break;

  case RT_OPACITY_HALPHA:
    //
    // This calculate the H-alpha intensity along a line-of-sight, assuming
    // a constant source function within a cell, and a grey opacity which
    // scales with the total gas density dtau=n(H)*5.0e-22*ds.
    // S=0.22*n(e)*n(H+)/n(H)/T^0.9.  Use cell_col for source function S.
    // Then I(out) = S +exp(-dtau)*(I(in)-S)
    //
    cell_col[0] = 0.22*(c->Ph[RO]*mpptr->get_X_H()/pconst.m_p())
              *1.1*c->Ph[source->s->opacity_var+first_tr]
                  *c->Ph[source->s->opacity_var+first_tr]
              *pow(MP->Temperature(c->Ph,gamma),-0.9);
    // now set cell_col[0] to be I(out)
    cell_col[0] = cell_col[0]+
        exp(-(c->Ph[RO]*mpptr->get_X_H()/pconst.m_p())*5.0e-22*ds)
        *(col2cell[0]-cell_col[0]);
    break;


  case RT_OPACITY_RR:
    // Get the radiative recombination rate *Vshell, for the
    // photoionisation equilibrium assumption (with caseB).
    // We assume the cell is fully ionised, so we set the
    // first tracer value to 1.0 in c->Ph[] when we call the rate
    // function, and reset it afterwards.
    // HARDCODED ASSUMPTION THAT 1ST TRACER IS ION FRACTION.
    temp = c->Ph[first_tr];
    c->Ph[first_tr] = 1.0;
    cell_col[0] = CI.get_cell_Vshell(c,source->s->id)/
                source->s->strength*
                MP->get_recombination_rate(0, c->Ph,gamma);
    c->Ph[first_tr] = temp;
    break;

  case RT_OPACITY_HHE:
    //
    // Here we have four column densities to set.
    // dTaui[i] = n[i]*sigma[i]*ds
    // n[i] = rho*X_i/m_i
    // 
    // Here 0 = H0,  1 = He0,  2 = He+,  3 = Dust
    // Element/Ion constants are defined in constants.h
    //
    // ---------- H0 -----------
    ix = source->s->opacity_var+first_tr;
    cell_col[0] = ds*MP->get_n_el(c->Ph, EL_H)*
                  MP->get_th_xsection(ION_H_N)*c->Ph[ix];
    // ---------- He0 -----------
    ix++;
    cell_col[1] = ds*MP->get_n_el(c->Ph, EL_HE)*
                  MP->get_th_xsection(ION_HE_N)*c->Ph[ix];
    // ---------- He+ -----------
    ix++;
    cell_col[2] = ds*MP->get_n_el(c->Ph, EL_HE)*
                  MP->get_th_xsection(ION_HE_P)*c->Ph[ix];
    // ---------- DUST -----------
    cell_col[3] = ds*MP->get_n_el(c->Ph, EL_H)*
                  MP->get_th_xsection(ION_DUST);
    break;

  default:
    rep.error("RT-ProcessCell opacity?",source->s->opacity_src);
    break;
  }
  CI.set_cell_col(c, source->s->id, cell_col);


  if (source->s->opacity_src == RT_OPACITY_HALPHA) {
    // set col2cell[0] to be I(out)-I(in)
    col2cell[0] = cell_col[0]-col2cell[0];
  }
  else {
    for (unsigned short int iT=0; iT<source->s->NTau; iT++)
      col2cell[iT] += cell_col[iT];
  }
  CI.set_col(c, source->s->id, col2cell);

#ifdef RT_TESTING
  cout <<"NG grid cell: "<<c->id <<" isdomain="<<c->isdomain;
  cout <<", col2cell="<<*col2cell<<", cell_col="<<*cell_col<<"  ";
  CI.get_col(c,source->s->id, col2cell);
  CI.get_cell_col(c,source->s->id, cell_col);
  cout <<", c2c="<<*col2cell<<", cc="<<*cell_col<<"  ";
  rep.printVec("pos",c->pos,ndim);
#endif
  
  return 0;
}



// ##################################################################
// ##################################################################



void raytracer_USC_infinity::set_Vshell_in_cell(
      cell *c, ///< current cell.
      double ds,           ///< Path Length through cell.
      const rad_source *rs ///< pointer to source struct.
      )
{
  //
  // ds *is* Vshell for parallel rays.
  //
  //cout <<"SRC-at-INF: vshell="<<ds<<", ds="<<ds<<"\n";
  CI.set_cell_Vshell(c, rs->s->id, ds);
  CI.set_cell_deltaS(c, rs->s->id, ds);
  return;
}



// ##################################################################
// ##################################################################


// *********************************************************************
// *********************************************************************
// ********************* END PARALLEL RAYS *****************************
// *********************************************************************
// *********************************************************************







// ##################################################################
// ##################################################################








// *********************************************************************
// *********************************************************************
// ****************** BEGIN SHORT CHARACTERISTICS **********************
// *********************************************************************
// *********************************************************************


// ##################################################################
// ##################################################################



raytracer_USC::raytracer_USC(
      class GridBaseClass *ggg,    ///< Pointer to grid
      class microphysics_base *mmm, ///< Pointer to MicroPhysics Class.
      int nd,     ///< number of dimensions of grid
      int csys,   ///< coordinate system
      int nv,     ///< number of variables in state vector
      int ftr,    ///< index of first tracer variable in state vector
      int Nsources  ///< Number of radiation sources
      )
  : raytracer_USC_infinity(ggg,mmm,nd,csys,nv,ftr)
{
#ifdef RT_TESTING
  cout <<"SC raytracer class constructor!\n";
  cout <<"  ndim = "<<ndim<<" and grid pointer = "<<gridptr<<"\n";
#endif
  if (ndim <1 || ndim>3) rep.error("bad ndim for raytracing!",ndim);
  if (coord_sys==COORD_CYL && ndim!=2)
    rep.error("only know 2d cylindrical coords",ndim);
  raytracer_USC::SrcDir = std::vector<enum direction> (ndim,NO);
   
  // These are the outward directions from the source for each octant.
  // XP,XN,XP,XN,XP,XN,XP,XN
  // YP,YP,YN,YN,YP,YP,YN,YN
  // ZP,ZP,ZP,ZP,ZN,ZN,ZN,ZN
  for (int i=0; i<8; i++) {
    if ( i%2     ==0) dir1[i]=XP; else dir1[i]=XN;
    if ( (i/2)%2 ==0) dir2[i]=YP; else dir2[i]=YN;
    if ( (i/4)%2 ==0) dir3[i]=ZP; else dir3[i]=ZN;
  }
  //  rep.printVec("dir1",dir1,8);
  //  rep.printVec("dir2",dir2,8);
  //  rep.printVec("dir3",dir3,8);

  //
  // Set the source Tau-min values vector to have Nsources elements.
  //
  TauMin.resize(Nsources);
  return;
}



// ##################################################################
// ##################################################################



raytracer_USC::~raytracer_USC()
{
#ifdef RT_TESTING
  cout <<"SC raytracer class destructor!\n";
#endif
  TauMin.clear();
}



// ##################################################################
// ##################################################################



int raytracer_USC::Add_Source(
      struct rad_src_info *src ///< source info.
      )
{
  //  cout <<"AddSource() adding source to list.\n";
#ifdef RT_TESTING
  cout <<"\n--BEGIN-----raytracer_USC::AddSource()------------\n";
#endif
  //
  // First see if the source is at infinity, and if so then call the parallel
  // rays setup function.
  //
  if (src->at_infinity) {
    raytracer_USC_infinity::add_source_to_list(src);
  }
  else {
    add_source_to_list(src);
  }

  //
  // Now we need to set Vshell in every grid point for this source.  
  //
  set_Vshell_for_source(&SourceList.back());

  return SourceList.back().s->id;
}



// ##################################################################
// ##################################################################



void raytracer_USC::add_source_to_list(
      struct rad_src_info *src ///< source info.
      )
{
  //
  // Create a new radiation source struct.
  //
  rad_source rs;

  //
  // Set pointer in rad_source to *src.  This sets the position, id,
  // and all the parameters necessary to set up the source.
  //
  rs.s = src;

#ifdef RT_TESTING
  cout <<"\t\t"; rep.printVec("Input Source Position",rs.s->pos,ndim);
#endif


  //
  // Set the source-dependent parts of the rt_source_data struct.
  //
  rs.data.id       = rs.s->id;
  rs.data.type     = rs.s->type;
  rs.data.strength = rs.s->strength;
  rs.data.Vshell   = 0.0;
  rs.data.dS       = 0.0;
  rs.data.NTau     = rs.s->NTau;
#ifdef RT_TESTING
  cout <<"***** NTAU = "<<rs.s->NTau<<"\n";
#endif
  for (unsigned short int iT=0; iT<rs.s->NTau; iT++) {
    rs.data.DelCol[iT] = 0.0;
    rs.data.Column[iT] = 0.0;
  }

  //
  // now find the source cell so we can set up a pointer to it.  If
  // source is not on grid, find edge cell closest to it.  This
  // function also centres the source on a cell-vertex.
  //
#ifdef RT_TESTING
  cout <<"RTxxx: AddSource() finding source.\n";
#endif
  rs.sc = find_source_cell(rs.s->pos);

  //
  // now we have set the source's position, get the integer position.
  //
  CI.get_ipos_vec(rs.s->pos,rs.ipos);
#ifdef RT_TESTING
  cout <<"\t\t"; rep.printVec("Assigned source ipos",rs.ipos,ndim);
  cout <<"\t\t"; rep.printVec("Assigned source dpos",rs.s->pos,ndim);
  cout <<"\tSERIAL: AddSource() source->sc = "<<rs.sc<<"\n";
#endif

  //
  // Now determine if source is on grid or not.
  //
  rs.src_on_grid = true;
  //
  // For non-cell-centred source, if it is at a cell corner then I
  // assume the source cell is the one with the source at its most
  // negative corner.  We set the source off-grid in a direction even
  // if it is on the boundary.
  //
  for (int i=0; i<ndim; i++) {
    if (rs.ipos[i]<=gridptr->iXmin_all(static_cast<axes>(i)))
      rs.src_on_grid=false;
    if (rs.ipos[i]>=gridptr->iXmax_all(static_cast<axes>(i)))
      rs.src_on_grid=false;
  }
  //
  // ipos can be overloaded if source is at infinity, so extra check:
  //
  if (rs.s->at_infinity) rs.src_on_grid=false;

  //
  // nearly all done, so now add the source to the list, update local
  // source lists, set TauMin for this source, and return its id.
  //
  SourceList.push_back(rs);
  update_local_variables_for_new_source(rs);
  set_TauMin_for_source(rs);

#ifdef RT_TESTING
  if (TauMin[rs.s->id]<0) rep.error("Duhhh",TauMin[rs.s->id]);
  cout <<"\t\tSource id:"<<rs.s->id<<" has TauMin[id]=";
  cout <<TauMin[rs.s->id]<<"\n";
  cout <<"\tSERIAL: AddSource() all done.\n";
  cout <<"Add_Source() source->sc = "<<rs.sc<<", id="<<rs.sc->id<<"\n";
  CI.print_cell(rs.sc);
  cout <<"--END-----raytracer_USC::AddSource()------------\n\n";
#endif

  return;
}



// ##################################################################
// ##################################################################



void raytracer_USC::set_TauMin_for_source(
      const struct rad_source rs
      )
{
  //
  // THIS CODE IS A MESS!  NEED A FUNCTION IN THE MICROPHYSICS CLASS
  // THAT TAKES AN INPUT PRIMITIVE STATE VECTOR AND A SOURCE AND THE
  // PATH LENGTH, AND RETURNS DELTA-TAU FOR EACH ENERGY RANGE.
  // THEN RAYTRACING CAN GO BACK TO HAVING TauMin=0.7 AS A CONSTANT.
  //
  // Now set TauMin for this source.
  //
  if (TauMin.size() <= static_cast<unsigned int>(rs.s->id)) {
    rep.error("Source id is larger than Nsources!",rs.s->id-TauMin.size());
  }
  TauMin[rs.s->id]=-1.0;

  if      ((rs.s->update==RT_UPDATE_EXPLICIT) &&
           (rs.s->opacity_src==RT_OPACITY_TOTAL) &&
           (rs.s->effect==RT_EFFECT_UV_HEATING)) {
    // This is for the UV heating, sigma=5e-22, Tau=1.9Av=1.9*1.086*NH*sigma
    //cout <<"UV heating source, TauMIN=0.0016\n";
    TauMin[rs.s->id] = 1.6e-3;
  }
  else if ((rs.s->update==RT_UPDATE_EXPLICIT) &&
           (rs.s->opacity_src==RT_OPACITY_MINUS) &&
           (rs.s->effect==RT_EFFECT_PION_MONO) ) {
    // Photoionisation, sigma=6.3e-18, Tau=sigma*NH*(1-x)
    //TauMin[rs.s->id] = 6.1e-7;
    TauMin[rs.s->id] = 1.86e-7;
    // --------------------------------------------------------------
  }
  else if ((rs.s->update==RT_UPDATE_EXPLICIT) &&
           (rs.s->opacity_src==RT_OPACITY_MINUS) &&
           (rs.s->effect==RT_EFFECT_PION_MULTI)) {
    // Photoionisation, sigma=6.3e-18, Tau=sigma*NH*(1-x)
    // This value seem to give best results for multifrequency sources, tested
    // for Teff of both 30000 and 39750K.
#ifdef PUREHYDROGEN
    // Iliev et al. (2009) test problems (Cosmology Code Comparison
    // Project).
    // Teff=1e5K means mean-photon-erg=29.62eV. This means cross section is
    // ~6.1e-19 cm2. mu=1.6738e-24g.
    // So TauMin=0.7 means SigmaMin=0.7*mu/sigma=1.92e-6.
    // 3.0e-6 seems to work better though.
    //
    TauMin[rs.s->id] = 3.0e-6;
#else
    // 35kK, mu=2.338e-24g, gives about 5.9e-7
    //TauMin[rs.id] = 6.0e-7; // original value
    TauMin[rs.s->id] = 7.0e-7;
#endif
    TauMin[rs.s->id] = 1.86e-7;

  }
  else if ((rs.s->update==RT_UPDATE_EXPLICIT) &&
           (rs.s->effect==RT_EFFECT_UV_HEATING) ) {
    // This is for Harpreet's column densities (value used is unimportant).
    rep.warning("Special case opacity for Harpreet's project!",
                rs.s->id,rs.s->id);
    TauMin[rs.s->id] = 5.0e-7;
  }
  else if ((rs.s->update==RT_UPDATE_EXPLICIT) &&
           (rs.s->effect==RT_EFFECT_PION_EQM) ) {
    //
    // This is for photoionisation equilibrium; the variable is the
    // fractional attenuation per cell.  So tau=0.7 corresponds to a
    // value of (1-exp(-0.7))=0.5034
    //
    TauMin[rs.s->id] = 0.5034;
  }

  else if (rs.s->effect==RT_EFFECT_HHE_MFQ) {
    //
    // This class uses actual optical depths (I'm going to change
    // everything so that it is default from now on).  Experiments
    // showed that TauMin=0.7 is the best compromise in terms of
    // photon conservation and shadowing and spherical HII regions.
    //
    TauMin[rs.s->id] = 0.7;
  }


  else rep.error("Unhandled case setting TauMin",rs.s->id);

  //
  // Seems that TauMin=0.7 is best in 2D, but =0.6 is better in 3D.
  //
  if (ndim==3) TauMin[rs.s->id] *= 6.0/7.0;

  return;
}



// ##################################################################
// ##################################################################



int raytracer_USC::NSources()
{
  return SourceList.size();
}



// ##################################################################
// ##################################################################



void raytracer_USC::Print_SourceList()
{
  vector<rad_source>::iterator i;
  for (i=SourceList.begin(); i!=SourceList.end(); ++i) {
    cout <<"Source "<<(*i).s->id<<" at position ["<<(*i).s->pos[XX];
    if (ndim>1) cout <<", "<<(*i).s->pos[YY];
    if (ndim>2) cout <<", "<<(*i).s->pos[ZZ];
    cout <<"] has strength "<<(*i).s->strength;
    if ( (*i).sc !=0 )
      cout <<" and cell id "<<(*i).sc->id;
    cout <<" on grid? (Y=1,N=0) "<<(*i).src_on_grid;
    cout <<".   at infinity? "<<(*i).s->at_infinity<<".  ";
    rep.printVec("integer posn",(*i).ipos,ndim);
  }
  return;
}



// ##################################################################
// ##################################################################


int raytracer_USC::RayTrace_SingleSource(
      const int s_id,  ///< Source id
      const double dt, ///< Timestep
      const double g   ///< eos gamma.
      )
{
#ifdef RT_TESTING
  cout <<"raytracer_USC::RayTrace_SingleSource() starting.\n";
  Print_SourceList();
#endif
  delt = dt;
  gamma = g;

#ifdef RT_TESTING
  //
  // Testing, to make sure we assign column densities to all cells.
  //
  cell *sc = gridptr->FirstPt();
  double tau[MAX_TAU];
  for (unsigned short int iT=0;iT<MAX_TAU;iT++) tau[iT]=-1.0;
  do {
    CI.set_col(sc, s_id, tau);
  } while ( (sc=gridptr->NextPt(sc)) !=0);
#endif // RT_TESTING

  rad_source *source=0;
  for (vector<rad_source>::iterator i=SourceList.begin();
                                    i!=SourceList.end(); ++i) {
    if ( (*i).s->id==s_id ) {
      source=&(*i);
      current_src = &(*i);
    }
  }
  if (!source) {
    cerr <<"Couldn't find source "<<s_id<<" in source list.\n";
    return 1;
  }
  //  cout <<"found source. moving to it.\n";

#ifdef RT_TESTING
  sc=source->sc;
#endif // RT_TESTING

  enum direction src_off_grid[ndim], src_at_infty[ndim];
  for (int i=0;i<ndim;i++) src_off_grid[i] = NO;
  for (int i=0;i<ndim;i++) src_at_infty[i] = NO;
  if (!source->src_on_grid) {
    for (int i=0;i<ndim;i++) {
      enum axes a = static_cast<axes>(i);
      enum direction posdir = static_cast<direction>(2*static_cast<int>(a)+1);
      enum direction negdir = gridptr->OppDir(posdir);
      
      //
      // Use integer positions (note these are not used if the source is
      // at infinity, so it doesn't matter if ipos[] is out_of_range.
      //
      if (source->ipos[a] <= gridptr->iXmin_all(a)) src_off_grid[a] = negdir;
      if (source->ipos[a] >= gridptr->iXmax_all(a)) src_off_grid[a] = posdir;
      
      if (source->s->pos[a] <=-1.e99) src_at_infty[a] = negdir;
      if (source->s->pos[a] >= 1.e99) src_at_infty[a] = posdir;
    }
  }
  int err=0;

  // check if the source is at infinity and, if so, call
  // parallel_rays routine.
  int ct=0;
  for (int i=0;i<ndim;i++) if (src_at_infty[i]!=NO) ct++;
  if (ct) {
    if (ct>1)
      rep.error("source at infinity not parallel to axis!",ct);
    enum direction dir=NO;
    for (int i=0;i<ndim;i++)
      if (src_at_infty[i]!=NO)
        dir = src_at_infty[i];
    err += trace_parallel_rays(source, dir);
    return err;
  }
  
  // source is not at infinity, so set list of start cells for each
  // quadrant.
  int ndirs=1; for (int i=0;i<ndim;i++) ndirs*=2;
  cell *startcell[ndirs];
  for (int v=0;v<ndirs;v++) startcell[v]=0;
#ifdef RT_TESTING
  cell *temp=sc; // for testing.
#endif // RT_TESTING
  set_startcells(source->sc,startcell,src_off_grid);
#ifdef RT_TESTING
  if (sc!=temp) rep.error("set_startcells() changed sc!",sc); // testing only
  cout <<"start cells!: \n";
  for (int v=0;v<ndirs;v++) CI.print_cell(startcell[v]);
#endif // RT_TESTING

  // Loop over all lines/quads/octs.
  for (int oct=0; oct<ndirs; oct++) {
    enum direction dirs[ndim];
    cell *c = startcell[oct];
    //    cout <<"oct "<<oct<<" and startcell = "<<c<<" dirs: ";
    //    cout <<dir1[oct]<<" "<<dir2[oct]<<" "<<dir3[oct]<<"\n";
    //
    // now, if c!=0, then we have cell(s) in the octant, so trace the octant.
    //
    if (c) {
      // set outward directions for this octant.
      dirs[XX] = dir1[oct];
      dirs[YY] = dir2[oct];
      dirs[ZZ] = dir3[oct];
      // octant directions are always away from source.
      for (int i=0;i<ndim;i++) SrcDir[i] = gridptr->OppDir(dirs[i]);

      if      (ndim==3)
        err += trace_octant(source, c, dirs[XX], dirs[YY], dirs[ZZ]);
      else if (ndim==2)
        err += trace_plane (source, c, dirs[XX], dirs[YY]);
      else if (ndim==1)
        err += trace_column(source, c, dirs[XX]);
      else rep.error("bad ndim in RayTrace_SingleSource()",ndim);
    }

  } // trace out 2 lines/4 quadrants/8 octants.
  
  //cout <<"raytracer_USC::RayTrace_SingleSource() finished.\n";
  return 0;
}



// ##################################################################
// ##################################################################



void raytracer_USC::set_startcells(
      cell *sc,      ///< source cell, or cell nearest to source if off grid.
      cell **startcell, ///< list of startcells for each line/quadrant/octant.
      enum direction *src_off_grid ///< list of dirs if source is off grid.
      )
{
  if (ndim==1) {
    bool in_my_half = false;
    // Loop over two directions along the 1D line.
    for (int i=0; i<2; i++) {
      cell *c = sc;
      enum direction dir = dir1[i];

      // First make sure this direction contains some cells to traverse.
      if      (dir==XP) {
	if (src_off_grid[XX]==XP) in_my_half = false;
	else in_my_half = true;
	if (in_my_half) startcell[i] = c; else startcell[i]=0;
      }
      else if (dir==XN) {
	if      (src_off_grid[XX]==XN) in_my_half = false;
	else if (src_off_grid[XX]==NO) {
	  c = gridptr->NextPt(c,dir);
	  if ( (c!=0) ) in_my_half = true;
	  else in_my_half=false;
	}
	else if (src_off_grid[XX]==XP) {
	  if ( (c!=0) ) in_my_half = true;
	  else rep.error("logic error in RayTraceSource()",c);
	}
	else rep.error("Bad direction",src_off_grid[XX]);
	if (in_my_half) startcell[i] = c; else startcell[i]=0;
      }
      else rep.error("Bad direction in 1D RayTraceSource()",dir);
    } // directions
  } //1D

  else if (ndim==2) {
    bool in_my_quad = false;
    // loop over quads
    //  cout <<"looping over quadrants. source cell ="<<sc->id<<"\n";
    for (int i=0; i<NQUADS; i++) {
      // Now have 4 quadrants to go around, instead of two sides of
      // a line.
      cell *c = sc;
      //  test if starting cell is in my quad
      if      (i==Q1) { // XP,YP quadrant.
	if (src_off_grid[XX]==XP || src_off_grid[YY]==YP) in_my_quad = false;
	else in_my_quad = true;
	if (in_my_quad) startcell[i] = c; else startcell[i]=0;
	// cout <<"in Q1, ? bool = "<<in_my_quad<<" ";
        // cout<<src_off_grid[YY]<<" "<<src_off_grid[XX]<<"\n";
      }
      
      else if (i==Q2) { // XN,YP quadrant
	if (src_off_grid[XX]==XN || src_off_grid[YY]==YP) {
          in_my_quad = false;
          // src off grid to top or left.
        }
	else if (src_off_grid[XX]!=XP) {
          // source is on grid in x-dir, so sc is in Q1, so move one
          // cell left.
	  c = gridptr->NextPt(c,XN);
	  if (c) in_my_quad = true;
	  else in_my_quad = false;
	}
	else // src off grid to right, so sc is first cell in Q2.
	  in_my_quad = true;
	if (in_my_quad) startcell[i] = c; else startcell[i]=0;
      }
      
      else if (i==Q4) { // XP,YN quadrant.
	if (src_off_grid[XX]==XP || src_off_grid[YY]==YN) {
          in_my_quad = false;
          // src off grid to bottom or right.
        }
	else if (src_off_grid[YY]!=YP) {
          // source is on grid in y-dir, so sc is in Q1, so move one cell YN.
	  c = gridptr->NextPt(c,YN);
	  if (c) in_my_quad = true;
	  else in_my_quad = false;
	}
	else {
          // src off grid to top but not to right, so sc is first cell in Q4.
	  in_my_quad = true;
        }
	if (in_my_quad) startcell[i] = c; else startcell[i]=0;
      }
      
      else if (i==Q3) { // XN,YN quadrant.
	if (src_off_grid[YY]!=YP) {
          // source is not off grid in +ve y-dir, so move one cell YN.
	  c = gridptr->NextPt(c,YN);
	}
	if (c!=0 && src_off_grid[XX]!=XP) {
          // source is not off-grid in +ve x-dir, so move one cell XN.
	  c = gridptr->NextPt(c,XN);
	}
	if (c!=0) in_my_quad = true;
	else in_my_quad = false;
	if (in_my_quad) startcell[i] = c; else startcell[i]=0;
      }
      else rep.error("logic error!!",i);
    }
  } // 2D

  else if (ndim==3) {
    enum direction dirs[3];
    for (int oct=0; oct<8; oct++) { // Now have 8 octants to go around
      bool in_my_oct = true;
      // set outward directions for this octant.
      dirs[XX] = dir1[oct];
      dirs[YY] = dir2[oct];
      dirs[ZZ] = dir3[oct];
      // first check the obvious... if the source is off grid in the
      // octant directions, then octant has no cells in it.
      for (int i=0;i<ndim;i++) {
	if (src_off_grid[static_cast<axes>(i)] == dirs[i]) {
          in_my_oct = false;
        }
      }
      
      if (in_my_oct==false) {
	startcell[oct]=0;
#ifdef RT_TESTING
	cout <<"octant "<<oct<<", source is off grid in this octant";
        cout <<"so start cell is NULL.\n";
#endif //RT_TESTING
      }
      else {
	//
	// octant may have cells, so check each one.
	//
#ifdef RT_TESTING
	cout <<"octant "<<oct<<", source cell:";
        CI.print_cell(sc);
	cout <<"octant "<<oct<<", startcell["<<oct<<"] = ";
        CI.print_cell(startcell[oct]);
#endif //RT_TESTING
	if      (oct==OCT1) startcell[oct] = sc;
	else if (oct==OCT2) {
	  if      (startcell[OCT1] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT1],XN);
            // may be zero!
          }
	  else if (src_off_grid[XX]==XP) {
            startcell[oct] = sc;
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT2",oct);
          }
	}
	else if (oct==OCT4) {
	  if      (startcell[OCT1] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT1],YN);
            // may be zero!
          }
	  else if (src_off_grid[YY]==YP) {
            startcell[oct] = sc;
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT4",oct);
          }
	}
	else if (oct==OCT3) {
#ifdef RT_TESTING
	  cout <<"OCT4:"<<oct<<" sc="<<sc<<" startcell[OCT2]=";
          cout <<startcell[OCT2];
	  cout <<" startcell[OCT4]="<<startcell[OCT4];
          rep.printVec(" src_off_grid",src_off_grid,ndim);
#endif // RT_TESTING
	  if      (startcell[OCT2] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT2],YN);
            // may be zero!
          }
	  else if (startcell[OCT4] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT4],XN);
            // may be zero!
          }
	  else if (src_off_grid[XX]==XP && src_off_grid[YY]==YP) {
            startcell[oct] = sc;
          }
	  // Extra condition: if 2 and 4 are both zero, then whole
          // plane must be (1,2,3,4)
	  else if (!startcell[OCT2] && !startcell[OCT4]) {
            startcell[oct] = 0;
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT3",oct);
          }
	}
	else if (oct==OCT5) { // XP,YP,ZN
	  if      (startcell[OCT1] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT1],ZN);
            // may be zero!
          }
	  else if (src_off_grid[ZZ]==ZP) {
            startcell[oct] = sc;
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT5",oct);
          }
	}
	else if (oct==OCT6) { // XN,YP,ZN
	  if      (startcell[OCT5] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT5],XN);
            // may be zero!
          }
	  else if (startcell[OCT2] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT2],ZN);
            // may be zero!
          }
	  else if (src_off_grid[XX]==XP && src_off_grid[ZZ]==ZP) {
            startcell[oct] = sc;
          }
	  // Extra condition: if 2 and 5 are both zero, then whole
          // plane must be (1,2,3,4)
	  else if (!startcell[OCT2] && !startcell[OCT5]) {
            startcell[oct] = 0; // if both zero, it must be zero.
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT6",oct);
          }
	}
	else if (oct==OCT8) { // XP,YN,ZN
	  if      (startcell[OCT5] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT5],YN);
            // may be zero!
          }
	  else if (startcell[OCT4] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT4],ZN);
            // may be zero!
          }
	  else if (src_off_grid[YY]==YP && src_off_grid[ZZ]==ZP) {
            startcell[oct] = sc;
          }
	  // Extra condition: if 4 and 5 are both zero, then whole
          // plane must be (1,2,3,4)
	  else if (!startcell[OCT4] && !startcell[OCT5]) {
            startcell[oct] = 0; // if both zero, it must be zero.
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT8",oct);
          }
	}
	else if (oct==OCT7) {
	  if      (startcell[OCT6] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT6],YN);
            // may be zero!
          }
	  else if (startcell[OCT8] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT8],XN);
            // may be zero!
          }
	  else if (startcell[OCT3] !=0) {
            startcell[oct] = gridptr->NextPt(startcell[OCT3],ZN);
            // may be zero!
          }
	  else if (src_off_grid[XX]==XP && src_off_grid[YY]==YP
                   && src_off_grid[ZZ]==ZP) {
            startcell[oct] = sc;
          }
	  else if (!startcell[OCT3] && !startcell[OCT8]) {
            // Extra conditions:
            // if 3 and 8 are zero, then all in that plane must be
            // (3,4,7,8)
            startcell[oct] = 0;
          }
	  else if (!startcell[OCT3] && !startcell[OCT6]) {
            // if 3 and 6 are zero, then all in that plane must be
            // (2,3,6,7)
            startcell[oct] = 0;
          }
	  else {
            rep.error("RayTracing3D::RayTraceSource logic error OCT7",oct);
          }
	}
	else rep.error("RayTracing3D::RayTraceSource bad octant",oct);
	
      } // octant may have cells, so set the starting cell.
#ifdef RT_TESTING
      cout <<"octant "<<oct<<", start_cell:"; CI.print_cell(startcell[oct]);
#endif //RT_TESTING
    } // loop over octants
  } // 3D
  else {
    rep.error("bad ndim in raytracer_USC::set_startcells()",ndim);
  }
  return;
}



// ##################################################################
// ##################################################################



cell * raytracer_USC::find_source_cell(
      double *pos ///< position of source.
      )
{
#ifdef RT_TESTING
  cout <<"find source cell N-Dim algorithm. ndim="<<ndim<<"\n";
#endif
  cell *sc=gridptr->FirstPt();
  for (int i=0;i<ndim;i++) {
    enum axes           a = static_cast<axes>     (i);
    enum direction posdir = static_cast<direction>(2*static_cast<int>(a)+1);
    enum direction negdir = gridptr->OppDir(posdir);

    // First move the source to a cell vertex.
    centre_source_on_cell(pos,a);

#ifdef RT_TESTING
    cout <<"have centred source on vertex, move on grid to it.\n";
#endif
    if      (pos[a]<gridptr->Xmin(a) ||
             pconst.equalD(pos[a],gridptr->Xmin(a))) {
      // If source is off grid in negative direction, then set
      // source cell to be the 1st grid cell in this dir.
#ifdef RT_TESTING
      cout <<"don't need to do anything as we are already at";
      cout <<" most negative cell in this axis.\n";
#endif
    }
    else if (pos[a]>gridptr->Xmax(a) ||
             pconst.equalD(pos[a],gridptr->Xmax(a))) {
      // If source is off grid in positive direction, then set
      // source cell to be the last grid cell in this dir.
#ifdef RT_TESTING
      cout <<"source off/at the positive end of grid, so go to";
      cout <<" last cell.\n";
#endif
      // source off the positive end of grid, so go to last cell.
      while (gridptr->NextPt(sc,posdir)->isgd) {
	sc=gridptr->NextPt(sc,posdir);
      }
    }
    else {
#ifdef RT_TESTING
      cout <<"source is on grid, so go in posdir to find it.\n";
#endif
      // source is on grid, so go in posdir until we find it.
      sc = find_src_on_grid(pos, sc, a);
    }
  } // loop over ndim axes.
    
#ifdef RT_TESTING
  cout <<"raytracer_USC::find_source_cell() returning.\n";
#endif
  return sc;
}



// ##################################################################
// ##################################################################



void raytracer_USC::centre_source_on_cell(
      double *pos,   ///< position of source (size ndim).
      enum axes axis ///< axis to find source along.
      )
{
  //
  // CORNER-CENTRED SOURCES: This function moves the source to a cell
  // corner, if it's not already there.  If the source is within a
  // cell it will get moved to the nearest cell boundary along the
  // current axis.
  //
  //bool changed_pos=false;

  double dist = pos[axis]-gridptr->Xmin(axis);
#ifdef RT_TESTING
  cout <<"axis:"<<axis<<"\tdist="<<dist<<" pos="<<pos[axis];
  cout <<" xmin="<<gridptr->Xmin(axis)<<"  ";
#endif
  double dx=gridptr->DX();
  int x = static_cast<int>(dist/dx);
  dist -= x*dx;    dist /= dx;
#ifdef RT_TESTING
  cout <<", new dist="<<dist<<" in units of cell size."<<"\n";
#endif
  if (!pconst.equalD(dist,0.0)) {
    if (fabs(dist)>1.0) {
      rep.error("centre_source_on_cell() wrong distance",dist);
    }
    //
    // Source is not at a cell boundary, so find the nearest one.  if
    // |dist|=0.5, then move it to the negative cell boundary, and do
    // an explicit check to make sure this happens.
    //
    if (pconst.equalD(dist,0.5)) {
      //cout <<"resetting position. dist="<<dist<<", dx="<<dx;
      //cout <<"  dist*dx="<<dist*dx<<"\n";
      if (dist>0)                    pos[axis] -=      dist *dx;
      else                           pos[axis] -= (1.0+dist)*dx;
    }
    else if (dist>0.0 && dist>  0.5) pos[axis] += (1.0-dist)*dx;
    else if (dist>0.0 && dist<= 0.5) pos[axis] -=      dist *dx;
    // 1+dist is positive, so we subtract to nearest edge.
    else if (dist<0.0 && dist<=-0.5) pos[axis] -= (1.0+dist)*dx;
    // dist is negative, so we add to position to get to edge
    else if (dist<0.0 && dist> -0.5) pos[axis] -=      dist *dx;
    else rep.error("logic failed in centre_source_on_cell()",dist);
#ifdef RT_TESTING
    cout << "New Source Location: x = "<<pos[axis]<<"\n";
#endif
    //changed_pos=true;
  }

  //
  // We don't change the global position here, but we do check back
  // in the Add_source() function if the position is changed, and
  // then propagate the change to the global data.
  //

#ifdef RT_TESTING
  cout <<"raytracer_USC::centre_source_on_cell() returning.\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



cell * raytracer_USC::find_src_on_grid(
      double *pos,   ///< position of source.
      cell *sc,      ///< starting cell.
      enum axes axis ///< axis to find source along.
      )
{
#ifdef RT_TESTING
  cout <<"finding source on grid!\n";
#endif
  if (!sc) rep.error("No starting cell for find_src_on_grid()",sc);
  enum direction posdir; //,negdir;
  posdir = static_cast<direction>(2*axis+1);
  //negdir = static_cast<direction>(2*axis);
  double dist = pos[axis]-CI.get_dpos(sc,axis);
  //    cout <<"dist="<<dist<<" and dx/2 = "<<halfdx<<"\n";

  //
  // Assume we always start at a X/Y/ZN boundary, so pos should be
  // greater than sc->pos.
  //
  if (dist <0.0) rep.error("Logic error in find_src_on_grid()",axis);

  while (dist>0.0 && gridptr->NextPt(sc,posdir)!=0) {
    sc=gridptr->NextPt(sc,posdir);
    dist = pos[axis]-CI.get_dpos(sc,axis);
  }
#ifdef RT_TESTING
  cout <<"dist="<<dist<<", finding source on grid! done!\n";
#endif

  dist = fabs(pos[axis]-CI.get_dpos(sc,axis));
  return sc;
}
 


// ##################################################################
// ##################################################################



int raytracer_USC::trace_column(
      const rad_source *source, ///< source we are tracing from.
      cell *c,                  ///< cell to start from.
      const enum direction dir  ///< direction we are looking.
      )
{
  double ds=0.0, Nc[MAX_TAU];
  for (unsigned short int iT=0; iT<source->s->NTau; iT++)
    Nc[iT] =0.0;

  int    err=0;
  
  if ( c!=0 ) {
    // need to check in case we have moved from source cell off the grid.
#ifdef RT_TESTING
    cout <<"raytracer_USC::trace_column() running.\n";
#endif 
    do {
#ifdef TESTING
      dp.c = c;
#endif
      err += get_cell_columns(source, c, Nc, &ds);
      err += ProcessCell(c,Nc,ds,source,delt);
    } while ( (c=gridptr->NextPt(c,dir))!=0 );
  }
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC::trace_plane(
      const rad_source *source,
      cell *cy, 
      const enum direction xdir,
      const enum direction ydir
      )
{
#ifdef RT_TESTING
      cout <<"tracing 2d plane, xdir="<<xdir<<", ydir="<<ydir;
      cout <<", cy="<<cy<<", id="<<cy->id<<"\n";
      //CI.print_cell(cy);
#endif 
  int err = 0;
  cell *cx = 0;
  if (cy!=0) {
    do {
#ifdef RT_TESTING
      cout <<"new column in 2d. cy="<<cy<<", id="<<cy->id<<"\n";
      //CI.print_cell(cy);
#endif 
      cx = cy;
      err += trace_column(source,cx,xdir);
    } while ((cy=gridptr->NextPt(cy,ydir))!=0);
  }
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC::trace_octant(
      const rad_source *source,  ///< source we are dealing with.
      cell *cz,                  ///< starting cell in octant
      const enum direction xdir, ///< x-direction from starting cell to go in.
      const enum direction ydir, ///< y-direction from starting cell to go in.
      const enum direction zdir  ///< z-direction from starting cell to go in.
      )
{
  int err = 0;
  cell *cy = 0;
  if (cz!=0) {
    do {
      //      cout <<"new plane in 3d.\n";
      cy = cz;
      err += trace_plane(source,cy,xdir,ydir);
    } while (gridptr->NextPt(cz,zdir)!=0);
  }
  return err;
}



// ##################################################################
// ##################################################################



int raytracer_USC::get_cell_columns(
      const rad_source *s,
      cell *c,
      double *Nc,
      double *ds
      )
{
  int err=0;
  if (ndim==3) {
    err += cell_cols_3d(s,c,Nc,ds);
  }
  else if (ndim==2) {
    err += cell_cols_2d(s,c,Nc,ds);
  }
  else if (ndim==1) {
    // SrcDir is set to point from the cell back to the source.
    err += cell_cols_1d(s,c,SrcDir[XX],Nc,ds);
  }
  else rep.error("bad ndim in get_cell_columns()",ndim);
  return err;
}



// ##################################################################
// ##################################################################


int raytracer_USC::cell_cols_2d(
      const rad_source *src,
      cell *c,
      double Nc[],
      double *ds
      )
{
  //  cout <<"raytracer_USC::cell_cols_2d() start\n";
  double delta=0.0;
  enum direction entryface=NO, perpface=NO;
  int diffx = abs(CI.get_ipos(c,XX) - src->ipos[XX]);
  int diffy = abs(CI.get_ipos(c,YY) - src->ipos[YY]);
  int mindiff=0, maxdiff=0;
  //cout <<"diffx="<<diffx<<" and diffy="<<diffy<<"\n";
  //cout <<"srcdir: xx="<<SrcDir[XX]<<" yy="<<SrcDir[YY]<<"\n";
  *ds = 0.0;
#ifdef RT_TESTING
  cout <<"\tGetting cell-cols 2D for cell id="<<c->id;
  cout <<": diffx="<<diffx<<" and diffy="<<diffy<<".  ";
  cout <<"srcdir: xx="<<SrcDir[XX]<<" yy="<<SrcDir[YY]<<".  ";
#endif // RT_TESTING

  //
  // See which side the cell enters through, and set delta accordingly.
  //
  if (diffx>=diffy) {
    entryface = SrcDir[XX];
    perpface  = SrcDir[YY];
    delta = static_cast<double>(diffy)/static_cast<double>(diffx);
    mindiff=diffy;
    maxdiff=diffx;
  }
  else {
    entryface = SrcDir[YY];
    perpface  = SrcDir[XX];
    delta = static_cast<double>(diffx)/static_cast<double>(diffy);
    mindiff=diffx;
    maxdiff=diffy;
  }

  //
  // ds is in physical units.
  //
  *ds = gridptr->DX() *sqrt(1.0+ delta*delta);
#ifdef RT_TESTING
  cout <<"ds="<<*ds<<" delta="<<delta<<"\n";
#endif // RT_TESTING
  double idx = gridptr->idx(); // cell size in integer units.
  double idxo2 = 0.5*idx;
  //double dist = sqrt(static_cast<double>(diffx*diffx+diffy*diffy));
  
  //
  // if the source is within the column (i.e. within x+-dx/2) then do a 1D column.
  // NOTE: for integer cell positions, the cell size is "idx".
  //
  if      (mindiff<idx && maxdiff<idx) {
#ifdef RT_TESTING
    cout <<"At source cell, setting col2cell to zero.\n";
#endif // RT_TESTING
    // at source cell, so column to cell is zero.
    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      Nc[iT] = 0.0;
    }
    //CI.set_col(c, src->s->id, Nc);
  }
  else if (mindiff<idx) {
#ifdef RT_TESTING
    cout <<"Within the source column of cells! mindiff="<<mindiff<<"\n";
#endif // RT_TESTING
    //
    // Not at source cell, but are within column, so instead of an averaging we
    // can just take the distance from source to cell, taking care to change the
    // distance because the angle to the source has changed.
    //
    cell *ngb = gridptr->NextPt(c,entryface);
    //
    // assume if neighbour doesn't exist, that the source is coming
    // in from off grid to cell c.
    //
    if (ngb) {
      CI.get_col(ngb, src->s->id, Nc);
    }
    else {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
        Nc[iT] = 0.0;
      }
    }

#ifdef RT_TESTING
    if (Nc[0] < 0.0) {
      cout <<"column is negative:"<<Nc[0]<<" coming from off grid???\n";
      CI.print_cell(c);
      CI.print_cell(ngb);
      if (gridptr->NextPt(c,XP))
        CI.print_cell(gridptr->NextPt(c,XP));
      Nc[0]=0.0;
      rep.error("Got negative column from a cell when we shouldn't have!",Nc[0]);
    }
#endif // RT_TESTING

    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      Nc[iT] = std::max(Nc[iT], 0.0);
    }
    
    //
    // Now we have the column, need to scale it by the changed
    // distance due to the different angle from source to cell
    // centre.
    //
    if (maxdiff<10*idx) {
      double maxd = static_cast<double>(maxdiff);
      double maxmin2 = maxd - idx;
      for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
        Nc[iT] *= sqrt((maxd*maxd+idxo2*idxo2)/
                       (maxmin2*maxmin2+idxo2*idxo2))*
                       maxmin2/maxd;
      }
#ifdef RT_TESTING
      cout <<"cell_cols_2d: close, Nc="<<Nc[0]<<"\n";
#endif
    }

  } // mindiff<2
  else {
    //
    // source is not in 1D column, so do the 2D averaging.
    //
    col2cell_2d(src, c, entryface, &perpface, &delta, Nc);
#ifdef RT_TESTING
    cout <<"cell_cols_2d: normal, Nc="<<Nc[0]<<"\n";
#endif
  }

#ifdef RT_TESTING
  if (maxdiff<idx) {
    cout <<"\t...\tcol = "<<*Nc<<"\t";
    cout <<"cell id: "<<c->id<<"  and col in cell = "<<*ds<<"\n";
    cout <<"raytracer_USC::cell_cols_2d() done\n";
  }
  //CI.print_cell(c);
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int raytracer_USC::cell_cols_3d(
      const rad_source *src,
      cell *c,
      double Nc[],
      double *ds
      )
{

  int dx[3]; // relative position vector.
  dx[XX] = abs(CI.get_ipos(c,XX) - src->ipos[XX]);
  dx[YY] = abs(CI.get_ipos(c,YY) - src->ipos[YY]);
  dx[ZZ] = abs(CI.get_ipos(c,ZZ) - src->ipos[ZZ]);

  int o[3]; // ordering vector.
  o[XX]=XX; o[YY]=YY; o[ZZ]=ZZ;

  //
  // sort the dx[] elements in decreasing order, using o[] to track which is which.
  //
  if (dx[o[2]]> dx[o[1]]) std::swap(o[2],o[1]);
  if (dx[o[0]]<=dx[o[1]]) std::swap(o[0],o[1]);
  if (dx[o[1]]<=dx[o[2]]) std::swap(o[1],o[2]);
  //   rep.printVec("dx",dx,ndim);
  //   rep.printVec("oo",o ,ndim);

  //
  // now dx[o[0]] >= dx[o[1]] >= dx[0[2]]
  // ray enters through SrcDir[o[0]], and points more in SrcDir[o[1]] than SrcDir[o[2]]
  //
  enum direction entryface, perpdirs[2];
  double deltas[2];
  entryface = SrcDir[o[0]];
  perpdirs[0] = SrcDir[o[1]];
  perpdirs[1] = SrcDir[o[2]];

  deltas[0] = static_cast<double>(dx[o[1]])/static_cast<double>(dx[o[0]]);
  deltas[1] = static_cast<double>(dx[o[2]])/static_cast<double>(dx[o[0]]);
  //
  // ds is in physical units
  //
  *ds = gridptr->DX()* sqrt(1.+ deltas[0]*deltas[0]+deltas[1]*deltas[1]);
  double idx = gridptr->idx(); // cell size in integer units.
  double idxo2 = 0.5*idx;



#ifdef RT_TESTING
  cout <<"******* Getting col2cell:\n"; //CI.print_cell(c);
#endif // RT_TESTING
  //
  // Need to do something more careful if dx[o[1]] and/or dx[o[2]]
  // are <2 since we don't want to take an average from a neigbour
  // which shouldn't contribute at all.
  //
  if (dx[o[0]]<idx) {
    //
    // We are at the source cell, since max-dist=1, so tau=0.
    //
#ifdef RT_TESTING
    cout <<"At the source Cell, so setting col2cell=0.0\n";
#endif // RT_TESTING
    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      Nc[iT]=0.0;
    }
  }
  else if (dx[o[1]]<idx) {
    //
    // We're within a cell distance in 2 of 3 directions, so we don't
    // need to do any averaging, just the geometric change.
    //
#ifdef RT_TESTING
    cout <<"In source column, so doing 1D column calculation.";
#endif // RT_TESTING
    cell *ngb = gridptr->NextPt(c,entryface);
    //
    // assume if neighbour doesn't exist, that the source is coming
    // in from off grid to cell c.
    //
    if (ngb) {
      CI.get_col(ngb, src->s->id, Nc);
    }
    else {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
        Nc[iT] = 0.0;
      }
    }

#ifdef RT_TESTING
    if (Nc[0] < 0.0) {
      if (ngb == src->sc) {
        cout <<"source cell has negative column density, resetting.\n";
        Nc[0]=0.0;
      }
      else if (pconst.equalD(Nc[0],0.0)) {
        cout <<"column is negative:"<<Nc[0]<<" but close to zero.";
        cout <<" ... resetting to zero.\n";
        Nc[0]=0;
      }
      else {
        cout <<"column is negative:"<<Nc[0]<<" coming from off grid?\n";
        CI.print_cell(c);
        CI.print_cell(ngb);
        if (gridptr->NextPt(c,XP))
          CI.print_cell(gridptr->NextPt(c,XP));
        Nc[0]=0.0;
        rep.error("wrongly got negative column from cell!",Nc[0]);
      }
    }
#endif // RT_TESTING
    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      Nc[iT] = std::max(Nc[iT], 0.0);
    }

    //
    // Now we have the column, need to scale it by the changed
    // distance due to the different angle from source to cell
    // centre.  Only do this near the source.
    //
    if (dx[o[0]]<15.0*idxo2) {
      double max  = dx[o[0]];
      double max2 = max - idx;
      //if (max<1.499999*idx) rep.error("maxdiff should be >=1.5*idx",max);
#ifdef RT_TESTING
      cout <<" tau="<<*Nc<<" scale factor=";
      cout <<sqrt((max*max+idxo2*idxo2)/(max2*max2+idxo2*idxo2))*max2/max<<"\n";
#endif // RT_TESTING
      for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
        Nc[iT] *= sqrt((max*max+idxo2*idxo2)/(max2*max2+idxo2*idxo2))*
                  max2/max;
      }
    }
  }
  else if (dx[o[2]]<idx) {
    //
    // Now only one distance is <idx, so we are in a plane with the
    // source bordering it.  So we need a 2D column calculation,
    // scaled by the distance from the source including the offset in
    // the third direction.  i.e. the intersection with the plane
    // dx2=1 happens closer to the source in the neighbouring cells,
    // so we need to make the path length shorter, at least for cells
    // close to the source.
    //
#ifdef RT_TESTING
    cout <<"In source plane; doing 2D average.";
#endif // RT_TESTING
    col2cell_2d(src, c, entryface, perpdirs, deltas, Nc);
#ifdef RT_TESTING
    cout <<" tau="<<Nc[0]<<"\n";
#endif // RT_TESTING
    //
    // For cells with dx0>=5*idx, the scaling factor is >= 0.9974, so we
    // won't bother scaling in this case.  For the cell with
    // dx=(3,3,1), the scaling is 0.83887, and for the other cells we
    // can do a approximation which for the (5,5,1) cell gives 0.982
    // instead of 0.983, so it's pretty good. For (5,3,1) --> 0.9755
    // vs 0.9733.
    //
    // The ratio is, for cell size 2 units, and assuming dx2=1, and r^2=dx0^2+dx1^2, 
    // and a=(dx0-2)/dx0:
    // ratio = a*sqrt([r^2+1]/[(ar)^2+1]) \simeq (1+1/(2r^2))(1-1/(2a^2r^2)) for ar>>1.
    //
    if (dx[o[0]]<5*idx) {
      if (dx[o[0]]==3*idxo2) {
        for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
          Nc[iT] *= 0.8388704928; // hard-coded value for sqrt((1+1/18)/(1+1/2))
        }
      }
      else {
        double one_over_r2 = 1.0/(dx[o[0]]*dx[o[0]] +dx[o[1]]*dx[o[1]]);
        double one_over_a2r2 = static_cast<double>(dx[o[0]]*dx[o[0]])
                        /((dx[o[0]]-idx)*(dx[o[0]]-idx))*one_over_r2;
#ifdef RT_TESTING
        cout <<"source plane averaging... 1/r2="<<one_over_r2;
        cout <<" 1/(ar)2="<<one_over_a2r2;
        cout <<" Scaling factor="<<(1.0+one_over_r2)*(1.0-one_over_a2r2)<<"\n";
#endif // RT_TESTING
        for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
          Nc[iT] *= (1.0+one_over_r2)*(1.0-one_over_a2r2);
        }
      }
    }
  } // if dx2<2 (so in source plane)
  else {
    //
    // Now we aren't in the source plane, so do the full 3D column calculation.
    //
#ifdef RT_TESTING
    cout <<"Not in source plane, so usual average.\n";
#endif // RT_TESTING
    col2cell_3d(src, c, entryface, perpdirs, deltas, Nc);
  }

  //   cout <<"RayTracing3D::GetCellColumns()\tcell id: "<<c->id<<"\t: col to cell = "<<*Nc<<"  ";
  //   cout <<"\tand ds="<<*ds<<"\n";   
  //   cout <<"\traytracer_USC::cell_cols_3d() done\n";
  return 0;
}



// ##################################################################
// ##################################################################



void raytracer_USC::col2cell_2d(
      const rad_source *src,        ///< source we are working on.
      const cell *c,                  ///< cell to get column to.
      const enum direction entryface, ///< face ray enters cell through.
      const enum direction *perpdir,///< array of perp dirs tosource (only 1 el)
      const double *delta,  ///< array of tan(theta) (1 el) (angle in [0,45]deg)
      double *Nc ///< Column densities.
      )
{
  double col1[MAX_TAU], col2[MAX_TAU];
  cell *c1 = gridptr->NextPt(c,  entryface);
  if (!c1) {
    for (short unsigned int iT=0; iT<src->s->NTau; iT++)
      col1[iT] = col2[iT] = 0.0;
  }
  else {
    cell *c2 = gridptr->NextPt(c1, (*perpdir) );
    if (!c2) {
      CI.get_col(c1, src->s->id, col1);
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col2[iT] = 0.0;
    }
    else {
      CI.get_col(c1, src->s->id, col1);
      CI.get_col(c2, src->s->id, col2);
    }
  }
  
  //
  // INTERPOLATION SCHEMES -- BASICALLY ALL WERE CRAP EXCEPT C2RAY...
  // 
  //
  for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
    Nc[iT] = interpolate_2D(src->s->id, *delta, col1[iT], col2[iT]);
  }
  return;
}


// ##################################################################
// ##################################################################



void raytracer_USC::col2cell_3d(
      const rad_source *src,          ///< source we are working on.
      const cell *c,                  ///< cell to get column to.
      const enum direction entryface, ///< face ray enters cell through.
      const enum direction *perpdir,  ///< array of perp dirs to source (1 el)
      const double *dx,               ///< array of tan(theta) (angle in [0,45]deg)
      double *Nc ///< Column densities.
      )
{
  // Algorithm is the same as that describe in Mellema et al.,2006, 
  // NewA, 11,374, appendix A.  Good for 3D cartesian geometry.
  //  cout <<"3D ShortChars:: entrydir = "<<entryface<<" and perps = ["<<perpdirs[0]<<", "<<perpdirs[1]<<"]\n";
  cell *c1=0, *c2=0, *c3=0, *c4=0;
  double col1[MAX_TAU], col2[MAX_TAU], col3[MAX_TAU], col4[MAX_TAU];

  c1 = gridptr->NextPt(c,  entryface);
  if (!c1) {
    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      col1[iT] = col2[iT] = col3[iT] = col4[iT] = 0.0;
    }
  }
  else {
    CI.get_col(c1,src->s->id, col1);

    c2 = gridptr->NextPt(c1,  perpdir[0]);
    if (!c2) {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col2[iT] = 0.0;
    }
    else {
      CI.get_col(c2,src->s->id, col2);
    }
    
    c3 = gridptr->NextPt(c1,  perpdir[1]);
    if (!c3) {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col3[iT] = 0.0;
    }
    else {
      CI.get_col(c3, src->s->id, col3);
    }
    
    if (c2 && c3) {
      c4 = gridptr->NextPt(c2, perpdir[1]);
      if (!c4)
        rep.error("lost on grid -- corner cell doesn't exist",c4);
      CI.get_col(c4, src->s->id, col4);
    }
    else {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col4[iT] = 0.0;
    }

  }
  //  cout <<"3D ShortChars:: col1="<<col1<<" col2="<<col2<<" col3="<<col3<<" col4="<<col4;
  //  cout <<"\t dx = ["<<dx[0]<<", "<<dx[1]<<"]"<<"\n";

  //
  //  0: C2Ray inverse Tau with minTau=0.7: (see Mellema et al.,2006, NewA, 11,374, eq.A.5)
  //
  for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
    Nc[iT] = interpolate_3D(src->s->id, dx[0], dx[1], 
                            col1[iT], col2[iT], col3[iT], col4[iT]);
  }
  return;

}


// ##################################################################
// ##################################################################




///
/// Apply the appropriate weighting scheme to get an interpolated optical depth for 2D
///
double raytracer_USC::interpolate_2D(
      const int src_id, ///< source id
      const double delta0, ///< delta = min(abs(dy/dx),abs(dx/dy));
      const double tau1, ///< first optical depth tau_1
      const double tau2  ///< second optical depth tau_2
      )
{
  //
  // This is the standard weighting used in C2Ray.
  // For other weighting schemes see raytracer_USC::col2cell_2d()
  //
  double w1,w2,mintau2d;
  mintau2d=TauMin[src_id];
  w1 = (1.-delta0)/max(mintau2d,tau1);
  w2 =     delta0 /max(mintau2d,tau2);
  mintau2d = (w1+w2);
  w1 /= mintau2d; w2 /= mintau2d;
  return w1*tau1 + w2*tau2;
}



// ##################################################################
// ##################################################################



///
/// Apply the appropriate weighting scheme to get an interpolated
/// optical depth for 3D
///
double raytracer_USC::interpolate_3D(
      const int src_id, ///< source id
      const double delta0, ///< delta0 = abs(dy/dx)
      const double delta1, ///< delta1 = abs(dz/dx)
      const double tau1, ///< first optical depth tau_1
      const double tau2, ///< second optical depth tau_2
      const double tau3, ///< third optical depth tau_3
      const double tau4 ///< fourth optical depth tau_4
      )
{
  //
  // This is the standard weighting used in C2Ray.
  // (see Mellema et al.,2006, NewA, 11,374, eq.A.5)
  //
  double w1,w2,w3,w4,mintau3d;
  mintau3d=TauMin[src_id];
  w1 = (1.-delta0)*(1.-delta1)/max(mintau3d,tau1);
  w2 =     delta0 *(1.-delta1)/max(mintau3d,tau2);
  w3 = (1.-delta0)*    delta1 /max(mintau3d,tau3);
  w4 =     delta0 *    delta1 /max(mintau3d,tau4);
  mintau3d = (w1+w2+w3+w4);
  w1/=mintau3d; w2/=mintau3d; w3/=mintau3d; w4/=mintau3d;
  return w1*tau1 + w2*tau2 + w3*tau3 + w4*tau4;

}


// ##################################################################
// ##################################################################



void raytracer_USC::set_Vshell_in_cell(
      cell *c, ///< current cell.
      double ds,           ///< Path Length through cell.
      const rad_source *source ///< pointer to source struct.
      )
{
  //
  // First set ds through the cell.
  //
  CI.set_cell_deltaS(c, source->s->id, ds);

  //
  // If the source is at infinity then we just set Vshell=ds=delta-x
  //
  if (source->s->at_infinity) {
#ifdef RT_TESTING
    cout <<"raytracer_USC::set_Vshell_in_cell() src at infinity!\n";
#endif // RT_TESTING
    CI.set_cell_Vshell(c,source->s->id,ds);
    return;
  }
  
  //
  // If source is at a finite distance, then we need the volume
  // of the shell bounded by the entry and exit points of the 
  // ray in the cell.
  //
  double Vshell=0.0;
  double c_pos[ndim]; // cell position vector.
  double rs; // distance from source to entry point of cell.
  
  //
  // Get the distance to the geometric cell-centre, *NOT* the
  // centre of volume since that makes the short characteristics
  // method horrendously complicated!!  The method is only 1st 
  // order anyway, and the position offset is 2nd order, so it
  // shouldn't make any difference.
  // Distance is from source to entry point of ray to cell.
  //
  CI.get_dpos(c,c_pos);
  rs = gridptr->distance(source->s->pos,c_pos) -0.5*ds;
#ifdef RT_TESTING
  cout <<"\tSetting Vshell for cell id="<<c->id<<": ";
  cout <<"rs="<<rs<<" ds="<<ds<<" dx="<<gridptr->DX()<<": ";
#endif // RT_TESTING

  //
  // 3D shell, volume is 4/3.Pi.((r+dr)^3-r^3) 
  //
  Vshell = 4.0*M_PI*((rs+ds)*(rs+ds)*(rs+ds) -rs*rs*rs)/3.0;

#ifdef RT_TESTING
  cout <<"Vshell="<<Vshell<<"\n";
#endif // RT_TESTING

  //
  // Now set Vshell in the cell-data.
  //
  CI.set_cell_Vshell(c, source->s->id, Vshell);
  return;
}



// ##################################################################
// ##################################################################



