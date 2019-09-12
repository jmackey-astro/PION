/// \file stellar_wind_boundaries.cpp
/// \brief Class definitions for stellar_wind boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/stellar_wind_boundaries.h"
#include "sim_params.h"
#include "grid/stellar_wind_BC.h"
#include "grid/stellar_wind_angle.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################


//
// Add internal stellar wind boundaries -- these are (possibly
// time-varying) winds defined by a mass-loss-rate and a terminal
// velocity.  A region within the domain is given fixed values
// corresponding to a freely expanding wind from a
// cell-vertex-located source.
//
int stellar_wind_bc::BC_assign_STWIND(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // Check that we have an internal boundary struct, and that we have
  // a stellar wind source to set up.
  //
  if (b->dir != NO)
    rep.error("STWIND not external boundary!",b->dir);
#ifdef TESTING
  cout <<"Assigning data to STWIND boundary. Nsrc=";
  cout <<SWP.Nsources<<"\n";
#endif
  if (SWP.Nsources<1) {
    rep.error("BC_assign_STWIND() No Sources!",SWP.Nsources);
  }

  //
  // Setup reference state vector and initialise to zero.
  //
  if (b->refval) {
    rep.error("Initialised STWIND boundary refval",b->refval);
  }
  b->refval = mem.myalloc(b->refval, par.nvar);
  if (!b->data.empty()) {
    rep.error("BC_assign_STWIND: Not empty boundary data",b->itype);
  }
  for (int v=0;v<par.nvar;v++) b->refval[v]=0.0;

  //
  // New structure: we need to initialise the stellar wind class with
  // all of the wind sources in the global parameter list (this was
  // formerly done in DataIOBase::read_simulation_parameters()).
  //
  // The type of class we set up is determined first.
  // stellar_wind_evolution is derived from stellar_wind, and 
  // stellar_wind_angle is derived from stellar_wind_evolution.
  //
  int err=0, wtype=0;
  int Ns = SWP.Nsources;
  for (int isw=0; isw<Ns; isw++) {
    if (SWP.params[isw]->type ==WINDTYPE_EVOLVING) wtype=1;
  }
  for (int isw=0; isw<Ns; isw++) {
    if (SWP.params[isw]->type ==WINDTYPE_ANGLE) wtype=2;
  }

  //
  // check values of xi.  At the moment we assume it is the same for
  // all wind sources, so it is not source-dependent.
  //
  double xi=0.0;
  for (int isw=0; isw<Ns; isw++) {
    if (isw==0) xi = SWP.params[isw]->xi;
    else {
      rep.errorTest("wind xi values don't match",xi,
                    SWP.params[isw]->xi);
    }
  }

  if (Ns>0) {
    cout <<"\n----------- SETTING UP STELLAR WIND CLASS ----------\n";
    if      (wtype==0) {
      grid->Wind = new stellar_wind(par.ndim, par.nvar, par.ntracer,
            par.ftr, par.coord_sys, par.eqntype,
            par.EP.MinTemperature);
    }
    else if (wtype==1) {
      grid->Wind = new stellar_wind_evolution(par.ndim, par.nvar,
            par.ntracer, par.ftr, par.coord_sys, par.eqntype,
            par.EP.MinTemperature, par.starttime, par.finishtime);
      err=0;
    }
    else if (wtype==2) {
      cout <<"Setting up stellar_wind_angle class\n";
      grid->Wind = new stellar_wind_angle(par.ndim, par.nvar,
            par.ntracer, par.ftr, par.coord_sys, par.eqntype,
            par.EP.MinTemperature, par.starttime, par.finishtime,
            xi);
    }
  }

  //
  // Run through sources and add sources.
  //
  for (int isw=0; isw<Ns; isw++) {
    cout <<"\tBC_assign_STWIND: Adding source "<<isw<<"\n";
    if (SWP.params[isw]->type==WINDTYPE_CONSTANT) {
      //
      // This is for spherically symmetric winds that are constant
      // in time.
      //
      err = grid->Wind->add_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->Mdot,
        SWP.params[isw]->Vinf,
        SWP.params[isw]->Vrot,
        SWP.params[isw]->Tstar,
        SWP.params[isw]->Rstar,
        SWP.params[isw]->Bstar,
        SWP.params[isw]->tr
        );
    }
    else  {
      //
      // This works for spherically symmetric winds and for
      // latitude-dependent winds that evolve over time.
      //
      cout <<"Adding source "<<isw<<" with filename ";
      cout <<SWP.params[isw]->evolving_wind_file<<"\n";
      err = grid->Wind->add_evolving_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->tr,
        SWP.params[isw]->evolving_wind_file,
        SWP.params[isw]->enhance_mdot,
        SWP.params[isw]->time_offset,
        par.simtime,
        SWP.params[isw]->update_freq,
        SWP.params[isw]->t_scalefactor
        );
    }
    if (err) rep.error("Error adding wind source",isw);
  }

  //
  // loop over sources, adding cells to boundary data list in order.
  //
  for (int id=0;id<Ns;id++) {
    cout <<"\tBC_assign_STWIND: Adding cells to source ";
    cout <<id<<"\n";
    BC_assign_STWIND_add_cells2src(par,grid, id);
  }
  //
  // Now we should have set everything up, so we assign the boundary
  // cells with their boundary values.
  //
  //err += BC_update_STWIND(par.simtime, b,0,0);
  cout <<"\tFinished setting up wind parameters\n";
  cout <<"------ DONE SETTING UP STELLAR WIND CLASS ----------\n\n";
  return err;
}



// ##################################################################
// ##################################################################



int stellar_wind_bc::BC_assign_STWIND_add_cells2src(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const int id ///< source id
      )
{
  //
  // We run through each cell, and if it is within the 
  // source's radius of influence, then we add it to the lists.
  //
  int err=0;
  int ncell=0;
  double srcpos[MAX_DIM];
  double srcrad;
  grid->Wind->get_src_posn(id,srcpos);
  grid->Wind->get_src_drad(id,&srcrad);

#ifdef TESTING
  cout <<"*** srcrad="<<srcrad<<"\n";
  rep.printVec("src", srcpos, par.ndim);
#endif

  cell *c = grid->FirstPt_All();
  do {
#ifdef TESTING
    cout <<"cell: "<<grid->distance_vertex2cell(srcpos,c)<<"\n";
#endif
    if (grid->distance_vertex2cell(srcpos,c) <= srcrad) {
      ncell++;
      err += grid->Wind->add_cell(grid, id,c);
    }
  } while ((c=grid->NextPt_All(c))!=0);
  
  err += grid->Wind->set_num_cells(id,ncell);

//#ifdef TESTING
  cout <<"BC_assign_STWIND_add_cells2src: Added "<<ncell;
  cout <<" cells to wind boundary for WS "<<id<<"\n";
//#endif
  return err;
}



// ##################################################################
// ##################################################################



//
// Update internal stellar wind boundaries -- these are (possibly
// time-varying) winds defined by a mass-loss-rate and a terminal
// velocity.  If fixed in time the wind is updated with b->refval,
// otherwise with a (slower) call to the  stellar wind class SW
//
int stellar_wind_bc::BC_update_STWIND(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      boundary_data *b, ///< Boundary to update.
      const int ,  ///< current fractional step being taken.
      const int    ///< final step (not needed b/c fixed BC).
      )
{
  //
  // The stellar_wind class already has a list of cells to update
  // for each source, together with pre-calculated state vectors,
  // so we just call the set_cell_values() function.
  //
#ifdef TESTING
  cout <<"stellar_wind_bc: updating wind boundary\n";
#endif
  int err=0;
  for (int id=0;id<grid->Wind->Nsources();id++) {
#ifdef TESTING
    cout <<"stellar_wind_bc: updating wind boundary for id="<<id<<"\n";
#endif
    err += grid->Wind->set_cell_values(grid, id,simtime);
  }

  return err;
}



// ##################################################################
// ##################################################################





