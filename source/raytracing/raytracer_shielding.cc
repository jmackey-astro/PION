///
/// \file raytracer_shielding.cc
/// \author Jonathan Mackey
///
/// \brief Definitions for ray-tracing class which only calculates optical depth
///        along grid axes to every cell.
///
/// Modifications:\n
/// - 2011.02.17 JM: File created.
/// - 2011.03.01 JM: Updated Process-cell with new column-density calculation.
///    It is still just testing code for now.
/// - 2011.03.02 JM: Added parallelised class.
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function.
/// - 2011.04.15 JM: Minor changes.
/// - 2011.04.22 JM: Fixed bug i generated by changing how Vshell is set.
///
//#include "../defines/functionality_flags.h"
//#include "../defines/testing_flags.h"
//#include "../global.h"
#include "raytracer_shielding.h" // short-characteristics ray-tracer.
using namespace std;

#ifndef HCORR
#error "Need the new H-correction--inspired cell interface for shielding calc!"
#endif

raytracer_shielding::raytracer_shielding(
                        class GridBaseClass    *ggg, ///< Pointer to grid
                        class MicroPhysicsBase *mmm  ///< Pointer to MicroPhysics Class.
                        )
: raytracer_USC_infinity(ggg,mmm)
{
  cout <<"raytracer_shielding:: constructor reporting for duty.\n";
  
  //
  // Set have_set_Vshell to have Nsources elements, and initialise to false.
  //
  have_set_Vshell.clear();
  have_set_Vshell.resize(SimPM.RS.Nsources,false);

  return;
}

raytracer_shielding::~raytracer_shielding()
{
  cout <<"raytracer_shielding:: destructor: mission complete.\n";
}

int raytracer_shielding::RayTrace_Column_Density(
                const int s_id,  ///< Source id
                const double dt, ///< Timestep
                const double g   ///< EOS Gamma.
                )
{
  if (!have_set_Vshell[s_id]) {
    //
    // First we set Vshell if it hasn't already been set. Only needs to be done once.
    // Trace through the grid setting Vshell in each cell (for parallel rays this is
    // identically grid->DX()).
    //
    cout <<"raytracer_shielding::RayTrace_Column_Density: setting Vshell for source "<<s_id<<".\n";
    cell *c=gridptr->FirstPt();
    double dx=gridptr->DX();
    do {
      CI.set_cell_Vshell(c, s_id, dx);
      CI.set_cell_deltaS(c, s_id, dx);
    } while ( (c=gridptr->NextPt(c)) !=0);
    have_set_Vshell[s_id] = true;
  }

  //
  // Here we have redefined the Process-cell function, so we can just call the
  // RayTrace_SingleSource() function and it will do what we want it to do.
  // This serial version of the shielding class will use raytracer_USC_infinity::
  // version of the update, whereas the later derived parallelised shielding class
  // will use its own version of the update.
  //
  int err = RayTrace_SingleSource(s_id,dt,g);
  return err;
}



int raytracer_shielding::ProcessCell(
                        class cell *c,       ///< Current cell.
                        double col2cell,     ///< Column to cell [N(H) per cm2].
                        double ds,           ///< Path Length through cell (physical units!).
                        const rad_source *source, ///< pointer to source struct.
                        const double dt       ///< Timestep
                        )
{

  //-----------------------------------------
  //-- PUT THIS IN HD'S MICROPHYSICS CLASS?--
  //-----------------------------------------
  //
  // Column density is stored for rays from the source to the far side of the
  // current cell.  As such, the stored value will be: col2cell +ds*n(H).
  // Sources have a unique id, numbered from zero, which determines which
  // variable in cell->extra_data their column density is stored in.
  //
  
  //
  // Column through cell is n(H)*ds.
  // 
  // Note that diffuse radiation has no geometric focussing or dilution.  It
  // is purely absorption along a ray.  So we don't need to worry about the 
  // centre-of-volume or any of that.
  //
  double local_col = c->Ph[RO]*ds;
  //
  // We may need to multiply the projected density by a tracer variable,
  // depending on what is providing the opacity.
  //
  switch (source->opacity_src) {
    case RT_OPACITY_TOTAL:
    // don't need to do anything here.
    CI.set_col(c, source->id, col2cell+local_col);
    CI.set_cell_col(c, source->id, local_col);
    break;

    case RT_OPACITY_MINUS:
    //
    // opacity provided by (1-y_i)*rho
    //
    local_col *= (1.0-c->Ph[source->opacity_var]);
    CI.set_col(c, source->id, col2cell+local_col);
    CI.set_cell_col(c, source->id, local_col);
    break;

    case RT_OPACITY_TRACER:
    //
    // Opacity provided by y_i*rho
    //
    local_col *= c->Ph[source->opacity_var];
    CI.set_col(c, source->id, col2cell+local_col);
    CI.set_cell_col(c, source->id, local_col);
    break;

    case RT_OPACITY_VSHELL:
    //
    // this is only applied in the AddSource() function to set Vshell in 
    // each grid cell at the start of the simulation.
    //
    set_Vshell_in_cell(c, ds, source);
    CI.set_cell_col(c, source->id, ds);
    CI.set_col     (c, source->id, col2cell+ds);
    break;

    default:
    cout <<"source id:"<<source->id<<", opactity-src="<<source->opacity_src<<" str="<<source->strength<<"\n";
    rep.error("Bad opacity_src idendifier in  raytracer_shielding",source->id);
    break;
  }
  //-----------------------------------------
  //-- PUT THIS IN HD'S MICROPHYSICS CLASS?--
  //-----------------------------------------

  return 0;
}


#ifdef PARALLEL
raytracer_shielding_pllel::raytracer_shielding_pllel(
          class GridBaseClass *ggg,   ///< Pointer to grid
          class MicroPhysicsBase *mmm ///< Pointer to MicroPhysics Class.
          )
:
  raytracer_USC_infinity(ggg,mmm),
  raytracer_shielding(ggg,mmm)
{
  cout <<"raytracer_shielding_pllel:: constructor reporting for duty.\n";
  return;
}

raytracer_shielding_pllel::~raytracer_shielding_pllel()
{
  cout <<"raytracer_shielding_pllel:: destructor doing some destructing.\n";
  return;
}

int raytracer_shielding_pllel::Add_Source(const struct rad_src_info *src ///< source info.
                                         )
{
  cout <<"\n--BEGIN-----raytracer_shielding_pllel::AddSource()------------\n";
  //
  // First call serial version.  This finds the source, and centres it
  // on a cell if needed.
  //
  cout <<"\t**** PARALLEL DIFFUSE Add_Source: calling serial version.\n";
  int id = raytracer_USC_infinity::Add_Source(src);
  cout <<"\t**** PARALLEL DIFFUSE Add_Source: serial version returned with id="<<id<<"\n";
  
  //
  // Now tell the parallel grid to decide which boundaries it needs to 
  // receive data from before processing this source, and which it needs
  // to send data to after processing.
  //
  cout <<"\t**** PARALLEL DIFFUSE Add_Source: Setting up extra RT boundaries on grid.\n";
  int err = gridptr->Setup_RT_Boundaries(id);
  if (err) rep.error("Failed to setup RT Boundaries",err);

  cout <<"\t**** PARALLEL DIFFUSE Add_Source: all done..\n";
  cout <<"\n--END-----raytracer_shielding_pllel::AddSource()------------\n";
  return id;
}

int raytracer_shielding_pllel::RayTrace_SingleSource(
                      const int s_id,  ///< Source id
                      const double dt, ///< Timestep
                      const double g   ///< eos gamma.
                      )
{
#ifdef RT_TESTING
  cout <<"Running raytracer_shielding_pllel::RayTrace_SingleSource().\n";
#endif // RT_TESTING

  int err=0;
  //cout <<"RT: Starting Raytracing for source: "<<s_id<<"\n";

  string t1="totalRT", t2="waitingRT", t3="doingRT", t4="tempRT";
  double total=0.0, wait=0.0, run=0.0;

  GS.start_timer(t1);
  //
  // First Receive RT boundaries from processors nearer source.
  //
  GS.start_timer(t2);
  //GS.start_timer(t4);
  err += gridptr->Receive_RT_Boundaries(s_id);
  //cout <<"RT: waiting to receive for "<<GS.stop_timer(t4)<<" secs.\n";
  GS.pause_timer(t2);

  //
  // Now we have the boundary conditions, so call the serial Raytracer.
  //
  GS.start_timer(t3);
  //GS.start_timer(t4);
  err += raytracer_USC_infinity::RayTrace_SingleSource(s_id, dt, g);
  //cout <<"RT: Tracing over domain took "<<GS.stop_timer(t4)<<" secs.\n";
  run = GS.pause_timer(t3);

  //
  // Finally, send the new column densities to processors further from source.
  //
  GS.start_timer(t2);
  //GS.start_timer(t4);
  err += gridptr->Send_RT_Boundaries(s_id);
  //cout <<"RT: Sending boundaries/Waiting for "<<GS.stop_timer(t4)<<" secs.\n";
  wait  = GS.pause_timer(t2);
  total = GS.pause_timer(t1);

  //cout <<"Diffuse RT: step:"<<SimPM.timestep<<" Total RT time="<<total;
  //cout <<" secs; processing="<<run<<" secs; waiting="<<wait<<"\n";
  return err;
}


#endif // PARALLEL

