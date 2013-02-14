///
/// \file microphysics_lowZ.cc
/// \author 
/// \date 2010.10.12
///
/// Microphysics class for low metallicity gas at high redshift.
///
/// Modifications:
///
/// - 2010.10.12 JM: Created template class.
///
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.03.21 JM: updated to work with Harpreet's code, and diffuse RT.
/// - 2011.11.09 JM: Added Set_Temp() function for correcting negative pressure
///                  in dynamics solver.
/// - 2012.01.19 JM: made timescales() function work.  Fixed bug in Temperature() function.
/// - 2012.03.12 JM: Modified interface to handle point source column
///     densities.  It is non-photon-conserving in the formulation (just uses
///     distance to source and not the shell volume from Mellema+2006), but
///     should be good enough for Harpreet.  Also updated timestep calculator
///     to deal with multiple column densities.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.

//
// First include all the defines, to make sure that we need to compile in 
// Harpreet's module to the final exe.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifndef EXCLUDE_HD_MODULE




#include "microphysics_lowZ.h"
using namespace std;




// ##################################################################
// ##################################################################



microphysics_lowz::microphysics_lowz(const int nv,
				     const int ntracer,
				     const std::string &trtype,
				     struct which_physics *ephys
				     )
  :
  kB(GS.kB()),
  m_p(GS.m_p()),
  nv_prim(nv)
{
  cout <<"Welcome to microphysics_lowz: the next generation chemistry solver!\n";
  //
  // Note this is only called once at the start of the simulation so there is 
  // no need to be efficient.
  //

  //
  // set flags for what processes we are and aren't doing.
  //
  ep.dynamics          = ephys->dynamics;
  ep.cooling           = ephys->cooling;
  ep.chemistry         = ephys->chemistry;
  ep.coll_ionisation   = ephys->coll_ionisation;
  ep.rad_recombination = ephys->rad_recombination;
  ep.phot_ionisation   = ephys->phot_ionisation;
  ep.raytracing        = ephys->raytracing;
  ep.update_erg        = ephys->update_erg;
  if (!ep.chemistry) {ep.coll_ionisation = ep.rad_recombination = ep.phot_ionisation = 0;}
  //  cout <<"\t\tExtra Physics flags set.\n";

  //
  // Set up tracer variables.  We will assume that ALL tracers are Harpreet's species, and
  // that there are no others.
  //
  cout <<"\t\tSetting up Tracer Variables.";
  microphysics_lowz::Nspecies = SimPM.ntracer;
  microphysics_lowz::Yvector_length = Nspecies+1;
  //
  // Local state vector.  We have density, internal energy,
  // shielding factor, the Y-values, and that's it.
  //
  Yi.resize(Yvector_length);
  Yf.resize(Yvector_length);
  microphysics_lowz::density = -1.0;  // mass density.
  microphysics_lowz::Eint = -1.0;     // internal energy per unit volume.
  microphysics_lowz::Column_density = -1.0;    // mass column density, int(rho*dx)
  microphysics_lowz::col_data_set=false;
  cout <<" Nspecies = "<<Nspecies<<", Vec-len="<<Yvector_length<<"\n";

//#define TEST_KVECTOR
#ifdef TEST_KVECTOR
  Column_density = 0.0;
  density=2.4e-24;
  double flux=0.0;
  double P[3];
  Y_init(Yi);
  conversion_JMcode( density, Column_density ,  P );
  calculate_X(Yi,X);
  double T_now=1.0e3;
  ofstream outf("Kvector.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"## Temperature(K) K[0] K[1] ... \n";
  outf.setf( ios_base::scientific );
  outf.precision(3);
  do {
    //calculate_T(Yi,X,P[0],T_now);
    calculate_K(K,flux,T_now,flux); 
    cout <<"current temperature is "<<T_now;
    outf << T_now<<" ";
    for (int v=0;v<56;v++) outf <<K[v]<<" ";
    outf << "\n";
    T_now *=1.1;
  } while (T_now<1.0e7);
  outf.close();
  rep.error("Bugging out",2);
#endif // TEST_KVECTOR

  //
  // All done
  //
  cout <<"microphysics_lowz: set up finished. Returning.\n"; 
  return;
}





// ##################################################################
// ##################################################################





microphysics_lowz::~microphysics_lowz()
{
  Yi.clear();
  Yf.clear();
  return;
}






// ##################################################################
// ##################################################################




double microphysics_lowz::Temperature(const double *pv, ///< primitive vector
				const double g   ///< eos gamma
				)
{
  //
  // put data from p_in into local vectors.
  //
  double ttt=0.0;
  convert_prim2local(pv,g);
  Column_density = 0.0;
  vector<double>col(5);
  double P[3];
  conversion_JMcode(//Yi,
                    density, col , P );
  calculate_X(Yi,X);
  calculate_T(Yi,X,P[0],ttt);
  return ttt;
}



// ##################################################################
// ##################################################################


///
/// Set the gas temperature to a specified value.
///
int microphysics_lowz::Set_Temp(
                double *p_in,     ///< primitive vector.
	        const double T_required, ///< temperature
	        const double g            ///< eos gamma.
	        )
{
  //
  // First calculate current temperature:
  //
  double T_now=0.0;
  convert_prim2local(p_in,g);
  Column_density = 0.0;
  vector<double>col(5);
  double P[3];
  conversion_JMcode(//Yi,
                    density, col , P );
  calculate_X(Yi,X);
  calculate_T(Yi,X,P[0],T_now);
  cout <<"current temperature is "<<T_now;
  //
  // Now multiply current internal energy by ratio of current to required
  // temperatures.  E = (g-1)*n*k*Ti*(Tf/Ti)
  //
  cout <<", corr. to E_int="<<Yi[Yvector_length-1];
  Yi[Yvector_length-1] *= T_required/T_now;
  cout <<", resetting to E_int="<<Yi[Yvector_length-1];
  cout <<" which corresponds to T="<<T_required;
  //
  // Now convert back to primitive variables (this uses Yf as the source
  // of microphysics quantities, so first copy Yi to Yf).
  //
  Yf = Yi;
  cout <<"; Yf[E]="<<Yf[Yvector_length-1]<<"\n";
  convert_local2prim(p_in, p_in, gamma);
  return 0;
}



// ##################################################################
// ##################################################################


int microphysics_lowz::convert_prim2local(
          const double *p_in,
          const double gam
          )
{
  density = p_in[RO];
  Eint    = p_in[PG]/(gam-1.0);

#ifdef MP_DEBUG
  if (p_in[PG]<=0. || !isfinite(p_in[PG])) {
    cout <<"neg.pres. input to MP: e="<<Eint<<endl;
    rep.error("Negative/infinite pressure input to RT solver!",p_in[PG]);  
  }
  if (density<=0.0) {
    rep.error("Negative density input to microphysics!",density);
  }
#endif // MP_DEBUG

  //
  // Now convert all the other variables.  Maybe just copy from pv_XXX to lv_XXX
  // for tracers.
  //
  for (int v=0; v<Nspecies; v++) {
    Yi[v] = p_in[SimPM.ftr+v];
#ifdef MP_DEBUG
    if (!isfinite(Yi[v])) {
      rep.error("INF/NAN input to microphysics Yi",v);
    }
#endif // MP_DEBUG
  }
  Yi[Yvector_length-1] = Eint;

  return 0;
}



// ##################################################################
// ##################################################################


int microphysics_lowz::convert_local2prim(
          const double *p_in,
          double *p_out,
          const double gam
          )
{
  //
  // This is so we don't forget the velocity, B-field, etc.
  //
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];
  //
  // convert output internal energy to gas pressure:
  //
  p_out[PG] = Yf[Yvector_length-1]*(gam-1.0);
#ifdef MP_DEBUG
  if (!isfinite(p_out[PG])) {
    rep.error("INF/NAN output from microphysics Eint",p_out[PG]);
  }
#endif // MP_DEBUG
 
  //
  // Now copy all the new tracer vales to p_out[]
  //
  for (int v=0; v<Nspecies; v++) {
    p_out[SimPM.ftr+v] = Yf[v];
#ifdef MP_DEBUG
    if (!isfinite(Yf[v])) {
      rep.error("INF/NAN output from microphysics Yf",v);
    }
#endif // MP_DEBUG
  }
  return 0;
}



// ##################################################################
// ##################################################################


int microphysics_lowz::TimeUpdateMP(
          const double *p_in,
          double *p_out,
          const double dt,
          const double g,
          const int sw_int,
          double *ttt
          )
{
  int err = 0;
  microphysics_lowz::gamma = g;

  //
  // put data from p_in into local vectors.
  //
  err += convert_prim2local(p_in,gamma);
  Column_density = 0.0;
  vector<double>cols(5,1.0e50);

#ifdef MP_DEBUG
  //
  // TESTING !!!
  //
  double P[3];
  conversion_JMcode( density, col , P );
  calculate_X(Yi,X);
  calculate_T(Yi,X,P[0],*ttt);
  if (*ttt<SimPM.EP.MinTemperature) {
    cout <<"Input temperature to MP is too low!  T_in="<<*ttt<<"\n";
  }
  //
  // TESTING !!!
  //
#endif // MP_DEBUG

  //
  // CHANGE THIS FUNCTION TO MATCH WHATEVER HARPREET WANTS
  //
  //Interface_with_JMs_code(density,Column_density,Yi, junk, Yf);
  // int interface_JMcode(
  //		     std::vector<double> &, // Yi, 
  //		     std::vector<double> & , //Yf, 
  //		     double, // rho , 
  //		     double, //rho_dx
  //		     double,//tf-dynamic timestep
  //		     double //flux
  //		     );
  
  err += interface_JMcode(Yi,Yf,density,cols,dt,0.0);
  // if (err) cout <<"\n\nerror\n\n";

  Eint = Yf[Yvector_length-1];
  //
  // Return gas temperature -- but we will just do the internal energy.
  // This isn't really needed for anything.
  //
  *ttt = Eint;

#ifdef MP_DEBUG
  //
  // TESTING !!!
  //
  Column_density = 0.0;
  conversion_JMcode( density, col, P );
  calculate_X(Yf,X);
  calculate_T(Yf,X,P[0],*ttt);
  if (*ttt<SimPM.EP.MinTemperature) {
    cout <<"Output temperature to MP is too low!  T_in="<<*ttt<<"\n";
  }
  //
  // TESTING !!!
  //
#endif // MP_DEBUG

  //
  // put updated state vector into p_out.
  //
  err += convert_local2prim(p_in, p_out, gamma);

  return err;
}



// ##################################################################
// ##################################################################


int microphysics_lowz::TimeUpdateMP_RTnew(
          const double *p_in,
          const int N_heating_srcs,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &heating_srcs,
          ///< list of UV-heating column densities and source properties.
          const int N_ionising_srcs,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &ionising_srcs,
          ///< list of ionising src column densities and source properties.
          // const std::vector<int> &sources, ///< (Ndiff-ION, Ndiff-UV, Npt-ION, Npt-UV)
          //const std::vector<double> &projected_density, ///< list of column densities.
          double *p_out,
          const double dt,
          const double g,
          const int sw_int,
          double *ttt    ///< final temperature (not strictly needed).
          )
{
  int err = 0;
  microphysics_lowz::gamma = g;

  //
  // put data from p_in into local vectors.
  //
  err += convert_prim2local(p_in,gamma);

  //
  // column densities:  Here we just consider UV radiation from diffuse sources,
  // so any other sources are just ignored.
  //
#ifdef MP_DEBUG
  if (heating_srcs.size() != static_cast<unsigned int>(N_heating_srcs)) {
    rep.error("Update: N_heating_srcs doesn't match vector size in Harpreet's MP integrator",
              heating_srcs.size());
  }
  if (ionising_srcs.size() != static_cast<unsigned int>(N_ionising_srcs)) {
    rep.error("Update: N_ionising_srcs doesn't match vector size in Harpreet's MP integrator",
              ionising_srcs.size());
  }
#endif // MP_DEBUG

  //
  // This array hold the column densities, organised as follows:
  //  -- the neutral H column density to point source,
  //  -- the total column density to point sourc, 
  //  -- and the H2 column density to point sourc,
  //  -- distance to point source
  //  -- minimum column density to diffuse radiation.
  // 
  vector<double> cols(5,1.0e100);
  get_column_densities(N_heating_srcs, heating_srcs,
                       N_ionising_srcs, ionising_srcs,
                       cols);
  

  //
  // CHANGE THIS FUNCTION TO MATCH WHATEVER HARPREET WANTS (last arg is UV-flux)
  //
  err = interface_JMcode(Yi,Yf,density,cols,dt,0.0);

  //
  // Return gas temperature -- but we will just do the internal energy.
  // This isn't really needed for anything.
  //
  Eint = Yf[Yvector_length-1];
  *ttt = Eint;
  //
  // put updated state vector into p_out.
  //
  err += convert_local2prim(p_in, p_out, gamma);

  return err;
}



// ##################################################################
// ##################################################################



void microphysics_lowz::get_column_densities(
          const int N_heating_srcs,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &heating_srcs,
          ///< list of UV-heating column densities and source properties.
          const int N_ionising_srcs,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &ionising_srcs,
          ///< list of ionising src column densities and source properties.
          std::vector<double> &cols
          )
{
  //
  // The "cols" array hold the column densities, organised as follows:
  //  -- the neutral H column density to point source,
  //  -- the total column density to point sourc, 
  //  -- and the H2 column density to point sourc,
  //  -- distance to point source
  //  -- minimum column density to diffuse radiation.
  // 
  
  //
  // First see if we have already setup column density data, and if not, set
  // it up.
  //
  if (!col_data_set) {
    have_pt_src=false;
    have_diff_r=false;
    ion_src_index=-1;
    uvh_src_index=-1;
    dis_src_index=-1;
    //
    // UV heating from point source; must at minimum have UV heating from a
    // point source.
    //
    for (int v=0; v<N_heating_srcs; v++) {
      if (heating_srcs[v].type==RT_SRC_SINGLE && 
          SimPM.RS.sources[heating_srcs[v].id].opacity_src==RT_OPACITY_TOTAL) {
        uvh_src_index=v;
        have_pt_src=true;
      }
      //
      // We have at most one ionising source:
      //
      if (N_ionising_srcs>1)
        rep.error("MP_LOWZ, too many ionising srcs.",N_ionising_srcs);
      if (N_ionising_srcs==1) {
        ion_src_index=0;
        if (!have_pt_src) rep.error("Need UVH src with PION src, MP_LOWZ",0);
      }
      //
      // And dissociating source.  This is packaged with UV heating sources
      // for now, but that is a hack, and I would like to update all of this at
      // some point so that the raytracer just has a point source and does all
      // of the column densities together.
      //
      for (int v=0; v<N_heating_srcs; v++) {
        if (heating_srcs[v].type==RT_SRC_SINGLE &&
            SimPM.RS.sources[heating_srcs[v].id].opacity_src==RT_OPACITY_TRACER) {
          dis_src_index=v;
        }
      }
      //
      // Now diffuse radiation:
      //
      ds_index.clear();
      for (int v=0; v<N_heating_srcs; v++) {
        if (heating_srcs[v].type == RT_SRC_DIFFUSE) {
          ds_index.push_back(v);
        }
      }
      if (ds_index.size()>=1) have_diff_r=true;
      //
      // finished, so set the bool indicator to avoid this step next time the
      // function is called.
      //
      col_data_set=true;
    }
  }

  //
  // Now actually set the column densities in the cols() array.
  // We assume the values are initialised to something large, so we don't have
  // to set them if there is no source.
  //
  if (have_pt_src) {
    if (ion_src_index>=0)
      cols[0] = ionising_srcs[ion_src_index].Column;
    if (dis_src_index>=0)
      cols[2] = heating_srcs[dis_src_index].Column;
    // UV-heating source certainly exists, so use it to set distance.
    cols[1] = heating_srcs[uvh_src_index].Column;
    //
    // Distance to source.  For R>>dR, Vshell=4.Pi.R^2.dR, and dS=dR, so R
    // is easy to obtain.  This is not exact for R~dR, but close enough.
    //
    cols[3] = sqrt(heating_srcs[uvh_src_index].Vshell/(4.0*M_PI*heating_srcs[uvh_src_index].dS));
  }
  if (have_diff_r) {
    for (size_t v=0; v<ds_index.size(); v++) {
      Column_density = min(
        Column_density,
        heating_srcs[ds_index[v]].Column +0.5*heating_srcs[ds_index[v]].DelCol
        );
    }
    cols[4] = Column_density;
  }
  
  //cout <<"microphysics_lowZ: cols = [";
  //for (int v=0;v<5;v++) cout <<cols[v]<<", ";
  //cout <<"\n";

  return;

}

// ##################################################################
// ##################################################################




int microphysics_lowz::TimeUpdate_RTsinglesrc(
          const double *p_in,   ///< Primitive Vector to be updated.
          double *p_out,        ///< Destination Vector for updated values.
          const double dt,      ///< Time Step to advance by.
          const double g,       ///< EOS gamma.
          const int sw_int,     ///< Switch for what type of integration to use.
          ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
          const double phot_in, ///< flux in per unit length along ray (F/ds or L/dV)
          const double ds,      ///< path length ds through cell.
          const double tau2cell, ///< Optical depth to entry point of ray into cell.
          double *deltau        ///< return optical depth through cell in this variable.
          )
{
  if (!ep.phot_ionisation) rep.error("RT requested, but phot_ionisation not set!",ep.phot_ionisation);

  //
  // This is more complicated, and should only be called by the ray-tracer.
  // See MP_Hydrogen to see how to code it!
  //
  rep.error("Don't call me!!!",999);
  return 1;
}




// ##################################################################
// ##################################################################




///
/// This returns the minimum timescale of the times flagged in the
/// arguments.  Time is returned in seconds.
///
double microphysics_lowz::timescales(const double *p_in,  ///< Current cell primitive vector.
				     const double gam,    ///< EOS gamma.
				     const bool f_cool,   ///< set to true if including cooling time.
				     const bool f_recomb, ///< set to true if including recombination time.
				     const bool f_photoion ///< set to true if including photo-ionsation time.
				     )
{
  double mintime = 1.0e99;
  std::vector<struct rt_source_data> temp;
  mintime = timescales_RT(p_in, 0, temp, 0, temp, gam);
  
#ifdef MP_DEBUG
  //cout <<"hd_step_size dt="<<mintime<<"\n";
  if (!isfinite(mintime))
    rep.error("mp_lowz::timescales() hd_step_size dt returned NAN/INF",mintime);
  if (mintime<1000.0)
    cout <<"small mintime! mintime="<<mintime<<"\n";
#endif // MP_DEBUG


  return mintime;
}




// ##################################################################
// ##################################################################




///
  /// This returns the minimum timescale of all microphysical processes, including
  /// reaction times for each species and the total heating/cooling time for the gas.
  /// It requires the radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.
  ///
double microphysics_lowz::timescales_RT(
          const double *p_in,  ///< Current cell primitive vector.
          const int N_heating_srcs,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &heating_srcs,
          ///< list of UV-heating column densities and source properties.
          const int N_ionising_srcs,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &ionising_srcs,
          ///< list of ionising src column densities and source properties.
          const double gam                ///< EOS gamma
          )
{
  int err = 0;
  microphysics_lowz::gamma = gam;

  //
  // put data from p_in into local vectors.
  //
  err += convert_prim2local(p_in,gamma);

  //
  // column densities:  Here we just consider UV radiation from diffuse sources,
  // so any other sources are just ignored.
  //
#ifdef MP_DEBUG
  if (heating_srcs.size() != static_cast<unsigned int>(N_heating_srcs)) {
    rep.error("Timescales: N_heating_srcs doesn't match vector size in Harpreet's MP integrator",
              heating_srcs.size());
  }
  if (ionising_srcs.size() != static_cast<unsigned int>(N_ionising_srcs)) {
    rep.error("Timescales: N_ionising_srcs doesn't match vector size in Harpreet's MP integrator",
              ionising_srcs.size());
  }
#endif // MP_DEBUG
  //
  // This array hold the column densities, organised as follows:
  //  -- the neutral H column density to point source,
  //  -- the total column density to point sourc, 
  //  -- and the H2 column density to point sourc,
  //  -- distance to point source
  //  -- minimum column density to diffuse radiation.
  // 
  vector<double> cols(5,1.0e100);
  get_column_densities(N_heating_srcs, heating_srcs,
                       N_ionising_srcs, ionising_srcs,
                       cols);

#ifdef RT_TESTING
  //cout <<"column density="<<Column_density<<"\n";
#endif // RT_TESTING
  //
  // Call Harpreet's get_timescales function for (Yi,density,Eint);
  //
  double mintime = 1.0e99;
  mintime = hd_step_size(Yi,density,cols);


  if (mintime<=0.0) {
    rep.error("Harpreet's timescales function. negative chem time.",mintime);
  }

  return mintime;
}




// ##################################################################
// ##################################################################





int microphysics_lowz::Init_ionfractions(double *p, ///< Primitive vector to be updated.
			const double g, ///< eos gamma.
			const double    ///< optional gas temperature to end up at. (-ve means use pressure)
			)
{
  //
  // We should have Nspecies tracers, so put them in a vector.
  //
  //std::vector<double> Yi;
  //for (int v=0;v<SimPM.ntracer; v++) {
  //  Yi.push_back(-1.0);
  //}
  // Yi should be the right size already.
  //
  Y_init(Yi);
  cout <<"Yi = [";
  for (int v=0;v<SimPM.ntracer; v++) {
    cout <<Yi[v]<<", ";
  }
  cout <<"]\n";

  for (int v=0;v<SimPM.ntracer; v++) {
    p[SimPM.ftr+v] = Yi[v];
  }

  return 0;
}





// ##################################################################
// ##################################################################

#endif // if not excluding Harpreet's module.


