///
/// \file MPv10.cpp
/// \author Maggie Celeste Goulden
/// \date 2018.10
///
/// Description:
/// - multi-species ionization/recombination non-equilibrium
///   chemistry solver.
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from 
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.
///
/// Modifications:
/// - 2018.10.09 JM: edited header
///
///

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ================================================================
// ================================================================

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include <iomanip>


// ##################################################################
// ##################################################################


#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include <set>
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/MPv10.h"

using namespace std;

//#define MPv10_DEBUG


//
// Timestep-limiting is important for making chemistry consistent
// with hydrodynamics.
// This is a good value for MPv3 (see Mackey,2012,A&A,539,A147)
// Will have to do some tests for MPv10...
//
#define DTFRAC 0.25

// ##################################################################
// ##################################################################


void MPv10::get_error_tolerances(
      double *reltol, ///< relative error tolerance.
      double atol[] ///< absolute error tolerances
      )
{
  *reltol = MPv10_RELTOL;
  for (int i=0;i<N_equations-1;i++) atol[i] = MPv10_ABSTOL; ///< minimum neutral fraction I care about.
  atol[N_equations-1] = MPv10_MINERG; ///< E_int: for n=1.0, T=1.0e4, ==> E=2.07e-12, so say 1e-17?
  return;
}


// ##################################################################
// ##################################################################


void MPv10::get_problem_size(
      int *ne, ///< number of equations
      int *np  ///< number of parameters in user_data vector.
      )
{
  *ne = N_equations;
  *np = N_extradata;
  return;
}

// ##################################################################
// ##################################################################

MPv10::MPv10(
      const int nd,   ///< grid dimensions
      const int csys,   ///< Coordinate System flag
      const int nv,   ///< Total number of variables in state vector
      const int ntracer,  ///< Number of tracer variables in state vector.
      const std::string *tracers,  ///< List of what the tracer variables mean.
      struct which_physics *ephys, ///< pointer to extra physics flags.
      struct rad_sources *rsrcs,   ///< radiation sources.
      const double g  ///< EOS Gamma
      )
: microphysics_base(ephys,rsrcs),
  ndim(nd), nv_prim(nv), eos_gamma(g), coord_sys(csys),
  T_min(1e0), T_max(1e9), Num_temps(100), photo_xsections()//photo_xsections(&Emin[0],&Emax[0],Nbins)
  {
  /// ===================================================================
  ///  Initialise Temperature, Recombination (radiative + dielectronic),
  ///       Ionisation, and (recomb and ionisation) Slopes Tables
  /// ===================================================================
  k_B = pconst.kB();  // Boltzmann constant.
  m_p = pconst.m_p(); // Proton mass.
  m_H = pconst.m_H(); // Hydrogen mass.
  m_He = pconst.m_He(); // Helium mass.
  m_C = pconst.m_C(); //Carbon mass
  m_N = pconst.m_N(); //Nitrogen mass
  m_O = pconst.m_O(); //Oxygen mass

  generate_lookup_tables();

  //***********************************************************
  

  cout <<"\n---------------------------------------------------------------------\n";
  cout <<"MPv10: a microphysics class.\n";

  // -----------------------------------------------------------------------------
  // --------- Set up tracer variables:                                  ---------
  // --------- (i) Identify elements present in tracer list              ---------
  // --------- (ii) Record X_mass_frac_index vector, N_elem              ---------
  // --------- (iii) Record: y_ion_index_prim (index in primitive vector),    ---------
  // ---------       y_ion_num_elec (# electrons of y_ion_index_prim species),---------
  // ---------       N_species_by_elem (ordered like X_mass_frac_index)  ---------
  // ---------                                                           ---------
  // -----------------------------------------------------------------------------
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  
  int len = ntracer;
  ftr = nv_prim -ntracer; // first tracer variable.
    
  string s; //pv_H1p=-1;
  
  N_elem = 0; N_species=0;
  
  // Set up corrector vector for fluxes, all entries = 1, i.e. unmodified
  
  /*pion_flt temp_vec[nv_prim] = {1};
  corrector = temp_vec;*/
  for (int i=0;i<nv_prim;i++) corrector[i] = 1;

  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    
    // copy tracer into map for Tr() function
    tracer_list[s] = i+ftr;
    
    // ================================================================
    //          (i) Identify elements present in tracer list
    //          (ii) Record X_mass_frac_index vector, N_elem,
    //               X_elem_atomic_mass, X_elem_number_density
    //          (iii) Record xsection arrays generated in photo_xsections
    // ================================================================
    if (s[0]=='X'){
      X_mass_frac_index.push_back(ftr + N_elem); ///<record primitive vector index of each element
      N_species_by_elem.push_back(0);
      
      //=======Helium========
      if (s.substr(2,2)=="He"){
        element_list["He"] = N_elem;
        X_elem_atomic_mass.push_back(m_He);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
      }
      //=======Hydrogen======
      else if (s[2]=='H'){
        element_list["H"]= N_elem;
        X_elem_atomic_mass.push_back(m_H);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
     }
      //=======Carbon======
      else if (s[2]=='C'){
        element_list["C"] = N_elem;
        X_elem_atomic_mass.push_back(m_C);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
      }
      //=======Nitrogen======
      else if (s[2]=='N'){
        element_list["N"] = N_elem;
        X_elem_atomic_mass.push_back(m_N);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
       }
      //=======Oxygen======
      else if (s[2]=='O'){
        element_list["O"] = N_elem;
        X_elem_atomic_mass.push_back(m_O);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
      }
      N_elem++;
    } 
  }

  //
  // ================================================================
  // (iii) Record:   N_species,    N_species_by_elem,  y_ion_num_elec
  //                 y_ion_index_prim,  y_ip1_index,        y_im1_index
  // ================================================================
  //
  
  //lv_y_ion_index_offset = ftr + N_elem; // gives the index at which ions first occur in primitive vector, maps to first index of local vector

  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    //=======Helium========
    if (s.substr(0,2)=="He"){
      cout << "\n\nTesting " << s << "\n";
      species_tracer_initialise(tracers, i, s, "He", 2, element_list["He"], len);
      if (s[2]=='1'){
        cout << s <<"\n";
        y_ip1_index_tables.push_back(4); //index of He2+ in tables
        y_ion_index_tables.push_back(3); //index of He1+ in tables
        y_im1_index_tables.push_back(2); //index of He0 in tables
      }
      else if (s[2]=='2'){
        y_ip1_index_tables.push_back(-1); //doesn't exist in tables
        y_ion_index_tables.push_back(4); //index of He2+ in tables
        y_im1_index_tables.push_back(3); //index of He1+ in tables
      }
    }
    //=======Hydrogen========
    else if (s[0] =='H'){
      cout << "\n\nTesting " << s << "\n";
      cout << "i= " << i << ", s=" << s <<", H_index = " << element_list["H"] <<"\n";
      species_tracer_initialise(tracers, i, s, "H", 1, element_list["H"], len);
      y_ip1_index_tables.push_back(-1); //doesn't exist in tables
      y_ion_index_tables.push_back(1); //index of H1+ in tables
      y_im1_index_tables.push_back(0); //index of H0 in tables
    }
    //=======Carbon========
    else if (s[0] =='C'){
      cout << "\n\nTesting " << s << "\n";
      cout << "i= " << i << ", s=" << s <<", C_index = " << element_list["C"] <<"\n";
      species_tracer_initialise(tracers, i, s, "C", 1, element_list["C"], len);
      cout << "Carbon tracer initialised";
      if (s[1]=='1'){
        y_ip1_index_tables.push_back(7); //index of C2+ in tables
        y_ion_index_tables.push_back(6); //index of C1+ in tables
        y_im1_index_tables.push_back(5); //index of C0 in tables
      }
      else if (s[1]=='2'){
        y_ip1_index_tables.push_back(8); //index of C3+
        y_ion_index_tables.push_back(7); //index of C2+
        y_im1_index_tables.push_back(6); //index of C1+
      }
      else if (s[1]=='3'){
        y_ip1_index_tables.push_back(9); //index of C4+
        y_ion_index_tables.push_back(8); //index of C3+
        y_im1_index_tables.push_back(7); //index of C2+
      }
      else if (s[1]=='4'){
        y_ip1_index_tables.push_back(10);
        y_ion_index_tables.push_back(9); //index of C4+
        y_im1_index_tables.push_back(8); 
      }
      else if (s[1]=='5'){
        y_ip1_index_tables.push_back(11); 
        y_ion_index_tables.push_back(10); //index of C5+
        y_im1_index_tables.push_back(9);
      }
      else if (s[1]=='6'){
        y_ip1_index_tables.push_back(-1); //doesn't exist
        y_ion_index_tables.push_back(11); //index of C6+
        y_im1_index_tables.push_back(10);
      }
    }
    else if (s[0] =='N'){
      cout << "\n\nTesting " << s << "\n";
      cout << "i= " << i << ", s=" << s <<", N_index = " << element_list["N"] <<"\n";
      species_tracer_initialise(tracers, i, s, "N", 1, element_list["N"], len);
      cout << "Nitrogen tracer initialised";
      if (s[1]=='1'){
        y_ip1_index_tables.push_back(14); //index of N2+ in tables
        y_ion_index_tables.push_back(13);  //index of N1+ in tables
        y_im1_index_tables.push_back(12); //index of N0 in tables
      }
      else if (s[1]=='2'){
        y_ip1_index_tables.push_back(15); //index of N3+
        y_ion_index_tables.push_back(14); //index of N2+
        y_im1_index_tables.push_back(13); //index of N1+
      }
      else if (s[1]=='3'){
        y_ip1_index_tables.push_back(16); //index of C4+
        y_ion_index_tables.push_back(15); //index of C3+
        y_im1_index_tables.push_back(14); //index of C2+
      }
      else if (s[1]=='4'){
        y_ip1_index_tables.push_back(17);
        y_ion_index_tables.push_back(16); //index of C4+
        y_im1_index_tables.push_back(15); 
      }
      else if (s[1]=='5'){
        y_ip1_index_tables.push_back(18); 
        y_ion_index_tables.push_back(17); //index of C5+
        y_im1_index_tables.push_back(16);
      }
      else if (s[1]=='6'){
        y_ip1_index_tables.push_back(19); //doesn't exist
        y_ion_index_tables.push_back(18); //index of C6+
        y_im1_index_tables.push_back(17);
      }
      else if (s[1]=='7'){
        y_ip1_index_tables.push_back(-1); //doesn't exist
        y_ion_index_tables.push_back(19); //index of C6+
        y_im1_index_tables.push_back(18);
      }
    }
    else if (s[0] =='O'){
      cout << "\n\nTesting " << s << "\n";
      cout << "i= " << i << ", s=" << s <<", O_index = " << element_list["O"] <<"\n";
      species_tracer_initialise(tracers, i, s, "O", 1, element_list["O"], len);
      cout << "Oxygen tracer initialised";
      if (s[1]=='1'){
        y_ip1_index_tables.push_back(22); //index of N2+ in tables
        y_ion_index_tables.push_back(21);  //index of N1+ in tables
        y_im1_index_tables.push_back(20); //index of N0 in tables
      }
      else if (s[1]=='2'){
        y_ip1_index_tables.push_back(23); //index of N3+
        y_ion_index_tables.push_back(22); //index of N2+
        y_im1_index_tables.push_back(21); //index of N1+
      }
      else if (s[1]=='3'){
        y_ip1_index_tables.push_back(24); //index of C4+
        y_ion_index_tables.push_back(23); //index of C3+
        y_im1_index_tables.push_back(22); //index of C2+
      }
      else if (s[1]=='4'){
        y_ip1_index_tables.push_back(25);
        y_ion_index_tables.push_back(24); //index of C4+
        y_im1_index_tables.push_back(23); 
      }
      else if (s[1]=='5'){
        y_ip1_index_tables.push_back(26); 
        y_ion_index_tables.push_back(25); //index of C5+
        y_im1_index_tables.push_back(24);
      }
      else if (s[1]=='6'){
        y_ip1_index_tables.push_back(27); //doesn't exist
        y_ion_index_tables.push_back(26); //index of C6+
        y_im1_index_tables.push_back(25);
      }
      else if (s[1]=='7'){
        y_ip1_index_tables.push_back(28); //doesn't exist
        y_ion_index_tables.push_back(27); //index of C6+
        y_im1_index_tables.push_back(26);
      }
      else if (s[1]=='8'){
        y_ip1_index_tables.push_back(-1); //doesn't exist
        y_ion_index_tables.push_back(28); //index of C6+
        y_im1_index_tables.push_back(27);
      }
    }
  }

  cout << "\nAfter reading in tracers, N_species=" << N_species;
  cout << ", N_elements=" << N_elem << "\n\n";  
  
  // ================================================================
  // ================================================================
#ifdef TESTING
  cout <<"MPv10:: EP and RS: "<<EP<<"\t"<<RS<<endl;
#endif

  // ----------------------------------------------------------------
  // --- Set up local variables: ion fraction and internal energy density.
  // ----------------------------------------------------------------

  // Get the mean mass per H atom from the He and Z mass fractions.
  // NOTE \Maggie{or let's not assume that...} Assume metal content is low enough to ignore it.
  double X = 1.0-EP->Helium_MassFrac;
  mean_mass_per_H = m_p/X;
  METALLICITY = EP->Metal_MassFrac/0.0142; // in units of solar.
  cout <<"Metallicity = "<<METALLICITY<<" of solar (0.0142)\n";

  setup_local_vectors();
  gamma_minus_one = eos_gamma -1.0;
  Min_NeutralFrac     = MPv10_ABSTOL;
  Max_NeutralFrac     = 1.0-MPv10_ABSTOL;

  // ----------------------------------------------------------------
  // We want to set up the CIE cooling function for metals-only from WSS09
  // (i.e. with H+He cooling subtracted out).
  // ----------------------------------------------------------------
  setup_WSS09_CIE_OnlyMetals();
  // ================================================================
  // ================================================================


  // ----------------------------------------------------------------
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  setup_cvode_solver_without_Jacobian();
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // Set flags for whether we have radiation sources.
  // ----------------------------------------------------------------
  for (int isrc=0; isrc<RS->Nsources; isrc++) {
    if (RS->sources[isrc].type==RT_SRC_SINGLE &&
        RS->sources[isrc].effect==RT_EFFECT_MFION) {
      N_rad_src++;
      rt_data.resize(N_rad_src);
      int err=set_multifreq_source_properties(&RS->sources[isrc],
                                      rt_data[N_rad_src-1].strength);
      if (err)
        rep.error("multifreq photoionisation setup MPv10",err);
    }
  }
  cout <<"\t\tMPv10: got "<<N_rad_src<<" radiation sources.\n";
  // ----------------------------------------------------------------
  
  cout <<"MPv10: Constructor finished and returning.\n";
  cout <<"---------------------------------------------------------------------\n\n";
  return;
}



// ##################################################################
// ##################################################################



int MPv10::set_multifreq_source_properties(
      const struct rad_src_info *rsi, ///< source data
      double *str  ///< O/P source luminosity per energy bin (erg/s/bin).
      )
{
  if (rsi->effect!=RT_EFFECT_MFION)
    rep.error("Wrong source type for id",rsi->id);

  // TODO:
  // Now we need to figure out how to get the luminosity of the star
  // in each frequency bin, in erg/s/bin.
  // * rsi->strength gives the luminosity of the star in erg/s
  // * rsi->Tstar    gives the effective temperature of the star.
  // * rsi->Rstar    gives the Radius of the star.
  // If the star were a blackbody, then this would be enough to
  // calculate the luminosity in each bin, if we have the bin ranges
  // set (which we do).  Unfortunately a BB is a bad approximation.
  // Maybe it is the best we can do for now.
  //
  // We want to add the luminosity in each bin to the array "str".
  return 1;
}



// ##################################################################
// ##################################################################



void MPv10::species_tracer_initialise(
    const std::string *tracers,  ///< List of what the tracer variables mean.
    int i, ///<index of current tracer in for loop
    string s, /// < current tracer in for loop
    string el_symbol, ///< element symbol, e.g. "He", "H"
    int el_symbol_length, ///< e.g. "H" is of length 1, "He" of length 2.
    int el_index, /// element index in N_species_by_elem, used in for loops for densities etc
    int length /// < length of tracers vector
    ){
  double *xsections = parameters[s].xsections;
  int len = get_nbins();
  y_ion_xsections.push_back( std::vector<double>(xsections,xsections+len) );
  
  N_species_by_elem[el_index]++;
  y_ion_index_prim.push_back(ftr + N_elem + N_species);//lv_y_ion_index_offset + N_species);
  y_ion_index_local.push_back(N_species);
     
  /// Use stringstream to convert electron number string to int.
  int num_elec_int;
  stringstream ss_e; 
  ss_e << s.substr(el_symbol_length,1); 
  ss_e >> num_elec_int; 
  /// Record the electron number and He_ion_index vector (used to identify tracer index from string)
  y_ion_num_elec.push_back(num_elec_int);
  
  /// Define the tracer for the next ion up / down in tracers list
  stringstream ss_ip1;        
  stringstream ss_im1;
  stringstream ss_neutral;
  ss_ip1 << el_symbol << num_elec_int +1 << "+";  
  ss_im1 << el_symbol << num_elec_int -1 << "+";
  ss_neutral << el_symbol << 0 << "+";
  string neutral = ss_neutral.str();
  string ip1 = ss_ip1.str();
  string im1 = ss_im1.str();
      
  /// Check if the next ion up / down exists in tracer list (or if it's neutral) to define y_ip1_index and y_im1_index.
  ///< due to ordering in tracer list, if the next ion exists, its index will be 1 above the current index.
  if (i+1 < length){
    if ( ip1==tracers[i+1] ){
      cout << "Next ion up = " << ip1<<"\n";
      y_ip1_index_local.push_back(N_species + 1); 
    }
    else{
      cout <<"Next ion up doesn't exist\n";
      y_ip1_index_local.push_back(-1); // -1 flag => species doesn't exist 
    }
  }
  else{
    cout <<"Next ion up doesn't exist\n";
    y_ip1_index_local.push_back(-1); // -1 flag => species doesn't exist 
  }
  ///< due to ordering in tracer list, if the previous ion exists, its index will be 1 below the current index.
  if ( im1 == neutral){
    cout <<"Next ion down is neutral \n";
    y_im1_index_local.push_back(-2); // -2 flag => species is neutral
  }
  else if (im1==tracers[i-1] ){//check if the next lowest ion is the neutral species
    cout << "Next ion down = " << im1<<"\n";
    y_im1_index_local.push_back(N_species - 1); 
  }
  else{
    y_im1_index_local.push_back(-1); // -1 flag => species doesn't exist
    cout <<"Next ion down doesn't exist.\n";
  }
  N_species++;
  
  return;
}




// ##################################################################
// ##################################################################



void MPv10::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  nvl     = N_species +1;    // two local variables to integrate
  N_extradata = 0;
  N_equations = nvl;
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
  lv_eint = N_species;
  //cout<<"!!!!!!!!!!!!!!!!!! nvl="<<nvl<<"\n";
  return;
  }

  

// ##################################################################
// ##################################################################



MPv10::~MPv10()
{
  //
  // Free vector memory
  //
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}



// ##################################################################
// ##################################################################



double MPv10::get_n_elec(
        const pion_flt *P ///< primitive state vector.
        )
{
  // get electron number density from P
  int species_counter=0;
  double ne;
  for (int elem=0;elem<N_elem;elem++) {//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    double X_elem_density = P[ftr + elem];
    
    for (int s=0;s<N_elem_species;s++) {//loop over every species in THIS element
      double y_ion_frac = P[ftr + N_elem + species_counter];
      
      //add to ne based on the number of electrons liberated to obtain this ion
      pion_flt number_density = X_elem_density*y_ion_frac;
      
      int num_elec = y_ion_num_elec[species_counter];
      
      ne += num_elec*number_density;
      species_counter ++;
    }
  }
  return 10.0;
}
      


// ##################################################################
// ##################################################################



double MPv10::get_n_ion(
        string name, ///< ion name
        const pion_flt *P ///< primitive state vector.
        )
{
  if ( Tr(name) == -1)
    return -1;
  else{
    //index of ION in primitive vector
    int ion_frac_index = Tr(name);
    
    //get index of element number density X_ELEM in primitive vector
    int elem_index;
    if (name.substr(0,2)=="He") { 
      elem_index = element_list["He"];
    }
    else{
      elem_index = element_list[name.substr(0,1)];
    }
    
    // element number density times ion fraction from P
    double elem_number_density = P[ftr+elem_index];
    double ion_fraction = P[ion_frac_index];
    
    return elem_number_density * ion_fraction;
  }
}
      

// ##################################################################
// ##################################################################



int MPv10::Tr(const string s)
{
  if (tracer_list.find(s)==tracer_list.end())
    return -1;
  else
    return tracer_list[s]; 
}



// ##################################################################
// ##################################################################



double MPv10::get_temperature(
      double *y_ion_frac, //< y_ion_fraction (by y_ion_index_prim)
      vector<pion_flt>& X_number_density,//const double nH, ///< nH
      const double E ///< E_int (per unit volume)
      )
{
  //
  // returns gas temperature according to E=nkT/(g-1), => T = E*(g-1)/(K*n)
  /*cout << "temp=" << gamma_minus_one*E/(k_B*get_ntot(y_ion_frac,X_number_density));
  cout << ", ntot = " <<get_ntot(y_ion_frac,X_number_density) ;
  cout << ", E = " << E << "\n";*/
  return gamma_minus_one*E/(k_B*get_ntot(y_ion_frac,X_number_density));
}



// ##################################################################
// ##################################################################



double MPv10::get_ntot(
      double *y_ion_frac,//<y_ion_fraction (by y_ion_index_prim)
      vector<pion_flt>& X_number_density//const double nH, ///< nH
      )
{
  int species_counter = 0;
  pion_flt n_tot=0;
  
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    pion_flt neutral_frac = 1; //neutral frac, got by subtracting off the ion fractions in the next loop.
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      pion_flt number_density = X_number_density[e]*y_ion_frac[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      
      n_tot += (1+num_elec)*number_density; //add on number of particles got by electrons + ionised atom
      
      neutral_frac -= y_ion_frac[species_counter];
      species_counter ++;
    }
    n_tot += neutral_frac*X_number_density[e];
  }
  return n_tot;
}



// ##################################################################
// ##################################################################



int MPv10::convert_prim2local(
      const pion_flt *p_in, ///< primitive vector from grid cell (length nv_prim)
      double *p_local,
      int function_flag /// < flag to say which function called this function
      )
{
  //
  // ==============================================================
  //                 Set NUMBER DENSITY VECTOR.
  // ==============================================================
  //
  for (int i=0;i<N_elem;i++){//loop over every element
    double n_X = p_in[RO]*( p_in[ X_mass_frac_index[i]] / X_elem_atomic_mass[i]);
    X_elem_number_density[i] = n_X;
  }
  
  //rep.printSTLVec("n_x",X_elem_number_density);
  //
  // ==============================================================
  //                 Set INTERNAL ENERGY in local vector.
  // ==============================================================
  //
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);
  //
  // ==============================================================
  //            Set ION FRACTIONS OF ELEMENT in LOCAL vector,
  //            converted from MASS FRACTIONS in PRIM vector.
  // ==============================================================
  //  
  //loop over every species and set the local vector, noting that loc_vec[i]=prim_vec[i]/X_e
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      p_local[ y_ion_index_local[species_counter]] = 
        p_in[y_ion_index_prim[species_counter]]/p_in[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }

  species_counter=0;
  int alter_flag = 0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
#ifdef MPv10_DEBUG
      cout <<"e="<<e<<", s="<<s;
      cout <<", species_counter="<<species_counter;
      cout <<", ion index=";
      cout <<y_ion_index_local[species_counter]<<"\n";
      cout <<"Min_NeutralFrac="<<Min_NeutralFrac;
      cout <<", 1-Max_NeutralFrac="<<1.0-Max_NeutralFrac<<"\n";
      if (p_local[y_ion_index_local[species_counter]]>1.01 || 
          p_local[y_ion_index_local[species_counter]]<-0.01) {
        cout <<"MPv10::convert_prim2local: bad ion fraction: ";
        cout <<"x(H0)="<<p_local[y_ion_index_local[species_counter]];
        cout <<", resetting to [0,1]\n";
      }
      cout << "function_flag = " << function_flag << "\n";
#endif

      // Introducing sense checks -- make sure value is positive, less than 1
        if ( p_local[y_ion_index_local[species_counter]] < (-2.0 * MPv10_ABSTOL) ) {
          alter_flag = 1;
        }
        else if ( p_local[y_ion_index_local[species_counter]] > (1.0 + MPv10_ABSTOL)) {
          alter_flag = 1;
        };
      species_counter ++;
    }
  }
  
  if (alter_flag == 1){
    sCMA(corrector, p_in);
    for (int t=0;t<nv_prim;t++){
      //cout << "correction = " << corrector[t] << "\n";
      p_local[t] = p_local[t]*corrector[t];}
  }

  //
  // Check for negative pressure (note this shouldn't happen, so we output a
  // warning) and set to 10K if we find it.
  //
  if (p_local[lv_eint]<=0.0) {
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_prim2local: negative pressure input: p=";
    cout <<p_local[lv_eint]<<", setting to "<<EP->MinTemperature<<"K.\n";
#endif
    //Define y_ion_frac
    pion_flt y_ion_frac[N_species];
    for (int s=0;s<N_species;s++){
      y_ion_frac[ s] = p_local[ s];
    }
    //reset the internal energy (requires using y_ion_frac in get_ntot)
    p_local[lv_eint] = get_ntot(y_ion_frac,X_elem_number_density)*k_B*EP->MinTemperature/(gamma_minus_one);
  }


#ifdef MPv10_DEBUG
  //
  // Check for NAN/INF
  //
  for (int v=0;v<2;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics",p_local[v]);
  }
  if (mpv_nH<0.0 || !isfinite(mpv_nH))
    rep.error("Bad density input to MPv10::convert_prim2local",mpv_nH);
#endif // MPv10_DEBUG

  //rep.printVec("prim2local local",p_local,nvl);
  //rep.printVec("prim2local prim ",p_in,nv_prim);

  return 0;
}



// ##################################################################
// ##################################################################



int MPv10::convert_local2prim(
      const double *p_local,
      const pion_flt *p_in, ///< input primitive vector
      pion_flt *p_out,      ///< updated primitive vector
      int function_flag /// < flag to say which function called me 
      )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];
  //
  //=================================================================
  //     Update output primitive vector according to local vector
  //     Also, define y_ion_frac while you're at it.
  //=================================================================
  //
  p_out[PG]    = p_local[lv_eint]*(gamma_minus_one);
  pion_flt y_ion_frac[N_species];

  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      p_out[ y_ion_index_prim[ species_counter]] = p_local[ y_ion_index_local[species_counter]] * p_in[ X_mass_frac_index[e]];
      y_ion_frac[ species_counter] = p_local[ y_ion_index_local[species_counter]];
      species_counter ++;
    }
  }

  // Set mass fraction tracers to be within the required range (not too close to zero or 1)
  
  species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    p_out[ X_mass_frac_index[ e]] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[ X_mass_frac_index[ e]])));

    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element

      #ifdef MPv10_DEBUG
      // Introducing sense checks -- make sure value is positive, less than 1
      if ( static_cast<double>(p_out[ y_ion_index_prim[ species_counter]]) < (-2 * MPv10_ABSTOL) ) {
        cout << "convert_local2prim: " << function_flag << " mass fraction goes negative here. \n [";
        for (int v=0; v<nv_prim;v++) {cout << p_out[v] << ", ";}
        cout << "] \n";
        print_flag = 1;
      }

      else if ( static_cast<double>(p_out[ y_ion_index_prim[ species_counter]]) > (1 + MPv10_ABSTOL) * static_cast<double>(p_out[ X_mass_frac_index[ e]]) ) {
        cout << "convert_local2prim: " << function_flag << " mass frac too large for species " << s << ": X = " << p_out[ y_ion_index_prim[ species_counter]] << "\n";
        cout << "Prim vector: \n [";
        for (int v=0; v<nv_prim;v++) {cout << p_out[v] << ", ";}
        cout << "] \n";
        print_flag = 1;
      }
      #endif

      p_out[ y_ion_index_prim[ species_counter]] = max(Min_NeutralFrac, min(static_cast<double>(p_out[ X_mass_frac_index[ e]])*Max_NeutralFrac, static_cast<double>(p_out[ y_ion_index_prim[ species_counter]])));
      species_counter ++;
    }
  }
  
  

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected x(H+) from p_out[]).
  //
  
  double T = get_temperature(y_ion_frac, X_elem_number_density, p_local[lv_eint]);
  //cout <<"current temp = " << T<<"\n";
  if (T>1.01*EP->MaxTemperature) {
    Set_Temp(p_out,EP->MaxTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() HIGH temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MaxTemperature<<"\n";
#endif // MPv10_DEBUG
  }
  if (T<0.99*EP->MinTemperature) {
    cout << "Low temperature!!!\n";
    Set_Temp(p_out,EP->MinTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() LOW  temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MinTemperature<<"\n";
#endif // MPv10_DEBUG
  }
  p_out[PG] = p_local[lv_eint] *(gamma_minus_one);
  //cout << "pressure=" <<p_out[PG]<<"\n";
  //rep.printVec("local2prim local",p_local,nvl);
  //rep.printVec("local2prim prim ",p_out,nv_prim);
  return 0;
}



// ##################################################################
// ##################################################################



double MPv10::Temperature(
      const pion_flt *pv, ///< primitive vector
      const double      ///< eos gamma
      )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    //cout <<"MPv10::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];

  convert_prim2local(pv,P, 4);
  
  pion_flt y_ion_frac[N_species];
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      y_ion_frac[ species_counter] = pv[ y_ion_index_prim[ species_counter]] / pv[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  return (get_temperature(y_ion_frac, X_elem_number_density, P[lv_eint]));
}



// ##################################################################
// ##################################################################



int MPv10::Set_Temp(
      pion_flt *p_pv,   ///< primitive vector.
      const double T, ///< temperature
      const double    ///< eos gamma.
      )
{
  //
  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.
  //
  //cout << "set temperature\n";
  if (p_pv[PG]<=0.0) {
    cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p_pv[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }
  double P[nvl];

  int err = convert_prim2local(p_pv,P, 1);
  
  //Determine y_ion_frac from the primitive vector
  pion_flt y_ion_frac[N_species];
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      y_ion_frac[ species_counter] = p_pv[ y_ion_index_prim[ species_counter]] / p_pv[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  
  //Determine internal energy using get_ntot
  P[lv_eint] = get_ntot(y_ion_frac,X_elem_number_density)*k_B*T/(gamma_minus_one);
  
  //Call convert_local2prim with the new local vector; this will
  // generate a new temperature value;
  err += convert_local2prim(P, p_pv, p_pv, 1);
  return err;
}



// ##################################################################
// ##################################################################



int MPv10::TimeUpdateMP(
      const pion_flt *p_in,   ///< Primitive Vector to be updated.
      pion_flt *p_out,        ///< Destination Vector for updated values.
      const double dt,      ///< Time Step to advance by.
      const double,         ///< EOS gamma.
      const int,            ///< Switch for what type of integration to use.
      double *random_stuff  ///< Vector of extra data (column densities, etc.).
      )
{
  //
  // Call the new update function, but with zero radiation sources.
  //
  std::vector<struct rt_source_data> temp;
  int err = TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}



// ##################################################################
// ##################################################################


  
int MPv10::TimeUpdateMP_RTnew(
      const pion_flt *p_in, ///< Primitive Vector to be updated.
      const int,                                  ///< unused.
      const std::vector<struct rt_source_data> &, ///< unused.
      const int,                                  ///< unused.
      std::vector<struct rt_source_data> &ion_src,
        ///< list of ionising src column densities and source properties.
      pion_flt *p_out,  ///< Destination Vector for updated values
                        ///< (can be same as first Vector.
      const double dt,  ///< Time Step to advance by.
      const double,     ///< EOS gamma.
      const int, ///< Switch (unused)
      double *random_stuff ///< final temperature (debugging).
      )
{
  int err=0;
  double P[nvl];

  err = convert_prim2local(p_in,P, 2);
  /*rep.printVec("p2l start prim ",p_in,nv_prim);
  cout <<"\n";
  rep.printVec("p2l start local",P,nvl);
  cout <<"\n\n";*/
  if (err) {
    rep.error("Bad input state to MPv10::TimeUpdateMP()",err);
  }
  setup_radiation_source_parameters(p_in, P, ion_src);


  // Populates CVODE vector with initial conditions (input)
  for (int v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];

  //
  // Calculate y-dot[] to see if anything is changing significantly over dt
  //
  double maxdelta=0.0;
  err = ydot(0, y_in, y_out, 0);
  /*cout <<"H+-ion-frac="<< NV_Ith_S(y_in,0) <<"\n";
  cout <<"He+-ion-frac="<< NV_Ith_S(y_in,1) <<"\n";
  cout <<"He2+-ion-frac="<< NV_Ith_S(y_in,2) <<"\n";
  cout <<"E_int="<< NV_Ith_S(y_in,3) <<"\n";*/

  if (err) 
    rep.error("dYdt() returned an error in MPv10::TimeUpdateMP_RTnew()",err);
  for (int v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out,v)*dt/NV_Ith_S(y_in,v)));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < EULER_CUTOFF) {
    //cout <<"maxdelta="<<maxdelta<<", Doing euler integration...\n";
    for (int v=0;v<nvl;v++) {
      NV_Ith_S(y_out,v) = NV_Ith_S(y_in,v) + dt*NV_Ith_S(y_out,v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    //cout <<"maxdelta="<<maxdelta<<", Doing CVODE integration...\n";
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      rep.error("integration failed: MPv10::TimeUpdateMP_RTnew()",err);
    }
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v=0;v<nvl;v++) P[v] = NV_Ith_S(y_out,v);
  err = convert_local2prim(P,p_in,p_out, 2);

#ifdef TEST_INF
  for (int v=0;v<nv_prim;v++) {
    if (!isfinite(p_in[v]) || !isfinite(p_out[v])) {
      cout <<"NAN in MPv3 update: "<<v<<"\n";
      rep.printVec("Pin ",p_in,nv_prim);
      rep.printVec("Pout",p_out,nv_prim);
      rep.printVec("Ploc ",P,nvl);
      //rep.error("NAN in MPv3",P[2]);
      return 1;
    }
  }
#endif

  //rep.printVec("l2p end prim ",p_out,nv_prim);
  //rep.printVec("l2p end local",P,nvl);

  return err;
}



// ##################################################################
// ##################################################################



double MPv10::timescales(
      const pion_flt *p_in, ///< Current cell state vector.
      const double,   ///< EOS gamma.
      const bool, ///< set to 'true' if including cooling time.
      const bool, ///< set to 'true' if including recombination time.
      const bool  ///< set to 'true' if including photo-ionsation time.
      )
{
#ifdef MPv10_DEBUG
  if (RS->Nsources!=0) {
    cout <<"WARNING: MPv10::timescales() using non-RT version!\n";
  }
#endif // MPv10_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin= timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}



// ##################################################################
// ##################################################################



///
/// This returns the minimum timescale of all microphysical processes, including
/// reaction times for each species and the total heating/cooling time for the gas.
/// It requires the radiation field as an input, so it has substantially greater
/// capability than the other timescales function.
/// Default setting is DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double MPv10::timescales_RT(
      const pion_flt *p_in, ///< Current cell state vector.
      const int N_heat,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &heat_src,
      ///< list of UV-heating column densities and source properties.
      const int N_ion,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &ion_src,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      )
{
  
  /*
   * DTAU TESTING
  double dtau[get_nbins()] = {0};
  get_dtau(1e10, p_in, dtau);
  cout << "dtau = " << dtau[5] << "\n";*/
  
  
  int err=0;
  //
  // First convert to local variables.
  //
  double P[nvl];

  //cout << "\n timescales_RT calling convert_prim2local \n";
  err = convert_prim2local(p_in,P,3);
  if (err) {
    rep.error("Bad input state to MPv10::timescales_RT()",err);
  }
  for (int v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];
  NV_Ith_S(y_in,lv_eint) = P[lv_eint];

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in MPv10::timescales_RT()",err);
  }

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t=HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the energy change 
  // (heating timescale) is not required for accurately tracking ionisation fronts 
  // (although it may be needed for cooling!).
  // For testing purposes there are ifdefs to allow the code to use a relative 
  // change in neutral fraction and/or the relative change in energy as the 
  // timestep criterion, rather than the default of absolute change in neutral
  // fraction.
  //
  //LOOP OVER ALL THE IONS HERE, NOT JUST LV_H0, AND THAT SHOULD FIX YOUR SEGFAULT PROBLEM.
  for (int v=0;v<N_equations;v++) t = min(t,DTFRAC/(fabs(NV_Ith_S(y_out, v))+TINYVALUE));
  //cout <<"limit by dx: dt="<<DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\n";


#ifdef MPv10_DEBUG
  if (t<3.16e9) {
  //cout <<"MP timescales: xdot="<<NV_Ith_S(y_out, lv_H0);
  cout <<", Edot="<<NV_Ith_S(y_out, lv_eint)<<" t_x="<<t;
  }
#endif // MPv10_DEBUG

#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  //
  // Now cooling/heating time to limit to X% change in energy).
  //
  t = min(t,DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE));
  //cout <<"using fractional energy change.\n";
  //cout <<"limit by dE: dt="<<DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE)<<"\n";
#endif

#ifdef MPv10_DEBUG
  if (t<3.16e9) {
  cout <<" and min(t_x,t_e)="<<t<<",  "; rep.printVec("P[1-x,E]",P,nvl);
  }
#endif // MPv10_DEBUG
  return t;
}



// ##################################################################
// ##################################################################



void MPv10::sCMA(
    pion_flt *corrector, ///< input corrector vector
    const pion_flt *p_in) ///< input primitive vector from grid cell (length nv_prim)
{
  //  Re-initialise corrector every step
  for (int i=0;i<nv_prim;i++) corrector[i] = 1.0;
  int print_flagg = 0;
  double total_mass_frac = 0;
  
  //loop over every species and get the sum
  int species_counter=0;
  
  // Calculate all-element correction
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    total_mass_frac += p_in[ X_mass_frac_index[e]];
  }
  double e_correction = 1.0 / total_mass_frac;
  species_counter = 0;
  // apply all-element correction, calculate species correction, apply species correction
  for (int e=0;e<N_elem;e++)  {  //loop over every element
    int N_elem_species=N_species_by_elem[e];  
    corrector[ X_mass_frac_index[e]] = e_correction;   // correct THIS element
    // Calculate all-species-pr-element correction, if needed, i.e.
    double s_frac = 0;
    
    for (int s=0;s<N_elem_species;s++) {
      s_frac += p_in[ y_ion_index_prim[species_counter]];
      species_counter ++;
    }
    
    if ( s_frac > ((p_in[ X_mass_frac_index[e]]* e_correction)  - Min_NeutralFrac)) {
      print_flagg = 1;
      double s_correction = ((p_in[ X_mass_frac_index[e]]* e_correction)  - Min_NeutralFrac) / s_frac;
      int inner_species_counter = (species_counter - N_elem_species);
      for (int s=0;s<N_elem_species;s++) {
        corrector[ y_ion_index_prim[inner_species_counter]] = s_correction;
      }
    }
  }
  return;
}



// ##################################################################
// ##################################################################



void MPv10::setup_radiation_source_parameters(
      const pion_flt *p_in, ///< primitive input state vector.
      double *P,  ///< local input state vector (x_in,E_int)
      std::vector<struct rt_source_data> &ion_src
      ///< list of ionising src column densities and source properties.
      )
{
  for (unsigned int v=0; v<ion_src.size(); v++)
    rt_data[v] = ion_src[v];

//struct rt_source_data {
//  double Vshell;   ///< Shell volume for discrete photo-ionisation/-heating rates.
//  double dS;       ///< Path length through cell.
//  double strength[MAX_TAU]; ///< Luminosity (or flux if source at infinity).
//  double Column[MAX_TAU];  ///< integral of quantities along LOS to near edge of cell.
//  double DelCol[MAX_TAU];  ///< integral of quantities along LOS through cell.
//  int id;   ///< source id.
//  int type; ///< diffuse-radiation or a real source.
//  short unsigned int NTau; ///< Number of LOS quantities traced for the source.
//};
  
  // we want to use these data in ydot to loop over each radiation
  // source, and calculate the ionization rate (using the photon-
  // conserving formula) for each species.
  // - Strength = erg/s/bin luminosity of source.
  // - Column = Tau to edge of cell
  // - DelCol = dTau through cell (we'll re-calculate this in ydot)
  // - Vshell and dS are obvious from Mellema et al. paper.

  return;
}



// ##################################################################
// ##################################################################



void MPv10::get_dtau(
        const pion_flt ds,    ///< ds, thickness of the cell
        const pion_flt *p_in, ///< input primitive vector
        pion_flt *dtau_vec    ///< output dtau vector
        )
{   
  for (int bin=0; bin<get_nbins(); bin++){
    double dtau = 0; // sum dtau across all species within this bin
    int species_counter=0;
    
    for (int e=0;e<N_elem;e++){//loop over every element
      int N_elem_species=N_species_by_elem[e];
      
      for (int s=0;s<N_elem_species;s++){//loop over every species
        double n_s = p_in[RO]*( p_in[ y_ion_index_prim[species_counter]] / X_elem_atomic_mass[e]);
        double xsec = y_ion_xsections[species_counter][bin];
        dtau += n_s*xsec*ds;
        species_counter ++;
      }
    }
    dtau_vec[bin] = dtau;
  }
  return;
}


  
  
// ##################################################################
// ##################################################################



int MPv10::ydot(
      double,               ///< current time (UNUSED)
      const N_Vector y_now, ///< current Y-value
      N_Vector y_dot,       ///< vector for Y-dot values
      const double *        ///< extra user-data vector (UNUSED)
      )
{
  //
  //  ========================================================
  //        Determine ne, y_ion_frac, and number density
  //        Also determine neutral fraction for later
  //  ========================================================
  //
  double ne=0;
  double y_ion_frac[N_species];
  double X_neutral_frac[N_elem];
  
  int species_counter=0;
  for (int elem=0;elem<N_elem;elem++) {//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++) {//loop over every species in THIS element
      //add to y_ion_frac
      y_ion_frac[species_counter] = NV_Ith_S(y_now,y_ion_index_local[species_counter]);
      //add to ne based on the number of electrons liberated to obtain this ion
      pion_flt number_density = X_elem_number_density[elem]*y_ion_frac[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      ne += num_elec*number_density;
      species_counter ++;
    }
  }

  // Find internal energy
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double dydt[N_equations];
//   for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
//   rep.printVec("ydot going in",dydt,N_equations);
    //initialise ydot to zero.
  for (int v=0;v<N_equations;v++) NV_Ith_S(y_dot, v) =0;

  // Use E_in, y_ion_frac and number density to determine temperature
  double T = get_temperature(y_ion_frac, X_elem_number_density, E_in);
  double Edot=0.0;
  // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH
  // at the end.

  // get neutral fractions
  species_counter=0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    X_neutral_frac[elem] =1.0;
    //cout << "\n neutral_frac=" << neutral_frac << "\n";
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      X_neutral_frac[elem] -=  y_ion_frac[species_counter];
      species_counter ++;
    }
  }
  //rep.printVec("ynow X_neutral_frac",X_neutral_frac,N_elem);
  //rep.printVec("ynow y_ion_frac",y_ion_frac,N_species);
  double yyy[N_equations];
  for (int v=0;v<N_equations;v++) yyy[v] = NV_Ith_S(y_now,v);
  //rep.printVec("ynow",yyy,N_equations);
  //cout << "T=" << T << "\n";
  //
  //  ========================================================
  //          Get Rate of Change of Each Species
  //  ========================================================
  //
  /// Start by getting the relevant temperature index:
  if (T > T_max){
    //NOTE RESET ENERGY HERE TOO!!!
    T = T_max;
  }
  else if ( T<T_min){
    T = T_min;
  }
  int temp_index = int (( log10(T) - log10(T_min) ) / delta_log_temp );
  double dT = T - Temp_Table[temp_index];


  
  /// ====== Collisional ionisation INTO this species, OUT of previous species ========
  /// this_y_dot(ion) += ionise_rate(im1)*n_e*y(im1) <<< add ionisation from less ionised species to current species
  /// =================================================================================
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    double neutral_frac = X_neutral_frac[elem];   
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists in tracer list
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_ion_rate = ionise_rate_table[y_im1_index_tables[species_counter]] [ temp_index];
        double upper_ion_rate_contrib = dT * ionise_slope_table[ y_im1_index_tables[species_counter]] [temp_index];
        double ionise_rate_im1 = lower_ion_rate + upper_ion_rate_contrib;
        
        //if the less ionised species ISN'T neutral 
        if (y_im1_index_local[species_counter] !=-2){
          double this_y_dot = ionise_rate_im1 *  NV_Ith_S(y_now, y_im1_index_local[species_counter]) *ne;   
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
          NV_Ith_S(y_dot, y_im1_index_local[species_counter]) -= this_y_dot;
          
          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot * X_elem_number_density[elem];
        }
        
        //if the less ionised species IS neutral
        else{
          double this_y_dot = ionise_rate_im1 * neutral_frac * ne;
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot *X_elem_number_density[elem];
        }
      }
      species_counter ++;
    }
  }
  

  /// ============== Radiative recombination OUT OF this species and INTO NEXT species =====================
  /// y_dot(ion) -= recomb_rate(ion)*n_e*y(ion) <<< subtract recombination to less ionised species
  /// ================================================================================  
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_recomb_rate = recomb_rate_table [y_ion_index_tables[species_counter] ][temp_index];
        double upper_recomb_rate_contrib = dT * ionise_slope_table[y_ion_index_tables[species_counter] ][temp_index];
        double recomb_rate_ion = lower_recomb_rate + upper_recomb_rate_contrib;

        double this_y_dot = recomb_rate_ion * NV_Ith_S(y_now, y_ion_index_local[species_counter])*ne;
        //Subtract this rate of change from the current species
        NV_Ith_S(y_dot, y_ion_index_local[species_counter]) -= this_y_dot;
        //add this rate of change to the recombination species (provided it isn't neutral)
        if (y_im1_index_local[species_counter] !=-2){
          NV_Ith_S(y_dot, y_im1_index_local[species_counter]) += this_y_dot;
        }
//         cout << "y_dot [" << species_counter << "] -=" << this_y_dot <<"\n";
//         for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
//         rep.printVec("ydot",dydt,N_equations);
        /// =========  HEATING DUE TO RECOMBINATION OUT OF THIS SPECIES ===========
        double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
        //
        Edot -= (3./2.) * T * k_B * this_y_dot *X_elem_number_density[elem];
     
      }
    species_counter ++;
    }
  }
  
  //Edot -= cooling_rate_SD93CIE(T) *ne;
  // electron impact excitation
  //Edot = 0;
  NV_Ith_S(y_dot,lv_eint) = Edot;

  
  for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
  //rep.printVec("ydot returned ",dydt,N_equations);
  
  return 0;
}




































// ####################################################################
// ####################################################################

  
void MPv10::set_atomic_data(){
  double eV_per_Erg = 1.602e-12; // number of eV in 1erg
  double amu_to_mp = 1.661/1.673; // number of amu in one proton mass.
  struct ion_struct i;
  
    /** \section numfrac Number Fractions
   * Number fractions are from Lodders,K., 2003, ApJ, 591, 1220-1247
   *  \section Masses
   * Masses are from Wikipedia, in Atomic Mass Units (amu), converted to units of proton mass.
   * */
  i.ion = "H0"; i.ip1 = "H1+"; i.im1 = "";
  i.el  = "H";
  i.charge = 0;
  i.ion_pot = eV_per_Erg *13.59844;
  i.i = H_0;
  i.g_stat = 2;
  ion_props["H0" ] = i;
  
  i.ion = "H1+"; i.ip1 = ""; i.im1 = "H0";
  i.el  = "H";
  i.charge = 1;
  i.ion_pot = -1.0e99;
  i.i = H_1p;
  i.g_stat = 1;
  ion_props["H1+" ] = i;

  i.ion = "He0";  i.ip1 = "He1+";  i.im1 = "";
  i.el  = "He";
  i.charge = 0;
  i.ion_pot = eV_per_Erg *24.58741;
  i.i = He0;
  i.g_stat = 1;
  ion_props["He0" ] = i;

  i.ion = "He1+"; i.ip1 = "He2+";  i.im1 = "He0";
  i.el  = "He";
  i.charge = 1;
  i.ion_pot = eV_per_Erg *54.41778;
  i.i = He1p;
  i.g_stat = 2;
  ion_props["He1+" ] = i;
  
  i.ion = "He2+";  i.ip1 = "";  i.im1 = "He1+";
  i.el  = "He";
  i.charge = 2;
  i.ion_pot = -1.0e99;
  i.i = He2p;
  i.g_stat = 1;
  ion_props["He2+" ] = i;

  // carbon: C, C1+, C2+, C3+, C4+, C5+, C6+
  i.ion = "C0"; i.ip1 = "C1+"; i.im1 = ""; i.el  = "C";
  i.charge = 0; i.ion_pot = eV_per_Erg *11.3;
  i.i = C0; i.g_stat = 1; ion_props["C0" ] = i;

  i.ion = "C1+"; i.ip1 = "C2+"; i.im1 = "C0"; i.el  = "C";
  i.charge = 1; i.ion_pot = eV_per_Erg *24.4;
  i.i = C1p;  i.g_stat = 1; ion_props["C1+" ] = i;

  i.ion = "C2+"; i.ip1 = "C3+"; i.im1 = "C1+"; i.el  = "C";
  i.charge = 2; i.ion_pot = eV_per_Erg *47.9;
  i.i = C2p;  i.g_stat = 1; ion_props["C2+" ] = i;
  
  i.ion = "C3+"; i.ip1 = "C4+"; i.im1 = "C2+"; i.el  = "C";
  i.charge = 3; i.ion_pot = eV_per_Erg *64.5;
  i.i = C3p; i.g_stat = 1; ion_props["C3+" ] = i;

  i.ion = "C4+"; i.ip1 = "C5+"; i.im1 = "C3+"; i.el  = "C";
  i.charge = 4; i.ion_pot = eV_per_Erg *392.1;
  i.i = C4p; i.g_stat = 1; ion_props["C4+" ] = i;

  i.ion = "C5+"; i.ip1 = "C6+"; i.im1 = "C4+"; i.el  = "C";
  i.charge = 5; i.ion_pot = eV_per_Erg *490.0;
  i.i = C5p;  i.g_stat = 1; ion_props["C5+" ] = i;

  i.ion = "C6+"; i.ip1 = ""; i.im1 = "C5+"; i.el  = "C";
  i.charge = 6; i.ion_pot = -1.e99;
  i.i = C6p;  i.g_stat = 1; ion_props["C6+" ] = i;
 
  // nitrogen: N, N1+, N2+, N3+, N4+, N5+, N6+, N7+
  i.ion = "N0"; i.ip1 = "N1+"; i.im1 = ""; i.el  = "N";
  i.charge = 0; i.ion_pot = eV_per_Erg *14.5;
  i.i = N0; i.g_stat = 1; ion_props["N0" ] = i;

  i.ion = "N1+"; i.ip1 = "N2+"; i.im1 = "N0"; i.el  = "N";
  i.charge = 1; i.ion_pot = eV_per_Erg *29.6;
  i.i = N1p; i.g_stat = 1; ion_props["N1+" ] = i;

  i.ion = "N2+"; i.ip1 = "N3+"; i.im1 = "N1+"; i.el  = "N";
  i.charge = 2; i.ion_pot = eV_per_Erg *47.5;
  i.i = N2p; i.g_stat = 1; ion_props["N2+" ] = i;

  i.ion = "N3+"; i.ip1 = "N4+"; i.im1 = "N2+"; i.el  = "N";
  i.charge = 3; i.ion_pot = eV_per_Erg *77.5;
  i.i = N3p; i.g_stat = 1; ion_props["N3+" ] = i;

  i.ion = "N4+"; i.ip1 = "N5+"; i.im1 = "N3+"; i.el  = "N";
  i.charge = 4; i.ion_pot = eV_per_Erg *97.9;
  i.i = N4p; i.g_stat = 1; ion_props["N4+" ] = i;

  i.ion = "N5+"; i.ip1 = "N6+"; i.im1 = "N+4"; i.el  = "N";
  i.charge = 5; i.ion_pot = eV_per_Erg *552.1;
  i.i = N5p; i.g_stat = 1; ion_props["N5+" ] = i;

  i.ion = "N6+"; i.ip1 = "N7+"; i.im1 = "N5+"; i.el  = "N";
  i.charge = 6; i.ion_pot = eV_per_Erg *667.0;
  i.i = N6p; i.g_stat = 1; ion_props["N6+" ] = i;
  
  i.ion = "N7+"; i.ip1 = ""; i.im1 = "N6+"; i.el  = "N";
  i.charge = 7; i.ion_pot = -1.e99;
  i.i = N7p;  i.g_stat = 1; ion_props["N7+" ] = i;
  
  // oxygen: O, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+, 
  i.ion = "O0"; i.ip1 = "O1+"; i.im1 = ""; i.el  = "O";
  i.charge = 0; i.ion_pot = eV_per_Erg *13.6;
  i.i = O0; i.g_stat = 1; ion_props["O0" ] = i;
  
  i.ion = "O1+"; i.ip1 = "O2+"; i.im1 = "O0"; i.el  = "O";
  i.charge = 1; i.ion_pot = eV_per_Erg *35.1;
  i.i = O1p; i.g_stat = 1; ion_props["O1+" ] = i;

  i.ion = "O2+"; i.ip1 = "O3+"; i.im1 = "O1+"; i.el  = "O";
  i.charge = 2; i.ion_pot = eV_per_Erg *54.9;
  i.i = O2p; i.g_stat = 1; ion_props["O2+" ] = i;

  i.ion = "O3+"; i.ip1 = "O4+"; i.im1 = "O2+"; i.el  = "O";
  i.charge = 3; i.ion_pot = eV_per_Erg *77.4;
  i.i = O3p;  i.g_stat = 1; ion_props["O3+" ] = i;
 
  i.ion = "O4+"; i.ip1 = "O5+"; i.im1 = "O3+"; i.el  = "O";
  i.charge = 4; i.ion_pot = eV_per_Erg *113.9;
  i.i = O4p; i.g_stat = 1; ion_props["O4+" ] = i;

  i.ion = "O5+"; i.ip1 = "O6+"; i.im1 = "O4+"; i.el  = "O";
  i.charge = 5; i.ion_pot = eV_per_Erg *138.1;
  i.i = O5p; i.g_stat = 1; ion_props["O5+" ] = i;

  i.ion = "O6+"; i.ip1 = "O7+"; i.im1 = "O5+"; i.el  = "O";
  i.charge = 6; i.ion_pot = eV_per_Erg *739.3;
  i.i = O6p; i.g_stat = 1; ion_props["O6+" ] = i;

  i.ion = "O7+"; i.ip1 = "O8+"; i.im1 = "O6+"; i.el  = "O";
  i.charge = 7; i.ion_pot = eV_per_Erg *871.4;
  i.i = O7p;  i.g_stat = 1; ion_props["O7+" ] = i;
 
  i.ion = "O8+"; i.ip1 = ""; i.im1 = "O7+"; i.el  = "O";
  i.charge = 8; i.ion_pot = -1.e99;
  i.i = O8p;  i.g_stat = 1; ion_props["O8+" ] = i;

  /*
    //ion
  i.ion = ""; i.ip1 = ""; i.im1 = ""; i.el  = "";
  i.charge = ; i.ion_pot = eV_per_Erg *;
  i.i = ;
  i.g_stat = 1; ion_props["" ] = i;
  */
  // etc. for other ions.
  return;
}

double MPv10::Coll_Ion_rate(
       double temperature,   ///< Precalculated Temperature.
       const struct ion_struct *ci ///< current ion.
       )
{
  // This uses fitting formulae from Voronov (1997) ADANDT, 65, 1.
  // (Atomic Data And Nuclear Data Tables)
  // rate returned in cm^3/s
  double A=0.0,X=0.0,K=0.0; int PP=0;
  if      (ci->i ==H_0)
    {if (temperature < 5.0e3) return 0.0; PP=0; A=2.91e-8; X=0.232; K=0.39;}
  else if (ci->i ==H_1p)
    {return 0.0;}
  else if (ci->i ==He0 )
    {if (temperature < 5.0e3) return 0.0; PP=0; A=1.75e-8; X=0.180; K=0.35;}
  else if (ci->i ==He1p)
    {if (temperature < 1.0e4) return 0.0; PP=1; A=2.05e-9; X=0.265; K=0.25;}
  else if (ci->i ==He2p)
    {return 0.0;}
  // carbon
  else if (ci->i ==C0) 
    {if (temperature < 3.0e2) return 0.0; PP=0; A=0.685e-7; X=0.193; K=0.25;}
  else if (ci->i ==C1p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.186e-7; X=0.286; K=0.24;}
  else if (ci->i ==C2p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.635e-8; X=0.427; K=0.21;}
  else if (ci->i ==C3p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.150e-8; X=0.416; K=0.13;}
  else if (ci->i ==C4p) 
    {if (temperature < 5.0e4) return 0.0; PP=1; A=0.299e-9; X=0.666; K=0.02;}
  else if (ci->i ==C5p) 
    {if (temperature < 5.0e4) return 0.0; PP=1; A=0.123e-9; X=0.620; K=0.16;}
  else if (ci->i ==C6p)
    {return 0.0;}
  // nitrogen
  else if (ci->i ==N0)
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.482e-7; X=.0652; K=0.42;}
  else if (ci->i ==N1p)
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.298e-7; X=0.310; K=0.30;}
  else if (ci->i ==N2p)
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.810e-8; X=0.350; K=0.24;}
  else if (ci->i ==N3p)
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.371e-8; X=0.549; K=0.18;}
  else if (ci->i ==N4p)
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.151e-8; X=.0167; K=0.74;}
  else if (ci->i ==N5p)
    {if (temperature < 1.0e5) return 0.0; PP=0; A=0.371e-9; X=0.546; K=0.29;}
  else if (ci->i ==N6p)
    {if (temperature < 1.0e5) return 0.0; PP=1; A=0.777e-10;X=0.624; K=0.16;}
  else if (ci->i ==N7p)
    {return 0.0;}
  // oxygen
  else if (ci->i ==O0) 
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.359e-7; X=0.073; K=0.34;}
  else if (ci->i ==O1p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.139e-7; X=0.212; K=0.22;}
  else if (ci->i ==O2p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.931e-8; X=0.270; K=0.27;}
  else if (ci->i ==O3p) 
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.102e-7; X=0.614; K=0.27;}
  else if (ci->i ==O4p) 
    {if (temperature < 1.0e3) return 0.0; PP=1; A=0.219e-8; X=0.630; K=0.17;}
  else if (ci->i ==O5p) 
    {if (temperature < 1.0e3) return 0.0; PP=0; A=0.195e-8; X=0.360; K=0.54;}
  else if (ci->i ==O6p) 
    {if (temperature < 1.0e5) return 0.0; PP=0; A=0.212e-9; X=0.396; K=0.35;}
  else if (ci->i ==O7p) 
    {if (temperature < 1.0e5) return 0.0; PP=1; A=0.521e-10;X=0.629; K=0.16;}
  else if (ci->i ==O8p)
    {return 0.0;}
  temperature = ci->ion_pot/k_B/temperature;
  return A*(1.+PP*sqrt(temperature))*exp(K*log(temperature) -temperature)/(X+temperature);

};


double MPv10::Rad_Recomb_rate(
      double temperature,   ///< Precalculated Temperature.
      const struct ion_struct *ci ///< current ion.
      )
{
  // rate is for recombination from current ion *ci to one stage lower.
  double r, a1,a2,a3,a4,a5,a6;
  if      (ci->i ==H_1p) {
    // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
    // For case B recombination coefficient.
    r = 3.41202e-10*exp(-0.782991*log(temperature));
  }
  else if (ci->i ==He1p) {
    // recombination from He(1+) (from Fabio... Raga, deColle, et al., 2007, A&A, 465, 879).
    //r  = 4.3e-13*exp(-0.672*log(temperature/1.e4));
    //r += 0.0019*exp(-1.5*log(temperature) -4.7e5/temperature)*(1.0+0.3*exp(-9.4e4/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1=9.356e-10; a2=0.7892; a3=4.266e-2; a4=4.677e6;
    r = a1/(sqrt(temperature/a3)*exp((1.-a2)*log(1.+sqrt(temperature/a3))+(1.+a2)*log(1.+sqrt(temperature/a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature,ci->i);
  }
  else if (ci->i ==He2p) {
    // recombination from He(2+) (from Fabio... Raga, deColle, et al., 2007, A&A, 465, 879).
    //r  = 2.21e-9*exp(-0.79*log(temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1=1.891e-10; a2=0.7524; a3=9.370; a4=2.774e6;
    r = a1/(sqrt(temperature/a3)*exp((1.-a2)*log(1.+sqrt(temperature/a3))+(1.+a2)*log(1.+sqrt(temperature/a4))));
  }
  // recombination for Carbon (from Raga, deColle, et al., 2007, A&A, 465, 879).
  else if (ci->i ==C1p) {
    a1=4.7e-13; a2=0.624; a3=6.9e-4; a4=1.1e5; a5=3.0; a6=4.9e4;
    //r = 1.0e-22; //a1*exp(a2*log(temperature) +a3/temperature);
    r = a1*exp(-a2*log(temperature/1.e4)) 
      +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
  }  
  else if (ci->i ==C2p) {
    a1=2.3e-12; a2=0.645; a3=0.007; a4=1.5e5; a5=0.5; a6=2.3e5;
    //    r = a1*exp(a2*log(temperature) +a3/temperature);
    r = a1*exp(-a2*log(temperature/1.e4)) 
      +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
  }  
  else if (ci->i ==C3p) {
    a1=3.2e-12; a2=0.770; a3=3.8e-3; a4=9.1e4; a5=2.0; a6=3.7e5;
    //    r = a1*exp(a2*log(temperature) +a3/temperature);
    r = a1*exp(-a2*log(temperature/1.e4)) 
      +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
  }  
  else if (ci->i ==C4p) {
    //a1=; a2=; a3=; a4=; a5=; a6=;
    //r = a1*exp(a2*log(temperature) +a3/temperature);
    //r = a1*exp(-a2*log(temperature/1.e4)) 
    //  +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1=8.540e-11; a2=0.5247; a3=5.014e2; a4=1.479e7;// a5=; a6=;
    r = a1/(sqrt(temperature/a3)*exp((1.-a2)*log(1.+sqrt(temperature/a3))+(1.+a2)*log(1.+sqrt(temperature/a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature,ci->i);    
  }  
  else if (ci->i ==C5p) {
    //a1=; a2=; a3=; a4=; a5=; a6=;
    //r = a1*exp(a2*log(temperature) +a3/temperature);
    //r = a1*exp(-a2*log(temperature/1.e4)) 
    //  +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1=2.765e-10; a2=0.6858; a3=1.535e2; a4=2.556e7;// a5=; a6=;
    r = a1/(sqrt(temperature/a3)*exp((1.-a2)*log(1.+sqrt(temperature/a3))+(1.+a2)*log(1.+sqrt(temperature/a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature,ci->i);
  }  
  else if (ci->i ==C6p) {
    //a1=; a2=; a3=; a4=; a5=; a6=;
    //r = a1*exp(a2*log(temperature) +a3/temperature);
    //r = a1*exp(-a2*log(temperature/1.e4)) 
    //  +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1=6.556e-10; a2=0.7567; a3=6.523e1; a4=2.446e7;// a5=; a6=;
    r = a1/(sqrt(temperature/a3)*exp((1.-a2)*log(1.+sqrt(temperature/a3))+(1.+a2)*log(1.+sqrt(temperature/a4))));
  }  
  // recombination for Nitrogen (from Raga, deColle, et al., 2007, A&A, 465, 879).
  else if (ci->i == N1p) {
    a1= 1.5e-12; a2=0.693; a3=0.0031; a4=2.9e5; a5=0.6; a6=1.6e5;
    r = a1*exp(-a2*log(temperature/1.e4)) 
      +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
  }
  else if (ci->i == N2p) {
    a1= 4.4e-12; a2=0.675; a3=0.0075; a4=2.6e5; a5=0.7; a6=4.5e5;
    r = a1*exp(-a2*log(temperature/1.e4)) 
      +a3*exp(-1.5*log(temperature) -a4/temperature)*(1.+a5*exp(-a6/temperature));
  }

  else {
    return 0.0;
  }
  //else rep.error("unknown ion in Rad_Recomb_rate()",ii);
  return r;
}


// ##################################################################
// ##################################################################

double MPv10::rad_recomb(
      double T,
      enum species i
      )
{
  // rates are from Verner & Ferland (1996) ApJS, 103, 467 where possible,
  // and from Pequignot, Petitjean, & Boisson, (1991) A&A, 251, 680 otherwise,
  // refitted with VF96 formula.
  //
  // Fortran subroutine with list of rates at:
  // http://hea-www.harvard.edu/~mazzotta/sub_mazz2/rrfit.f
  // ftp://gradj.pa.uky.edu//dima//rec//rrfit.f via http://www.pa.uky.edu/~verner/fortran.html
  double a1=0.0,a2=0.0,a3=0.0,a4=0.0,r=0.0;
  if      (i ==H_1p) {
    // For hydrogen, we calculate the VF96 case A recomb. coeff., and
    // then use the minimum of this and the Storey & Hummer (1995) coeff.
    //a1=7.982e-11; a2=0.7480; a3=3.148; a4=7.036e5;
    //r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
    //    if (T>1.e3)
    //cout <<"T="<<T<<"\t r(VF)="<<r<<"\t r(SH)="<<3.41202e-10*exp(-0.782991*log(T))<<"\n";
    //r = min(r,3.41202e-10*exp(-0.782991*log(T)));

    // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
    // For case B recombination coefficient.
    r = 3.41202e-10*exp(-0.782991*log(T));
    //r = 2.59e-13; // basic approx.
  }
  else if (i ==He1p) {
    a1=9.356e-10; a2=0.7892; a3=4.266e-2; a4=4.677e6;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }
  else if (i ==He2p) {
    a1=1.891e-10; a2=0.7524; a3=9.370; a4=2.774e6;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }
  else if (i ==C1p) {
    // 7.651E-09,0.8027,1.193E-03,9.334E+12
    a1=7.651E-09; a2=0.8027; a3=1.193E-03; a4=9.334E+12;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  else if (i ==C2p) {
    // 8.577E-10,0.7837,7.286E-01,1.140E+07
    a1=8.577E-10; a2=0.7837; a3=7.286E-01; a4=1.140E+07;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  else if (i ==C3p) {
    // 2.020E-09,0.7798,6.690E-01,2.425E+06
    a1=2.020E-09; a2=0.7798; a3=6.690E-01; a4=2.425E+06;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  else if (i ==C4p) {
    // 8.540e-11,0.5247,5.014e+02,1.479e+07
    a1=8.540e-11; a2=0.5247; a3=5.014e2; a4=1.479e7;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  else if (i ==C5p) {
    // 2.765e-10,0.6858,1.535e+02,2.556e+07
    a1=2.765e-10; a2=0.6858; a3=1.535e2; a4=2.556e7;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  else if (i ==C6p) {
    // 6.556e-10,0.7567,6.523e+01,2.446e+07
    a1=6.556e-10; a2=0.7567; a3=6.523e1; a4=2.446e7;
    r = a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
  }  
  //else rep.error("unknown ion in Rad_Recomb_rate()",ii);
  return r;
}


// ##################################################################
// ##################################################################



double MPv10::dielec_recomb(
      double T,
      enum species i
      )
{
  // i is the ion we are recombining *from*
  // rates are from Mazzotta et al. (1998)
  double c[4]; c[0]=c[1]=c[2]=c[3]=0.0;
  double E[4]; E[0]=E[1]=E[2]=E[3]=0.0;
  if      (i==H_1p) {
    return 0.0;
  }
  else if (i==He1p) {
    c[0] = 1.12e-9; E[0] = 39.70;
  }
  else if (i==He2p) {
    return 0.0;
  }
  else if (i==C1p) {
    c[0]=1.0422e-09; c[1]=5.8484e-10; c[2]=5.6306e-11;
    E[0]=12.57;      E[1]=162.90;     E[2]=6.30;
  }
  else if (i==C2p) {
    c[0]=4.6178e-12; c[1]=2.8234e-09; E[0]=0.49; E[1]=11.78;
  }
  else if (i==C3p) {
    c[0]=5.3858e-10; c[1]=1.5056e-09; c[2]=6.0332e-10;
    E[0]=185.96;     E[1]=17.99;      E[2]=2.41;
  }
  else if (i==C4p) {
    c[0]=1.4008e-08; E[0]=287.34;
  }
  else if (i==C5p) {
    c[0]=3.3558e-08; E[0]=356.46;
  }
  else if (i==C6p) {
    return 0.0;
  }
  //else if (i=="") {
  //}
  else {
    cout <<"unknown ion: "<<i<<"\n";
    return -1.0;
  }
  T /= 1.16e4; // convert from K to eV.
  //  if (T<=1.e-100) {cout <<"tiny temperature\n";return -1.0;}
  
  double rate=0.0;
  for (int i=0;i<4;i++) rate += c[i]*exp(-E[i]/T);
  rate *= exp(-1.5*log(T));
  return rate;
}



//#############################################################################
//#############################################################################


void MPv10::generate_lookup_tables(){ 
  set_atomic_data();   
  
  //  Start by generating the logarithmic temperature scale:
  const int Num_temps = 1e2; //NB this needs to be const so can initialise arrays with it. If this is >=1e4, get stack overflow errors.
  delta_log_temp = (log10(T_max) - log10(T_min))/(Num_temps-1);
  double Temp_arr[Num_temps] = {};
  
  for (int i=0; i < Num_temps; i++) {
    float T_i = pow (10.0,  (log10f(T_min) + i* delta_log_temp));
    Temp_arr[i] = T_i;
  }
  Temp_Table.insert(Temp_Table.end(), &Temp_arr[0], &Temp_arr[Num_temps]);
    
  //  Now make lookup tables for the recombination rates and ionisation rates respectively.
  const int number_of_species = 29;
  
  species species_list[number_of_species] = {H_0, H_1p, 
                                             He0, He1p, He2p, 
                                             C0, C1p, C2p, C3p, C4p, C5p, C6p, 
                                             N0, N1p, N2p, N3p, N4p, N5p, N6p, N7p, 
                                             O0, O1p, O2p, O3p, O4p, O5p, O6p, O7p, O8p};
  string ion_names[number_of_species] = {"H0", "H1+", 
                                         "He0", "He1+", "He2+", 
                                         "C0", "C1+", "C2+", "C3+", "C4+", "C5+", "C6+", 
                                         "N0", "N1+", "N2+", "N3+", "N4+", "N5+", "N6+", "N7+", 
                                         "O0", "O1+", "O2+", "O3+", "O4+", "O5+", "O6+", "O7+", "O8+"};
  
  // ======================================================================================================
  //  Resize header-defined vectors to store data
  // ======================================================================================================
  recomb_rate_table.resize(number_of_species);
  for (int i = 0; i < number_of_species; ++i) recomb_rate_table[i].resize(Num_temps);
  
  recomb_slope_table.resize(number_of_species);
  for (int i = 0; i < number_of_species; ++i) recomb_slope_table[i].resize(Num_temps);
  
  ionise_rate_table.resize( number_of_species);
  for (int i = 0; i < number_of_species; ++i) ionise_rate_table[i].resize(Num_temps);
  
  ionise_slope_table.resize(number_of_species);
  for (int i = 0; i < number_of_species; ++i) ionise_slope_table[i].resize(Num_temps);
  
  
  //  Calculate recombination and ionisation rates
  for ( int s=0; s<number_of_species; s++ ) {
    for ( int T_i=0; T_i<Num_temps; T_i++){
      float this_recombination_rate;
      float this_ionisation_rate;
      float this_temperature = Temp_Table[T_i];
    
      ion_struct this_ion;
            
      this_ion.i = species_list[s];
      this_ion.ion_pot =  ion_props[ion_names[s]].ion_pot;
            
      this_recombination_rate = Rad_Recomb_rate(this_temperature, &this_ion);
      this_ionisation_rate = Coll_Ion_rate(this_temperature, &this_ion);
            
      recomb_rate_table[s][T_i] = this_recombination_rate;
      ionise_rate_table[s][ T_i] = this_ionisation_rate;
    }
  }
    
  //  Now, generate the slopes tables for interpolating temperatures that don't exactly fit here.
  for ( int s=0; s<number_of_species; s++ ) {
    for ( int T_i=0; T_i<Num_temps-1; T_i++){ //note this is <Num_temps -1, b/c i+1
      float this_recombination_slope;
      float this_ionisation_slope;
            
      this_recombination_slope =  ( recomb_rate_table[s][T_i+1] - recomb_rate_table[s][T_i] )
                                  / ( Temp_Table[T_i+1] - Temp_Table[T_i] );
                                      
      this_ionisation_slope =     ( ionise_rate_table[s][T_i+1] - ionise_rate_table[s][T_i] )
                                  / ( Temp_Table[T_i+1] - Temp_Table[T_i] );
            
      recomb_slope_table[s][T_i] = this_recombination_slope;
      ionise_slope_table[s][T_i] = this_ionisation_slope;
    }
  }
  
  
  double erg_per_eV = 1.60218e-12;
  /// ===================================================================
  /// Initialise ionisation potentials vector
  /// ===================================================================
  double ionisation_pot_arr[number_of_species] = {13.59844, -1.0e99,
                                                  24.58741, 54.41778, -1.0e99,
                                                  11.3, 24.4, 47.9, 64.5, 392.1, 490.0, -1.0e99,
                                                  14.5, 29.6, 47.5, 77.5, 97.9, 552.1, 667.0, -1.e99,
                                                  13.6, 35.1, 54.9, 77.4, 113.9, 138.1, 739.3, 871.4, 1.e99}; //energy (eV) required to raise ion from species i to species i+1
  for (int i=0; i<number_of_species; i++) ionisation_pot_arr[i]*=erg_per_eV; //convert eV to erg
  ionisation_potentials.insert(ionisation_potentials.end(), &ionisation_pot_arr[0], &ionisation_pot_arr[number_of_species]);
 
  return;
};



