#include <iostream>
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include "Riemann_solvers/HLLD_MHD.h"
#include "Riemann_solvers/riemannMHD.h"

 
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING
 
#include "microphysics/microphysics_base.h"
//#include "eqns_mhd_adiabatic.h"
#include <iostream>
#include <cmath>
using namespace std;



// ##################################################################
// ##################################################################


// Constructor
HLLD_MHD::HLLD_MHD(
      const int nv,      ///< Length of State Vectors, nvar
      const double g    ///< Gamma for state vector.
      )
 :  eqns_base(nv), eqns_mhd_ideal(nv)
{
  cout <<"(HLLD_MHD::HLLD_MHD) Initialising HLL Multi-State Solver Class.\n";
  if(eq_nvar<8){
    rep.error("#elements!=8, QUIT.");
  }
  HD_nvar  = 7; // Case for Bx=const
  eq_gamma = g;

  HDl_lambda     = mem.myalloc(HDl_lambda,     HD_nvar);
  HDr_lambda     = mem.myalloc(HDr_lambda,     HD_nvar);
  //HD_evalues     = mem.myalloc(HD_evalues,     HD_nvar);
  //HD_strengths   = mem.myalloc(HD_strengths,   HD_nvar);
  //HD_right_evecs = mem.myalloc(HD_right_evecs, HD_nvar);
  //for (int v=0;v<HD_Nvar;v++) 
  //  HD_right_evecs[v] = mem.myalloc(HD_right_evecs[v], HD_nvar);

  //HD_meanp = mem.myalloc(HD_meanp, eq_nvar);
  HD_UL   = mem.myalloc(HD_UL,   eq_nvar); // conserved
  HD_UR   = mem.myalloc(HD_UR,   eq_nvar);

  HD_FL   = mem.myalloc(HD_FL,   eq_nvar); // flux
  HD_FR   = mem.myalloc(HD_FR,   eq_nvar);
  
  HD_ULs  = mem.myalloc(HD_ULs,  eq_nvar); // U*
  HD_URs  = mem.myalloc(HD_URs,  eq_nvar);

  HD_FLs  = mem.myalloc(HD_FLs,  eq_nvar); // F*
  HD_FRs  = mem.myalloc(HD_FRs,  eq_nvar); 

  HD_ULss = mem.myalloc(HD_ULss, eq_nvar); // U**
  HD_URss = mem.myalloc(HD_URss, eq_nvar);

  HD_FLss = mem.myalloc(HD_FLss, eq_nvar); // F**
  HD_FRss = mem.myalloc(HD_FRss, eq_nvar); 
  

  cout << "(HLLD_MHD::HLLD_MHD) All set.\n"
  return;
}



// ##################################################################
// ##################################################################

// Destructor
HLLD_MHD::~HLLD_MHD()
{

  cout << "(riemann_MHD::riemann_MHD) Commencing Destruction." << "\n";
  HDl_lambda     = mem.myfree(HDl_lambda);
  HDr_lambda     = mem.myfree(HDr_lambda);
  //HD_evalue   = mem.myfree(HD_evalue);
  //HD_strength = mem.myfree(HD_strength);
  //for (int v=0;v<7;v++) 
  //  HD_right_evecs[v] = mem.myfree(HD_right_evecs[v]);

  //HD_meanp = mem.myfree(HD_meanp);
  HD_UL  = mem.myfree(HD_left); 	// conserved
  HD_UR  = mem.myfree(HD_right);

  HD_FL   = mem.myfree(HD_FL); 		// flux
  HD_FR   = mem.myfree(HD_FR);

  HD_ULs  = mem.myfree(HD_ULs); 	// U*
  HD_URs  = mem.myfree(HD_URs);

  HD_FLs  = mem.myfree(HD_FLs); 	// F*
  HD_FRs  = mem.myfree(HD_FRs); 

  HD_ULss  = mem.myfree(HD_ULss); 	// U**
  HD_URss  = mem.myfree(HD_URss);

  HD_FLss = mem.myfree(HD_FLss);	// F**
  HD_FRss = mem.myfree(HD_FRss); 

  cout << "(riemann_MHD::riemann_MHD) Destruction Accomplished." << "\n";
  return;
}



// ##################################################################
// ##################################################################



int HLLD_MHD::MHD_HLLD_flux_solver(
      const pion_flt *left,  ///< input left state
      const pion_flt *right, ///< input right state
      const double gamma,    ///< input gamma

      //const pion_flt etamax, ///< H-correction eta-max value.
      //pion_flt *out_ps,       ///< output p*
      //pion_flt *out_pss,       ///< output p**
      pion_flt *out_flux         ///< output flux
      )
{



  PtoU(left ,HD_UL,eq_gamma);
  PtoU(right,HD_UR,eq_gamma);


  PUtoFlux(left,HD_UL,HD_FL);
  PUtoFlux(right,HD_UR,HD_FR);
  //
  // compute right and left wave speeds
  // 
  double gp_l = gamma * left[eqPG];
  double BB_l = pow(left[eqBX],2)   + pow(left[eqBY],2)         + pow(left[eqBZ],2);
  double cf_l = sqrt((gp_l + BB_l   + sqrt(pow((gp_l + BB_l),2) - 4*gp_l * pow(left[eqBX],2)))/(2*left[eqRO]));
  double cs_l = sqrt((gp_l + BB_l   - sqrt(pow((gp_l + BB_l),2) - 4*gp_l * pow(left[eqBX],2)))/(2*left[eqRO]));
  double ca_l = abs(left[eqBX])/sqrt(left[eqRO]);


  HDl_lambda[0] = left[eqVX] - cf_l;
  HDl_lambda[1] = left[eqVX] - ca_l;
  HDl_lambda[2] = left[eqVX] - cs_l;
  HDl_lambda[3] = left[eqVX];
  HDl_lambda[0] = left[eqVX] + cf_l;
  HDl_lambda[1] = left[eqVX] + ca_l;
  HDl_lambda[2] = left[eqVX] + cs_l;

  double gp_r = gamma * right[eqPG];
  double BB_r = pow(right[eqBX],2)  + pow(right[eqBY],2)        + pow(right[eqBZ],2);
  double cf_r = sqrt((gp_r + BB_r   + sqrt(pow((gp_r + BB_r),2) - 4*gp_r * pow(right[eqBX],2)))/(2*right[eqRO]));
  double cs_r = sqrt((gp_r + BB_r   - sqrt(pow((gp_r + BB_r),2) - 4*gp_r * pow(right[eqBX],2)))/(2*right[eqRO]));
  double ca_r = abs(right[eqBX])/sqrt(right[eqRO]);


  HDr_lambda[0] = right[eqVX] - cf_r;
  HDr_lambda[1] = right[eqVX] - ca_r;
  HDr_lambda[2] = right[eqVX] - cs_r;
  HDr_lambda[3] = right[eqVX];
  HDr_lambda[0] = right[eqVX] + cf_r;
  HDr_lambda[1] = right[eqVX] + ca_r;
  HDr_lambda[2] = right[eqVX] + cs_r;

  double cf_max = max(cf_l,cf_r);


  double tp_r = right[eqPG] + 0.5 * BB_r; // total pressure
  double tp_l = left[eqPG] + 0.5 * BB_l; // total pressure

  double S_l = min(left[eqVX],right[eqVX]) - cf_max; // m05 eq 67
  double S_r = max(left[eqVX],right[eqVX]) - cf_max;

  double temp = (S_r - right[eqVX]) * right[eqRO]  - (S_l - left[eqVX]) * left[eqRO];
  double temp_l = S_l - left[eqVX];
  double temp_r = S_r - right[eqVX];

  double S_m = (temp_r  * HD_UR[eqMMx] - temp_l * HD_UL[eqMMx] - tp_r + tp_l) / temp; // Sm (m05 eq 38)
  double tp_s = (temp_r * right[eqRO]  * tp_l   - temp_l * left[eqRO] * tp_r  + left[eqRO] 
			* right[eqRO]  * temp_r * temp_l)/temp// total pressure* (m05 eq 41)
 

  //
  // calculate intermediate states
  //
  HD_ULs[eqRHO] = left[eqRO]  * (S_l - left[eqVX])/(S_l - S_m); // d* (m05 eq 43)
  HD_URs[eqRHO] = right[eqRO] * (S_r - left[eqVX])/(S_r - S_m);

  HD_ULs[eqMMX] = S_m * HD_ULs[eqRHO];    // Mx* (m05 eq 39)
  HD_URs[eqMMx] = S_m * right[eqRO];

  HD_ULs[eqBBX] = left[eqBBX];
  HD_URs[eqBBX] = right[eqBBX];


  double temp_l1 = left[eqRO]  * (S_l - left[eqVX])  * (S_l - S_m) - pow(left[eqBX],2);
  double temp_r1 = right[eqRO] * (S_r - right[eqVX]) * (S_r - S_m) - pow(right[eqBX],2);
  double temp_l2 = S_m - left[eqVX];
  double temp_r2 = S_m - right[eqVX];
  temp_l  = temp_l2/temp_l1;
  temp_r  = temp_r2/temp_r1;

  double vy_sl  = left[eqVY]  - left[eqBX]  * left[eqBY]  * temp_l;
  double vy_sr  = right[eqVY] - right[eqBX] * right[eqBY] * temp_r;
  HD_ULs[eqMMY] = vy_sl * HD_ULs[eqRHO]; // Vy* (m05 eq. 44)
  HD_URs[eqMMY] = vy_sr * HD_URs[eqRHO]; 

  double vz_sl  = left[eqVZ]  - left[eqBX]  * left[eqBZ]  * temp_l;
  double vz_sr  = right[eqVZ] - right[eqBX] * right[eqBZ] * temp_r;
  HD_ULs[eqMMZ] = vz_sl * HD_ULs[eqRHO]; // Vz* (m05 eq. 46)
  HD_URs[eqMMZ] = vz_sr * HD_URs[eqRHO]; 

  temp_l2 = left[eqRO]  * pow((S_l - left[eqVX]),2)  - pow(left[eqBx],2);
  temp_l2 = right[eqRO] * pow((S_r - right[eqVX]),2) - pow(right[eqBx],2);
	
  HD_ULs[eqBBY] = left[eqBY]  * temp_l; // By* (m05 eq. 45)
  HD_URs[eqBBY] = right[eqBY] * temp_r; 

  HD_ULs[eqBBZ] = left[eqBZ]  * temp_l; // Bz* (m05 eq. 47)
  HD_URs[eqBBZ] = right[eqBZ] * temp_r; 


  temp_l1 = left[eqVX]  * left[eqBX]  + left[eqVY]  * left[eqBY]  + left[eqVZ]  * left[eqBZ]; 
  temp_r1 = right[eqVX] * right[eqBX] + right[eqVY] * right[eqBY] + right[eqVZ] * right[eqBZ]; 
  temp_l2 = S_m * HD_ULs[eqBBX] + vy_sl * HD_ULs[eqBBY] + vz_sl * HD_ULs[eqBBZ];
  temp_r2 = S_m * HD_URs[eqBBX] + vy_sr * HD_URs[eqBBY] + vz_sr * HD_URs[eqBBZ];
  temp_l  = temp_l1 - temp_l2;
  temp_r  = temp_r1 - temp_r2;

  HD_ULs[eqERG] = ((S_l - left[eqVX])  * HD_UL[eqERG] - tp_l * left[eqVX]  + tp_s * S_m + left[eqBx]  * temp_l)/(S_l - S_m);
  HD_URs[eqERG] = ((S_r - right[eqVX]) * HD_UR[eqERG] - tp_r * right[eqVX] + tp_s * S_m + right[eqBx] * temp_r)/(S_r - S_m);
  // e* (m05 eq 48)


  double SS_l = S_m - abs(Bx)/np.sqrt(HD_ULs[eqRHO]);  // m05 eq 51
  double SS_r = S_m + abs(Bx)/np.sqrt(HD_URs[eqRHO]);

  double sgn =  (left[eqBX] > 0) - (left[eqBX] < 0);  // signum
  temp_l = sqrt(HD_ULs[eqRHO]);
  temp_r = sqrt(HD_URs[eqRHO]);
  temp   = temp_l + temp_r;

  HD_ULss[eqRHO] = HD_ULs[eqRHO];
  HD_URss[eqRHO] = HD_URs[eqRHO];

  HD_ULss[eqMMX] = S_m * HD_ULss[eqRHO] ; // Vx** eq 39
  HD_URss[eqMMX] = S_m * HD_URss[eqRHO] ; 

  double vy_ss   = temp_l * vy_sl + temp_r * vy_sr + (HD_URs[eqBBY] - left[eqBY]) * sgn)/temp; 
  HD_ULss[eqMMY] = vy_ss * HD_ULss[eqRHO]; // Vy** eq 59
  HD_URss[eqMMY] = vy_ss * HD_URss[eqRHO]; 

  double vz_ss   = temp_l * vz_sl + temp_r * vz_sr + (HD_URs[eqBBz] - left[eqBZ]) * sgn)/temp;
  HD_ULss[eqMMZ] = vz_ss  * HD_ULss[eqRHO]; // Vz** eq 60
  HD_URss[eqMMZ] = vz_ss  * HD_URss[eqRHO]; 

  HD_ULss[eqBBX] = HD_URss[eqBBX] = left[eqBX];
  HD_ULss[eqBBY] = HD_URss[eqBBY] = (temp_l * HD_URs[eqBBY] + temp_r * left[eqBY] + temp_l * temp_r * (vy_sr - vy_sl) * sgn)/temp; //By** eq 61
  HD_ULss[eqBBZ] = HD_URss[eqBBZ] = (temp_l * HD_URs[eqBBz] + temp_r * HD_ULs[eqBBz] + temp_l * temp_r * (vz_sr - vz_sl) * sgn)/temp; //Bz** eq 62

  temp = S_m * HD_ULss[eqBBX] + vy_ss * HD_ULss[eqBBY] + vz_ss * HD_ULss[eqBBZ];

  HD_ULss[eqERG] = HD_ULs[eqERG] - temp_l * (temp_l2 - temp) * sgn // E** eq 63
  HD_URss[eqERG] = HD_URs[eqERG] + temp_r * (temp_r2 - temp) * sgn



  // Fa*, Fa** from eq. 64,65
  for (int v=0; v<eq_nvar; v++){
  HD_FLs[v] = HD_FL[v] + S_l + HD_ULs[v] - S_l * HD_UL[v];
  HD_FRs[v] = HD_FR[v] + S_r + HD_URs[v] - S_r * HD_UR[v];

  HD_FLss[v] = HD_FLs[v] + SS_l * HD_ULss[v] - (SS_l - S_l) * HD_ULss[v] - S_l * HD_UL[v];
  HD_FRss[v] = HD_FRs[v] + SS_r * HD_URss[v] - (SS_r - S_r) * HD_URss[v] - S_r * HD_UR[v];
  }	

  for (int v=0; v<eq_nvar; v++){
    if     (S_l>0)    out_flux[v] = HD_FL[v]; 
    else if (SS_l>=0) out_flux[v] = HD_FLs[v]; 
    else if  (S_m>=0) out_flux[v] = HD_FLss[v];
    else if (SS_r>=0) out_flux[v] = HD_FRss[v];
    else if  (S_r>=0) out_flux[v] = HD_FRs[v];
    else:             out_flux[v] = HD_FR[v];


  // include sanity checks?




  return err;



  return 0;
}


// ###################################################################
// ###################################################################




