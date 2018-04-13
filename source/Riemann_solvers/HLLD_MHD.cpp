//    	references:
//	* Stone et al. 2008 	(s08)
//	* Miyoshi and Kusano 2005 (m05)
//	* Migone and Bodo 2008 	(mb08)


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
#include "equations/eqns_mhd_adiabatic.h"
#include <iostream>
#include <cmath>
using namespace std;



// ##################################################################
// ##################################################################


// Constructor
HLLD_MHD::HLLD_MHD(
      const int nv,     ///< Length of State Vectors, nvar
      const double g    ///< Gamma for state vector.
      // const double Tmin
      )
 :  eqns_base(nv), eqns_mhd_ideal(nv)
{
  cout <<"(HLLD_MHD::HLLD_MHD) Initialising HLL Multi-State Solver Class.\n";
  if(eq_nvar<8){
    rep.error("#elements!=8, QUIT.",eq_nvar);
  }
  HD_nvar  = 7; 
  eq_gamma = g;

  HD_lambda = mem.myalloc(HD_lambda, 5);   // wave speeds (one entropy, two Alfvén)

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
  

  cout << "(HLLD_MHD::HLLD_MHD) All set.\n";
  return;
}



// ##################################################################
// ##################################################################

// Destructor
HLLD_MHD::~HLLD_MHD()
{

  cout << "(riemann_MHD::riemann_MHD) Commencing Destruction." << "\n";

  HD_lambda = mem.myfree(HD_lambda); // wave speeds (one entropy, two Alfvén)

  HD_UL  = mem.myfree(HD_UL); 	     // conserved
  HD_UR  = mem.myfree(HD_UR);

  HD_FL   = mem.myfree(HD_FL); 	     // flux
  HD_FR   = mem.myfree(HD_FR);

  HD_ULs  = mem.myfree(HD_ULs);      // U*
  HD_URs  = mem.myfree(HD_URs);

  HD_FLs  = mem.myfree(HD_FLs);      // F*
  HD_FRs  = mem.myfree(HD_FRs); 

  HD_ULss  = mem.myfree(HD_ULss);    // U**
  HD_URss  = mem.myfree(HD_URss);

  HD_FLss = mem.myfree(HD_FLss);     // F**
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
      //pion_flt *out_p,         ///< output p*
      pion_flt *out_flux       ///< output flux
      )
{

  //
  // compute conserved U and Flux (m05 eq 3)
  // 
  double BX = 0.5*(left[eqBX]+right[eqBX]); // Bx is constant (Should be mean of left and right state)

  PtoU(left ,HD_UL,eq_gamma);
  PtoU(right,HD_UR,eq_gamma);

  PUtoFlux(left,HD_UL,HD_FL);
  PUtoFlux(right,HD_UR,HD_FR);

  //
  // compute wave speeds (m05 eq 3)
  // 
  double gp_l = gamma * left[eqPG];
  double  BB_l = pow(BX,2) + pow(left[eqBY],2)  + pow(left[eqBZ],2);
  double  temp = gp_l + BB_l;
  double  cf_l = cfast(left,gamma); //sqrt((temp + sqrt(pow(temp,2) - 4*gp_l * pow(BX,2)))/(2*left[eqRO]));

  double gp_r = gamma * right[eqPG];
  double  BB_r = pow(BX,2) + pow(right[eqBY],2) + pow(right[eqBZ],2);
  temp = gp_r + BB_r;
  double  cf_r = cfast(right,gamma);  //sqrt((temp + sqrt(pow(temp,2) - 4*gp_r * pow(BX,2)))/(2*right[eqRO]));

  double cf_max = max(cf_l,cf_r);
 
  HD_lambda[0] = min(left[eqVX],right[eqVX]) - cf_max; // S_l (m05 eq 67)
  HD_lambda[4] = max(left[eqVX],right[eqVX]) + cf_max; // S_r

  double sl_vl = HD_lambda[0] - left[eqVX];
  double sr_vr = HD_lambda[4] - right[eqVX];
  double tp_r  = Ptot(right,gamma); //right[eqPG]  + 0.5 * BB_r; // total pressure
  double tp_l  = Ptot(left, gamma); //left[eqPG]   + 0.5 * BB_l;
  temp = sr_vr  * right[eqRO]  - sl_vl * left[eqRO];

  HD_lambda[2] = (sr_vr * HD_UR[eqMMX] - sl_vl * HD_UL[eqMMX] - tp_r + tp_l)/temp; // S_m (m05 eq 38)
  double tp_s  = (sr_vr * right[eqRO]  * tp_l  - sl_vl * left[eqRO] * tp_r  
		 + left[eqRO] * right[eqRO] * sr_vr * sl_vl * (right[eqVX] - left[eqVX]))/temp; // total pressure* (m05 eq 41)
  
  //
  // calculate intermediate states
  //
  double sl_sm = HD_lambda[0] - HD_lambda[2];
  double sr_sm = HD_lambda[4] - HD_lambda[2];


  HD_ULs[eqRHO] = left[eqRO]  * sl_vl/sl_sm; // d* (m05 eq 43)
  HD_URs[eqRHO] = right[eqRO] * sr_vr/sr_sm;

  HD_ULs[eqMMX] = HD_lambda[2] * HD_ULs[eqRHO]; // Vx* (m05 eq 39)
  HD_URs[eqMMX] = HD_lambda[2] * HD_URs[eqRHO];


  double temp_l1 = HD_lambda[2] - left[eqVX];
  double temp_l2 = left[eqRO]   * sl_vl * sl_sm - pow(BX,2);
  double temp_r1 = HD_lambda[2] - right[eqVX];
  double temp_r2 = right[eqRO] * sr_vr * sr_sm - pow(BX,2);

  double vys_l  = left[eqVY]  - BX * left[eqBY]  * temp_l1/temp_l2; // Vy* (m05 eq. 44)
  double  vys_r  = right[eqVY] - BX * right[eqBY] * temp_r1/temp_r2;
  
  HD_ULs[eqMMY] = vys_l * HD_ULs[eqRHO]; 
  HD_URs[eqMMY] = vys_r * HD_URs[eqRHO]; 

  double vzs_l  = left[eqVZ]  - BX  * left[eqBZ]  * temp_l1/temp_l2; // Vz* (m05 eq. 46)
  double  vzs_r  = right[eqVZ] - BX  * right[eqBZ] * temp_r1/temp_r2;
  
  HD_ULs[eqMMZ] = vzs_l * HD_ULs[eqRHO]; 
  HD_URs[eqMMZ] = vzs_r * HD_URs[eqRHO]; 


  HD_ULs[eqBBX] = HD_URs[eqBBX] = BX;

  temp_l1 = left[eqRO]  * pow(sl_vl,2) - pow(BX,2);
  temp_r1 = right[eqRO] * pow(sr_vr,2) - pow(BX,2);
	
  HD_ULs[eqBBY] = left[eqBY]  * temp_l1/temp_l2; // By* (m05 eq. 45)
  HD_URs[eqBBY] = right[eqBY] * temp_r1/temp_r2; 

  HD_ULs[eqBBZ] = left[eqBZ]  * temp_l1/temp_l2; // Bz* (m05 eq. 47)
  HD_URs[eqBBZ] = right[eqBZ] * temp_r1/temp_r2; 


  temp_l1 = left[eqVX]  * BX + left[eqVY]  * left[eqBY]  + left[eqVZ]  * left[eqBZ]; // B dot v
  temp_r1 = right[eqVX] * BX + right[eqVY] * right[eqBY] + right[eqVZ] * right[eqBZ]; 
  temp_l2 = HD_lambda[2] * HD_ULs[eqBBX] + vys_l * HD_ULs[eqBBY] + vzs_l * HD_ULs[eqBBZ]; // B* dot v*
  temp_r2 = HD_lambda[2] * HD_URs[eqBBX] + vys_r * HD_URs[eqBBY] + vzs_r * HD_URs[eqBBZ];


  HD_ULs[eqERG] = (sl_vl * HD_UL[eqERG] - tp_l * left[eqVX]  + tp_s * HD_lambda[2] + BX  
			 * (temp_l1 - temp_l2))/sl_sm;  // e* (m05 eq 48)
  HD_URs[eqERG] = (sr_vr * HD_UR[eqERG] - tp_r * right[eqVX] + tp_s * HD_lambda[2] + BX 
			 * (temp_r1 - temp_r2))/sr_sm;

  HD_lambda[1] = HD_lambda[2] - abs(BX)/sqrt(HD_ULs[eqRHO]);  // S_l* (m05 eq 51)
  HD_lambda[3] = HD_lambda[2] + abs(BX)/sqrt(HD_URs[eqRHO]);  // S_r*

  //rep.printVec("lambda",HD_lambda,5);


  // ADD CASE BX=0
  if (BX==0)
    for (int v=0; v<eq_nvar; v++) {
      HD_ULss[v] = HD_ULs[v];
      HD_URss[v] = HD_URs[v];
    }
  else {
    HD_ULss[eqRHO] = HD_ULs[eqRHO]; // rho** (m05 eq 49)
    HD_URss[eqRHO] = HD_URs[eqRHO];

    double sgn =  (BX > 0) - (BX < 0);  // signum
    temp_l1 = sqrt(HD_ULs[eqRHO]);
    temp_r1 = sqrt(HD_URs[eqRHO]);
    temp    = temp_l1 + temp_r1;

    HD_ULss[eqMMX] = HD_lambda[2] * HD_ULss[eqRHO] ; // Vx** (m05 eq 39)
    HD_URss[eqMMX] = HD_lambda[2] * HD_URss[eqRHO] ; 

    double vy_ss   = (temp_l1 * vys_l + temp_r1 * vys_r + (HD_URs[eqBBY] - HD_ULs[eqBBY]) * sgn)/temp; 
    HD_ULss[eqMMY] = vy_ss * HD_ULss[eqRHO]; // Vy** (m05 eq 59)
    HD_URss[eqMMY] = vy_ss * HD_URss[eqRHO]; 

    double vz_ss   = (temp_l1 * vzs_l + temp_r1 * vzs_r + (HD_URs[eqBBZ] - HD_ULs[eqBBZ]) * sgn)/temp;
    HD_ULss[eqMMZ] = vz_ss  * HD_ULss[eqRHO]; // Vz** (m05 eq 60)
    HD_URss[eqMMZ] = vz_ss  * HD_URss[eqRHO]; 

    HD_ULss[eqBBX] = HD_URss[eqBBX] = BX;
    HD_ULss[eqBBY] = HD_URss[eqBBY] = (temp_l1 * HD_URs[eqBBY] + temp_r1 * HD_ULs[eqBBY] 
					+ temp_l1 * temp_r1 * (vys_r - vys_l) * sgn)/temp; //By** eq 61
    HD_ULss[eqBBZ] = HD_URss[eqBBZ] = (temp_l1 * HD_URs[eqBBZ] + temp_r1 * HD_ULs[eqBBZ] 
					+ temp_l1 * temp_r1 * (vzs_r - vzs_l) * sgn)/temp; //Bz** eq 62

    temp = HD_lambda[2] * HD_ULss[eqBBX] + vy_ss * HD_ULss[eqBBY] + vz_ss * HD_ULss[eqBBZ]; // B** dot v**

    HD_ULss[eqERG] = HD_ULs[eqERG] - temp_l1 * (temp_l2 - temp) * sgn; // E** (m05 eq 63)
    HD_URss[eqERG] = HD_URs[eqERG] + temp_r1 * (temp_r2 - temp) * sgn;
  }

  //rep.printVec("HD_UL",HD_UL,eq_nvar);
  //rep.printVec("HD_ULs",HD_ULs,eq_nvar);
  //rep.printVec("HD_ULss",HD_ULss,eq_nvar);
  //rep.printVec("HD_URss",HD_URss,eq_nvar);
  //rep.printVec("HD_URs",HD_URs,eq_nvar);
  //rep.printVec("HD_UR",HD_UR,eq_nvar);

  // change out_p with UtoP(HD_UL[v],out_ps,gamma); // Tmin ?

  // fluxes (m05 eq 66)
  if      (HD_lambda[0]>0)
    for (int v=0; v<eq_nvar; v++) 
      out_flux[v] = HD_FL[v];
  //    out_p[v]    = HD_UL[v];
  else if (HD_lambda[1]>=0)
    for (int v=0; v<eq_nvar; v++) // Fl* (m05 eq. 64) 
      out_flux[v] = HD_FL[v] + HD_lambda[0] * (HD_ULs[v] - HD_UL[v]); 
  //    out_p[v]    = HD_ULs[v];
  else if (HD_lambda[2]>=0)
    for (int v=0; v<eq_nvar; v++)   // Fl** (m05 eq. 65)
      out_flux[v] = HD_FL[v] + HD_lambda[1] * HD_ULss[v] - (HD_lambda[1] - HD_lambda[0]) * HD_ULs[v] - HD_lambda[0] * HD_UL[v];
  //    out_p[v]    = HD_ULss[v];
  else if (HD_lambda[3]>=0)
    for (int v=0; v<eq_nvar; v++) 
      out_flux[v] = HD_FR[v] + HD_lambda[3] * HD_URss[v] - (HD_lambda[3] - HD_lambda[4]) * HD_URs[v] - HD_lambda[4] * HD_UR[v];
  //    out_p[v]    = HD_URss[v];
  else if (HD_lambda[4]>=0)
    for (int v=0; v<eq_nvar; v++) 
      out_flux[v] = HD_FR[v] + HD_lambda[4] * (HD_URs[v] - HD_UR[v]);
  //    out_p[v]    = HD_URs[v];
  else for (int v=0; v<eq_nvar; v++) out_flux[v] = HD_FR[v];
  //    out_p[v]    = HD_UR[v];

  
  //rep.printVec("out_Flux",out_flux,eq_nvar);


  return 0;
}


// ###################################################################
// ###################################################################




