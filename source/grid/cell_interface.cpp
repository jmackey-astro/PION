///
/// \file   cell_interface.cc
/// \author Jonathan Mackey
/// \date   12.11.2010
///
/// Purpose: Define a set of routines for accessing cell data; these
/// would be too memory intensive if each cell had to have a copy of
/// the functions.
/// 
/// History: Used to be in global.cc (up to SVN rev. 236).
///
/// Modifications:
/// - 2010.11.12 JM: Added support for H-correction speeds.  Changed
///   'col' to be an access function 'monochromatic_tau(cell *c)'.
/// - 2010.11.15/19 JM: Debugged.
/// - 2010.12.30 JM: Added DivV() set/get functions.
/// - 2011.02.17 JM: Enhanced optical depth storage for multiple sources.
/// - 2011.02.24 JM: Simplified optical depth storage a little. (tidied up 02.25).
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2011.03.21 JM: Updated optical-depth info, multiple variables per source.
/// - 2011.04.18 JM: Added storage for dS, path length through a cell for raytracing.
/// - 2011.10.17 JM: Updated RT storage.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.08.20 JM: Moved raytracing set/get functions to header and
///    made them inline.
/// - 2013.09.20 JM: Changed initialisation of unsigned ints to zero.
/// - 2015.01.10 JM: New include statements for new file structure.

#include "cell_interface.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
using namespace std;


#ifdef COUNT_ENERGETICS
struct energetics *GLOBAL_CE=0;
///< for tracking rates in microphysics/raytracing.
#endif

class cell_interface CI;
///< global class for accessing cell data, positions, neigbours.


/************************* CELL INTERFACE ***********************/

// ##################################################################
// ##################################################################

cell_interface::cell_interface()
{
  minimal_cell = false;
  dxo2 = -HUGEVALUE;
  xmin = 0;
  //
  // this means I can have grids with up to 5e5 zones before it fails...
  //
  int_converter = 1.0+2.0e-6; // ONE_PLUS_EPS;
  cell_size_int_units=2;

  /// This must be set to true to create a cell.
  have_setup_extra_data = false;
  /// Flag: 0=not doing RT,  N>=1 = N tau values.
  using_RT = 0;
  /// Flag: 0=no Hcorr, N=need N variables (Lapidus=1,Hcorr=Ndim).
  using_Hcorr = 0;
  /// Size of extra_data array (can be zero).
  N_extra_data = 0;
  //
  // index arrays initialise to zero.
  //
  NTau  = 0;
  iTau  = 0;
  iDTau = 0;
  iVsh  = 0;
  idS   = 0;
  /// indices of Hcorrection values in extra_data [XX,YY,ZZ].
  for (int v=0; v<MAX_DIM; v++)
    iHcorr[v] = 0;

}


// ##################################################################
// ##################################################################

cell_interface::~cell_interface()
{
  if (xmin) delete [] xmin;
  xmin=0;
  if (using_RT>0) {
    NTau  = mem.myfree(NTau);
    iTau  = mem.myfree(iTau);
    iDTau = mem.myfree(iDTau);
    iVsh  = mem.myfree(iVsh);
    idS   = mem.myfree(idS);
  }
}


// ##################################################################
// ##################################################################

void cell_interface::set_minimal_cell_data()
{
  minimal_cell = true;
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::unset_minimal_cell_data()
{
  minimal_cell = false;
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::set_dx(const double dx)
{
  cell_diameter = dx;
  dxo2 = 0.5*dx;
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::set_xmin(const double *xm)
{
  if (!xmin) {
    xmin = new double [SimPM.ndim];
  }
  for (int v=0; v<SimPM.ndim; v++) {
    xmin[v] = xm[v];
  }
  return;
}



// ##################################################################
// ##################################################################

//
// Set variables for extra_data based on what we need per cell.
// Currently monochromatic radiation needs one double, the
// H-correction needs Ndim doubles, and Lapidus viscosity one double.
//
void cell_interface::setup_extra_data(
        const struct rad_sources &rsi, ///< Flag for ray-tracing
        const int hc_flag,  ///< Flag for H-correction
        const int dv_flag   ///< Flag for Div(V).
        )
{
#ifdef TESTING
  cout <<"\ncell_interface::setup_extra_data():\n";
#endif
  //
  // Set up a 1D array for all of the extra data that a cell needs,
  // and set indices to access the required elements.
  //
  // Start with no extra data:
  //
  N_extra_data=0;

  //
  // Now add ray-tracing optical depth variable(s).
  //
  if    (rsi.Nsources <= 0) {
    // do nothing because there is no raytracing.
  }
  else {
    //
    // A number of sources, each with a number of variables determined by the source type.
    //
    using_RT = 1;
    //
    // Each source, s, numbered from 0 to Nsources-1, needs a number
    // of variables:
    // - NTau[s]:  the number of Tau and DTau variables for source s.
    // - iTau[s]:  index of first Tau variable in s (others follow).
    // - iDTau[s]: index of first DTau variable in s (others follow).
    // - iVsh[s]:  index of Vshell variable in extra_data[] for s.
    // - idS[s]:   index of dS variable in extra_data[] for s.
    //
    NTau  = mem.myalloc(NTau,  rsi.Nsources);
    iTau  = mem.myalloc(iTau,  rsi.Nsources);
    iDTau = mem.myalloc(iDTau, rsi.Nsources);
    iVsh  = mem.myalloc(iVsh,  rsi.Nsources);
    idS   = mem.myalloc(idS,   rsi.Nsources);

    for (int s=0; s<rsi.Nsources; s++) {
      //
      // Number of quantities traced from source:
      //
      NTau[s] = rsi.sources[s].NTau;
      //
      // New update with rates: need NTau[s] vars for Tau and DTau,
      // one var for Vshell and dS.
      //
      iTau[s]  = N_extra_data; N_extra_data += NTau[s];
      iDTau[s] = N_extra_data; N_extra_data += NTau[s];
      iVsh[s]  = N_extra_data; N_extra_data++;
      idS[s]   = N_extra_data; N_extra_data++;
    } // loop over radiation sources.
#ifdef TESTING
    cout <<"\t\t Adding RT: N="<<N_extra_data<<"\n";
#endif
  }

  //
  // Now add H-correction/Lapidus viscosity variables
  //
  if (hc_flag) {
    using_Hcorr = hc_flag;
    if (hc_flag > MAX_DIM)
      rep.error("Can't ask for more than MAX_DIM H-corr variables",hc_flag);

    for (int v=0; v<hc_flag; v++) {
      iHcorr[v] = N_extra_data;
      N_extra_data += 1;
    }
#ifdef TESTING
    cout <<"\t\t Adding HCORR: N="<<N_extra_data<<"\n";
#endif
  }

  //
  // DivV just needs one extra variable.
  //
  if (dv_flag) {
    using_DivV = dv_flag;
    iDivV = N_extra_data;
    N_extra_data += 1;
#ifdef TESTING
    cout <<"\t\t Adding DIVV: N="<<N_extra_data<<"\n";
#endif
  }


  have_setup_extra_data = true;
#ifdef TESTING
  cout <<"\n";
#endif
  return;
}


// ##################################################################
// ##################################################################


// returns true if using minimal cells.
bool cell_interface::query_minimal_cells()
{return minimal_cell;}


// ##################################################################
// ##################################################################

cell * cell_interface::new_cell()
{
  if (!have_setup_extra_data)
    rep.error("Setup extra data before calling new_cell!",using_RT);

  //
  // If this is the first cell we are assigning, make sure we have set
  // dx/2 and the xmin pointer correctly.
  //
  if (dxo2<0.0) {
    rep.error("Set dx before creating cells!",dxo2);
    dxo2 = SimPM.dx/2.0;
    xmin = SimPM.Xmin;
#ifdef TESTING
    cout <<"dx/2="<<dxo2<<", dx="<<SimPM.dx<<"   ";
    rep.printVec("CELL INTERFACE, xmin",xmin,SimPM.ndim);
#endif // TESTING
  }
  cell *c=0;
  c = mem.myalloc(c,1);
  c->pos = 0;


  //
  // Allocate memory and initialise to zero.
  //
  c->P   = mem.myalloc(c->P,  SimPM.nvar);
  c->ngb = mem.myalloc(c->ngb, 2*SimPM.ndim);
  c->pos = mem.myalloc(c->pos, SimPM.ndim);

  for (int v=0;v<SimPM.ndim;v++) c->pos[v]=0;
  for (int v=0;v<SimPM.nvar;v++) c->P[v] = 0.0;
  for (int v=0; v<2*SimPM.ndim; v++) c->ngb[v] = 0;
  c->npt = 0;
  c->id  = -9999; c->isedge=-999; c->isbd=c->isgd=false;

  //
  // If we need all the [dU,Ph] arrays, initialise them, but if we have
  // set "minimal_cell" to true, then skip them to save memory in 
  // analysis code.
  //
  if (minimal_cell) {
    //cout <<"Minimal cells!\n";
    c->Ph = 0;
    c->dU = 0;
  }
  else {
    c->Ph  = mem.myalloc(c->Ph, SimPM.nvar);
    c->dU  = mem.myalloc(c->dU, SimPM.nvar);
    for (int v=0;v<SimPM.nvar;v++) c->Ph[v] = c->dU[v] = 0.0;
  }

  //cout <<"Nxd="<<N_extra_data<<"\n";
  if (N_extra_data>=1) {
    c->extra_data = mem.myalloc(c->extra_data, N_extra_data);
    for (short unsigned int v=0;v<N_extra_data;v++)
      c->extra_data[v] = 0.0;
  }

  return c;
}


// ##################################################################
// ##################################################################

void cell_interface::delete_cell(cell *c)
{
  c->pos= mem.myfree(c->pos);
  c->P  = mem.myfree(c->P);
  c->Ph = mem.myfree(c->Ph);
  c->dU = mem.myfree(c->dU);
  c->ngb = mem.myfree(c->ngb);

  if (N_extra_data>=1)
    c->extra_data = mem.myfree(c->extra_data);

  c = mem.myfree(c);
  return;
}


// ##################################################################
// ##################################################################

void cell_interface::set_pos(
      cell *c, ///< pointer to cell
      const double *p_in ///< double array of size ndim, containing cell position.
      )
{
  //
  // Set position integer according to Xmin+i*DX/2=x
  //
  for (int v=0;v<SimPM.ndim;v++) {
    c->pos[v] = static_cast<int>(int_converter*((p_in[v]-xmin[v])/dxo2));
  }
  //  rep.printVec("int-pos",c->pos,SimPM.ndim);
  return;
}


// ##################################################################
// ##################################################################

void cell_interface::set_pos(
      cell *c, ///< pointer to cell
      const int *p_in ///< integer array of size ndim, containing cell position.
      )
{
  //
  // Set position integer according to Xmin+i*DX/2=x
  // This function assumes a clever person has set p_in to have these values!
  // If not, the code may fail catastrophically.
  //
  for (int v=0;v<SimPM.ndim;v++) {
    c->pos[v] = p_in[v];
  }
  //  rep.printVec("int-pos",c->pos,SimPM.ndim);
  return;
}


// ##################################################################
// ##################################################################

void cell_interface::get_dpos(
      const cell *c, ///< pointer to cell
      double *p_out ///< array to write position into.
      )
{
  for (int v=0;v<SimPM.ndim;v++)
    p_out[v] = xmin[v] +c->pos[v]*dxo2;
  return;
}


// ##################################################################
// ##################################################################

double cell_interface::get_dpos(
      const cell *c, ///< pointer to cell
      const int v ///< element of position vector we want
      )
{
  return xmin[v] +c->pos[v]*dxo2;
}


// ##################################################################
// ##################################################################

void cell_interface::get_ipos(
      const cell *c, ///< pointer to cell
      int *ipos_out  ///< array to write integer position into.
      )
{
  for (int v=0;v<SimPM.ndim;v++)
    ipos_out[v] = c->pos[v];
  return;
}  
  

// ##################################################################
// ##################################################################

int cell_interface::get_ipos(
      const cell *c, ///< pointer to cell
      const int v    ///< element of position we want.
      )
{
  return c->pos[v];
}


// ##################################################################
// ##################################################################

void cell_interface::get_ipos_vec(
      const double *p_in, ///< physical position (input)
      int *p_out          ///< integer position (output)
      )
{
  if (dxo2<0.0)
    rep.error("set up grid before trying to get integer positions!!!",dxo2);
  for (int v=0;v<SimPM.ndim;v++) {
    if (fabs(p_in[v])>VERY_LARGE_VALUE)
      p_out[v] = -1234567;
    else
      p_out[v] = static_cast<int>(int_converter*((p_in[v]-xmin[v])/dxo2));
  }
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::get_ipos_as_double(
      const double *p_in, ///< physical position (input)
      double *p_out       ///< integer position (output)
      )
{
  if (dxo2<0.0)
    rep.error("set up grid before trying to get integer positions!!!",dxo2);
  for (int v=0;v<SimPM.ndim;v++) {
    if (fabs(p_in[v])>VERY_LARGE_VALUE)
      p_out[v] = -VERY_LARGE_VALUE;
    else
      p_out[v] = (p_in[v]-xmin[v])/dxo2;
  }
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::get_dpos_vec(
      const int *p_in, ///< integer position (output)
      double *p_out    ///< physical position (input)
      )
{
  for (int v=0;v<SimPM.ndim;v++)
        p_out[v] = xmin[v] +(p_in[v])*dxo2;
  return;
}


// ##################################################################
// ##################################################################


void cell_interface::copy_cell(
    const cell *c1,
    cell *c2
    )
{
  for (int i=0;i<SimPM.ndim;i++) c2->pos[i] = c1->pos[i];
  for (int v=0;v<SimPM.nvar;v++) {
    c2->P[v]  = c1->P[v];
    c2->Ph[v] = c1->Ph[v];
    c2->dU[v] = c1->dU[v];
  }
  for (int i=0;i<2*SimPM.ndim;i++) c2->ngb[i]=c1->ngb[i];
  c2->npt = c1->npt;
  c2->id = c1->id;
  c2->isedge = c1->isedge;
  c2->isbd = c1->isbd;
  c2->isgd = c1->isgd;
  return;
}
 

// ##################################################################
// ##################################################################

void cell_interface::print_cell(const cell *c)
{
  if(c==0) {cout <<"Null Pointer!\n"; return;}
  cout <<"cell:\t id = "<<c->id<<"\n";
  cout <<"\tcell pointer= "<<c<<"\n";
  cout <<"\tisedge:"<<c->isedge<<"\tisbd:"<<c->isbd<<"\tisgd:"<<c->isgd<<"\n";
  cout<<"\tnpt: "<<c->npt;
  if (c->npt!=0) cout <<"\tnpt[id]: "<<c->npt->id<<"\n";
  else cout <<"\tnpt is not addressed (last point?).\n";
  if (N_extra_data>0) {
    cout <<"\t";
    rep.printVec("extra_data[]",c->extra_data,N_extra_data);
  }
  cout <<"\t"; rep.printVec("pos[]",c->pos,SimPM.ndim);
  cout <<"\t"; double p[SimPM.ndim]; get_dpos(c,p); rep.printVec("dpos[]",p,SimPM.ndim);
  cout <<"\t"; rep.printVec("P[]  ",c->P,SimPM.nvar);
  if (!minimal_cell) {
    cout <<"\t"; rep.printVec("Ph[] ",c->Ph,SimPM.nvar);
    cout <<"\t"; rep.printVec("dU[] ",c->dU,SimPM.nvar);
  }
  cout <<"\t"; rep.printVec("ngb[]",c->ngb,2*SimPM.ndim);
  return;
}




/************************* CELL INTERFACE ***********************/

