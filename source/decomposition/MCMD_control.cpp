/// \file MCMD_control.cpp
/// 
/// \brief Defines class for controlling Multi-Core-Multi-Domain
///        simulations.
/// 
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.27 JM: moved from sim_control_MPI.cpp
/// - 2016.02.02 JM: Added option to decompose only along one axis.
/// - 2018.01.24 JM: worked on making SimPM non-global


#include "decomposition/MCMD_control.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include "constants.h"
#include "tools/command_line_interface.h"
#include <algorithm>
using namespace std;

//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------

// ##################################################################
// ##################################################################


MCMDcontrol::MCMDcontrol()
{
  nproc = -1;
  myrank = -1;
  LocalNcell = -1;
  for (int i=0; i<MAX_DIM; i++) {
    LocalNG[i] = offsets[i] = ix[i] = nx[i] = -1;
    LocalXmin[i] = LocalXmax[i] = LocalRange[i] = -1.e99;
  }
  ngbprocs=0;
  parent_proc=-1;

  ReadSingleFile  =true; ///< If the ICs are in a single file, set this to true.
  WriteSingleFile =false; ///< If all processes to write to one file, set this
  WriteFullImage  =false; ///< If multiple fits files, each is the full domain size.
}

// ##################################################################
// ##################################################################


MCMDcontrol::~MCMDcontrol() {
  if (ngbprocs!=0)
    ngbprocs = mem.myfree(ngbprocs);
  return;
}


// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
      class SimParams &SimPM,  ///< simulation parameters
      class level &level     ///< parameters for NG grid level
      )
{
#ifdef TESTING
  cout << "---MCMDcontrol::decomposeDomain() decomposing domain. ";
  cout << " Nproc="<<nproc<<", myrank="<<myrank<<"\n";
#endif
  
  //
  // We subdivide the domain in half recursively, where the axis we cut along
  // is always the one in which the domain is longest.  In this way we minimize
  // domain interfaces, which should minimize communication.
  //
  if (SimPM.ndim==1) {
    // 1D is so simple it's worth putting in it's own section.
    LocalRange[XX] = level.Range[XX]/nproc;
    LocalXmin[XX] = level.Xmin[XX] + myrank*LocalRange[XX];
    LocalXmax[XX] = LocalXmin[XX] + LocalRange[XX];
    LocalNG[XX] = level.NG[XX]/nproc;
    offsets[XX] = 0;
    LocalNcell = LocalNG[XX];
    // Set up ngbprocs array to point to neighbouring processors
    nx[XX] = nproc;
    ix[XX] = myrank;
    pointToNeighbours(SimPM,level);
  } // If 1D
  
  else if (SimPM.ndim==2 || SimPM.ndim==3) {
    for (int i=0;i<SimPM.ndim;i++) {
      LocalRange[i] = level.Range[i];
      LocalXmin[i]  = level.Xmin[i];
      LocalXmax[i]  = level.Xmax[i];
      LocalNG[i]    = level.NG[i];
      LocalNcell    = level.Ncell;
    }
    double maxrange=0.;
    enum axes maxdir=XX;
    int npcounter=1;
    
    // Loop, dividing range in half n times, until npcounter=nproc.
    if (nproc==1) { 
      //cout <<"Only one processor!!! Use serial version...\n";
      //rep.error("Only one processor!!! Use serial version",nproc);
    }
    while (npcounter < nproc) {
      // find longest axis of subrange.
      int i=0, dsrc=-1;
      
      // --- Check if we are doing raytracing with a source at infinity. ---
      if (SimPM.EP.raytracing && SimPM.RS.Nsources>0) {
        //
	// if one source at infinity, b/c then we decompose to keep
	// rays on one processor all the time.
        //
        bool at_infinity=true;
        for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
          if (SimPM.RS.sources[isrc].type==RT_SRC_DIFFUSE) at_infinity=false;
          if ((SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) &&
              (!SimPM.RS.sources[isrc].at_infinity))       at_infinity=false;
        }
        //
        // Now if at_infinity is still set, we want to get the direction along
        // which we don't want to decompose the domain.
        // Now also checks if there are multiple sources at infinity in
        // different directions, in which case we can't make all of them fully
        // parallel.
        //
        if (at_infinity) {
          for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
            if ((SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) && 
                (SimPM.RS.sources[isrc].at_infinity) ) {
              for (int ii=0;ii<SimPM.ndim;ii++) {
                if (fabs(SimPM.RS.sources[isrc].pos[ii])>1.e99) {
                  if (dsrc==-1) dsrc=ii;
                  // srcs at infinity in multiple directions
                  else if (dsrc!=ii) at_infinity=false;
                }
              } // loop over directions
            }   // if we are at the ionising source.
          }     // loop over sources.
        }       // if at infinity
        
        if (at_infinity) {
          //cout <<"\t\tFound a source at infinity in direction "<<dsrc;
          //cout <<", so not decomposing domain in this direction.\n";
          if (dsrc<0) rep.error("no direction to source at infinity",dsrc);
        }
        else {
          //cout <<"\t\tEither multiple or no sources at infinity.\n";
          dsrc=-1;
        }
      }
      dsrc = -1;
      // --- end of RT source at infinity bit ---
      
      maxrange=0.; i=0;
      while (i<SimPM.ndim) {
        if (LocalRange[i] > maxrange && dsrc!=i) {
          maxrange=LocalRange[i];
          maxdir = static_cast<axes>(i);
        }
        i++;
      }
      // Half that range and multiply nproc by 2.
      LocalRange[maxdir] /= 2.0;
      npcounter *=2;
    }
    if (npcounter != nproc) rep.error("1: nproc not a power of 2!",nproc);
    
    //
    // Now we know the range of each subcell, so calculate where I fit into
    // the hierarchy, defined by
    // \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
    // This requires myrank to count from zero!
    //
    for (int i=0;i<SimPM.ndim;i++)
      nx[i] =static_cast<int>(ONE_PLUS_EPS*level.Range[i]/LocalRange[i]);
    int temp=myrank;
    if (SimPM.ndim==3) {
      ix[ZZ] = temp/nx[XX]/nx[YY];
      temp -= ix[ZZ]*nx[XX]*nx[YY];
    }
    ix[YY] = temp/nx[XX];
    temp -= ix[YY]*nx[XX];
    ix[XX] = temp;
    
    LocalNcell = level.Ncell;
    for (int i=0;i<SimPM.ndim;i++) {
      LocalXmin[i]  = level.Xmin[i] +ix[i]*LocalRange[i];
      LocalXmax[i]  = level.Xmin[i] +(ix[i]+1)*LocalRange[i];
      LocalNG[i]    = level.NG[i]/nx[i];
      LocalNcell   /= nx[i];
      offsets[i]    = ix[i]*LocalNG[i];
    }
    // Set up ngbprocs array to point to neighbouring processors
    pointToNeighbours(SimPM,level);
  } // if 2D or 3D
  
  else rep.error("Bad NDIM in DecomposeDomain",SimPM.ndim);

#ifdef TESTING
  // Display some debugging info.
  if (myrank==0) {
    for (int i=0;i<SimPM.ndim;i++) {
      cout <<"Sim: idim="<<i<<"  \tRange="<<level.Range[i];
      cout <<",\t  xmin,xmax = "<<level.Xmin[i]<<", "<<level.Xmax[i];
      cout <<"    \t Ncell = "<<level.Ncell<<"\n";
    }
  }
  for (int i=0;i<SimPM.ndim;i++) {
    cout <<"Proc "<<myrank<<": idim="<<i<<"  \tRange="<<LocalRange[i];
    cout <<",\t  xmin,xmax = "<<LocalXmin[i]<<", "<<LocalXmax[i];
    cout <<"    \t Ncell = "<<LocalNcell;
    cout <<"\t neighbours : "<<ngbprocs[2*i]<<", "<<ngbprocs[2*i+1]<<"\n";
  }
  cout << "---MCMDcontrol::decomposeDomain() Domain decomposition done.\n\n";
#endif
  return(0);
}



// ##################################################################
// ##################################################################



int MCMDcontrol::decomposeDomain(
      const enum axes daxis,   ///< Axis to decompose domain along.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class level &level     ///< pointer to domain parameters for NG grid level
      )
{
  cout << "---MCMDcontrol::decomposeDomain() decomposing domain";
  cout << " along axis " << static_cast<int>(daxis) <<"\n";

  //
  // We subdivide the domain in half recursively, along one axis
  //
  if (SimPM.ndim==1) {
    rep.error("Use nomal domain decomposition for 1D",daxis);
  } // If 1D
  
  else if (SimPM.ndim==2 || SimPM.ndim==3) {
    for (int i=0;i<SimPM.ndim;i++) {
      LocalRange[i] = level.Range[i];
      LocalXmin[i]  = level.Xmin[i];
      LocalXmax[i]  = level.Xmax[i];
      LocalNG[i]    = level.NG[i];
      LocalNcell    = level.Ncell;
    }
    //
    // npcounter counts how many sub-domains I have.
    //
    int npcounter=1;
    
    while (npcounter < nproc) {
      // find longest axis of subrange.
      // Half that range and multiply nproc by 2.
      LocalRange[daxis] /= 2.0;
      npcounter *=2;
    }
    if (npcounter != nproc) {
      rep.error("2: nproc not a power of 2!",nproc);
    }
    
    //
    // Now we know the range of each subcell, so calculate where I fit into
    // the hierarchy, defined by
    // \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
    // This requires myrank to count from zero!
    //
    for (int i=0;i<SimPM.ndim;i++)
      nx[i] =static_cast<int>(ONE_PLUS_EPS*level.Range[i]/LocalRange[i]);
    int temp=myrank;
    if (SimPM.ndim==3) {
      ix[ZZ] = temp/nx[XX]/nx[YY];
      temp -= ix[ZZ]*nx[XX]*nx[YY];
    }
    ix[YY] = temp/nx[XX];
    temp -= ix[YY]*nx[XX];
    ix[XX] = temp;
    
    LocalNcell = level.Ncell;
    for (int i=0;i<SimPM.ndim;i++) {
      LocalXmin[i]  = level.Xmin[i] +ix[i]*LocalRange[i];
      LocalXmax[i]  = level.Xmin[i] +(ix[i]+1)*LocalRange[i];
      LocalNG[i]    = level.NG[i]/nx[i];
      LocalNcell   /= nx[i];
      offsets[i]    = ix[i]*LocalNG[i];
    }
    // Set up ngbprocs array to point to neighbouring processors
    pointToNeighbours(SimPM,level);    
  } // if 2D or 3D
  
  else rep.error("Bad NDIM in DecomposeDomain",SimPM.ndim);

  rep.printVec("NXYZ",LocalNG,SimPM.ndim);
  rep.printVec("XMIN",LocalXmin,SimPM.ndim);
  rep.printVec("XMAX",LocalXmax,SimPM.ndim);
  rep.printVec("RANGE",LocalRange,SimPM.ndim);


  cout << "---MCMDcontrol::decomposeDomain() finished.\n\n";
  return 0;
}


// ##################################################################
// ##################################################################


int MCMDcontrol::pointToNeighbours(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class level &level     ///< pointer to domain parameters for NG grid level
      )
{
#ifdef TESTING
  cout <<"Setting up neighbour pointers in MCMD.\n";
#endif
  ///
  /// \section PBC Periodic Boundaries
  /// Note that if there are periodic boundaries, the pointers to the 
  /// wrapped around neighbouring processors are not set here.  They are
  /// set when the grid boundaries are set up, in the function
  /// UniformGridParallel::BC_setBCtypes
  /// I'm not completely happy with this as it makes the 'modular-ness' of 
  /// the code weaker, but it is something that will work.  The neighbouring
  /// processor list is a global variable, and I don't think it's the worst
  /// thing for it to be set in two places.
  /// Here any simulation boundaries have neighbouring processor id set to 
  /// -999.
  ///
  int nx[SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++) {
    nx[i] =static_cast<int>(ONE_PLUS_EPS*level.Range[i]/LocalRange[i]);
  }
  // Point to neigbours
  if (ngbprocs) {
    //rep.error("neigbours already set up!",ngbprocs);
    delete [] ngbprocs; ngbprocs=0;
  }
  ngbprocs = new int [2*SimPM.ndim];
  if (!ngbprocs) rep.error("Memory Allocation of neighbours.",ngbprocs);
  ngbprocs[XN] = myrank -1;
  ngbprocs[XP] = myrank +1;
  if (pconst.equalD(LocalXmin[XX],level.Xmin[XX])) ngbprocs[XN] = -999;
  if (pconst.equalD(LocalXmax[XX],level.Xmax[XX])) ngbprocs[XP] = -999;
  if (SimPM.ndim >1) {
    ngbprocs[YN] = myrank -nx[XX];
    ngbprocs[YP] = myrank +nx[XX];
    if (pconst.equalD(LocalXmin[YY],level.Xmin[YY])) ngbprocs[YN] = -999;
    if (pconst.equalD(LocalXmax[YY],level.Xmax[YY])) ngbprocs[YP] = -999;
  }
  if (SimPM.ndim >2) {
    ngbprocs[ZN] = myrank -nx[XX]*nx[YY];
    ngbprocs[ZP] = myrank +nx[XX]*nx[YY];
    if (pconst.equalD(LocalXmin[ZZ],level.Xmin[ZZ])) ngbprocs[ZN] = -999;
    if (pconst.equalD(LocalXmax[ZZ],level.Xmax[ZZ])) ngbprocs[ZP] = -999;
  }
#ifdef TESTING
  rep.printVec("ngbprocs",ngbprocs,2*SimPM.ndim);
#endif

  return(0);
}


// ##################################################################
// ##################################################################


///
/// Get a list of all abutting domains, including corner/edge intersections.
///
void MCMDcontrol::get_abutting_domains(
      const int ndim,  ///< grid dimensions
      std::vector<int> &dl ///< write list to this vector.
      )
{
  //
  // If we previously generated the list, then just copy elements.
  //
  //if (!full_ngb_list.empty()) {
  //  for (unsigned int v=0; v<full_ngb_list.size(); v++)
  //    dl.push_back(full_ngb_list[v]);
  //  return;
  //}
  full_ngb_list.clear();

  //
  // Otherwise we generate the list and then copy.
  //
  //
  // List ordering is as follows:
  //  XN,XP,YN,YP,ZN,ZP
  //  YNXN, YNXP, YPXN, YPXP,
  //  ZNXN, ZNXP, ZNYN, ZNYP,
  //  ZNYNXN, ZNYNXP, ZNYPXN, ZNYPXP
  //  ZPXN, ZPXP, ZPYN, ZPYP,
  //  ZPYNXN, ZPYNXP, ZPYPXN, ZPYPXP
  //

  //
  // Find which coordinate directions have neighbours and store in a
  // bool array.
  //
  bool d[2*MAX_DIM];
  for (int v=0;v<2*MAX_DIM;v++) d[v] = false;
  enum direction posdir, negdir;
  for (int v=0;v<ndim;v++) {
    negdir = static_cast<direction>(2*v);
    posdir = static_cast<direction>(2*v+1);
    if (ngbprocs[negdir]>=0) d[negdir] = true;
    if (ngbprocs[posdir]>=0) d[posdir] = true;
  }

  //
  // Add coordinate directions to list.
  //
  for (int v=0; v<2*ndim; v++) {
    if (d[v]) full_ngb_list.push_back(ngbprocs[v]);
  }

  //
  // X-Y plane -- can read off neighbour and add/subtract 1 from it.
  //
  if (d[YN] && d[XN]) full_ngb_list.push_back(ngbprocs[YN]-1);
  if (d[YN] && d[XP]) full_ngb_list.push_back(ngbprocs[YN]+1);
  if (d[YP] && d[XN]) full_ngb_list.push_back(ngbprocs[YP]-1);
  if (d[YP] && d[XP]) full_ngb_list.push_back(ngbprocs[YP]+1);

  //
  // X-Y-ZN plane -- need to calculate rank since offsets are not trivial.
  //
  if (d[ZN]) {
    if (d[XN]) 
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]) +(ix[XX]-1) );
    if (d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]) +(ix[XX]+1) );
    if (d[YN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]) );
    if (d[YP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]) );
    if (d[YN] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]-1) );
    if (d[YN] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]+1) );
    if (d[YP] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]-1) );
    if (d[YP] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]+1) );
  }
  //
  // X-Y-ZP plane
  //
  if (d[ZP]) {
    if (d[XN]) 
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]  ) +(ix[XX]-1) );
    if (d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]  ) +(ix[XX]+1) );
    if (d[YN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]  ) );
    if (d[YP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]  ) );
    if (d[YN] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]-1) );
    if (d[YN] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]+1) );
    if (d[YP] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]-1) );
    if (d[YP] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]+1) );
  }

  //
  // Copy to return vector and return.
  //
  for (unsigned int v=0; v<full_ngb_list.size(); v++)
    dl.push_back(full_ngb_list[v]);

  return;
}



// ##################################################################
// ##################################################################



//
// Returns the ix array for any requested rank.
//
void MCMDcontrol::get_domain_ix(
      const int ndim,  ///< grid dimensions
      const int r, ///< rank ix requested for.
      int *arr     ///< array to put ix into.
      )
{
  int temp=r;
  if (ndim==3) {
    arr[ZZ] = temp/nx[XX]/nx[YY];
    temp -= arr[ZZ]*nx[XX]*nx[YY];
  }
  if (ndim>1) {
    arr[YY] = temp/nx[XX];
    temp -= arr[YY]*nx[XX];
  }
  arr[XX] = temp;
  return;
}



// ##################################################################
// ##################################################################



int MCMDcontrol::get_grid_rank(
      class SimParams &par,  ///< simulation parameters
      const double *loc, ///< location sought
      const int l      ///< grid level
      )
{
  int proc=-1;
  // ir[] is index of grid in each direction on level l.
  int ir[par.ndim];
  //rep.printVec("NX",nx,par.ndim);
  //rep.printVec("Loc",loc,par.ndim);
  //rep.printVec("Xmin",par.levels[l].Xmin,par.ndim);
  //rep.printVec("Range",par.levels[l].Range,par.ndim);
  for (int i=0;i<par.ndim;i++) {
    ir[i] = static_cast<int>(floor(nx[i]*(loc[i]-par.levels[l].Xmin[i])
                                        / par.levels[l].Range[i]));
  }
  //rep.printVec("IR",ir,par.ndim);
  // check validity (>0 and <nx[])
  for (int i=0;i<par.ndim;i++) {
    if (ir[i]<0  || ir[i]>=nx[i]) {
      return -1;
    }
  }
  proc = ir[XX];
  if (par.ndim>1) proc += nx[XX]*ir[YY];
  if (par.ndim>2) proc += nx[XX]*nx[YY]*ir[ZZ];
  return proc;
}



// ##################################################################
// ##################################################################



void MCMDcontrol::set_NG_hierarchy(
      class SimParams &par,  ///< simulation parameters
      const int l  ///< level to work on
      )
{
#ifdef TEST_BC89FLUX
  cout <<"Setting up NG hierarchy (MCMDcontrol)\n";
#endif
  double centre[MAX_DIM];
  child_procs.clear();
  int ix[MAX_DIM];

  // set rank of parent for each grid except root level 0
  // Every child grid must have a parent because the nested grid is
  // entirely within the coarser level.
  if (l>0) {
    // centre[] is centre of local grid on level l.
    for (int i=0;i<par.ndim;i++)
      centre[i] = 0.5*(LocalXmin[i]+LocalXmax[i]);
    // get rank of parent grid in each direction.
    parent_proc = get_grid_rank(par, centre,l-1);
    //cout <<"level "<<l<<", parent process is "<<parent_proc<<"\n";

    pgrid.rank = parent_proc;
    get_domain_ix(par.ndim, parent_proc, ix);
    for (int i=0;i<par.ndim;i++) {
      // parent level has range 2x my range in each direction
      pgrid.Xmin[i] = par.levels[l-1].Xmin[i] + ix[i]*2.0*LocalRange[i];
      pgrid.Xmax[i] = pgrid.Xmin[i] + 2.0*LocalRange[i];
    }
    
    // now find neighbouring grids, if they exist.  First set ranks.
    pgrid_ngb.resize(2*par.ndim);
    pgrid_ngb[XN].rank = parent_proc -1;
    pgrid_ngb[XP].rank = parent_proc +1;
    if (pconst.equalD(pgrid.Xmin[XX],par.levels[l-1].Xmin[XX])) pgrid_ngb[XN].rank = -999;
    if (pconst.equalD(pgrid.Xmax[XX],par.levels[l-1].Xmax[XX])) pgrid_ngb[XP].rank = -999;
    if (par.ndim >1) {
      pgrid_ngb[YN].rank = parent_proc -nx[XX];
      pgrid_ngb[YP].rank = parent_proc +nx[XX];
      if (pconst.equalD(pgrid.Xmin[YY],par.levels[l-1].Xmin[YY])) pgrid_ngb[YN].rank = -999;
      if (pconst.equalD(pgrid.Xmax[YY],par.levels[l-1].Xmax[YY])) pgrid_ngb[YP].rank = -999;
    }
    if (par.ndim >2) {
      pgrid_ngb[ZN].rank = parent_proc -nx[XX]*nx[YY];
      pgrid_ngb[ZP].rank = parent_proc +nx[XX]*nx[YY];
      if (pconst.equalD(pgrid.Xmin[ZZ],par.levels[l-1].Xmin[ZZ])) pgrid_ngb[ZN].rank = -999;
      if (pconst.equalD(pgrid.Xmax[ZZ],par.levels[l-1].Xmax[ZZ])) pgrid_ngb[ZP].rank = -999;
    }
    // Now set Xmin, Xmax:
    for (int d=0;d<2*par.ndim;d++) {
      if (pgrid_ngb[d].rank>=0) {
        get_domain_ix(par.ndim, pgrid_ngb[d].rank, ix);
        for (int i=0;i<par.ndim;i++)
          pgrid_ngb[d].Xmin[i] = par.levels[l-1].Xmin[i] + ix[i]*2.0*LocalRange[i];
        for (int i=0;i<par.ndim;i++)
          pgrid_ngb[d].Xmax[i] = pgrid_ngb[d].Xmin[i] + 2.0*LocalRange[i];
      }
    }

#ifdef TEST_BC89FLUX
    // *** debugging info ***
    cout <<"level "<<l<<", parent proc = "<<parent_proc<<", parent grid xmin/xmax:\n";
    rep.printVec("Xmin",pgrid.Xmin,par.ndim);
    rep.printVec("Xmax",pgrid.Xmax,par.ndim);
    for (int d=0;d<2*par.ndim;d++) {
      if (pgrid_ngb[d].rank<0)
        cout <<"pproc ngb in direction "<<d<<" has no neighbour.\n";
      else {
        cout <<"pproc ngb in direction "<<d<<" has neighbour proc ";
        cout <<pgrid_ngb[d].rank<<"\n";
        rep.printVec("Xmin",pgrid_ngb[d].Xmin,par.ndim);
        rep.printVec("Xmax",pgrid_ngb[d].Xmax,par.ndim);
      }
    }
    cout.flush();
    // *** debugging info ***
#endif
  }

  // set rank of child grids, if they exist.
  // Domain is intersection of child full grid and this local grid.
  // Must be split in half/quadrant/octant of this grid, unless
  // nproc==1, for which there is only one child so it is trivial.
  struct cgrid child;
  int child_rank=-1;
  vector<int> children;

  if (l<par.grid_nlevels-1 && nproc==1) {
    child.rank = myrank;
    for (int v=0;v<par.ndim;v++)
      child.Xmin[v] = par.levels[l+1].Xmin[v];
    for (int v=0;v<par.ndim;v++)
      child.Xmax[v] = par.levels[l+1].Xmax[v];
    child_procs.push_back(child);
#ifdef TEST_BC89FLUX
    cout <<"v="<<0<<": only child has rank="<<child_rank;
    cout <<", child_procs_size="<<child_procs.size()<<"\n";
#endif
  }
  else if (l<par.grid_nlevels-1) {
    // split domain in 4 in each dimension, and get ranks of all 
    // grids in these sub-patches.  Eliminate duplicates, and then
    // set up child vector.
    if (par.ndim==1) {
      for (int i=0;i<4;i++) {
        centre[XX] = LocalXmin[XX] + 0.25*i*LocalRange[XX];
        child_rank = get_grid_rank(par, centre,l+1);
        if (child_rank>=0) children.push_back(child_rank);
      }
    }
    else if (par.ndim==2) {
      // find all ranks in 4x4 grid of sub-patches.
      for (int i=0;i<4;i++) {
        for (int j=0; j<4; j++) {
          centre[XX] = LocalXmin[XX] + 0.25*i*LocalRange[XX];
          centre[YY] = LocalXmin[YY] + 0.25*j*LocalRange[YY];
          child_rank = get_grid_rank(par, centre,l+1);
          if (child_rank>=0) children.push_back(child_rank);
        }
      }
    }
    else if (par.ndim==3) {
      // find all ranks in 4x4x4 grid of sub-patches.
      for (int i=0;i<4;i++) {
        for (int j=0; j<4; j++) {
          for (int k=0; k<4; k++) {
            centre[XX] = LocalXmin[XX] + 0.25*i*LocalRange[XX];
            centre[YY] = LocalXmin[YY] + 0.25*j*LocalRange[YY];
            centre[ZZ] = LocalXmin[ZZ] + 0.25*k*LocalRange[ZZ];
            child_rank = get_grid_rank(par, centre,l+1);
#ifdef TEST_BC89FLUX
    // *** debugging info ***
            cout <<"child_rank = "<<child_rank<<", pos=";
            rep.printVec("",centre,par.ndim);
            cout.flush();
    // *** debugging info ***
#endif
            if (child_rank>=0) children.push_back(child_rank);
          }
        }
      }
    }
    
    if (children.size()>0) {
      sort(children.begin(),children.end());
      children.erase( unique( children.begin(), children.end() ), children.end() );
    }
#ifdef TEST_BC89FLUX
    // *** debugging info ***
    for (size_t v=0;v<children.size();v++)
      cout <<"children["<<v<<"] = "<<children[v]<<"\n";
    // *** debugging info ***
#endif
    // set rank and xmin/xmax of child grid
    for (size_t v=0;v<children.size();v++) {
      child.rank = children[v];
      get_domain_ix(par.ndim, child.rank, ix);
      for (int i=0;i<par.ndim;i++) {
        // child level has range 0.5x my range in each direction
        child.Xmin[i] = par.levels[l+1].Xmin[i] + ix[i]*0.5*LocalRange[i];
        child.Xmax[i] = child.Xmin[i] + 0.5*LocalRange[i];
      }
      child_procs.push_back(child);
    }

#ifdef TEST_BC89FLUX
    // *** debugging info ***
    for (size_t v=0;v<child_procs.size();v++) {
      cout <<"child procs on level "<<l+1<<": rank="<<child_procs[v].rank<<"\n";
      rep.printVec("Xmin",child_procs[v].Xmin,par.ndim);
      rep.printVec("Xmax",child_procs[v].Xmax,par.ndim);
      cout.flush();
    }
    // *** debugging info ***
#endif

  } // if not on finest level grid (set children)
  
  // set rank of grids sharing a boundary with child grids, if they
  // exist.
  cgrid_ngb.resize(2*par.ndim);
  for (int d=0;d<2*par.ndim;d++) cgrid_ngb[d].resize(0);
  if (l==par.grid_nlevels-1 || nproc==1) {
    // do nothing
  }
  else {
    // for each outward normal direction on this grid, find the
    // level l+1 grids (if any) that do not intersect the grid, but
    // share a face.
    double dx=LocalRange[XX]/LocalNG[XX];
    int pdir[2];

    //
    // First get ranks of level l+1 grids into the cgrid_ngb[] vectors
    //
    if (par.ndim==1) {
      // max. one l+1 grid in each direction
      // negative direction
      centre[XX] = LocalXmin[XX] -dx;
      child_rank = get_grid_rank(par, centre,l+1);
      child.rank = child_rank;
      cgrid_ngb[XN].push_back(child);
      // positive direction
      centre[XX] = LocalXmax[XX] +dx;
      child_rank = get_grid_rank(par, centre,l+1);
      child.rank = child_rank;
      cgrid_ngb[XP].push_back(child);
    }
    else if (par.ndim==2) {
      for (int i=0;i<par.ndim;i++) {
        int nd = 2*i;
        int pd = 2*i+1;
        pdir[0] = i+1%par.ndim;
        // negative direction
        children.clear();
        centre[i] = LocalXmin[i] -dx;
        for (int p=0;p<4;p++) {
          centre[pdir[0]] = LocalXmin[pdir[0]] + 0.25*p*LocalRange[pdir[0]];
          child_rank = get_grid_rank(par, centre,l+1);
          if (child_rank>0) children.push_back(child_rank);
        }
        sort(children.begin(),children.end());
        children.erase( unique( children.begin(), children.end() ), children.end() );
        for (size_t v=0;v<children.size();v++) {
          child.rank = children[v];
          cgrid_ngb[nd].push_back(child);
        }
        // positive direction
        children.clear();
        centre[i] = LocalXmax[i] +dx;
        for (int p=0;p<4;p++) {
          centre[pdir[0]] = LocalXmin[pdir[0]] + 0.25*p*LocalRange[pdir[0]];
          child_rank = get_grid_rank(par, centre,l+1);
          if (child_rank>0) children.push_back(child_rank);
        }
        sort(children.begin(),children.end());
        children.erase( unique( children.begin(), children.end() ), children.end() );
        for (size_t v=0;v<children.size();v++) {
          child.rank = children[v];
          cgrid_ngb[pd].push_back(child);
        }
      } // dims
    } // 2D
    else {
#ifdef TEST_BC89FLUX
    // *** debugging info ***
      cout <<"3D: getting ranks for l+1 grids that share a face with my grid.\n";
      cout.flush();
    // *** debugging info ***
#endif
      for (int i=0;i<par.ndim;i++) {
        int nd = 2*i;
        int pd = 2*i+1;
        pdir[0] = (i+1)%par.ndim;
        pdir[1] = (i+2)%par.ndim;
#ifdef TEST_BC89FLUX
      // *** debugging info ***
        cout <<"3D: axis = "<<i<<", perp dirs = "<<pdir[0]<<", "<<pdir[1]<<"\n";
        cout.flush();
      // *** debugging info ***
#endif
        // negative direction
        children.clear();
        centre[i] = LocalXmin[i] -dx;
        for (int p=0;p<4;p++) {
          for (int q=0;q<4;q++) {
            centre[pdir[0]] = LocalXmin[pdir[0]] + 0.25*p*LocalRange[pdir[0]];
            centre[pdir[1]] = LocalXmin[pdir[1]] + 0.25*q*LocalRange[pdir[1]];
            child_rank = get_grid_rank(par, centre,l+1);
            if (child_rank>=0) children.push_back(child_rank);
          }
        }
        if (children.size()>0) {
          sort(children.begin(),children.end());
          children.erase( unique( children.begin(), children.end() ), children.end() );
        }
        for (size_t v=0;v<children.size();v++) {
          child.rank = children[v];
          cgrid_ngb[nd].push_back(child);
        }
        // positive direction
        children.clear();
        centre[i] = LocalXmax[i] +dx;
        for (int p=0;p<4;p++) {
          for (int q=0;q<4;q++) {
            centre[pdir[0]] = LocalXmin[pdir[0]] + 0.25*p*LocalRange[pdir[0]];
            centre[pdir[1]] = LocalXmin[pdir[1]] + 0.25*q*LocalRange[pdir[1]];
            child_rank = get_grid_rank(par, centre,l+1);
#ifdef TEST_BC89FLUX
    // *** debugging info ***
            //cout <<"+ve direction, facing child rank = "<<child_rank<<" : ";
            //rep.printVec("+ve dir, centre",centre,par.ndim);
    // *** debugging info ***
#endif
            if (child_rank>=0) children.push_back(child_rank);
          }
        }
        if (children.size()>0) {
          sort(children.begin(),children.end());
          children.erase( unique( children.begin(), children.end() ), children.end() );
        }
        for (size_t v=0;v<children.size();v++) {
          child.rank = children[v];
          cgrid_ngb[pd].push_back(child);
        }
      } // dims
    } // 3D

    // Now we have ranks in each direction, so we set xmin/xmax too
    for (int d=0;d<2*par.ndim;d++) {
      for (size_t cg=0; cg<cgrid_ngb[d].size(); cg++) {
        get_domain_ix(par.ndim, cgrid_ngb[d][cg].rank, ix);
        for (int i=0;i<par.ndim;i++) {
          // child level has range 0.5x my range in each direction
          cgrid_ngb[d][cg].Xmin[i] = par.levels[l+1].Xmin[i] + ix[i]*0.5*LocalRange[i];
          cgrid_ngb[d][cg].Xmax[i] = cgrid_ngb[d][cg].Xmin[i] + 0.5*LocalRange[i];
        }
      } // loop over grids
    } // loop over dims

#ifdef TEST_BC89FLUX
    // *** debugging info ***
    for (int d=0;d<2*par.ndim;d++) {
      for (size_t cg=0; cg<cgrid_ngb[d].size(); cg++) {
        cout <<"dir="<<d<<", cg="<<cg<<", l+1 ngb grid rank="<<cgrid_ngb[d][cg].rank<<"\n";
        rep.printVec("Xmin",cgrid_ngb[d][cg].Xmin,par.ndim);
        rep.printVec("Xmax",cgrid_ngb[d][cg].Xmax,par.ndim);
        cout.flush();
      }
    }
    // *** debugging info ***
#endif


  } // if there are child grids

  return;
}



// ##################################################################
// ##################################################################



void MCMDcontrol::get_parent_grid_info(
      struct cgrid *cg
      )
{
  cg->rank = pgrid.rank;
  for (int i=0;i<MAX_DIM;i++) cg->Xmin[i] = pgrid.Xmin[i];
  for (int i=0;i<MAX_DIM;i++) cg->Xmax[i] = pgrid.Xmax[i];
  //std::cout <<"pgrid rank = "<<pgrid.rank<<"  "<<cg->rank<<std::endl;
  return;
}



// ##################################################################
// ##################################################################



void MCMDcontrol::get_parent_ngb_grid_info(
      vector<struct cgrid>  &pgngb
      )
{
  pgngb.resize(pgrid_ngb.size());
  for (size_t iter=0; iter<pgrid_ngb.size(); iter++) {
    pgngb[iter].rank = pgrid_ngb[iter].rank;
    for (int i=0;i<MAX_DIM;i++)
      pgngb[iter].Xmin[i] = pgrid_ngb[iter].Xmin[i];
    for (int i=0;i<MAX_DIM;i++)
      pgngb[iter].Xmax[i] = pgrid_ngb[iter].Xmax[i];
  }
  return;
}



// ##################################################################
// ##################################################################



void MCMDcontrol::get_child_grid_info(
      vector<struct cgrid>  &cg
      )
{
  cg.resize(child_procs.size());
  for (size_t iter=0; iter<child_procs.size(); iter++) {
    cg[iter].rank = child_procs[iter].rank;
    for (int i=0;i<MAX_DIM;i++)
      cg[iter].Xmin[i] = child_procs[iter].Xmin[i];
    for (int i=0;i<MAX_DIM;i++)
      cg[iter].Xmax[i] = child_procs[iter].Xmax[i];
  }
  return;
}



// ##################################################################
// ##################################################################



void MCMDcontrol::get_level_lp1_ngb_info(
    vector< vector<struct cgrid> >  &cgngb
    )
{
  cgngb.resize(cgrid_ngb.size());
  for (size_t p=0; p<cgrid_ngb.size(); p++) {
    cgngb[p].resize(cgrid_ngb[p].size());
    for (size_t q=0; q<cgrid_ngb[p].size(); q++) {
      cgngb[p][q].rank = cgrid_ngb[p][q].rank;
      for (int i=0;i<MAX_DIM;i++) cgngb[p][q].Xmin[i] = cgrid_ngb[p][q].Xmin[i];
      for (int i=0;i<MAX_DIM;i++) cgngb[p][q].Xmax[i] = cgrid_ngb[p][q].Xmax[i];
    }
  }

  return;
}


// ##################################################################
// ##################################################################



//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------





