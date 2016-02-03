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


#include "MCMD_control.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include "constants.h"
#include "tools/command_line_interface.h"
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

  ReadSingleFile  =true; ///< If the ICs are in a single file, set this to true.
  WriteSingleFile =false; ///< If you want all the processors to write to one file, set this (BUGGY!)
  WriteFullImage  =false; ///< If you want multiple fits files, but each one is the full domain size (bad!), set this.
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


int MCMDcontrol::decomposeDomain()
{
  //  cout << "---MCMDcontrol::decomposeDomain() decomposing domain.\n";

  //
  // We subdivide the domain in half recursively, where the axis we cut along
  // is always the one in which the domain is longest.  In this way we minimize
  // domain interfaces, which should minimize communication.
  //
  if (SimPM.ndim==1) {
    // 1D is so simple it's worth putting in it's own section.
    LocalRange[XX] = SimPM.Range[XX]/nproc;
    LocalXmin[XX] = SimPM.Xmin[XX] + myrank*LocalRange[XX];
    LocalXmax[XX] = LocalXmin[XX] + LocalRange[XX];
    LocalNG[XX] = SimPM.NG[XX]/nproc;
    offsets[XX] = 0;
    LocalNcell = LocalNG[XX];
    // Set up ngbprocs array to point to neighbouring processors
    nx[XX] = nproc;
    ix[XX] = myrank;
    pointToNeighbours();
  } // If 1D
  
  else if (SimPM.ndim==2 || SimPM.ndim==3) {
    for (int i=0;i<SimPM.ndim;i++) {
      LocalRange[i] = SimPM.Range[i];
      LocalXmin[i]  = SimPM.Xmin[i];
      LocalXmax[i]  = SimPM.Xmax[i];
      LocalNG[i]    = SimPM.NG[i];
      LocalNcell    = SimPM.Ncell;
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
	// check if we have only one source at infinity, b/c then we decompose to keep
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
      // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
      // DOMAIN DECOMPOSITION.
      dsrc = -1;
      // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
      // DOMAIN DECOMPOSITION.

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
    if (npcounter != nproc) rep.error("nproc not a power of 2!",nproc);
    
    //
    // Now we know the range of each subcell, so calculate where I fit into
    // the hierarchy, defined by
    // \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
    // This requires myrank to count from zero!
    //
    for (int i=0;i<SimPM.ndim;i++)
      nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/LocalRange[i]);
    int temp=myrank;
    if (SimPM.ndim==3) {
      ix[ZZ] = temp/nx[XX]/nx[YY];
      temp -= ix[ZZ]*nx[XX]*nx[YY];
    }
    ix[YY] = temp/nx[XX];
    temp -= ix[YY]*nx[XX];
    ix[XX] = temp;
    
    LocalNcell = SimPM.Ncell;
    for (int i=0;i<SimPM.ndim;i++) {
      LocalXmin[i]  = SimPM.Xmin[i] +ix[i]*LocalRange[i];
      LocalXmax[i]  = SimPM.Xmin[i] +(ix[i]+1)*LocalRange[i];
      LocalNG[i]    = SimPM.NG[i]/nx[i];
      LocalNcell   /= nx[i];
      offsets[i]    = ix[i]*LocalNG[i];
    }
    // Set up ngbprocs array to point to neighbouring processors
    pointToNeighbours();    
  } // if 2D or 3D
  
  else rep.error("Bad NDIM in DecomposeDomain",SimPM.ndim);

  // Display some debugging info.
  //  if (myrank==0) {
  //    for (int i=0;i<SimPM.ndim;i++) {
  //     cout <<"Sim: idim="<<i<<"  \tRange="<<SimPM.Range[i];
  //     cout <<",\t  xmin,xmax = "<<SimPM.Xmin[i]<<", "<<SimPM.Xmax[i];
  //     cout <<"    \t Ncell = "<<SimPM.Ncell<<"\n";
  //   }
  // }
  // for (int i=0;i<SimPM.ndim;i++) {
  //   cout <<"Proc "<<myrank<<": idim="<<i<<"  \tRange="<<LocalRange[i];
  //   cout <<",\t  xmin,xmax = "<<LocalXmin[i]<<", "<<LocalXmax[i];
  //   cout <<"    \t Ncell = "<<LocalNcell;
  //   cout <<"\t neighbours : "<<ngbprocs[2*i]<<", "<<ngbprocs[2*i+1]<<"\n";
  // }
  // cout << "---MCMDcontrol::decomposeDomain() Domain decomposition done.\n\n";
  return(0);
}


// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
    const enum axes daxis   ///< Axis to decompose domain along.
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
      LocalRange[i] = SimPM.Range[i];
      LocalXmin[i]  = SimPM.Xmin[i];
      LocalXmax[i]  = SimPM.Xmax[i];
      LocalNG[i]    = SimPM.NG[i];
      LocalNcell    = SimPM.Ncell;
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
      rep.error("nproc not a power of 2!",nproc);
    }
    
    //
    // Now we know the range of each subcell, so calculate where I fit into
    // the hierarchy, defined by
    // \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
    // This requires myrank to count from zero!
    //
    for (int i=0;i<SimPM.ndim;i++)
      nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/LocalRange[i]);
    int temp=myrank;
    if (SimPM.ndim==3) {
      ix[ZZ] = temp/nx[XX]/nx[YY];
      temp -= ix[ZZ]*nx[XX]*nx[YY];
    }
    ix[YY] = temp/nx[XX];
    temp -= ix[YY]*nx[XX];
    ix[XX] = temp;
    
    LocalNcell = SimPM.Ncell;
    for (int i=0;i<SimPM.ndim;i++) {
      LocalXmin[i]  = SimPM.Xmin[i] +ix[i]*LocalRange[i];
      LocalXmax[i]  = SimPM.Xmin[i] +(ix[i]+1)*LocalRange[i];
      LocalNG[i]    = SimPM.NG[i]/nx[i];
      LocalNcell   /= nx[i];
      offsets[i]    = ix[i]*LocalNG[i];
    }
    // Set up ngbprocs array to point to neighbouring processors
    pointToNeighbours();    
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

int MCMDcontrol::pointToNeighbours()
{
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
    nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/LocalRange[i]);
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
  if (pconst.equalD(LocalXmin[XX],SimPM.Xmin[XX])) ngbprocs[XN] = -999;
  if (pconst.equalD(LocalXmax[XX],SimPM.Xmax[XX])) ngbprocs[XP] = -999;
  if (SimPM.ndim >1) {
    ngbprocs[YN] = myrank -nx[XX];
    ngbprocs[YP] = myrank +nx[XX];
    if (pconst.equalD(LocalXmin[YY],SimPM.Xmin[YY])) ngbprocs[YN] = -999;
    if (pconst.equalD(LocalXmax[YY],SimPM.Xmax[YY])) ngbprocs[YP] = -999;
  }
  if (SimPM.ndim >2) {
    ngbprocs[ZN] = myrank -nx[XX]*nx[YY];
    ngbprocs[ZP] = myrank +nx[XX]*nx[YY];
    if (pconst.equalD(LocalXmin[ZZ],SimPM.Xmin[ZZ])) ngbprocs[ZN] = -999;
    if (pconst.equalD(LocalXmax[ZZ],SimPM.Xmax[ZZ])) ngbprocs[ZP] = -999;
  }
  return(0);
}

// ##################################################################
// ##################################################################

///
/// Get a list of all abutting domains, including corner/edge intersections.
///
void MCMDcontrol::get_abutting_domains(
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
  for (int v=0;v<SimPM.ndim;v++) {
    negdir = static_cast<direction>(2*v);
    posdir = static_cast<direction>(2*v+1);
    if (ngbprocs[negdir]>=0) d[negdir] = true;
    if (ngbprocs[posdir]>=0) d[posdir] = true;
  }

  //
  // Add coordinate directions to list.
  //
  for (int v=0; v<2*SimPM.ndim; v++) {
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

//
// Returns the ix array for any requested rank.
//
void MCMDcontrol::get_domain_ix(const int r, ///< rank ix requested for.
				   int *arr     ///< array to put ix into.
				   )
{
  int temp=r;
  if (SimPM.ndim==3) {
    arr[ZZ] = temp/nx[XX]/nx[YY];
    temp -= arr[ZZ]*nx[XX]*nx[YY];
  }
  if (SimPM.ndim>1) {
    arr[YY] = temp/nx[XX];
    temp -= arr[YY]*nx[XX];
  }
  arr[XX] = temp;
  return;
}

// ##################################################################
// ##################################################################
//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------





