/// \file uniform_grid_pllel.cc
/// \author Jonathan Mackey
/// 
/// Definitions for parallel uniform grid.
/// 
/// Modifications:\n
/// - 2007-11-06 started to get it going.
/// - 2010-01-22 JM: Changed RT source finding if we have non-cell-centred sources.
/// - 2010-01-26 JM: Debugging 3D RT for vertex-centred sources.
/// - 2010-02-03 JM: changed variable names in destructor ('i' was defined twice!)
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; Added oneway-outflow BC
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.  Switched 'endl' statements for '\n'.
/// - 2010.11.19 JM: Got rid of testing myalloc() myfree() functions.
/// - 2010.12.07 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  Will probably keep it since it is the way
///   things will go eventually.  The new grid classes have extra
///   function for the distance between two points, two cells, and
///   between a vertex and a cell.
/// - 2011.02.24 JM: Worked on generalising multi-source radiative transfer.
///    Still some way to go to getting it functional.  Wraps old setup() 
///    function for point sources into a new function.  The setup function now
///    calls either this new function, or another new one for sources at infinity.
/// - 2011.02.25 JM: removed HCORR ifdef around new code. 
///     Added possibility of multiple sources in send/recv.
/// - 2011.02.28 JM: Got rid of RSP radiation-sources-parameters class.
/// - 2011.03.02 JM: Fixed bugs in diffuse radiation boundary comms.
/// - 2011.03.21 JM: Got rid of zero-ing of column densities.  It is done in
///    the cell constructor, so we don't need to worry about it.
/// - 2011.03.22 JM: small bugfixes.
/// - 2011.04.22 JM: Fixed setup_recv_boundaries() functions so that new boundary cells
///    are only created if they don't already exist.  It now works for multiple point
///    sources, and I don't think they need to be at the same place.
/// - 2013.09.05 JM: Debugged for new get/set col functions.
/// - 2015.01.28 JM: Removed GEOMETRIC_GRID (it is default now), and
///    updated for the new code structure.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"

#include "global.h"
#include "grid/uniform_grid.h"
#include "sim_control_MPI.h"
#include "microphysics/microphysics_base.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef PARALLEL



// ##################################################################
// ##################################################################



UniformGridParallel::UniformGridParallel(
      int nd,
      int nv,
      int eqt,
      double *xn,
      double *xp,
      int *nc,
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,xn,xp,nc)
{
  //cout <<"UniformGridParallel constructor.\n";
  //rep.printVec("Local Xmin",xn,nd);
  //rep.printVec("Local Xmax",xp,nd);
  //rep.printVec("Local Npt ",nc,nd);

#ifdef PLLEL_RT
  //RT_nbd=0;
  //RT_nbc=0;
  //  RT_bd=0;
#endif // PLLEL_RT

  //
  // If we need to know about the global size of the simulation:
  //
  SIM_ixmin  = mem.myalloc(SIM_ixmin, G_ndim); 
  SIM_ixmax  = mem.myalloc(SIM_ixmax, G_ndim); 
  SIM_irange = mem.myalloc(SIM_irange,G_ndim); 
  //
  // Set integer dimensions/location of grid.
  //
  CI.get_ipos_vec(SimPM.Xmin, SIM_ixmin );
  CI.get_ipos_vec(SimPM.Xmax, SIM_ixmax );
  for (int v=0;v<G_ndim;v++)
    SIM_irange[v] = SIM_ixmax[v]-SIM_ixmin[v];
  //rep.printVec("SIM iXmin ", SIM_ixmin, G_ndim);
  //rep.printVec("SIM iXmax ", SIM_ixmax, G_ndim);
  //rep.printVec("SIM iRange", SIM_irange,G_ndim);
  
  //
  // Set pointer to muli-core class.
  //
  mpiPM = p;

  return;
}



// ##################################################################
// ##################################################################



UniformGridParallel::~UniformGridParallel()
{
#ifdef TESTING
  //cout <<"UniformGridParallel destructor.\n";
#endif
#ifdef PLLEL_RT
  //
  // Delete cells created for RT receive boundaries.
  // Have to delete lists of cells, and boundary_data, and maybe cells.
  //
  for (unsigned int isrc=0;isrc<RT_source_list.size(); isrc++) {
    struct RT_source_comms_info *rc = &(RT_source_list[isrc]);
    if (!rc->RT_recv_list.empty()) {
      struct boundary_data *b;
      for (unsigned int ii=0; ii<rc->RT_recv_list.size(); ii++) {
        int n=rc->RT_recv_list[ii].dir;
        b = rc->RT_recv_list[ii].RT_bd;
        //
        // If a pre-existing boundary, just delete the list.
        //
        if (n==dir_XN || n==dir_XP || n==dir_YN || n==dir_YP || n==dir_ZN || n==dir_ZP) {
          //cout <<"\tNot deleting standard boundary data for dir: "<<n<<"\n";
          if      (!b) {
            //cout <<"\t    Zero boundary data pointer for dir: "<<n<<"\n";
          }
          else if (!b->data.empty()) {
            //cout <<"\t    Emptying list.\n";
            b->data.clear(); // No dynamically allocated memory in these boundaries.
          }
          else {
            //cout <<"\t    Empty boundary list (!)(?) Probably shouldn't happen.\n";
          }
        } // if pre-existing.
        
        //
        // Else we created the cells too, so we need to delete them.
        //
        else {
          //cout <<"\t    Deleting  created boundary data for dir: "<<n<<"\n";
          //cout <<"\t    Deleting created  boundary data for dir: "<<n<<"\n";
          if      (!b) {
            //cout <<"\t\tZero boundary data pointer for dir: "<<n<<"\n";
          }
          else if (b->data.empty()) {
            //cout <<"\t\tRT_BC destructor: No boundary cells to delete.\n";
          }
          else {
            //cout <<"\t    Deleting standard boundary data for dir: "<<n<<"\n";
            list<cell *>::iterator i=b->data.begin();
            do {
              //cout <<"\t\tDeleting an RT boundary cell...\n";
              if (isrc==0) {
                //
                // only delete the actual cell in the first pass through.  So the first 
                // source should always be the ionising source!
                //
                deleteCell(*i);  // This seems to work in terms of actually freeing the memory.
              }
              b->data.erase(i);
              i=b->data.begin();
            }  while(i!=b->data.end());
            if(b->data.empty()) {
              //cout <<"\t\tDeleting RT Boundary data done.\n";
            }
            else {
              //cout <<"\t not empty list! FIX ME!!!\n";
            }
          }
        }

        //cout <<"\t\tNow deleting boundary data...\n";
        if (b) {
          b = mem.myfree(b);
          //RT_recv_list[ii].RT_bd = mem.myfree(RT_recv_list[ii].RT_bd);
        }
        //cout <<"\t\tFinished with this Receive boundary\n";
      }
      rc->RT_recv_list.clear();
    }
    
    //
    // Delete send lists.  No cells to delete, just the lists, and boundary_data
    //
    if (!rc->RT_send_list.empty()) {
      //cout <<"\tDeleting RT Send List.\n";
      struct boundary_data *b;
      for (unsigned int ii=0; ii<rc->RT_send_list.size(); ii++) {
        b = rc->RT_send_list[ii].RT_bd;
        if      (!b)  {
          //cout <<"\t\tZero send boundary data pointer for dir: "<<RT_send_list[i].dir<<"\n";
        }
        else if (!b->data.empty()) {
          //cout <<"\t\tEmptying send list.\n";
          b->data.clear(); // No dynamically allocated memory in these boundaries.
        }
        else {
          //cout <<"\t\tEmpty send boundary list (!)(?) Probably shouldn't happen.\n";
        }
        b = mem.myfree(b);
      }
      rc->RT_send_list.clear();
      //cout <<"\tDeleting RT Send List... done!\n";
    }
    
  }
#endif // PLLEL_RT

  SIM_ixmin  = mem.myfree(SIM_ixmin);
  SIM_ixmax  = mem.myfree(SIM_ixmax);
  SIM_irange = mem.myfree(SIM_irange);

  return;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_setBCtypes(string bctype)
{
  string fname="UniformGridParallel::BC_setBCtypes";
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif 
  if(bctype=="FIXED" || bctype=="PERIODIC" || bctype=="ABSORBING") {
#ifdef TESTING
    cout <<"using old-style boundary condition specifier: "<<bctype<<" on all sides.\n";
    cout <<"Converting to new style specifier.\n";
#endif 
    
    if (bctype=="FIXED") {
      bctype = "XNfix_XPfix_";
      if (G_ndim>1) bctype += "YNfix_YPfix_";
      if (G_ndim>2) bctype += "ZNfix_ZPfix_";
    }
    else if (bctype=="PERIODIC") {
      bctype = "XNper_XPper_";
      if (G_ndim>1) bctype += "YNper_YPper_";
      if (G_ndim>2) bctype += "ZNper_ZPper_";
    }
    else if (bctype=="ABSORBING") {
      bctype = "XNout_XPout_";
      if (G_ndim>1) bctype += "YNout_YPout_";
      if (G_ndim>2) bctype += "ZNout_ZPout_";
    }
#ifdef TESTING
    cout <<"New bctype string = "<<bctype<<"\n";
#endif 
  }
  // Set up boundary data struct.  Assumes bctype is in format "XPper XNper " etc.,
  // so that each boundary is defined by 6 characters, and number of boundaries is
  // given by length/6.
  int len = bctype.length(); len = (len+5)/6;
  if (len < 2*G_ndim) rep.error("Need boundaries on all sides!",len);
#ifdef TESTING
  cout <<"Got "<<len<<" boundaries to set up.\n";
#endif 

  int i=0; string::size_type pos; 
  if (len > 2*G_ndim) {
#ifdef TESTING
    cout <<"\t checking if I need extra boundaries in this subdomain.\n";
#endif 
    // Go through all the internal boundaries I know about, and see if the
    // subdomain needs to know about them.
    for (i=2*G_ndim; i<=len; i++) {
      if ( (pos=bctype.find("INjet")) !=string::npos) {
        // Jet always comes in from the XN boundary, at the YN corner in 2d,
        // and at the centre of the boundary in 3d.
        // So in 2D I will always assume that proc 0 only needs to know about the jet
        // i.e. a processor always spans the jet width.
        // In 3D it's more complicated, and perhaps all processors on the XN boundary
        // should keep the jet BC, even though some might be empty.
        if (G_ndim==2) {
          if (mpiPM->get_myrank() !=0) {
            cout <<"Removing jet bc; bctype = "<<bctype<<"\n";
      bctype.erase(pos,6);
#ifdef TESTING
      cout <<"Removed jet bc;  bctype = "<<bctype<<"\n";
#endif 
      len -= 1;
    }
        }
  else if (G_ndim==3) {
    if (!pconst.equalD(SimPM.Xmin[XX], mpiPM->LocalXmin[XX])) {
#ifdef TESTING
      cout <<"Removing jet bc; bctype = "<<bctype<<"\n";
#endif 
      bctype.erase(pos,6);
#ifdef TESTING
      cout <<"Removed jet bc;  bctype = "<<bctype<<"\n";
#endif 
      len -= 1;
    }
        }
        else rep.error("Bad ndim for jet simulation in BC_setBCtypes()",G_ndim);
      }
      else if ( (pos=bctype.find("INdm2")) !=string::npos) {
  // This is a test problem (double mach reflection) with hard-coded
  // boundaries, and this boundary is along the YN boundary between
  // x=[0,1/6].
  if ( (!pconst.equalD(SimPM.Xmin[YY],mpiPM->LocalXmin[YY])) 
       || (mpiPM->LocalXmin[XX]>1./6.) ) {
#ifdef TESTING
    cout <<"Removing dmr2 bc; bctype = "<<bctype<<"\n";
#endif 
    bctype.erase(pos,6);
#ifdef TESTING
    cout <<"Removed dmr2 bc;  bctype = "<<bctype<<"\n";
#endif 
    len -= 1;
  }
      }
      else if ( (pos=bctype.find("INwnd")) !=string::npos) {
  //
  // This is a stellar wind problem -- it only sets cells which
  // are on the domain to be part of the wind, so we can let
  // every proc set it up without any trouble.
  //
      }
      else rep.error("Couldn't find boundary condition for extra boundary condition","IN");
    } // loop through extra boundaries.
  } // checking if extra boundaries are needed.
  
  i=0; enum axes dir;
  string temp;
  // Loop through boundaries, and if local boundary is not sim boundary,
  // set it to be a parallel boundary.
  for (i=0; i<G_ndim; i++) {
    dir = static_cast<axes>(i);
    if (!pconst.equalD(G_xmin[i], SimPM.Xmin[i])) {
      // local xmin is not Sim xmin, so it's an mpi boundary
      if      (dir==XX) temp = "XNmpi_";
      else if (dir==YY) temp = "YNmpi_";
      else if (dir==ZZ) temp = "ZNmpi_";
      else rep.error("Bad axis!",dir);
      bctype.replace(2*i*6,6,temp); cout <<"new bctype="<<bctype<<"\n";
    }
    if (!pconst.equalD(G_xmax[i], SimPM.Xmax[i])) {
      // local xmax is not Sim xmin, so it's an mpi boundary
      if      (dir==XX) temp = "XPmpi_";
      else if (dir==YY) temp = "YPmpi_";
      else if (dir==ZZ) temp = "ZPmpi_";
      else rep.error("Bad axis!",dir);
      bctype.replace((2*i+1)*6,6,temp); cout <<"new bctype="<<bctype<<"\n";
    }
  }
  
  // Now set up boundary data structs.
  BC_bd=0; 
  BC_bd = mem.myalloc(BC_bd,len);
  UniformGrid::BC_nbd = len;

  // Now go through the domain edge boundaries and assign them.
  string d[6] = {"XN","XP","YN","YP","ZN","ZP"};
  for (i=0; i<2*G_ndim; i++) {
    BC_bd[i].dir = static_cast<direction>(i); //XN=0,XP=1,YN=2,YP=3,ZN=4,ZP=5
    if ( (pos=bctype.find(d[i])) == string::npos)
      rep.error("Couldn't find boundary condition for ",d[i]);
    BC_bd[i].type = bctype.substr(pos+2,3);
    if      (BC_bd[i].type=="per") {BC_bd[i].itype=PERIODIC; BC_bd[i].type="PERIODIC";}
    else if (BC_bd[i].type=="out" ||
       BC_bd[i].type=="abs") {BC_bd[i].itype=OUTFLOW; BC_bd[i].type="OUTFLOW";}
    else if (BC_bd[i].type=="owo") {BC_bd[i].itype=ONEWAY_OUT; BC_bd[i].type="ONEWAY_OUT";}
    else if (BC_bd[i].type=="inf") {BC_bd[i].itype=INFLOW ; BC_bd[i].type="INFLOW";}
    else if (BC_bd[i].type=="ref") {BC_bd[i].itype=REFLECTING; BC_bd[i].type="REFLECTING";}
    else if (BC_bd[i].type=="jrf") {BC_bd[i].itype=JETREFLECT; BC_bd[i].type="JETREFLECT";}
    else if (BC_bd[i].type=="fix") {BC_bd[i].itype=FIXED; BC_bd[i].type="FIXED";}
    else if (BC_bd[i].type=="dmr") {BC_bd[i].itype=DMACH; BC_bd[i].type="DMACH";}
    else if (BC_bd[i].type=="mpi") {BC_bd[i].itype=BCMPI; BC_bd[i].type="BCMPI";}
    else if (BC_bd[i].type=="sb1") {
      BC_bd[i].itype=STARBENCH1;
      BC_bd[i].type="STARBENCH1";  // Wall for Tremblin mixing test.
    }
    else rep.error("Don't know this BC type",BC_bd[i].type);
    
    if(!BC_bd[i].data.empty()) rep.error("Boundary data not empty in constructor!",BC_bd[i].data.size());
    BC_bd[i].refval=0;
#ifdef TESTING
    cout <<"\tBoundary type "<<i<<" is "<<BC_bd[i].type<<"\n";
#endif 
  }

  // If there are any extra boundaries, deal with them too.
  if (i<BC_nbd) {
#ifdef TESTING
    cout <<"Got "<<i<<" boundaries, but have "<<BC_nbd<<" boundaries.\n";
    cout <<"Must have extra BCs... checking if I know what they are.\n";
#endif 
    do {
      BC_bd[i].dir = NO;
      // search string, starting at position 6i.
      // This ensures if we have more than one internal boundary, that we
      // don't keep finding the first one...
      if ( (pos=bctype.find("IN",i*6)) ==string::npos)
  rep.error("Couldn't find boundary condition for extra boundary condition","IN");
      BC_bd[i].type = bctype.substr(pos+2,3);
      if      (BC_bd[i].type=="jet") {BC_bd[i].itype=JETBC; BC_bd[i].type="JETBC";}
      else if (BC_bd[i].type=="dm2") {BC_bd[i].itype=DMACH2; BC_bd[i].type="DMACH2";}
      else if (BC_bd[i].type=="wnd") {BC_bd[i].itype=STWIND;   BC_bd[i].type="STWIND";}
      else rep.error("Don't know this BC type",BC_bd[i].type);
      if(!BC_bd[i].data.empty()) rep.error("Boundary data not empty in constructor!",BC_bd[i].data.size());
      BC_bd[i].refval=0;
#ifdef TESTING
      cout <<"\tBoundary type "<<i<<" is "<<BC_bd[i].type<<"\n";
#endif 
      i++;
    } while (i<BC_nbd);
  }
#ifdef TESTING
  cout <<"BC types and data set up.\n";
#endif 

  //
  // If we have periodic boundaries, need to set neighbouring processors to
  // wrap around.  So set the number of procs in each direction.
  //
  int nx[SimPM.ndim];
  for (i=0;i<G_ndim;i++)
    nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/mpiPM->LocalRange[i]);
  for (i=0; i<2*G_ndim; i++) {
    if (BC_bd[i].itype == PERIODIC) {
      switch (i) {
       case XN:
  mpiPM->ngbprocs[XN] = mpiPM->get_myrank() +nx[XX] -1; break;
       case XP:
  mpiPM->ngbprocs[XP] = mpiPM->get_myrank() -nx[XX] +1; break;
       case YN:
  mpiPM->ngbprocs[YN] = mpiPM->get_myrank() +(nx[YY]-1)*nx[XX]; break;
       case YP:
  mpiPM->ngbprocs[YP] = mpiPM->get_myrank() -(nx[YY]-1)*nx[XX]; break;
       case ZN:
  mpiPM->ngbprocs[ZN] = mpiPM->get_myrank() +(nx[ZZ]-1)*nx[YY]*nx[XX]; break;
       case ZP:
  mpiPM->ngbprocs[ZP] = mpiPM->get_myrank() -(nx[ZZ]-1)*nx[YY]*nx[XX]; break;
       default:
  rep.error("UniformGridParallel::BC_setBCtypes: Bad direction",i); break;
      } // set neighbour according to direction.
      if ( (mpiPM->ngbprocs[i]<0) || (mpiPM->ngbprocs[i]>=mpiPM->get_nproc()) )
  rep.error("UniformGridParallel::BC_setBCtypes: Bad periodic neighbour",mpiPM->ngbprocs[i]);
      if (mpiPM->ngbprocs[i] == mpiPM->get_myrank()) {
  //  cout <<"UniformGridParallel::BC_setBCtypes: only one proc in dir [i]: "<<i<<"\n";
  //  cout <<"UniformGridParallel::BC_setBCtypes: periodic on single proc, so setting ngb to -999.\n";
  mpiPM->ngbprocs[i] = -999;
      }
    } // if periodic  
#ifdef TESTING
    cout<<"Neighbouring processor in dir "<<i<<" = "<<mpiPM->ngbprocs[i]<<"\n";
#endif // TESTING
  } // loop over directions.
  return(0);
}



// ##################################################################
// ##################################################################



int UniformGridParallel::SetupBCs(int Nbc, string typeofbc)
{
  BC_nbc = Nbc;
  if (BC_setBCtypes(typeofbc) !=0)
    rep.error("UniformGrid::setupBCs() Failed to set types of BC",1);
  
  // This function loops through all points on the grid, and if it is an
  // edge point, it finds which direction(s) are off the grid and sets up
  // boundary cell(s) in that direction.
  class cell *c, *temppt;
  enum direction dir; int err=0;
  c = FirstPt();
  //  cout <<"Assigning edge data.\n";
  do {
    if (c->isedge !=0) {
      for (int i=0;i<2*G_ndim;i++) {
        dir = static_cast<enum direction>(i);
         temppt=NextPt(c,dir);
         if (temppt==0) {
           c->ngb[dir] = BC_setupBCcells(c, dir, G_dx, 1);
           if (c->ngb[dir]==0)
             rep.error("Failed to set up boundary cell",c->ngb[dir]);
         } // if dir points off the grid.
      } // Loop over all neighbours
    } // if an edge point.
    // If I use internal boundaries in the future:
    // if (c->isbd) bc->assignIntBC(c);
  } while ( (c=NextPt(c)) !=0);
  if (err!=0) rep.error("Failed to set up Boundary cells",err);
  // cout <<"Done.\n";
  // Loop through all boundaries, and assign data to them.
  for (int i=0; i<BC_nbd; i++) {
    switch (BC_bd[i].itype) {
     case PERIODIC:   err += BC_assign_PERIODIC(  &BC_bd[i]); break;
     case OUTFLOW:    err += BC_assign_OUTFLOW(   &BC_bd[i]); break;
     case ONEWAY_OUT: err += BC_assign_ONEWAY_OUT(&BC_bd[i]); break;
     case INFLOW:     err += BC_assign_INFLOW(    &BC_bd[i]); break;
     case REFLECTING: err += BC_assign_REFLECTING(&BC_bd[i]); break;
     case FIXED:      err += BC_assign_FIXED(     &BC_bd[i]); break;
     case JETBC:      err += BC_assign_JETBC(     &BC_bd[i]); break;
     case JETREFLECT: err += BC_assign_JETREFLECT(&BC_bd[i]); break;
     case DMACH:      err += BC_assign_DMACH(     &BC_bd[i]); break;
     case DMACH2:     err += BC_assign_DMACH2(    &BC_bd[i]); break;
     case BCMPI:      err += BC_assign_BCMPI(&BC_bd[i],BC_MPItag); break;
     case STWIND:     err += BC_assign_STWIND(    &BC_bd[i]); break;
     case STARBENCH1: err += BC_assign_STARBENCH1(&BC_bd[i]); break;
     default:
      rep.warning("Unhandled BC",BC_bd[i].itype,-1); err+=1; break;
    }
    if (i==XP || (i==YP && G_ndim>=2) || (i==ZP && G_ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      COMM->barrier("UniformGridParallel__SetupBCs");
    }
  }
  return(err);
}



// ##################################################################
// ##################################################################



int UniformGridParallel::TimeUpdateExternalBCs(const int cstep, const int maxstep)
{
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<BC_nbd;i++) {
    b = &BC_bd[i];
    switch (b->itype) {
    case PERIODIC:   err += BC_update_PERIODIC(   b, cstep, maxstep); break;
    case OUTFLOW:    err += BC_update_OUTFLOW(    b, cstep, maxstep); break;
    case ONEWAY_OUT: err += BC_update_ONEWAY_OUT( b, cstep, maxstep); break;
    case INFLOW:     err += BC_update_INFLOW(     b, cstep, maxstep); break;
    case REFLECTING: err += BC_update_REFLECTING( b, cstep, maxstep); break;
    case FIXED:      err += BC_update_FIXED(      b, cstep, maxstep); break;
    case JETBC:      err += BC_update_JETBC(      b, cstep, maxstep); break;
    case JETREFLECT: err += BC_update_JETREFLECT( b, cstep, maxstep); break;
    case DMACH:      err += BC_update_DMACH(      b, cstep, maxstep); break;
    case DMACH2:     err += BC_update_DMACH2(     b, cstep, maxstep); break;
    case BCMPI:      err += BC_update_BCMPI(      b, cstep, maxstep,BC_MPItag); break;
    case STARBENCH1: err += BC_update_STARBENCH1( b, cstep, maxstep); break;
    case STWIND: case RADSHOCK: case RADSH2:
      //
      // These are updated in TimeUpdateInternalBCs
      //
      break;
    default:
       rep.warning("Unhandled BC: parallel update external",b->itype,-1); err+=1; break;
    }
    if (i==XP || (i==YP && G_ndim>=2) || (i==ZP && G_ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      //      cout <<"Barrier synching for dir="<<i<<" with 1=XP,2=YP,3=ZP.\n";
      ostringstream tmp; tmp <<"UniformGridParallel_TimeUpdateExternalBCs_"<<i;
      COMM->barrier(tmp.str());
      tmp.str("");
    }
  }
  return 0;
} // time update BCs.



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_assign_PERIODIC(  boundary_data *b)
{
  // For parallel grid, periodic data is on a different processor,
  // which is already pointed to by mpiPM->ngbprocs[b->dir]
  // So I just have to call BC_assign_BCMPI and it will do the job.
  int err=0;
  if (mpiPM->ngbprocs[b->dir] <0) {
    //    cout <<"BC_assign_PERIODIC: non comm periodic in direction "<<b->dir<<"\n";
    err = UniformGrid::BC_assign_PERIODIC(b);
  }
  else {
    //    cout<<"BC_assign_PERIODIC: communicating periodic bc in direction "<<b->dir<<"\n";
    //    cout<<"BC_assign_PERIODIC: calling mpi assign BC function\n";
    err = UniformGridParallel::BC_assign_BCMPI(b,BC_PERtag);
  }
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_update_PERIODIC(   struct boundary_data *b,
                 const int cstep,
                 const int maxstep
                 )
{
  // For parallel grid, periodic data is on a different processor,
  // which is already pointed to by mpiPM->ngbprocs[b->dir]
  // So I just have to call BC_update_BCMPI and it will do the job.
  int err=0;
  if (mpiPM->ngbprocs[b->dir] <0) {
    //    cout <<"BC_update_PERIODIC: non-communicating periodic BC in direction "<<b->dir<<"\n";
    err = UniformGrid::BC_update_PERIODIC(b,cstep,maxstep);
  }
  else {
    //    cout<<"BC_update_PERIODIC: communicating periodic BC in direction "<<b->dir<<"\n";
    //    cout<<"BC_update_PERIODIC: calling mpi update BC function\n";
    err = UniformGridParallel::BC_update_BCMPI(b,cstep,maxstep,BC_PERtag);
  }
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_assign_BCMPI(boundary_data *b,
           int comm_tag
           )
{
  // first choose cells to send and send this boundary to appropriate
  // processor.
  list<cell *> cells;
  int ncell=0;
  int err  =0;
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: sending data in dir: "<<b->dir<<"\n";
#endif 
  err += comm_select_data2send(&cells, &ncell, b->dir);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: sending "<<ncell<<" cells.  Boundary data contains "<<b->data.size()<<" cells.\n";
#endif 

  //
  // New stuff: send the data.
  //
  string send_id;
#ifdef TESTING
  cout <<"BC_assign_BCMPI: sending data...\n";
#endif 
  err += COMM->send_cell_data(mpiPM->ngbprocs[b->dir], // to_rank
            &cells,        // cells list.
            ncell,        // number of cells.
            send_id,     // identifier for send.
            comm_tag
            );
  if (err) rep.error("sending data failed",err);

  //
  // New stuff: receive data.
  //
#ifdef TESTING
  cout <<"BC_assign_BCMPI: looking for data to receive...\n";
#endif 
  string recv_id; int recv_tag=-1; int from_rank=-1;
  err = COMM->look_for_data_to_receive(&from_rank, ///< rank of sender
               recv_id,    ///< identifier for receive.
               &recv_tag,  ///< comm_tag associated with data.
               COMM_CELLDATA ///< type of data we want.
               );
  if (err) rep.error("look for cell data failed",err);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: got data to receive from rank: "<<from_rank<<"\n";
#endif 

  //
  // See what direction we are getting data from:
  // This is complicated for the case where we get an internal and a periodic boundary
  // from the same process -- we need to decide which one we are getting.
  //
#ifdef TESTING
  cout <<"BC_assign_BCMPI: associating direction with rank.\n";
#endif 
  enum direction dir=NO;
  for (int i=0; i<2*G_ndim; i++) {
    if (from_rank == mpiPM->ngbprocs[i]) {
      if (dir==NO) dir = static_cast<direction>(i);
      else {
  //
  // must be receiving periodic and internal boundary from same proc.
  // In this case we need to decide which it is, and whether the source
  // direction is the positive or negative direction.
  // N.B. if none of the conditions below are satisfied, then the earlier 
  // assignment of dir was the correct one.  This is just for if the first
  // assignment was the wrong one.
  //
      if (recv_tag == BC_PERtag) {
    // Periodic boundary, so in neg.dir, from_rank should be > myrank, and vice versa.
    if      (from_rank>mpiPM->get_myrank() && (i==XN || i==YN || i==ZN)) dir = static_cast<direction>(i);
    else if (from_rank<mpiPM->get_myrank() && (i==XP || i==YP || i==ZP)) dir = static_cast<direction>(i);
  }
  else if (recv_tag == BC_MPItag) {
    // Internal boundary, so in neg.dir. from_rank should be < myrank.
    if      (from_rank<mpiPM->get_myrank() && (i==XN || i==YN || i==ZN)) dir = static_cast<direction>(i);
    else if (from_rank>mpiPM->get_myrank() && (i==XP || i==YP || i==ZP)) dir = static_cast<direction>(i);
  }
  else rep.error("bad tag",recv_tag);
      }
    }
  }
  if (dir==NO) rep.error("Message is not from a neighbour!",from_rank);
#ifdef TESTING
  cout <<mpiPM->get_myrank()<<"\tBC_assign_BCMPI: Receiving Data type ";
  cout <<recv_tag<<" from rank: "<<from_rank<<" from direction "<<dir<<"\n";
#endif 

  //
  // Now choose the boundary data associated with boundary we are receiving:
  //
  struct boundary_data *recv_b = 0;
  recv_b = &BC_bd[static_cast<int>(dir)];

  //
  // Receive the data:
  //
#ifdef TESTING
  cout <<"BC_assign_BCMPI: calling COMM->receive_cell_data().\n";
#endif 
  err = COMM->receive_cell_data(from_rank,  ///< rank of process we are receiving from.
        &(recv_b->data),     ///< list of cells to get data for. 
        recv_b->data.size(), ///< number of cells in list (extra checking!)
        recv_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
        recv_id      ///< identifier for receive, for any book-keeping that might be needed.
        );
  if (err) rep.error("COMM->receive_cell_data() returned error",err);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: returned from COMM->receive_cell_data().\n";
#endif 

  //
  // Set P=Ph for received boundary.
  //
  struct boundary_data *b2 = &BC_bd[static_cast<int>(dir)];
  list<cell*>::iterator c=b2->data.begin();
  for (c=b2->data.begin(); c!=b2->data.end(); ++c) {
    for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
  } // all cells.
    
  //
  // Wait until send is finished, and then return.
  //
  //  cout <<"BC_update_BCMPI: waiting for finish "<<b->dir<<"\n";
  err = COMM->wait_for_send_to_finish(send_id);
  if (err) rep.error("Waiting for send to complete!", err);
  //cout <<"BC_update_BCMPI: send and receive finished.\n";
  //cout <<"BC_update_BCMPI: Sent BC "<<b->dir<<", received BC "<<dir<<"\n";
  //cout <<"*******************************************\n\n";

  return 0;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_update_BCMPI(boundary_data *b,
           const int cstep,
           const int maxstep,
           int comm_tag
           )
{
  // first choose cells to send and send this boundary to appropriate
  // processor.
  list<cell *> cells;
  int ncell=0;
  int err  =0;
  //cout <<"*******************************************\n";
  //cout <<"BC_update_BCMPI: sending data in dir: "<<b->dir<<"\n";
  err += comm_select_data2send(&cells, &ncell, b->dir);
  //  cout <<"BC_update_BCMPI: sending "<<ncell<<" cells.  Boundary data contains "<<b->data.size()<<" cells.\n";

  //
  // New stuff: send the data.
  //
  string send_id;
  err += COMM->send_cell_data(mpiPM->ngbprocs[b->dir], // to_rank
            &cells,        // cells list.
            ncell,        // number of cells.
            send_id,     // identifier for send.
            comm_tag
            );
  if (err) rep.error("sending data failed",err);

  //
  // New stuff: receive data.
  //
  string recv_id; int recv_tag=-1; int from_rank=-1;
  err = COMM->look_for_data_to_receive(&from_rank, ///< rank of sender
               recv_id,    ///< identifier for receive.
               &recv_tag,  ///< comm_tag associated with data.
               COMM_CELLDATA ///< type of data we want.
               );
  if (err) rep.error("look for cell data failed",err);

  //
  // See what direction we are getting data from:
  // This is complicated for the case where we get an internal and a periodic boundary
  // from the same process -- we need to decide which one we are getting.
  //
  enum direction dir=NO;
  for (int i=0; i<2*G_ndim; i++) {
    if (from_rank == mpiPM->ngbprocs[i]) {
      if (dir==NO) dir = static_cast<direction>(i);
      else {
  //
  // must be receiving periodic and internal boundary from same proc.
  // In this case we need to decide which it is, and whether the source
  // direction is the positive or negative direction.
  // N.B. if none of the conditions below are satisfied, then the earlier 
  // assignment of dir was the correct one.  This is just for if the first
  // assignment was the wrong one.
  //
      if (recv_tag == BC_PERtag) {
    // Periodic boundary, so in neg.dir, from_rank should be > myrank, and vice versa.
    if      (from_rank>mpiPM->get_myrank() && (i==XN || i==YN || i==ZN)) dir = static_cast<direction>(i);
    else if (from_rank<mpiPM->get_myrank() && (i==XP || i==YP || i==ZP)) dir = static_cast<direction>(i);
  }
  else if (recv_tag == BC_MPItag) {
    // Internal boundary, so in neg.dir. from_rank should be < myrank.
    if      (from_rank<mpiPM->get_myrank() && (i==XN || i==YN || i==ZN)) dir = static_cast<direction>(i);
    else if (from_rank>mpiPM->get_myrank() && (i==XP || i==YP || i==ZP)) dir = static_cast<direction>(i);
  }
  else rep.error("bad tag",recv_tag);
      }
    }
  }
  if (dir==NO) rep.error("Message is not from a neighbour!",from_rank);
  //  cout <<mpiPM->get_myrank()<<"\tcomm_receive_any_data: Receiving Data type ";
  //  cout <<recv_tag<<" from rank: "<<from_rank<<" from direction "<<dir<<"\n";

  //
  // Now choose the boundary data associated with boundary we are receiving:
  //
  struct boundary_data *recv_b = 0;
  recv_b = &BC_bd[static_cast<int>(dir)];

  //
  // Receive the data:
  //
  err = COMM->receive_cell_data(from_rank,  ///< rank of process we are receiving from.
        &(recv_b->data),     ///< list of cells to get data for. 
        recv_b->data.size(), ///< number of cells in list (extra checking!)
        recv_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
        recv_id      ///< identifier for receive, for any book-keeping that might be needed.
        );
  if (err) rep.error("COMM->receive_cell_data() returned error",err);

  //
  // If on a full timestep update, set P=Ph for received boundary.
  //
  if (cstep==maxstep) {
    struct boundary_data *b2 = &BC_bd[static_cast<int>(dir)];
    list<cell*>::iterator c=b2->data.begin();
    for (c=b2->data.begin(); c!=b2->data.end(); ++c) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    } // all cells.
  } // if full timestep.
  
  //
  // Wait until send is finished, and then return.
  //
  //  cout <<"BC_update_BCMPI: waiting for finish "<<b->dir<<"\n";
  err = COMM->wait_for_send_to_finish(send_id);
  if (err) rep.error("Waiting for send to complete!", err);
  //cout <<"BC_update_BCMPI: send and receive finished.\n";
  //cout <<"BC_update_BCMPI: Sent BC "<<b->dir<<", received BC "<<dir<<"\n";
  //cout <<"*******************************************\n\n";

  return 0;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::comm_select_data2send(list<cell *> *l,
                 int *nc,
                 enum direction d)
{
  // Check inputs.
  if ( !(*l).empty() ) {
    rep.warning("comm_select_data2send: List not empty! Emptying it now.",0,(*l).size());
    (*l).clear(); if(!(*l).empty()) rep.error("emptying list.",(*l).empty());
  }
  if ( *nc !=0 ) {
    rep.warning("comm_select_data2send: uninitialized counter in input. setting it to zero.",0,*nc);
    *nc=0;
  }
  
  // Calculate perpendicular directions.
  enum direction perp2d=NO, perp3d=NO, oppdir=NO;
  oppdir = OppDir(d);
  switch (d) {
   case XN: case XP:
    perp2d = YP; perp3d = ZP; break;
   case YN: case YP:
    perp2d = XP; perp3d = ZP; break;
   case ZN: case ZP:
    perp2d = XP; perp3d = YP; break;
   default:
    rep.error("comm_select_data2send: Bad direction",d); break;
  }
  
  // First get to boundary from first point on grid [XN,YN,ZN].
  cell *c= FirstPt(); if (!c) rep.error("lost on grid!",c);
  while ( (NextPt(c,d)) && (NextPt(c,d)->isgd)) {c = NextPt(c,d);}
  //  rep.printVec("select_data2send: first point to send",c->x,G_ndim);
  
  // Now go through and add all cells on the boundary to the list.
  // Safe to assume grid is bigger than depth of boundary cells!!
  // (so i don't check that cells exist, etc.)
  if (G_ndim==1) {
    for (int n=0; n<BC_nbc; n++) {
      //      rep.printVec("1D select_data2send: point to send",c->x,G_ndim);
      (*l).push_back(c); (*nc)++;
      c=NextPt(c,oppdir);
    }
  } // if 1D
  
  else if (G_ndim==2) {
    cell *temp=0;
    do {
      // Get column of cells at current point.
      temp=c;
      for (int n=0; n<BC_nbc; n++) {
  (*l).push_back(temp); (*nc)++; // have at least one cell.
  temp = NextPt(temp,oppdir);
      } // get column of depth nbc
      c = NextPt(c,perp2d);
    } while ( (c) && (c->isgd) );
  } // if 2D
  
  else if (G_ndim==3) {
    cell *temp=0; cell *t2=0;
    do {
      // Get column of cells at current point.
      t2 = c;
      for (int n=0; n<BC_nbc; n++) {
  // rep.printVec("select_data2send: point to send",t2->x,G_ndim);
  (*l).push_back(t2); (*nc)++; // have at least one cell.
  t2 = NextPt(t2,oppdir);
      }
      temp = c;
      while ( (temp=NextPt(temp,perp2d)) && (temp->isgd) ) {
  t2 = temp;
  for (int n=0; n<BC_nbc; n++) {
    (*l).push_back(t2); (*nc)++; // have at least one cell.
    t2 = NextPt(t2,oppdir);
  } // get column of depth nbc
      } // loop over 2d dir.
      c = NextPt(c,perp3d);
    } while ( (c) && (c->isgd) );
  } // if 3D
  else rep.error("Bad NDIM",G_ndim);
  return 0;
}



// ##################################################################
// ##################################################################




/***************************************************************************************/
/***************************** RAYTRACING BOUNDARY STUFF *******************************/
/***************************************************************************************/
#ifdef PLLEL_RT


// ##################################################################
// ##################################################################



int UniformGridParallel::Setup_RT_Boundaries(const int src_id)
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::Setup_RT_Boundaries() starting!\n";
#endif 

  // ----------------------- SANITY CHECKS --------------------------
  // processor ids are numbered from zero; if they're not (by MPI implementation)
  // then my code will bug out before it gets to here...
  //if (!RT_recv_list.empty() || !RT_send_list.empty())
  //  rep.error("send/recv lists not empty!!! more than one source???",RT_recv_list.size());

  //
  // Check that send/recv lists for this source have not already been set up!
  // Element i of the list should NOT correspond to source i, unless by chance.
  // Look through the ids in RT_source_list, just to make sure it is not already there.
  //
  int sle=-1; // source list element.
  for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
    if (RT_source_list[v].source_id == src_id) sle = v;
  }
  if (sle>=0)
    rep.error("source is already present in RT_source_list!",sle);

  //
  // Check that we are not past the correct number of sources!
  //
  if (src_id >= SimPM.RS.Nsources) {
      rep.error("Requested RT source which does not exist",
                src_id-SimPM.RS.Nsources);
  }
  // ----------------------- SANITY CHECKS --------------------------

  // ------------------------ SETUP BOUNDARIES ----------------------
  //
  // Now we need to add a new send list and a new recv list for this source.
  //
  struct RT_source_comms_info  this_src_comms;
  this_src_comms.source_id = src_id;
  this_src_comms.RT_recv_list.clear();
  this_src_comms.RT_send_list.clear();
  
  //std::vector<struct RT_boundary_list_element>  send_list_i;
  //RT_send_listarray.push_back(send_list_i);
  //std::vector<struct RT_boundary_list_element>  recv_list_i;
  //RT_recv_listarray.push_back(recv_list_i);
  //
  // Now we call the appropriate function, depending on what type of source it is
  // The only two relevant types of source, as far as boundary communication is
  // concerned, are sources at infinity and those at a finite position, because the
  // boundary communication differs between them.
  //
  int err=0;
  if (!SimPM.RS.sources[src_id].at_infinity) {
    err += setup_RT_finite_ptsrc_BD(this_src_comms.source_id,
                  this_src_comms.RT_recv_list,
                  this_src_comms.RT_send_list
                  );
  }
  else {
    err += setup_RT_infinite_src_BD(this_src_comms.source_id,
                  this_src_comms.RT_recv_list,
                  this_src_comms.RT_send_list
                  );
  }

  //
  // Now we add this_src_comms to the list.
  //
  RT_source_list.push_back(this_src_comms);

  // ------------------------ SETUP BOUNDARIES ----------------------

#ifdef RT_TESTING
  cout <<"UniformGridParallel::Setup_RT_Boundaries() finished!\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::setup_RT_infinite_src_BD(
                const int src_id,
                std::vector<struct RT_boundary_list_element>  &RT_recv_list,
                std::vector<struct RT_boundary_list_element>  &RT_send_list
                )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_diffuse_radn_src_BD() starting.\n";
#endif 
  int err=0;
  //
  // Make sure src is at infinity, and has the correct type.
  //
  if (!SimPM.RS.sources[src_id].at_infinity) {
    rep.error("Source for diffuse radiation is not at infinity",src_id);
  }
  //
  // If the previous function returned at all, then the source exists, so no
  // need to check that.  For a source at infinity there can be only one send
  // and one receive domain (or none, if current domain an edge domain).
  //
  
  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
    rep.error("diffuse radn src recv-list not empty!",
              RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
    rep.error("diffuse radn src send-list not empty!",
              RT_send_list.size());
  }

  //
  // Get source direction:
  //
  enum direction srcdir = RT_src_at_infty_direction(src_id);
  if (srcdir==NO) {
    rep.error("src position is not at infinity",srcdir);
  }
  enum direction recv_dir = srcdir;
  enum direction send_dir = OppDir(srcdir);
  //
  // Recv boundary in the source direction, send boundary in opp.dir.
  // Can have only one send and one recv source, so the lists should
  // get exactly one element.
  //
  int recv_proc = mpiPM->ngbprocs[recv_dir];
  int send_proc = mpiPM->ngbprocs[send_dir];

  //
  // Set up the lone recv boundary list element.  Note we leave an
  // empty vector if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempR;
  tempR.RT_bd=0;
  switch (recv_dir) {
    case XN: tempR.dir = static_cast<int>(dir_XN); break;
    case XP: tempR.dir = static_cast<int>(dir_XP); break;
    case YN: tempR.dir = static_cast<int>(dir_YN); break;
    case YP: tempR.dir = static_cast<int>(dir_YP); break;
    case ZN: tempR.dir = static_cast<int>(dir_ZN); break;
    case ZP: tempR.dir = static_cast<int>(dir_ZP); break;
    default:
    cout <<"\t No processor in receive direction d="<<recv_dir;
    cout <<": proc="<<recv_proc<<".\n";
    tempR.dir=-1;
    break;
  }
  tempR.rank = recv_proc;
  if (tempR.rank >=0) {
#ifdef RT_TESTING
    cout <<"\t\tFound diffuse-recv-proc in dir="<<recv_dir;
    cout <<" (R.dir="<<tempR.dir<<")";
    cout <<", rank="<<tempR.rank<<", setting up boundary data.\n";
#endif 
    RT_recv_list.push_back(tempR);
    err = setup_RT_recv_boundary(RT_recv_list.front());
    if (err) 
      rep.error("failed to set up diffuse radn recv boundary",err);
  }
  else {
#ifdef RT_TESTING
    cout <<"\t\tDiffuse-recv-proc doesn't exist in dir="<<recv_dir;
    cout <<"; I must be at edge of domain.\n";
#endif 
  }

  //
  // Set up the lone send boundary list element.  Note we leave an empty vector
  // if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempS;
  tempS.RT_bd=0;
  switch (send_dir) {
    case XN: tempS.dir = static_cast<int>(dir_XN); break;
    case XP: tempS.dir = static_cast<int>(dir_XP); break;
    case YN: tempS.dir = static_cast<int>(dir_YN); break;
    case YP: tempS.dir = static_cast<int>(dir_YP); break;
    case ZN: tempS.dir = static_cast<int>(dir_ZN); break;
    case ZP: tempS.dir = static_cast<int>(dir_ZP); break;
    default:
#ifdef RT_TESTING
    cout <<"\t No processor in send direction d="<<send_dir<<": proc="<<send_proc<<".\n";
#endif 
    tempS.dir=-1;
    break;
  }
  tempS.rank = send_proc;
  if (tempS.rank >=0) {
#ifdef RT_TESTING
    cout <<"\t\tFound diffuse-send-proc in dir="<<send_dir;
    cout <<", rank="<<tempS.rank<<", setting up boundary data.\n";
#endif 
    RT_send_list.push_back(tempS);
    err = setup_RT_send_boundary(RT_send_list.front());
    if (err) rep.error("failed to set up diffuse radn send boundary",err);
  }
  else {
#ifdef RT_TESTING
    cout <<"\t\tDiffuse-send-proc doesn't exist in dir="<<send_dir;
    cout <<"; I must be at edge of domain.\n";
#endif 
  }

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_diffuse_radn_src_BD() returning.\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



enum direction UniformGridParallel::RT_src_at_infty_direction(const int src_id)
{
  //
  // Get source position vector; compare values to find the unreasonably
  // value which identifies the direction of the infinite source.
  //
  double srcpos[MAX_DIM];
  for (int v=0;v<G_ndim;v++) {
    srcpos[v] = SimPM.RS.sources[src_id].pos[v];
  }
  enum direction srcdir=NO;
  for (int v=0;v<G_ndim;v++) {
    if (srcpos[v] < -1.e99) srcdir = static_cast<direction>(2*v);
    if (srcpos[v] >  1.e99) srcdir = static_cast<direction>(2*v+1);
  }
  if (srcdir==NO) {
    cout <<"WARNING: RT source is not at infinity! RT_src_at_infty_direction()\n";
    rep.printVec("Srcpos",srcpos,G_ndim);
  }
  return srcdir;
}



// ##################################################################
// ##################################################################




int UniformGridParallel::setup_RT_finite_ptsrc_BD(
                              const int src_id,
                              std::vector<struct RT_boundary_list_element>  &RT_recv_list,
                              std::vector<struct RT_boundary_list_element>  &RT_send_list
                              )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_monochromatic_ptsrc_BD() starting.\n";
#endif
  //
  // NOTE THIS FUNCTION IS OLD CODE, AND PROBABLY NOT THE MOST EFFICIENT WAY OF 
  // DOING THIS!  IN PARTICULAR, IT IS DEFINITELY NOT MODULAR IN ANY WAY, AND 
  // IS VERY COMPLICATED, AND HAD LOTS OF BUGS -- TOOK ME A MONTH AT LEAST TO 
  // FIX ALL OF THEM.  SO I DIDN'T WANT TO TOUCH IT WHEN UPGRADING -- I JUST 
  // WRAPPED IT IN THIS FUNCTION.
  //

  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
          rep.error("Monochromatic point src recv-list not empty!",RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
          rep.error("Monochromatic point src send-list not empty!",RT_send_list.size());
  }

  //
  // First we find the source:
  //
#ifdef CELL_CENTRED_SRC
  //
  // Find Source.
  //
  enum direction srcdir[G_ndim];
  double srcpos[G_ndim];
  for (int v=0;v<G_ndim;v++) {
    srcpos[v] = SimPM.RS.sources[src_id].pos[v];
  }

  for (int i=0;i<G_ndim;i++) {
    if      (srcpos[i]<=G_xmin[i]) srcdir[i] = static_cast<direction>(2*i);
    else if (srcpos[i]>=G_xmax[i]) srcdir[i] = static_cast<direction>(2*i+1);
    else                           srcdir[i] = static_cast<direction>(-1);
  }
#endif // CELL_CENTRED_SRC

#ifdef NON_CELL_CENTRED_SRC
  //
  // Source should be at a cell corner, so we need to be more careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processors!), so
  // this should work fine.
  //
  enum direction srcdir[G_ndim];
  double srcpos[MAX_DIM];
  for (int v=0;v<G_ndim;v++) {
    srcpos[v] = SimPM.RS.sources[src_id].pos[v];
  }


  int i_srcpos[G_ndim];
  CI.get_ipos_vec(srcpos,i_srcpos);

  for (int i=0;i<G_ndim;i++) {
    if      (i_srcpos[i]<G_ixmin[i]) {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin="<<G_xmin[i]<<" src is off grid in neg.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i);
    }
    else if (i_srcpos[i]>G_ixmax[i])  {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmax="<<G_xmax[i]<<" src is off grid in pos.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i+1);
    }
    else {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin="<<G_xmin[i]<<" xmax="<<G_xmax[i];
      cout <<" src is on local domain.\n";
#endif
      srcdir[i] = static_cast<direction>(-1);
    }
  }
#endif // NON_CELL_CENTRED_SRC

  //
  //
  // Choose processors I need to receive data from and send data to.
  //
  //
  int nx[SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++)
    nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/mpiPM->LocalRange[i]);

  //
  // First see if direction of source off grid has a neigbour processor and
  // if the opposite direction has a neighbour processor.
  //
  // recv procs:
  bool recv_proc_exists[G_ndim];
  for (int i=0;i<G_ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if      ((srcdir[i]==negdir) && (!pconst.equalD(G_xmin[i],SimPM.Xmin[i])))
      recv_proc_exists[i]=true;
    else if ((srcdir[i]==posdir) && (!pconst.equalD(G_xmax[i],SimPM.Xmax[i])))
      recv_proc_exists[i]=true;
    else
      recv_proc_exists[i]=false;
  }

  //
  // send procs:
  //
#ifdef CELL_CENTRED_SRC
  bool send_proc_exists[2*G_ndim];
  for (int i=0;i<G_ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if (srcpos[i]>=G_xmin[i] && (!pconst.equalD(mpiPM->LocalXmin[i],SimPM.Xmin[i])) )
      send_proc_exists[negdir] = true;
    else 
      send_proc_exists[negdir] = false; // either doesn't exist, or we don't need it.

    if (srcpos[i]<=G_xmax[i] && (!pconst.equalD(mpiPM->LocalXmax[i],SimPM.Xmax[i])) )
      send_proc_exists[posdir] = true;
    else 
      send_proc_exists[posdir] = false; // either doesn't exist, or we don't need it.
  }
#endif // CELL_CENTRED_SRC
#ifdef NON_CELL_CENTRED_SRC
  //
  // Source should be at a cell corner, so we need to be more careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processors!), so
  // this should work fine.
  //
  // So we only send data to a neighbour if no part of it is in a
  // source plane.
  //
  bool send_proc_exists[2*G_ndim];
  for (int i=0;i<G_ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if ( (i_srcpos[i]>G_ixmin[i]) && (!pconst.equalD(G_xmin[i],SimPM.Xmin[i])) )
      send_proc_exists[negdir] = true;
    else 
      send_proc_exists[negdir] = false; // either doesn't exist, or we don't need it.

    if ( (i_srcpos[i]<G_ixmax[i]) && (!pconst.equalD(G_xmax[i],SimPM.Xmax[i])) )
      send_proc_exists[posdir] = true;
    else 
      send_proc_exists[posdir] = false; // either doesn't exist, or we don't need it.
  }
#endif // NON_CELL_CENTRED_SRC


  //
  // Now get id's and directions of receive processors.
  //
  struct RT_boundary_list_element temp;  temp.RT_bd=0;
  // x-dir
  if      (srcdir[XX]==XN && recv_proc_exists[XX]==true) {
    temp.rank=mpiPM->get_myrank()-1; temp.dir = dir_XN;
    RT_recv_list.push_back(temp);
  }
  else if (srcdir[XX]==XP && recv_proc_exists[XX]==true) {
    temp.rank=mpiPM->get_myrank()+1; temp.dir = dir_XP;
    RT_recv_list.push_back(temp);
  }

  // y-dir
  if (G_ndim>1) {
    if (srcdir[YY]==YN && recv_proc_exists[YY]==true) {
      temp.rank=mpiPM->get_myrank()-nx[XX]; temp.dir = dir_YN;
      RT_recv_list.push_back(temp);
      if      (srcdir[XX]==XN && recv_proc_exists[XX]==true) {
  temp.rank=mpiPM->get_myrank()-nx[XX]-1; temp.dir += dir_XN; //srcdir[YY]+srcdir[XX];
  RT_recv_list.push_back(temp);
      }
      else if (srcdir[XX]==XP && recv_proc_exists[XX]==true) {
  temp.rank=mpiPM->get_myrank()-nx[XX]+1; temp.dir += dir_XP; //srcdir[YY]+srcdir[XX];
  RT_recv_list.push_back(temp);
      }
    }
    else if (srcdir[YY]==YP && recv_proc_exists[YY]==true) {
      temp.rank=mpiPM->get_myrank()+nx[XX]; temp.dir = dir_YP;
      RT_recv_list.push_back(temp);
      if      (srcdir[XX]==XN && recv_proc_exists[XX]==true) {
  temp.rank=mpiPM->get_myrank()+nx[XX]-1; temp.dir += dir_XN; //srcdir[YY]+srcdir[XX];
  RT_recv_list.push_back(temp);
      }
      else if (srcdir[XX]==XP && recv_proc_exists[XX]==true) {
  temp.rank=mpiPM->get_myrank()+nx[XX]+1; temp.dir += dir_XP; //srcdir[YY]+srcdir[XX];
  RT_recv_list.push_back(temp);
      }
    }
  } // y-dir

  // z-dir
  if (G_ndim>2 && recv_proc_exists[ZZ]==true && srcdir[ZZ]!=NO) {
    int rank; int dir;
    if       (srcdir[ZZ]==ZN) {
      rank = mpiPM->get_myrank()-nx[XX]*nx[YY]; dir = dir_ZN;
    }
    else if (srcdir[ZZ]==ZP) {
      rank = mpiPM->get_myrank()+nx[XX]*nx[YY]; dir = dir_ZP;
    }
    else {
      rank=999; dir=999;
      rep.error("Bad Logic -- z-dir RT setup",srcdir[ZZ]);
    }
    temp.rank=rank; temp.dir=dir; temp.RT_bd=0;

    // Have up to 3 extra boundaries if there is a Z neighbour.
    RT_recv_list.push_back(temp);
    if (recv_proc_exists[YY] && srcdir[YY]!=NO) {
      if      (srcdir[YY]==YN) {temp.rank = rank-nx[XX]; temp.dir=dir +dir_YN;} //srcdir[ZZ]+srcdir[YY];}
      else if (srcdir[YY]==YP) {temp.rank = rank+nx[XX]; temp.dir=dir +dir_YP;} //srcdir[ZZ]+srcdir[YY];}
      else rep.error("Bad Logic -- yz-dir RT setup",srcdir[ZZ]);
      RT_recv_list.push_back(temp);
    }
    if (recv_proc_exists[XX] && srcdir[XX]!=NO) {
      if      (srcdir[XX]==XN) {temp.rank = rank-1; temp.dir=dir +dir_XN;}//srcdir[ZZ]+srcdir[XX];}
      else if (srcdir[XX]==XP) {temp.rank = rank+1; temp.dir=dir +dir_XP;}//srcdir[ZZ]+srcdir[XX];}
      else rep.error("Bad Logic -- xz-dir RT setup",srcdir[ZZ]);
      RT_recv_list.push_back(temp);
    }
    if (recv_proc_exists[YY] && srcdir[YY]!=NO &&
  recv_proc_exists[XX] && srcdir[XX]!=NO) {
      // dir=srcdir[ZZ]; rank=mpiPM->get_myrank()(+-)nx[XX]*nx[YY];
      //dir = srcdir[ZZ]+srcdir[YY]+srcdir[XX];
      if      (srcdir[YY]==YN) {rank -= nx[XX]; dir += dir_YN;}
      else if (srcdir[YY]==YP) {rank += nx[XX]; dir += dir_YP;}
      if      (srcdir[XX]==XN) {rank -= 1;      dir += dir_XN;}
      else if (srcdir[XX]==XP) {rank += 1;      dir += dir_XP;}
      temp.rank=rank; temp.dir=dir;
      RT_recv_list.push_back(temp);
    }
  } //z-dir
  // Up to Seven processors to receive from, or as few as None.

  //
  // Choose processors I need to send data to.
  //
  // x-dir
  int rank=mpiPM->get_myrank(); int dir;
  if (send_proc_exists[XN]) {
    temp.rank = rank-1; temp.dir = dir_XN; RT_send_list.push_back(temp);
  }
  if (send_proc_exists[XP]) {
    temp.rank = rank+1; temp.dir = dir_XP; RT_send_list.push_back(temp);
  }

  // y-dir
  if (G_ndim>1) {
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[YN]) {
      rank -= nx[XX]; dir=dir_YN;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
      if (send_proc_exists[XN]) {
  temp.rank = rank-1; temp.dir=dir +dir_XN; RT_send_list.push_back(temp);
      }
      if (send_proc_exists[XP]) {
  temp.rank = rank+1; temp.dir=dir +dir_XP; RT_send_list.push_back(temp);
      }
    }
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[YP]) {
      rank += nx[XX]; dir=dir_YP;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
      if (send_proc_exists[XN]) {
  temp.rank = rank-1; temp.dir=dir +dir_XN; RT_send_list.push_back(temp);
      }
      if (send_proc_exists[XP]) {
  temp.rank = rank+1; temp.dir=dir +dir_XP; RT_send_list.push_back(temp);
      }
    }
  }
  // z-dir
  if (G_ndim>2) {
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[ZN]) {
      rank -= nx[XX]*nx[YY]; dir = dir_ZN;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
      if (send_proc_exists[XN])
  {temp.rank = rank-1; temp.dir=dir +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[XP])
  {temp.rank = rank+1; temp.dir=dir +dir_XP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN])
  {temp.rank = rank-nx[XX]; temp.dir=dir +dir_YN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP])
  {temp.rank = rank+nx[XX]; temp.dir=dir +dir_YP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN] && send_proc_exists[XN])
  {temp.rank = rank-nx[XX]-1; temp.dir=dir +dir_YN +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN] && send_proc_exists[XP])
  {temp.rank = rank-nx[XX]+1; temp.dir=dir +dir_YN +dir_XP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP] && send_proc_exists[XN])
  {temp.rank = rank+nx[XX]-1; temp.dir=dir +dir_YP +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP] && send_proc_exists[XP])
  {temp.rank = rank+nx[XX]+1; temp.dir=dir +dir_YP +dir_XP; RT_send_list.push_back(temp);}
    }
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[ZP]) {
      rank += nx[XX]*nx[YY]; dir = dir_ZP;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
      if (send_proc_exists[XN])
  {temp.rank = rank-1; temp.dir=dir +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[XP])
  {temp.rank = rank+1; temp.dir=dir +dir_XP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN])
  {temp.rank = rank-nx[XX]; temp.dir=dir +dir_YN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP])
  {temp.rank = rank+nx[XX]; temp.dir=dir +dir_YP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN] && send_proc_exists[XN])
  {temp.rank = rank-nx[XX]-1; temp.dir=dir +dir_YN +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YN] && send_proc_exists[XP])
  {temp.rank = rank-nx[XX]+1; temp.dir=dir +dir_YN +dir_XP; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP] && send_proc_exists[XN])
  {temp.rank = rank+nx[XX]-1; temp.dir=dir +dir_YP +dir_XN; RT_send_list.push_back(temp);}
      if (send_proc_exists[YP] && send_proc_exists[XP])
  {temp.rank = rank+nx[XX]+1; temp.dir=dir +dir_YP +dir_XP; RT_send_list.push_back(temp);}
    }
  } // z-dir
  // Up to 26 neighbours in total for sending (if source is on this processor's domain).  

  //
  // Figure out how many RT boundaries I need and set them up.
  //
#ifdef RT_TESTING
  cout <<"send list: "<<RT_send_list.size()<<" boundaries.\n";
  cout <<"recv list: "<<RT_recv_list.size()<<" boundaries.\n";
  cout <<"SEND LIST DIRECTIONS: ";
  for (unsigned int i=0;i<RT_send_list.size();i++) cout <<RT_send_list[i].dir<<" ";
  cout <<"\n";
  cout <<"RECV LIST DIRECTIONS: ";
  for (unsigned int i=0;i<RT_recv_list.size();i++) cout <<RT_recv_list[i].dir<<" ";
  cout <<"\n";
#endif

  //RT_nbd = static_cast<int>(RT_send_list.size()+RT_recv_list.size());
  //if (RT_bd) rep.error("RT_bd already initialised",RT_bd);
  //try {
  //  RT_bd = new struct boundary_data [RT_nbd];
  //}
  //catch (std::bad_alloc) {
  //  rep.error("out of memory: RT_bd (UniformGridParallel::Setup_RT_Boundaries)",RT_bd);
  //}
  

  //
  // Set up receive boundaries.  First initialise all boundary_data pointers
  // to zero, then setup faces, edges, and corners, in that order (so ngb 
  // pointers are set up to get to corners.
  //
  int err=0;
  for (unsigned int i=0;i<RT_recv_list.size();i++)
    RT_recv_list[i].RT_bd=0;
  //
  // First set up face boundaries.
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    int d=RT_recv_list[i].dir;
    if (d==dir_XN || d==dir_XP || d==dir_YN || d==dir_YP ||
  d==dir_ZN || d==dir_ZP) {
      err += setup_RT_recv_boundary(RT_recv_list[i]);
    }
  }
  //
  // Now set up edge boundaries.
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    int d=RT_recv_list[i].dir;
    if (d==dir_YNXN || d==dir_YNXP || d==dir_YPXN || d==dir_YPXP ||
  d==dir_ZNXN || d==dir_ZNXP || d==dir_ZPXN || d==dir_ZPXP ||
  d==dir_ZNYN || d==dir_ZNYP || d==dir_ZPYN || d==dir_ZPYP) {
      err += setup_RT_recv_boundary(RT_recv_list[i]);
    }
  }
  //
  // Now set up corner boundaries
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    int d=RT_recv_list[i].dir;
    if (d==dir_ZNYNXN || d==dir_ZNYNXP || d==dir_ZNYPXN || d==dir_ZNYPXP ||
  d==dir_ZPYNXN || d==dir_ZPYNXP || d==dir_ZPYPXN || d==dir_ZPYPXP) {
      err += setup_RT_recv_boundary(RT_recv_list[i]);
    }
  }
  //
  // Make sure we initialised all the boundaries correctly!
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    if (!(RT_recv_list[i].RT_bd)) rep.error("Missed out on receive boundary!",RT_recv_list[i].dir);
  }
  if (err) rep.error("failed to setup RT recv boundaries",err);

  for (unsigned int i=0;i<RT_send_list.size();i++) {
    RT_send_list[i].RT_bd=0;
    err += setup_RT_send_boundary(RT_send_list[i]);
  }
  for (unsigned int i=0;i<RT_send_list.size();i++) {
    if (!(RT_send_list[i].RT_bd)) rep.error("Missed out on send boundary!",i);
  }
  if (err) rep.error("failed to setup RT send boundaries",err);

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_monochromatic_ptsrc_BD() finished!\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################





class cell * UniformGridParallel::get2startcell(const enum rt_dirs d ///< direction of boundary.
            )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::get2startcell starting.\n";
#endif 
  class cell *c = FirstPt();
  if (!c) rep.error("Called get2startcell before grid initialised!",c);

  if (d%3==2) {
    // XP dir
    do {c=NextPt(c,XP);} while (NextPt(c,XP)!=0 && NextPt(c,XP)->isgd);
  }
  if ((d/3)%3==2) {
    // YP dir
    do {c=NextPt(c,YP);} while (NextPt(c,YP)!=0 && NextPt(c,YP)->isgd);
  }
  if ((d/9)%3==2) {
    // ZP dir
    do {c=NextPt(c,ZP);} while (NextPt(c,ZP)!=0 && NextPt(c,ZP)->isgd);
  }
#ifdef RT_TESTING
  cout <<"UniformGridParallel::get2startcell returning.\n";
#endif 
  return c;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::setup_RT_recv_boundary(
        struct RT_boundary_list_element &b ///< boundary info.
        )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_recv_boundary() starting (dir="<<b.dir<<").\n";
#endif 
  int err=0;
  cell *c = get2startcell(static_cast<enum rt_dirs>(b.dir));

  //
  // Go through the cases one by one...
  //
  enum direction d1,d2,d3;
  int d;
  cell *t=0;
  //
  // First set up boundary data, since it starts as an uninitialised pointer.
  //
  if (b.RT_bd) rep.error("Already set up RT boudary data!",b.RT_bd);
  b.RT_bd = mem.myalloc(b.RT_bd,1);

  switch (b.dir) {
    //
    // First do the faces, which already have boundaries set up, so we only need
    // to connect the boundary cells to each other laterally.
    //
  case dir_XN:
    d1=XN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_XP:
    d1=XP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_YN:
    d1=YN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_YP:
    d1=YP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_ZN:
    d1=ZN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_ZP:
    d1=ZP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_connect_face_cells(b.RT_bd, &(BC_bd[d1]), d1);
    break;

    //
    // now more difficult ones, where we have to create new cells.
    //
  case dir_ZNYNXN:
  case dir_ZNYNXP:
  case dir_ZNYPXN:
  case dir_ZNYPXP:
  case dir_ZPYNXN:
  case dir_ZPYNXP:
  case dir_ZPYPXN:
  case dir_ZPYPXP:
    //
    // in this case, should have created edge boundaries first, so i can navigate across them.
    //
    if (G_ndim!=3) rep.error("Need 3d for corner cells!!! logic error",b.dir);
    d=(b.dir %3);
    switch (d) {
      //case 0: d1=NO; break; // there must be a direction!
    case 1: d1=XN; break;
    case 2: d1=XP; break;
    default: d1=NO; rep.error("logic!!!",d);
    }
    d=((b.dir/3) %3);
    switch (d) {
      //case 0: d2=NO; break; // there must be a direction!
    case 1: d2=YN; break;
    case 2: d2=YP; break;
    default: d2=NO; rep.error("logic!!!",d);
    }
    d=((b.dir/9) %3);
    switch (d) {
      //case 0: d3=NO; break; // there must be a direction!
    case 1: d3=ZN; break;
    case 2: d3=ZP; break;
    default: d3=NO; rep.error("logic!!!",d);
    }
    //
    // Now create a new cell for the corner, and add it to the list.
    //
#ifdef RT_TESTING
    cout <<"\t\t\tAdding corner cell!!\n";
#endif 
    t = RT_new_corner_cell(c,d1,d2,d3);
#ifdef RT_TESTING
    cout <<"Recv grid cell at corner: "; PrintCell(c);
    cout <<"Recv Corner Cell: "; PrintCell(t);
#endif 
    b.RT_bd->data.push_back(t);
    break;

  case dir_YNXN:
  case dir_YNXP:
  case dir_YPXN:
  case dir_YPXP:
    if (G_ndim<2) rep.error("Need 2d/3d for edge cells!!! logic error",b.dir);
    d1=static_cast<direction>((b.dir %3)-1);     // [0,1] for XN,XP
    d2=static_cast<direction>(((b.dir/3) %3)+1); // [2,3] for YN,YP
#ifdef RT_TESTING
    cout <<"dir="<<b.dir<<" and d1="<<d1<<" and d2="<<d2<<"\n";
#endif 
    d3=ZP;
    //
    // Now pass boundary data to setup function, which will create cells and
    // add them to the boundary data list.
    //
    err += RT_edge_cells_setup(c, b.RT_bd, d1, d2, d3);
    break;

  case dir_ZNXN:
  case dir_ZNXP:
  case dir_ZPXN:
  case dir_ZPXP:
    if (G_ndim<2) rep.error("Need 2d/3d for edge cells!!! logic error",b.dir);
    d1=static_cast<direction>((b.dir %3)-1);     // [0,1] for XN,XP
    d2=static_cast<direction>(((b.dir/9) %3)+3); // [4,5] for ZN,ZP
#ifdef RT_TESTING
    cout <<"dir="<<b.dir<<" and d1="<<d1<<" and d2="<<d2<<"\n";
#endif 
    d3=YP;
    //
    // Now pass boundary data to setup function, which will create cells and
    // add them to the boundary data list.
    //
    err += RT_edge_cells_setup(c, b.RT_bd, d1, d2, d3);
    break;

  case dir_ZNYN:
  case dir_ZNYP:
  case dir_ZPYN:
  case dir_ZPYP:
    if (G_ndim<2) rep.error("Need 2d/3d for edge cells!!! logic error",b.dir);
    d1=static_cast<direction>(((b.dir/3) %3)+1); // [2,3] for YN,YP
    d2=static_cast<direction>(((b.dir/9) %3)+3); // [4,5] for ZN,ZP
#ifdef RT_TESTING
    cout <<"dir="<<b.dir<<" and d1="<<d1<<" and d2="<<d2<<"\n";
#endif 
    d3=XP;
    //
    // Now pass boundary data to setup function, which will create cells and
    // add them to the boundary data list.
    //
    err += RT_edge_cells_setup(c, b.RT_bd, d1, d2, d3);
    break;



  default:
    rep.error("bad dir!",b.dir);
  }

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_recv_boundary() returning (dir="<<b.dir<<").\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################





int UniformGridParallel::RT_edge_cells_setup(cell *c, ///< grid corner cell.
               struct boundary_data *b, ///< pointer to boundary data.
               const enum direction d1, ///< 1st dir
               const enum direction d2, ///< 2nd dir
               const enum direction d3 ///< direction to trace along.
               )
{
  //
  // THIS CODE ASSUMES THAT EDGE BOUNDARY DATA MAY NOT EXIST (AS OPPOSED TO FACE
  // BOUNDARY DATA WHICH MUST EXIST).  EDGES ARE AT DOMAIN CORNERS IN 2D, AND IN
  // 3D ARE THE EXTENSION OF THESE CORNERS ALONG THE 3RD DIMENSION.
  // THE FUNCTION RT_new_edge_cell() RETURNS A POINTER TO THE CELL IF IT HAS ALREADY
  // BEEN CREATED, OR CREATES A NEW CELL, CONNECTS IT TO THE GRID, AND RETURNS ITS
  // POINTER IF IT DOESN'T PRE-EXIST.
  //
#ifdef RT_TESTING
  cout <<"RT_edge_cells_setup starting\n";
#endif 
  if (!b) rep.error("Empty RT Boundary data passed to RT_edge_cells_setup()",b);
  if (G_ndim==1) rep.error("need 2d/3d for edge cells setup",d1);

  cell *t=0, *t2=0;
  unsigned int count=0;
  if (G_ndim==2) {
#ifdef RT_TESTING
    cout <<"calling RT_new_edge_cell(c,d1,d2)\n";
#endif 
    t = RT_new_edge_cell(c,d1,d2);
#ifdef RT_TESTING
    cout <<"Back from RT_new_edge_cell(c,d1,d2)\n";
#endif 
    if (b->data.empty()) {
#ifdef RT_TESTING
      cout <<"data is empty, adding cell to it.\n";
#endif 
    }
    b->data.push_back(t);
#ifdef RT_TESTING
    cout <<"Added cell to boundary data.\n";
#endif 
    count++;
  }
  else {
    //
    // set t2 to be the previous cell in the d3 column (if it exists).
    //
    t2=NextPt(c,OppDir(d3));
    if (!t2) rep.error("no face boundary (or got lost), RT_edge_cells_setup",t2);
    t2=NextPt(t2,d1); if (t2) t2=NextPt(t2,d2); // may or may not exist.
    //
    // now trace through the d3 column, setting the d3 ngb pointers as we go.
    // RT_new_edge_cell() sets the d1,d2 pointers.
    //
    do {
      t = RT_new_edge_cell(c,d1,d2);
      if (t2) t2->ngb[d3]= t;
      t->ngb[OppDir(d3)] = t2;
      t2=t;
      b->data.push_back(t);
      count++;
    } while ( (c=NextPt(c,d3))!=0 && c->isgd);
  }

  if (count != b->data.size()) {
    cerr <<"UniformGridParallel::RT_edge_cells_setup() wrong number of cells!\n";
    return (count-b->data.size());
  }
  else {
#ifdef RT_TESTING
    cout <<"RT_edge_cells_setup returning.\n";
#endif 
    return 0;
  }
}




// ##################################################################
// ##################################################################




class cell * UniformGridParallel::RT_new_edge_cell(const cell *c, ///< grid edge cell.
               const enum direction d1, ///< first dir
               const enum direction d2  ///< second dir
               )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::RT_new_edge_cell() starting.\n";
#endif 
  cell *t=0;
  cell *move=0;
  move=NextPt(c,d1); if (!move || move->isgd) rep.error("lost on grid: edge",d1);
  //
  // So now we are on the boundary and the new cell needs to be in dir d2.  If
  // we have multiple sources, however, it may already have been created when
  // setting up boundary data for a previous source.  So if it exists we just 
  // return a pointer to the cell.  If not, then we create a new cell and set
  // up its connectivity to the rest of the grid, and then return it.
  //
  if (NextPt(move,d2)!=0 && !NextPt(move,d2)->isgd) {
    t = NextPt(move,d2);
  }
  else if (NextPt(move,d2)!=0) {
    rep.error("edge cell already exists! RT_new_edge_cell",NextPt(move,d2));
  }
  else {
    t = newCell();
    t->isbd = true; 
    t->isgd = false;
    t->isedge = -1;
    // Boundary ids are irrelevant -10 means not set, I set to something meaningful later.
    t->id = -98;
    double tau[MAX_TAU];
    for (short unsigned int v=0; v<MAX_TAU; v++) tau[v]=0.0;
    for (int v=0;v<SimPM.RS.Nsources;v++) {
      CI.set_col(t,v,tau);
    } // no harm initialising it!!!

    //
    // Assign position of cell, based on positions of neighbours
    //
    int ix[MAX_DIM], delta=2;
    CI.get_ipos(move,ix);
    enum axes a = static_cast<axes>(d2/2);
    switch (d2) {
    case XN: case YN: case ZN:
      ix[a] -= delta; break;
    case XP: case YP: case ZP:
      ix[a] += delta; break;
    default: rep.error("Bad dir",d2);
    }
    CI.set_pos(t,ix);

    // assign neighbour pointers in the plane d1,d2 -- so two each-way pointers
    enum direction opp;
    move->ngb[d2] = t;
    opp = OppDir(d2); t->ngb[opp] = move;
    move=NextPt(c,d2); if (!move || move->isgd) rep.error("lost on grid: edge2",d2);
    move->ngb[d1] = t;
    opp = OppDir(d1); t->ngb[opp] = move;

#ifdef RT_TESTING
    if (mpiPM->get_myrank()==60 && c->id==31) {
      cout <<"RT_new_edge_cell() grid cell, c, new cell, t.\n";
      cout <<"d1="<<d1<<" and d2="<<d2<<"\n";
      PrintCell(c); PrintCell(t);
    }
#endif
  } // if cell doesn't already exist.

#ifdef RT_TESTING
  cout <<"UniformGridParallel::RT_new_edge_cell() returning.\n";
#endif 
  return t;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::RT_connect_face_cells(
        struct boundary_data *b,        ///< pointer to RT boundary data.
        const struct boundary_data *b2, ///< pointer to BC boundary data.
        const enum direction offdir     ///< face dir
        )
{
  //if (G_ndim==1) return 0; // no connections to make!
  if (!b || !b2)
    rep.error("UniformGridParallel::RT_connect_face_cells() Null b/b2 pointer",b);

  //
  // First set all the directions we need:
  //
  enum direction ondir=OppDir(offdir);
#ifdef RT_TESTING
    cout <<"RT_connect_face_cells() offdir="<< offdir;
    cout <<",  ondir="<<ondir<<"\n";
#endif // RT_TESTING
  enum direction dir1=NO, opp1=NO, dir2=NO, opp2=NO;
  switch (offdir) {
  case XN: case XP:
    dir1=YN; dir2=ZN; opp1=YP; opp2=ZP; break;
  case YN: case YP:
    dir1=XN; dir2=ZN; opp1=XP; opp2=ZP; break;
  case ZN: case ZP:
    dir1=XN; dir2=YN; opp1=XP; opp2=YP; break;
  default:
    rep.error("BAD DIR! UniformGridParallel::RT_connect_face_cells",offdir);
  }

  //
  // Run through all the BC boundary points (may be 2 deep, so skip 2nd ones).
  // Connect all the 1st boundary points to each other.
  // At the same time, add 1st boundary cells to the RT boundary list.
  //
  list<cell*>::const_iterator bpt=b2->data.begin();
  cell *t=0;
  do{
#ifdef RT_TESTING
    cout <<"RT_connect_face_cells() cpos="<< (*bpt)->pos[0]<<"\n";
#endif // RT_TESTING
    if (NextPt(*bpt,ondir)->isgd) {
      //
      // Add to RT list.
      //
      b->data.push_back(*bpt);

      if (G_ndim>1) {
        //
        // Make connections to neighbouring boundary cells.
        //
        t=NextPt(*bpt,ondir);
        t=NextPt(t,dir1); 
        (*bpt)->ngb[dir1] = NextPt(t,offdir);
        t=NextPt(NextPt(t,opp1),opp1);
        (*bpt)->ngb[opp1] = NextPt(t,offdir);
        if (G_ndim==3) {
          t=NextPt(t,dir1);
          t=NextPt(t,dir2); 
          (*bpt)->ngb[dir2] = NextPt(t,offdir);
          t=NextPt(NextPt(t,opp2),opp2);
          (*bpt)->ngb[opp2] = NextPt(t,offdir);
        }
      }
    }
    //
    // Move to next cell.
    //
    ++bpt;
  } while (bpt !=b2->data.end());
  return 0;
}




// ##################################################################
// ##################################################################




class cell * UniformGridParallel::RT_new_corner_cell(const cell *c, ///< grid corner cell.
                 const enum direction d1, ///< x-dir
                 const enum direction d2, ///< y-dir
                 const enum direction d3  ///< z-dir
                 )
{
  //
  // CORNER OF DOMAIN IN 3D -- THE APEX OF THE THREE COORDINATE DIRECTIONS HAS 
  // A BOUNDARY CELL WITH NO FACES CONNECTED TO THE LOCAL DOMAIN, BUT IT IS 
  // STILL NEEDED FOR RAY-TRACING TO GET OPTICAL DEPTH ALONG RAYS ENTERING THE
  // DOMAIN.  THIS FUNCTION SETS UP THE CORNER CELL.
  //
  cell *t=0;
  cell *move=0;
  move=NextPt(c,d1); if (!move || move->isgd) rep.error("lost on grid: corner",d1);
  move=NextPt(move,d2); if (!move || move->isgd) rep.error("lost on grid: corner",d2);

  //
  // If we have multiple sources, the cell may already have been created when
  // setting up boundary data for a previous source.  So if it exists we just 
  // return a pointer to the cell.  If not, then we create a new cell and set
  // up its connectivity to the rest of the grid, and then return it.
  //
  if (NextPt(move,d3)!=0 && !NextPt(move,d3)->isgd) {
    t = NextPt(move,d3);
  }
  else if (NextPt(move,d3)!=0) {
    rep.error("corner cell already exists! RT_new_corner_cell",NextPt(move,d3));
  }
  else {
    t = newCell();
    t->isbd = true; 
    t->isgd = false;
    t->isedge = -1;
#ifdef RT_TESTING
    cout <<"Move cell c+d1="<<d1<<" +d2="<<d2;PrintCell(move);
    t->id = -500;
#else
    t->id = -99;  // Boundary ids are irrelevant
#endif

    //
    // Assign position of cell, based on positions of neighbours
    //
    int ix[MAX_DIM], delta=2;
    CI.get_ipos(move,ix);
    enum axes a = static_cast<axes>(d3/2);
    switch (d3) {
    case ZN: ix[a] -= delta; break;
    case ZP: ix[a] += delta; break;
    default: rep.error("Bad dir",d3);
    }
    CI.set_pos(t,ix);

    //
    // assign neighbour pointers both ways in the 3 directions:
    //
    enum direction opp;
    opp = OppDir(d3); t->ngb[opp]=move;  move->ngb[d3]=t;

    opp = OppDir(d2); move=NextPt(move,opp);
    t->ngb[opp] = NextPt(move,d3);
    if (!t->ngb[opp]) rep.error("setup edges before corners!",opp);
    NextPt(move,d3)->ngb[d2] = t;

    opp = OppDir(d1); move=NextPt(move,d2); move=NextPt(move,opp);  
    t->ngb[opp] = NextPt(move,d3);
    if (!t->ngb[opp]) rep.error("setup edges before corners!",opp);
    NextPt(move,d3)->ngb[d1] = t;
  } // if cell didn't already exist.

  return t;
}



// ##################################################################
// ##################################################################





int UniformGridParallel::setup_RT_send_boundary(struct RT_boundary_list_element &b ///< boundary info.
            )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_send_boundary() starting (dir="<<b.dir<<").\n";
#endif 
  int err=0;
  cell *c = get2startcell(static_cast<enum rt_dirs>(b.dir));
  cell *t = 0;

  //
  // Set the directions associated with the boundary.
  //
  enum direction d1,d2,d3;
  int d;
  d=(b.dir %3);
  switch (d) {
  case 0: d1=NO; break;
  case 1: d1=XN; break;
  case 2: d1=XP; break;
  default: d1=NO; rep.error("logic!!!",d);
  }
  d=((b.dir/3) %3);
  switch (d) {
  case 0: d2=NO; break;
  case 1: d2=YN; break;
  case 2: d2=YP; break;
  default: d2=NO; rep.error("logic!!!",d);
  }
  d=((b.dir/9) %3);
  switch (d) {
  case 0: d3=NO; break;
  case 1: d3=ZN; break;
  case 2: d3=ZP; break;
  default: d3=NO; rep.error("logic!!!",d);
  }

  //
  // Set up boundary data, since I need to add cells to the list.
  //
  if (b.RT_bd) rep.error("Already set up RT boudary data!",b.RT_bd);
  b.RT_bd = mem.myalloc(b.RT_bd,1);
  
  switch (b.dir) {
    //
    // First do the corners, which are just one cell, and I should be there already!
    //
  case dir_ZNYNXN:
  case dir_ZNYNXP:
  case dir_ZNYPXN:
  case dir_ZNYPXP:
  case dir_ZPYNXN:
  case dir_ZPYNXP:
  case dir_ZPYPXN:
  case dir_ZPYPXP:
#ifdef RT_TESTING
    cout <<"Send Corner Cell: "; PrintCell(t);
#endif 
    b.RT_bd->data.push_back(c);
    break;

    //
    // Now do the edges, where I need to add a column of cells.
    // Two of the dirs label the edge, and the other is the direction
    // along the column.
    //
  case dir_YNXN:
  case dir_YNXP:
  case dir_YPXN:
  case dir_YPXP:
  case dir_ZNXN:
  case dir_ZNXP:
  case dir_ZPXN:
  case dir_ZPXP:
  case dir_ZNYN:
  case dir_ZNYP:
  case dir_ZPYN:
  case dir_ZPYP:
    if (G_ndim==1) rep.error("Bad ndim!",G_ndim);
    //
    // if 2D, only add one cell...
    //
    if (G_ndim==2) {
      b.RT_bd->data.push_back(c);
    }
    else {
      //
      // Get normal direction:
      //
      if      (d1==NO) d1=XP;
      else if (d2==NO) d1=YP;
      else if (d3==NO) d1=ZP;
      else rep.error("Need a direction for tracing a column!!!",d1);
      
      //
      // trace along that column, starting from cell c.
      //
      do {
        b.RT_bd->data.push_back(c);
      } while ((c=NextPt(c,d1))->isgd);
    }
    break;

  case dir_XN:
  case dir_XP:
  case dir_YN:
  case dir_YP:
  case dir_ZN:
  case dir_ZP:
    if      (G_ndim==1) {
      b.RT_bd->data.push_back(c);
    }
    else if (G_ndim==2) {
      //
      // Find Perpendicular Direction
      //
      if      (d1==NO) d1=XP;
      else if (d2==NO) d1=YP;
      else rep.error("Need a direction for tracing a column[face]!!!",d1);
      //
      // trace along that column, starting from cell c.
      //
      do {
        b.RT_bd->data.push_back(c);
      } while ((c=NextPt(c,d1))->isgd);
    }
    else {
      //
      // Find Perpendicular Directions
      //
      if      (d1!=NO) d=XP;
      else if (d2!=NO) d=YP;
      else if (d3!=NO) d=ZP;
      else rep.error("Need a direction for tracing a column[face]!!!",d1);
      if      (d==XP) {d1=YP; d2=ZP;}
      else if (d==YP) {d1=XP; d2=ZP;}
      else if (d==ZP) {d1=XP; d2=YP;}
      else rep.error("BAD DIRECTION!!!",d);
      //
      // trace along that plane, starting from cell c.
      //
      do {
        t=c;
        do {
          b.RT_bd->data.push_back(t);
        } while ((t=NextPt(t,d1))->isgd);
      } while ((c=NextPt(c,d2))->isgd);
    }
    break;
    
  default:
    rep.error("bad dir!",b.dir);
  }

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_send_boundary() returning (dir="<<b.dir<<").\n";
#endif 
  return err;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::Receive_RT_Boundaries(const int src_id ///< source id
                 )
{
#ifdef RT_TESTING
  cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": Starting\n";
#endif 
  int err=0;
  
  //
  // We need the recv-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle=-1; // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
   cout <<"\tRecv_RT_Boundaries() src="<<src_id<<": no procs to recv from.\n";
#endif 
   return 0;
  }
  else {
    for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle<0) rep.error("source not found in RT_source_list",sle);
#ifdef RT_TESTING
  else
   cout <<"\tRevcv_RT_Boundaries() src="<<src_id<<": receiving data.\n";
#endif 

//  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
//    rep.error("Can't receive boundary for src_id b/c boundary was never set up",src_id);
//  }
//  if (RT_source_list[src_id].source_id != src_id) {
//    rep.error("source id for boundary indexed by src_id is wrong",src_id);
//  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);

  //
  // See how many boundaries to receive, and create a list of them
  //
  int n_recv = i_src->RT_recv_list.size();
  if (n_recv==0) {
#ifdef RT_TESTING
    cout <<"UniformGridParallel::Receive_RT_Boundaries() nothing to receive... returning.\n";
#endif 
  }
  else {
    int  r_dirs[n_recv]; // List of directions.
    bool r_done[n_recv]; // Will set these to true as I receive them.
#ifdef RT_TESTING
    cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": to recv from: ";
#endif 
    for (int i=0;i<n_recv;i++) {
      r_dirs[i] = i_src->RT_recv_list[i].dir;
      r_done[i] = false;
#ifdef RT_TESTING
      cout <<" "<<r_dirs[i]<<"[rank="<<i_src->RT_recv_list[i].rank<<"]";
#endif 
    }
#ifdef RT_TESTING
    cout <<"\n";
#endif 
    
    //
    // Now receive the n_recv boundaries, in any order:
    //
    double *buf=0;
    struct boundary_data *b;
    int from_rank, comm_tag, ct, i_recv;
    string id;
    for (int i=0;i<n_recv;i++) {
      id.erase();
      //
      // Find a message that has been sent to us.
      //
      err = COMM->look_for_data_to_receive(&from_rank, id, &comm_tag, COMM_DOUBLEDATA);
      if (err) rep.error("Receive_RT_Boundaries() look4data returned abnormally",err);
      if (comm_tag!=BC_RTtag) {
        rep.error("Getting data, but not RT and waiting for RT data!",comm_tag);
      }

      // 
      // Set pointer to boundary we are receiving, and note which one it is in i_recv.
      //
      b=0; bool test=false; i_recv=-1;
      for (int v=0; v<n_recv; v++) {
        if (i_src->RT_recv_list[v].rank==from_rank) {
          b=i_src->RT_recv_list[v].RT_bd;
          i_recv = v;
          if (!test) test=true;
          else rep.error("Got two boundaries with same fromrank!!!",from_rank);
        }
      }
      if (!b) rep.error("No Receive Boundary corresponding to data from this proc!",from_rank);
      if (r_dirs[i_recv] != i_src->RT_recv_list[i_recv].dir) {
        rep.error("direction mismatch",i_recv); // this is paranoid checking!
      }
      if (r_done[i_recv]) rep.error("This boundary already received!!!",i_recv);
      
      //
      // See how many doubles we are expecting, and allocate memory.
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      ct = b->data.size()*SimPM.RS.sources[src_id].NTau;
      if (ct<1) {
        cerr <<"data size = "<< b->data.size();
        cerr <<", NTau = "<<SimPM.RS.sources[src_id].NTau<<"\n";
        rep.error("Empty boundary!",ct);
      }
      buf=0;
      buf = mem.myalloc(buf,ct);
#ifdef RT_TESTING
      cout <<"\tReceive_RT_Boundaries() rad_src="<<src_id<<": getting "<<ct<<" cells from proc. "<<from_rank<<"\n";
#endif
      
      //
      // Receive data into buffer.
      //
      err = COMM->receive_double_data(from_rank, ///< rank of process we are receiving from.
              comm_tag,  ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
              id, ///< identifier for receive, for any book-keeping that might be needed.
              ct, ///< number of doubles to receive
              buf ///< Pointer to array to write to (must be already initialised).
              );
      if (err) rep.error("Receive_RT_Boundaries() getdata failed",err);

      //
      // Loop over cells in boundary and copy received col values into c->col.
      // Hopefully they are ordered the same way in packing as unpacking!!!
      //
      list<cell*>::iterator c=b->data.begin();
      int count=0;
      double tau[MAX_TAU];
      for (short unsigned int v=0; v<MAX_TAU; v++) tau[v]=0.0;

      for (c=b->data.begin(); c!=b->data.end(); ++c) {
        if (count>=ct) rep.error("too many cells!!!",count-ct);
#ifdef RT_TESTING
        if (count<32) {
          cout <<"col = "<<buf[count]<<" for cell "<<count<<": cell-id="<<(*c)->id;
                cout <<", pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (SimPM.ndim>2) cout<<","<<(*c)->pos[ZZ]<<"]"<<"\n";
          else cout <<"]\n";
        }
#endif  
        for (short unsigned int v=0;
              v<SimPM.RS.sources[src_id].NTau; v++) {
          tau[v] = buf[count];
          count++;
        }
        CI.set_col(*c, src_id, tau);
      }
      if (count != ct) rep.error("BIG ERROR!",count-ct);
      
      //
      // Note that this boundary has been received, so we know if we get the same data twice.
      //
      r_done[i_recv] = true;
      buf = mem.myfree(buf);
      if (err) rep.error("Error receiving boundary data from rank: ",from_rank);
    } // Loop over boundaries.
  } // else have boundaries.

#ifdef RT_TESTING
  cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": Finished\n";
#endif 
  return err;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::Send_RT_Boundaries(const int src_id ///< source id
              )
{
  int err=0;
#ifdef RT_TESTING
  cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Starting\n";
#endif 
  //
  // We need the send-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle=-1; // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
   cout <<"\tSend_RT_Boundaries() src="<<src_id<<": no procs to send to.\n";
#endif 
   return 0;
  }
  else {
    for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle<0) rep.error("source not found in RT_source_list",sle);
#ifdef RT_TESTING
  else
   cout <<"\tSend_RT_Boundaries() src="<<src_id<<": sending data.\n";
#endif 

  //  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
  //    rep.error("Can't send boundary for src_id b/c boundary was never set up",src_id);
  //  }
  //  if (RT_source_list[src_id].source_id != src_id) {
  //    rep.error("source id for boundary indexed by src_id is wrong",src_id);
  //  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);


  //
  // See how many boundaries to send:
  // If none, just skip on to the COMM->barrier() call at the end.
  // If we have some, do the sends.
  //
  int n_send = i_src->RT_send_list.size();
  if (n_send==0) {
#ifdef RT_TESTING
    cout <<"UniformGridParallel::Send_RT_Boundaries() nothing to send... returning.\n";
#endif 
  }
  else {
    //
    // Allocate memory for send buffers.
    //
#ifdef RT_TESTING
    cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Allocating "<<n_send<<" data buffers.\n";
#endif 
    double *data=0;
    std::string *id=0;
    id = mem.myalloc(id, n_send);
    //
    // Loop over Boundaries, put c->col into data, send it.
    // I'm going to use non-blocking sends, so they all send at once.
    //
    for (int i=0; i<n_send; i++) {
      struct boundary_data *b = i_src->RT_send_list[i].RT_bd;
      data=0;

#ifdef RT_TESTING
      cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Sending boundary in dir "<<i_src->RT_send_list[i].dir<<"\n";
#endif 
      //
      // How many cell column densities to send:
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      int nc = b->data.size()*SimPM.RS.sources[src_id].NTau;
      data = mem.myalloc(data, nc);

      //
      // Put columns into data
      //
      list<cell*>::iterator c=b->data.begin();
      int count=0;
      double tau[MAX_TAU];
      for (short unsigned int v=0; v<MAX_TAU; v++) tau[v]=0.0;
      
      for (c=b->data.begin(); c!=b->data.end(); ++c) {
        if (count>=nc) rep.error("too many cells!!!",count-nc);

        CI.get_col(*c, src_id, tau);
        for (short unsigned int v=0;
             v<SimPM.RS.sources[src_id].NTau; v++) {
          data[count] = tau[v];
          count++;
        }
#ifdef RT_TESTING
        if (count<32) {
          cout <<"send data ["<<i<<"]: col[0] = ";
          CI.get_col(*c,src_id,tau);
          cout << tau[0] <<" for cell ";
          cout <<count<<": pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (SimPM.ndim>2)cout<<","<<(*c)->pos[ZZ]<<"]\n";
          else cout <<"]\n";
        }
        //if (count<32)
        //  cout <<"send data ["<<i<<"]: col = "<<(*c)->col<<" for cell "<<count<<": pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY]<<","<<(*c)->pos[ZZ]<<"]\n";
#endif 
      }

      //
      // Send data
      //
#ifdef RT_TESTING
      cout <<"Send_BC["<<i<<"]: Sending "<<nc<<" doubles to proc: "<<i_src->RT_send_list[i].rank;
      cout <<" in direction "<<i_src->RT_send_list[i].dir<<"\n";
#endif
      err += COMM->send_double_data(i_src->RT_send_list[i].rank, ///< rank to send to.
            nc,      ///< size of buffer, in number of doubles.
            data,    ///< pointer to double array.
            id[i],   ///< identifier for send, for tracking delivery later.
            BC_RTtag ///< comm_tag, to say what kind of send this is.
            );
      if (err) rep.error("Send_BC[] send_data failed.",err);
#ifdef RT_TESTING
      cout <<"Send_BC["<<i<<"]: send returned with id="<<id[i]<<"\n";
#endif
    //
    // Delete data array for reuse in next send.
    //
    data = mem.myfree(data);
    } // loop over all sends to send data.

    //
    // Loop over all sends, and wait until they complete!
    //
#ifdef RT_TESTING
    cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Waiting for sends to complete.\n";
#endif 
    for (int i=0; i<n_send; i++) {
      err += COMM->wait_for_send_to_finish(id[i]);
    }

    //
    // delete array of ids
    //
    id = mem.myfree(id);
  } // else (have boundaries to send).

  //
  // Wait for all processors to finish RT for this source, so no messages get mixed up!
  // Not sure if this is strictly necessary... so should try to do without it, once I 
  // know the algorithm is sound.
  //
#ifdef RT_TESTING
  cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Waiting for all procs to finish.\n";
#endif
  COMM->barrier("Send_RT_Boundaries_finished");
  return err;
}


// ##################################################################
// ##################################################################



#endif // PLLEL_RT
//-------------------------------------------------------------
//------------------- CARTESIAN GRID END ----------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID START ------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_cyl_parallel::uniform_grid_cyl_parallel(
      int nd,
      int nv,
      int eqt,
      double *xn,
      double *xp,
      int *nc,
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,xn,xp,nc),
  UniformGridParallel(nd,nv,eqt,xn,xp,nc,p),
  VectorOps_Cyl(nd),
  uniform_grid_cyl(nd,nv,eqt,xn,xp,nc)
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel:: cyl. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim>2)
    rep.error("Need to write code for 3 dimensions",G_ndim);

#ifdef TESTING
  cout <<"cylindrical grid: dR="<<G_dx<<"\n";
#endif 
  set_dx(G_dx);
}



// ##################################################################
// ##################################################################


uniform_grid_cyl_parallel::~uniform_grid_cyl_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel destructor. Present and correct!\n";
#endif 
}



// ##################################################################
// ##################################################################


double uniform_grid_cyl_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  //cout <<" Cell radius: "<< R_com(c)/CI.phys_per_int() +G_ixmin[Rcyl];
  //rep.printVec("  cell centre",c->pos,G_ndim);
  return (R_com(c)-SimPM.Xmin[Rcyl])/CI.phys_per_int() +SIM_ixmin[Rcyl];
}



// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID END --------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- SPHERICAL GRID START --------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_sph_parallel::uniform_grid_sph_parallel(
      int nd,
      int nv,
      int eqt,
      double *xn,
      double *xp,
      int *nc,
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,xn,xp,nc),
  UniformGridParallel(nd,nv,eqt,xn,xp,nc,p),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd),
  uniform_grid_sph(nd,nv,eqt,xn,xp,nc)
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel:: sph. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim!=1)
    rep.error("Need to write code for >1 dimension",G_ndim);

#ifdef TESTING
  cout <<"spherical grid: dr="<<G_dx<<"\n";
#endif 
  set_dx(G_dx);
}


// ##################################################################
// ##################################################################



uniform_grid_sph_parallel::~uniform_grid_sph_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel destructor. Present and correct!\n";
#endif 
}


// ##################################################################
// ##################################################################



double uniform_grid_sph_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  return (R_com(c)-SimPM.Xmin[Rsph])/CI.phys_per_int() +SIM_ixmin[Rsph];
}


// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------
#endif //PARALLEL
