/// \file raytracer_SC_pllel.cc
/// 
/// \author Jonathan Mackey
/// 
/// Definitions for parallel uniform grid.
/// 
/// Modifications:\n
///  - 2010-01-19 JM: call new weighting scheme functions instead of duplicating serial code.
///  - 2010-01-22 JM: worked on sources at cell-corners to try and get it working.
///  - 2010-01-23 JM: made sure source cell geometry is only for cell-centred sources.
///  - 2010-01-26 JM: Debugging sources at cell-vertices in 3D.
///
///  - 2010.07.23 JM: New RSP source position class interface.
///
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
///    Changed source indexing so IDs are set by SimPM.RS and not assigned locally.
///    Changed Add_Source() interface so that all info gets passed to class.
///    Should handle multiple sources now, and sources at infinity.  Of course I 
///    need a "rates-only" update function to actually use multiple emitting sources.
///
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function.  Just a
///    dummy function for now.  (2011.04.15 JM: it points to the old update fn now! The
///    deciding of whether to update or calculate column-density is in ProcessCell()).
///
/// - 2011.04.22 JM: Updated Add_Source() to set Vshell in cells if Column-density
///    update is required.
///
/// - 2011.04.23 JM: Added TauMin[] array of values of TauMin for each source, ordered
///    by source id.  This removes the need for interpolate_2D_RHO() so I got rid of it.
///    Now the interpolate functions read TauMin[src_id] and apply that to the interpolation.
/// - 2012.02.27 JM: Only output timings for source with id=0 (reduces filesize).
/// - 2012.03.31 JM: Added separators between functions.  Updated Add_Source()
///    function to call add_source_to_list() and set_Vshell_for_source() at the
///    appropriate places.
/// - 2013.09.05 JM: Debugged for new get/set col functions.

#include "../global.h"
#include "raytracer_SC.h"
#include <fstream>
#include <iostream>
using namespace std;

#ifdef PARALLEL
#ifdef PLLEL_RT




// ##################################################################
// ##################################################################



raytracer_USC_pllel::raytracer_USC_pllel(class GridBaseClass *ggg,   ///< Pointer to grid
					 class MicroPhysicsBase *mmm ///< Pointer to MicroPhysics Class.
					 )
  : raytracer_USC(ggg,mmm)
{
#ifdef RT_TESTING
  cout <<"SC PARALLEL raytracer class constructor!\n";
#endif
  return;
}

raytracer_USC_pllel::~raytracer_USC_pllel()
{
#ifdef RT_TESTING
  cout <<"SC PARALLEL raytracer class destructor!\n";
#endif
}




// ##################################################################
// ##################################################################




int raytracer_USC_pllel::Add_Source(struct rad_src_info *src ///< source info.
                                   )
{
  cout <<"\n--BEGIN-----raytracer_USC_pllel::AddSource()------------\n";
  //
  // First call serial version.  This finds the source, and centres it
  // on a cell if needed.
  //
#ifdef RT_TESTING
  cout <<"\t**** PARALLEL Add_Source: calling serial version.\n";
#endif
  if (src->at_infinity) {
    raytracer_USC_infinity::add_source_to_list(src);
  }
  else {
    add_source_to_list(src);
  }
  int id=SourceList.back().s->id;
#ifdef RT_TESTING
  cout <<"\t**** PARALLEL Add_Source: serial version returned with id="<<id<<"\n";
#endif
  
  //
  // Now tell the parallel grid to decide which boundaries it needs to 
  // receive data from before processing this source, and which it needs
  // to send data to after processing.
  //
#ifdef RT_TESTING
  cout <<"\t**** PARALLEL Add_Source: Setting up extra RT boundaries on grid.\n";
#endif
  int err = gridptr->Setup_RT_Boundaries(id);
  if (err) rep.error("Failed to setup RT Boundaries",err);

  //
  // Set Vshell for every cell on the grid.
  //
  set_Vshell_for_source(&SourceList.back());

#ifdef RT_TESTING
  cout <<"\t**** PARALLEL Add_Source: all done..\n";
#endif
  cout <<"--END-----raytracer_USC_pllel::AddSource()------------\n";
  return id;
}




// ##################################################################
// ##################################################################




int raytracer_USC_pllel::RayTrace_SingleSource(const int s_id,  ///< Source id
					       const double dt, ///< Timestep
					       const double g   ///< eos gamma.
					       )
{
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
  err += raytracer_USC::RayTrace_SingleSource(s_id, dt, g);
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
	
  if (mpiPM.myrank==0 && (SimPM.timestep%100==0) && s_id==0) {
    cout <<"RT: step:"<<SimPM.timestep<<" Total RT time="<<total;
    cout <<" secs; processing="<<run<<" secs; waiting="<<wait<<"\n";
  }
  return err;
}




// ##################################################################
// ##################################################################



void raytracer_USC_pllel::col2cell_2d(
        const rad_source *src,          ///< source we are working on.
        const cell *c,                  ///< cell to get column to.
        const enum direction entryface, ///< face ray enters cell through.
        const enum direction *perpdir,  ///< array of perp directions towards source (only 1 el in 2D)
        const double *delta,            ///< array of tan(theta) (1 el in 2D) (angle in [0,45]deg)
        double *Nc                      ///< Column densities.
        )
{
  //
  // Get column densities of the two cells closer to the source.
  // If there are no cells, we have zero column to the source.
  // It is assumed that internal boundaries have cells with appropriately updated
  // column data already set up in the boundaries where they need to be.
  //
  double col1[MAX_TAU], col2[MAX_TAU];
  
  cell *c1 = gridptr->NextPt(c,  entryface);
  if (!c1) {
    for (short unsigned int iT=0; iT<src->s->NTau; iT++)
      col1[iT] = col2[iT] = 0.0;
  }
  else {
    cell *c2 = gridptr->NextPt(c1, (*perpdir) );
    if (!c2) {
      CI.get_col(c1, src->s->id, col1);
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col2[iT] = 0.0;
    }
#ifndef NO_SOURCE_CELL_GEOMETRY
#ifdef CELL_CENTRED_SRC
    else if (c2 == src->sc && src->src_on_grid ) {
      // Need to check if c2 is source cell, b/c if it is, the column is wrong by root2.
      CI.get_col(c2, src->s->id, col2);
      for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
        col1[iT]  = 0.0;
        col2[iT] *= sqrt(2.0);
    }
#endif // CELL_CENTRED_SRC
#endif // NO_SOURCE_CELL_GEOMETRY
    else {
      CI.get_col(c1, src->s->id, col1);
      CI.get_col(c2, src->s->id, col2);
    }
  }

  //
  //  INTERPOLATION SCHEMES -- BASICALLY ALL WERE CRAP EXCEPT C2RAY...
  // The column-density variables are normalised differently for each source,
  // now handled by TauMin[s->id].  For C2Ray it is Tau, and for the New udpate
  // it is the mass density integrated along the line of sight (weighted by e.g. 1-x).
  //
  //  0: C2Ray inverse Tau with minTau=0.7: (see Mellema et al.,2006, NewA, 11,374, eq.A.5)
  //
  for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
    Nc[iT] = interpolate_2D(src->s->id, *delta, col1[iT], col2[iT]);
  }
  return;
}



// ##################################################################
// ##################################################################




void raytracer_USC_pllel::col2cell_3d(
        const rad_source *src,            ///< source we are working on
        const cell *c,                  ///< cell to get column to.
        const enum direction entryface, ///< face ray enters cell through.
        const enum direction *perpdir,  ///< array of perp directions towards source (2 els in 3d)
        const double *dx,               ///< array of tan(theta) (angle in [0,45]deg)
        double *Nc                      ///< Column densities.
        )
{
  //
  // Algorithm is the same as that describe in Mellema et al.,2006, NewA, 11,374,
  // appendix A.  Good for 3D cartesian geometry.

  //
  // First get column densities/opacities of the two cells closer to the source.
  // If there are no cells, we have zero column/opacity to the source.
  // It is assumed that internal boundaries have cells with appropriately updated
  // column data already set up in the boundaries where they need to be.
  //
#ifdef RT_TESTING
  //if (mpiPM.myrank==60 && c->id==31) {
  //cout <<"3D ShortChars:: entrydir = "<<entryface<<" and perps = ["<<perpdir[0]<<", "<<perpdir[1]<<"]\n";
  //}
#endif
  cell *c1=0, *c2=0, *c3=0, *c4=0;
  double col1[MAX_TAU], col2[MAX_TAU], col3[MAX_TAU], col4[MAX_TAU];

  c1 = gridptr->NextPt(c,  entryface);
  if (!c1) {
    for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
      col1[iT] = col2[iT] = col3[iT] = col4[iT] = 0.0;
    }
  }
  else {
    CI.get_col(c1,src->s->id, col1);

    c2 = gridptr->NextPt(c1,  perpdir[0]);
    if (!c2) {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col2[iT] = 0.0;
    }
    else {
      CI.get_col(c2,src->s->id, col2);
      
#ifndef NO_SOURCE_CELL_GEOMETRY
#ifdef CELL_CENTRED_SRC
      if (c2==s->sc && s->src_on_grid ) {
        for (short unsigned int iT=0; iT<src->s->NTau; iT++)
          col2[iT] *= sqrt(2.0);
      }
#endif // CELL_CENTRED_SRC
#endif
    }
    
    c3 = gridptr->NextPt(c1,  perpdir[1]);
    if (!c3) {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col3[iT] = 0.0;
    }
    else {
      CI.get_col(c3, src->s->id, col3);

#ifndef NO_SOURCE_CELL_GEOMETRY
#ifdef CELL_CENTRED_SRC
      if (c3==s->sc && s->src_on_grid ) {
        for (short unsigned int iT=0; iT<src->s->NTau; iT++)
          col3[iT] *= sqrt(2.);
      }
#endif // CELL_CENTRED_SRC
#endif
    }
    
    if (c2 && c3) {
      c4 = gridptr->NextPt(c2, perpdir[1]);
      if (!c4) {
	//
	// This can happen if beside a boundary that isn't sent/received with MPI, so 
	// the cells are not connected to each other, but the columns in the boundary
	// are zero in this case, so we can set it to zero here.
	//
#ifdef RT_TESTING
	//cout <<"...corner cell doesn't exist...\n";
	//gridptr->PrintCell(c);
	//gridptr->PrintCell(c1);
	//gridptr->PrintCell(c2);
	//gridptr->PrintCell(c3);
	//  else cout <<"cell c3 doesn't exist either\n.";
	//rep.error("lost on grid -- corner cell doesn't exist",c4);
#endif
        for (short unsigned int iT=0; iT<src->s->NTau; iT++)
          col4[iT] = 0.0;
	//cout <<"...............................\n";
#ifdef RT_TESTING
	//	if (mpiPM.myrank==60) {
	//  cout <<"RT_BD corner cell doesn't exist... cells c,c1,c2,c3\n";
	// gridptr->PrintCell(c);
	//  gridptr->PrintCell(c1);
	//  gridptr->PrintCell(c2);
	//  if (c3) gridptr->PrintCell(c3);
	//  else cout <<"cell c3 doesn't exist either\n.";
	//}
#endif
      }
      else {
        CI.get_col(c4, src->s->id, col4);

#ifndef NO_SOURCE_CELL_GEOMETRY
#ifdef CELL_CENTRED_SRC
	if (c4==s->sc && s->src_on_grid ) {
        for (short unsigned int iT=0; iT<src->s->NTau; iT++)
          col4[iT] *= sqrt(3.0);
        }
#endif // CELL_CENTRED_SRC
#endif // NO_SOURCE_CELL_GEOMETRY
#ifdef RT_TESTING
	//if (c4->id==-500 && mpiPM.myrank==191)
	//  cout <<"RT_BD corner cell found! col="<<col4<<"\n";
#endif
      }
    }
    else {
      for (short unsigned int iT=0; iT<src->s->NTau; iT++)
        col4[iT] = 0.0;
    }
  }
  //  cout <<"3D ShortChars:: col1="<<col1<<" col2="<<col2<<" col3="<<col3<<" col4="<<col4;
  //  cout <<"\t dx = ["<<dx[0]<<", "<<dx[1]<<"]"<<"\n";

#ifdef RT_TESTING
  //cout <<"3D ShortChars:: cols for cell id="<<c->id<<"\n";
  if (!GS.equalD(dx[0],0.0) && !GS.equalD(dx[1],0.0) &&
      (col1[0]<0.0 || col2[0]<0.0 || col3[0]<0.0 || col4[0]<0.0)) {
    cout <<"3D ShortChars:: col1="<<col1[0]<<" col2="<<col2[0]<<" col3="<<col3[0]<<" col4="<<col4[0];
    cout <<"\t dx = ["<<dx[0]<<", "<<dx[1]<<"]"<<"\n";
    gridptr->PrintCell(c);
    gridptr->PrintCell(c1);
    gridptr->PrintCell(c2);
    gridptr->PrintCell(c3);
    gridptr->PrintCell(c4);
    //if (col1<0.0) {cout <<"col1="<<col1<<", setting to zero.\n"; col1=0.0;}
    //if (col2<0.0) {cout <<"col2="<<col2<<", setting to zero.\n"; col2=0.0;}
    //if (col3<0.0) {cout <<"col3="<<col3<<", setting to zero.\n"; col3=0.0;}
    //if (col4<0.0) {cout <<"col4="<<col4<<", setting to zero.\n"; col4=0.0;}
  }
  //new max, cell posn=Vector pos : [1.54412e+17, 2.49265e+18, 2.44853e+18 ]
  //2.20588e+16, 1.30147e+18, 1.03676e+18
  double dd[3] = {1.54412e+17, 2.49265e+18, 2.44853e+18};
  //
  // geometry: This should be to centre--of--volume of cell for
  // non-cartesian geometries!  But it's testing--only code, so it
  // doesn't matter.
  //
  double dpos[3]; CI.get_dpos(c,dpos);
  if (GS.distance(dd,dpos,SimPM.ndim) < gridptr->DX()/2.) {
    cout.setf(ios_base::scientific); cout.precision(15);
    cout <<"cell with max.diff in step 1: id="<<c->id<<"\n";
    cout <<"pnt  pos="; rep.printVec("pnt",dd,SimPM.ndim);
    cout <<"cell pos="; rep.printVec("pos",dpos,SimPM.ndim);
    cout <<"3D ShortChars:: col1="<<col1[0]<<" col2="<<col2[0]<<" col3="<<col3[0]<<" col4="<<col4[0];
    cout <<"\t dx = ["<<dx[0]<<", "<<dx[1]<<"]"<<"\n";
    double w1,w2,w3,w4,mintau3d;
    mintau3d=0.7;
    w1 = (1.-dx[0])*(1.-dx[1])/max(mintau3d,col1[0]);
    w2 =     dx[0] *(1.-dx[1])/max(mintau3d,col2[0]);
    w3 = (1.-dx[0])*    dx[1] /max(mintau3d,col3[0]);
    w4 =     dx[0] *    dx[1] /max(mintau3d,col4[0]);
    mintau3d = (w1+w2+w3+w4);
    w1/=mintau3d; w2/=mintau3d; w3/=mintau3d; w4/=mintau3d;
    cout <<"col2cell = ";
    cout<<w1*col1[0] + w2*col2[0] + w3*col3[0] + w4*col4[0]<<"\n";
    cout.precision(6); cout.setf(ios_base::fixed);
  }
#endif

#ifdef RT_TESTING
  for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
    if (col1[iT]<0.0) {cout <<"col1="<<col1[iT]<<", setting to zero.\n"; col1[iT]=0.0;}
    if (col2[iT]<0.0 && !GS.equalD(dx[0],0.0))
      {cout <<"col2="<<col2[iT]<<", setting to zero.\n"; col2[iT]=0.0;}
    //else if (col2<0.) cout <<"dx[0]="<<dx[0]<<"\n";
    if (col3[iT]<0.0 && !GS.equalD(dx[1],0.0))
      {cout <<"col3="<<col3[iT]<<", setting to zero.\n"; col3[iT]=0.0;}
    //else if (col3<0.) cout <<"dx[1]="<<dx[1]<<"\n";
    if (col4[iT]<0.0 && !GS.equalD(dx[0],0.0)  && !GS.equalD(dx[1],0.0))
      {cout <<"col4="<<col4[iT]<<", setting to zero.\n"; col4[iT]=0.0;}
    //else if (col4<0.) cout <<"dx[0]="<<dx[0]<<"\tdx[1]="<<dx[1]<<"\n";
  }
#endif


  //
  //  0: C2Ray inverse Tau with minTau=0.7: (see Mellema et al.,2006, NewA, 11,374, eq.A.5)
  //
  for (short unsigned int iT=0; iT<src->s->NTau; iT++) {
    Nc[iT] = interpolate_3D(src->s->id, dx[0], dx[1], 
                            col1[iT], col2[iT], col3[iT], col4[iT]);
  }
  return;

}



// ##################################################################
// ##################################################################




#endif // PLLEL_RT
#endif // PARALLEL

