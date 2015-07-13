///
/// \file sim_projection.cc
/// \author Jonathan Mackey
///
/// - 2009-06-03 Adapted from trunks.cc
/// - 2009-06-11 Adapted from map_etrunks.cc
/// - 2009-06-23 Added in Doppler Broadening for line-of-sight velocity profiles.
/// - 2009-06-25 Split into other files; main() moved to main_projection.cc
/// - 2010-03-22 JM: Moved to MHD_ET_2010, modified to also get
///   projected magnetic field Q and U components
/// - 2010-04-22 JM: Fixed a bug where in get_point_v_los() the call
///   to get_point_neutraldensity() had hardcoded the ionisation
///   fraction element of the state vector to be 5, rather than the
///   value passed to the function.  This is only a problem for MHD
///   runs where 5 is an incorrect value.
/// - 2010.12.13 JM: Added  NEW_STOKES_CALC ifdef to Makefile; the new
///    code in the ifdef does a different Stokes Q,U calculation and
///    replaces the projected Bx,By with values calculated from Q,U.
/// - 2013.05.24 JM: Modified H-alpha emission calculation to output
///    intensity in units of erg/cm2/s/arcsec (set absorption to 0).
/// - 2013.10.15 JM: Updated to use microphysics classes to get the
///    gas temperature and n(H), and added [N II] forbidden line
///    emission.
/// - 2015.07.03 JM: Got rid of NEW_STOKES_CALC (because old code was
///    broken, so I just deleted the old code).
/// - 2015.07.03 JM: updated for pion_dev: uses MCMD, SimSetup,
///    constants.h
/// - 2015.07.13 JM: Multithreaded add_integration_pts_to_pixels

//
// File to analyse a sequence of files from a photo-evaporating random clumps
// simulation.  First we get the directory listing, then for each file we load
// it onto the grid, run some analysis on it, output results to a file, and
// continue to the next file.
//
// The main thing this code does it creates a 2D image of either projected,
// mean, or integrated quantities along the line of sight in each pixel.
// The lines of sight can be aligned with the computational grid, in which 
// case it is very easy; or else at an angle to the grid, in which case we
// need to do a lot more work.
//
// Run with e.g.: 
// ./sim_projection /mnt/local/jm/mysims/multi3d_v3/ R3d_n192_M110_std_0000 M110std_vprof_const 1 1 YP ZP 40 2 64 -8.0e5 8.0e5 1
// ./sim_projection /mnt/local/jm/mysims/multi3d_v3/ R3d_n192_M110_std_0000 M110std_vprof_const 1 1 YP ZP 40 2 64 -8.0e5 8.0e5 1
// ./sim_projection /mnt/local/jm/mysims/multi3d_v3/ R3d_n192_M110_c15_0000 vprof_c15 1 1 YP ZP 40 2 64 -8.0e5 8.0e5 [1/2]
// ./sim_projection /mnt/local/jm/mysims/multi3d_v3/ R3d_n192_M110_std_0000 vprof_std 1 1 YP ZP 40 2 64 -8.0e5 8.0e5 [1/2]
// ./sim_projection /mnt/projects/astrophysics/jmackey/EagleNebula/multi3d_v2/ R3dnew_n192_M9_0000.0 test2 1 1 YP ZP 40 2 20 -10.0e5 10.0e5
//



#include "sim_projection.h"




// void print_array(string name, double *arr, int nel)
// {
//   cout <<name<<": [";
//   for (int v=0;v<nel-1;v++) cout <<arr[v]<<", ";
//   cout <<arr[nel-1]<<"]\n";
//   return;
// }

// ------------------------------------------------------------
// ************************************************************
// ***************  END OF CLASS DECLARATIONS *****************
// ------------------------------------------------------------



// -------------------------------------------------------------
// *************************************************************
// ********************** AXES_DIRECTIONS **********************
// *************************************************************
// -------------------------------------------------------------

enum direction axes_directions::get_posdir(const enum axes a)
{
  int i=static_cast<int>(a);
  if (i<0 || i>2) {
    cout <<"BAD AXIS passed to axes_directions::get_posdir(): "<<i<<"\n";
    return NO;
  }
  return static_cast<direction>(2*i+1);
}



// ##################################################################
// ##################################################################


enum direction axes_directions::get_negdir(const enum axes a)
{
  int i=static_cast<int>(a);
  if (i<0 || i>2) {
    cout <<"BAD AXIS passed to axes_directions::get_negdir(): "<<i<<"\n";
    return NO;
  }
  return static_cast<direction>(2*i);
}



// ##################################################################
// ##################################################################


enum axes axes_directions::get_axis_from_dir(const enum direction dir)
{
  enum axes a;
  if      (dir==XP || dir==XN) a=XX;
  else if (dir==YP || dir==YN) a=YY;
  else if (dir==ZP || dir==ZN) a=ZZ;
  else {
    a=XX; rep.error("Bad direction to get axis from",dir);
  }
  return a;
}



// ##################################################################
// ##################################################################


enum direction axes_directions::cross_product(const enum direction d1,
					      const enum direction d2
					      )
{
  if (d1==d2) return NO;
  if (d1==NO || d2==NO) return NO;
  if      (d1==XP && d2==XN) return NO;
  else if (d1==XP && d2==YN) return ZN;
  else if (d1==XP && d2==YP) return ZP;
  else if (d1==XP && d2==ZN) return YP;
  else if (d1==XP && d2==ZP) return YN;

  else if (d1==XN && d2==XP) return NO;
  else if (d1==XN && d2==YN) return ZP;
  else if (d1==XN && d2==YP) return ZN;
  else if (d1==XN && d2==ZN) return YN;
  else if (d1==XN && d2==ZP) return YP;

  else if (d1==YN && d2==XN) return ZN;
  else if (d1==YN && d2==XP) return ZP;
  else if (d1==YN && d2==YP) return NO;
  else if (d1==YN && d2==ZN) return XP;
  else if (d1==YN && d2==ZP) return XN;

  else if (d1==YP && d2==XN) return ZP;
  else if (d1==YP && d2==XP) return ZN;
  else if (d1==YP && d2==YN) return NO;
  else if (d1==YP && d2==ZN) return XN;
  else if (d1==YP && d2==ZP) return XP;

  else if (d1==ZN && d2==XN) return YP;
  else if (d1==ZN && d2==XP) return YN;
  else if (d1==ZN && d2==YN) return XN;
  else if (d1==ZN && d2==YP) return XP;
  else if (d1==ZN && d2==ZP) return NO;

  else if (d1==ZP && d2==XN) return YN;
  else if (d1==ZP && d2==XP) return YP;
  else if (d1==ZP && d2==YN) return XP;
  else if (d1==ZP && d2==YP) return XN;
  else if (d1==ZP && d2==ZN) return NO;

  else rep.error("Unhandled combination of directions in cross_product!!!!",d1);
  return NO;
}
// -------------------------------------------------------------
// *************************************************************
// ********************** AXES_DIRECTIONS **********************
// *************************************************************
// -------------------------------------------------------------



// -------------------------------------------------------------
// *************************************************************
// ******************* COORDINATE_CONVERSION *******************
// *************************************************************
// -------------------------------------------------------------

coordinate_conversion::coordinate_conversion(const enum direction los, ///< Line of sight direction
					     const int angle,          ///< Angle of LOS w.r.t. los direction.
					     const enum direction up_dir, ///< vertical direction, which stays in image plane.
					     class GridBaseClass *ggg ///< pointer to grid of data.
					     )
{
  gptr = ggg;
  if (!gptr) rep.error("Need a valid grid pointer to set up image!!!",gptr);
  //
  // First get grid domain boundaries and cell sizes:
  //
  for (int v=0; v<3; v++) {
    sim_xminP[v]  = gptr->Xmin(static_cast<axes>(v));
    sim_xmaxP[v]  = gptr->Xmax(static_cast<axes>(v));
    sim_rangeP[v] = gptr->Range(static_cast<axes>(v));
  }
  sim_dxP = gptr->DX();
  if (gptr->Ndim() !=3) rep.error("Need 3D sim for projections",gptr->Ndim());

  sim_dxI = CI.get_integer_cell_size();
  CI.get_ipos_vec(sim_xminP, sim_xminI);
  CI.get_ipos_vec(sim_xmaxP, sim_xmaxI);
  for (int v=0; v<3; v++) {
    sim_rangeI[v] = sim_xmaxI[v]-sim_xminI[v];
    sim_ncell[v]  = sim_rangeI[v]/sim_dxI;
    // SANITY CHECK!
    if (sim_ncell[v] != SimPM.NG[v])
      rep.error("Cells dont match at all!!!", sim_ncell[v]-SimPM.NG[v]);
  }

  //
  // Set image and sim. orientations.  It is assumed the image plane is 
  // X-Y, and the normal is at an angle Theta to Z, and in the XZ plane.
  //
  // img_dir[0]=XP; img_dir[1]=YP; img_dir[2]=ZP;
  // img_axs[0]=XX; img_axs[1]=YY; img_axs[2]=ZZ;
  // NOTE THE IMAGE COORDINATE SYSTEM IS A LEFT-HANDED SYSTEM!!!
  //
  sd[2] = los;
  sd[1] = up_dir;
  sd[0] = cross_product(los, up_dir);
  for (int v=0;v<3;v++) {
    sa[v] = get_axis_from_dir(sd[v]);

    if (sd[v]==NO) rep.error("one direction is NO in image normals",v);
    else if (sd[v]==XN || sd[v]==YN || sd[v]==ZN)
      ss[v] = -1;
    else
      ss[v] =  1;
  }


  //
  // Angles
  //
  theta_deg = angle;
  if (theta_deg==0)
    zero_angle = true;
  else
    zero_angle = false;

  if (abs(theta_deg)>45)
    rep.error("Angle must be in range [-45,45].  For larger angle project along a different axis",theta_deg);

  theta = M_PI*theta_deg/180.0;
  sintheta = sin(theta);
  costheta = cos(theta);
  tantheta = tan(theta);


  //
  // Set image pixel size equal to unity (square pixels).
  //
  im_dx=1;

  //
  // Set the domain in image coordinates, and get number of pixels from that.
  //
  set_sim_extents_in_image_coords();
  set_npix();

  //
  // Output some information:
  //
  cout <<"image directions: x="<<sd[0]<<", y="<<sd[1];
  cout <<", los is "<<theta_deg<<" degrees from direction "<<sd[2]<<endl;

  return;
}



// ##################################################################
// ##################################################################


coordinate_conversion::~coordinate_conversion()
{}



// ##################################################################
// ##################################################################


void coordinate_conversion::set_sim_extents_in_image_coords()
{
  //
  // Bottom, left, near corner of image domain is [0,0,0]
  //
  for (int v=0;v<3;v++) s_xmin_img[v] = 0.0;
  //
  // Height of domain is number of cells in y-dir.
  // IMPLICITLY ASSUME A PIXEL IS UNIT WIDTH/AREA.
  //
  s_xmax_img[YY] = sim_ncell[sa[YY]];
  //
  // For XX and ZZ, we have a rotated box, so use cosines/sines
  // If theta<0, then -sin(t)>0, which is what we want.
  // Could have just used fabs() I guess.
  //
  // We set the simulation origin in the X-Z plane as the most negative
  // point relative to our sd[XX],sd[ZZ] directions.
  //
  if (zero_angle) {
    s_xmax_img[XX] = sim_ncell[sa[XX]];
    s_xmax_img[ZZ] = sim_ncell[sa[ZZ]];
    s_origin_img[XX] = 0.0;
    s_origin_img[ZZ] = 0.0;
  }
  else if (theta>0.0) {
    s_xmax_img[XX] = sim_ncell[sa[XX]]*costheta +sim_ncell[sa[ZZ]]*sintheta;
    s_xmax_img[ZZ] = sim_ncell[sa[XX]]*sintheta +sim_ncell[sa[ZZ]]*costheta;
    s_origin_img[XX] = sim_ncell[sa[ZZ]]*sintheta;
    s_origin_img[ZZ] = 0.0;;
    //cout <<"HIJONO:::"; rep.printVec("ncell",sim_ncell,3);
    //cout <<"HIJONO:::"; rep.printVec("xmaxI",s_xmax_img,3);
    //cout <<"HIJONO:::"; rep.printVec("orign",s_origin_img,3);
  }
  else {
    s_xmax_img[XX] =  sim_ncell[sa[XX]]*costheta -sim_ncell[sa[ZZ]]*sintheta;
    s_xmax_img[ZZ] = -sim_ncell[sa[XX]]*sintheta +sim_ncell[sa[ZZ]]*costheta;
    s_origin_img[XX] = 0.0;
    s_origin_img[ZZ] = -sim_ncell[sa[XX]]*sintheta; // >0 b/c sin(t)<0
    //cout <<"HIJONO:::"; rep.printVec("ncell",sim_ncell,3);
    //cout <<"HIJONO:::"; rep.printVec("xmaxI",s_xmax_img,3);
    //cout <<"HIJONO:::"; rep.printVec("orign",s_origin_img,3);
  }

  s_origin_img[YY]=0.0; // This value is never used.

  return;
}



// ##################################################################
// ##################################################################



bool coordinate_conversion::point_in_Isim_domain(const double *x /// Point in sim coords (integer)
						 )
{
  bool inside=true;
  for (int v=0;v<3;v++) {
    if (x[v]<=sim_xminI[v] || x[v]>=sim_xmaxI[v]) inside=false;
  }
  return inside;
}
    


// ##################################################################
// ##################################################################


void coordinate_conversion::set_npix()
{
  //
  // Adding 1-0.5*(cos(t)+sin(t)) ensures that no cell-centres can be outside
  // the image, and also that no pixel has no cell-centres in it!
  // (cos(t)+sin(t))/2 is the projected distance from the edge cell to the corner of 
  // the domain, so adding 1-(this) means we are exactly 1 cell further (in projected
  // coords) than the last one.  So taking the integer part of this means the last pixel
  // always has at least one cell in it.
  // Note we need to take the abs value of sintheta, whereas costheta>0 on [-45,45]
  //
  im_npixels=1;
  for (int v=0;v<2;v++) {
    im_npix[v] = static_cast<int>(s_xmax_img[v] -s_xmin_img[v] +(1.0-0.5*(costheta+fabs(sintheta))));
    im_npixels *= im_npix[v];
    cout <<"\t\tdirection: "<<v<<"\t got "<<im_npix[v]<<" pixels, so "<<im_npixels<<" in total.\n";
  }

  return;
}  



// ##################################################################
// ##################################################################



void coordinate_conversion::get_npix(int *n ///< 2D array to put number of pixels in each direction.
				     )
{
  for (int v=0;v<2;v++) n[v] = im_npix[v];
  return;
}



// ##################################################################
// ##################################################################



void coordinate_conversion::get_image_Ipos(const int *spos, ///< input position, sim coords, integer units.
					   double *im_pos   ///< converted position in image coords.
					   )
{
  //
  // first get delta, the distance between the left hand edge of the sim and 
  // the point in question.  Then divide by dx to get it in units of number of cells, 
  // which is the image unit.
  //
  double delta[3];
  for (int v=0;v<3;v++) {
    if (ss[v]>0)
      delta[v] = static_cast<double>(spos[sa[v]] - sim_xminI[sa[v]])/sim_dxI;
    else 
      delta[v] = static_cast<double>(sim_xmaxI[sa[v]] - spos[sa[v]])/sim_dxI;
  }
  //
  // Now we have this, I can use simple geometry to get from the sim 'origin' to
  // the point in question.
  //
  im_pos[XX] = s_origin_img[XX] + delta[XX]*costheta - delta[ZZ]*sintheta;
  im_pos[YY] = delta[YY];
  im_pos[ZZ] = s_origin_img[ZZ] + delta[XX]*sintheta + delta[ZZ]*costheta;

  return;
}



// ##################################################################
// ##################################################################



void coordinate_conversion::get_image_Dpos(const double *spos, ///< input position, sim coords, integer units.
					   double *im_pos      ///< converted position in image coords.
					   )
{
  //
  // first get delta, the distance between the left hand edge of the sim and 
  // the point in question.  Then divide by dx to get it in units of number of cells, 
  // which is the image unit.
  //
  double delta[3];
  for (int v=0;v<3;v++) {
    if (ss[v]>0)
      delta[v] = (spos[sa[v]] - sim_xminI[sa[v]])/sim_dxI;
    else 
      delta[v] = (sim_xmaxI[sa[v]] - spos[sa[v]])/sim_dxI;
  }
  //
  // Now we have this, I can use simple geometry to get from the sim 'origin' to
  // the point in question.
  //
  im_pos[XX] = s_origin_img[XX] + delta[XX]*costheta - delta[ZZ]*sintheta;
  im_pos[YY] = delta[YY];
  im_pos[ZZ] = s_origin_img[ZZ] + delta[XX]*sintheta + delta[ZZ]*costheta;


  return;
}



// ##################################################################
// ##################################################################



void coordinate_conversion::get_sim_Dpos(const double *im_pos, ///< position in image coordinates.
					 double *spos        ///< position in sim coords (dx=2).
					 )
{
  //
  // First get the deltas, which are perpendicular distances from the
  // negative (in image coords) edges of the simulation box.
  //
  double delta[3];
  delta[XX] = (im_pos[XX]-s_origin_img[XX])*costheta +
    (im_pos[ZZ]-s_origin_img[ZZ])*sintheta;
  delta[YY] = im_pos[YY];
  delta[ZZ] = (im_pos[ZZ]-s_origin_img[ZZ])*costheta -
    (im_pos[XX]-s_origin_img[XX])*sintheta;

  //
  // Now based on whether s[v] is positive or negative we set the 
  // actual position relative to xmax or xmin.
  //
  for (int v=0;v<3;v++) {
    if (ss[v]>0)
      spos[sa[v]] = delta[v]*sim_dxI +sim_xminI[sa[v]];
    else 
      spos[sa[v]] = sim_xmaxI[sa[v]] -delta[v]*sim_dxI;
  }

  return;
}

// -------------------------------------------------------------
// *************************************************************
// ******************* COORDINATE_CONVERSION *******************
// *************************************************************
// -------------------------------------------------------------



// -------------------------------------------------------------
// *************************************************************
// ************************ IMAGE CLASS ************************
// *************************************************************
// -------------------------------------------------------------



// ##################################################################
// ##################################################################



image::image(const enum direction los, ///< Line of sight direction
	const int t,            ///< Angle of LOS w.r.t. los direction.
	const enum direction perp, ///< vertical direction, which stays in image plane.
	class GridBaseClass *ggg ///< pointer to grid of data.
	)
  :
  coordinate_conversion(los,t,perp,ggg)
{
  cell_positions_set=false;

  //
  // coordinate_conversion constructor sets npix[2], and im_npixels,
  // so I can initialise the image array of pixels now.
  //
  pix=0;
  pix = mem.myalloc(pix,im_npixels);

  initialise_pixels();

  return;
}



// ##################################################################
// ##################################################################



image::~image()
{
  if (cell_positions_set)
    delete_cell_positions();

  if (pix) {
    for (int i=0;i<im_npixels;i++) delete_pixel_data(&(pix[i]));
    mem.myfree(pix);
  }

  return;
}



// ##################################################################
// ##################################################################

#ifdef THREADS

struct pix_int_args {
  class image *IMG; ///< pointer to image class.
  struct pixel *px; ///< pointer to pixel
  size_t i;      ///< pixel id
  double sim_dxP;
  double s_xmax_img[3];
};

//
// void function for threading with Andy's threadpool library
//
void calculate_pix_integration_pts(void *arg)
{
  struct pix_int_args *pia = reinterpret_cast<struct pix_int_args *>(arg);
  size_t i = pia->i;
  pixel *p = pia->px;
  class image *img = pia->IMG;
  double sim_dxP = pia->sim_dxP;
  double *s_xmax_img = pia->s_xmax_img;
  //
  // Set integration points dx = half the cell size.
  //
  double hh=0.5;
  p->int_pts.dx      = hh;
  p->int_pts.dx_phys = sim_dxP*hh;
  p->int_pts.npt     = static_cast<int>((s_xmax_img[ZZ]+0.5)/hh) +1;
  p->int_pts.p = mem.myalloc(p->int_pts.p, p->int_pts.npt);
  //
  // Set positions of each point along the line of sight, and 
  // assign neighbouring cells and their associated weights.
  //
  struct point_4cellavg *pt;
  double ppos_isim[3];
  //    cell *c = (*(p->cells.begin()));
  cell *c = p->inpixel;

  for (int ipt=0;ipt<p->int_pts.npt; ipt++) {
    pt = &(p->int_pts.p[ipt]);
      
    pt->pos[XX] = p->ix[XX] +0.5;
    pt->pos[YY] = p->ix[YY] +0.5;
    pt->pos[ZZ] = ipt*hh;
    //
    // Convert image position to a simulation position.
    //
    img->get_sim_Dpos(pt->pos, ppos_isim);
    //
    // pass in position, starting cell, and get out the four surrounding
    // cells (or some nulls if it's not surrounded), and the bilinear
    // interpolation weights associated with each cell.
    //
    img->find_surrounding_cells(ppos_isim, c, pt->ngb, pt->wt);
    //
    // Now for each point, we have set its position, neighbours, weights.
    //
  }

  pia = mem.myfree(pia);
  return;
}
#endif //THREADS



// ##################################################################
// ##################################################################


void image::add_integration_pts_to_pixels()
{
  if (!pix)
    rep.error("Can't add integration points to uninitialised pixels!",pix);

  pixel *p=0;

  for (int i=0;i<im_npixels;i++) {

#ifdef THREADS
    struct pix_int_args *pia = mem.myalloc(pia,1);
    pia->i = static_cast<size_t>(i);
    pia->px = &(pix[i]);
    pia->IMG = this;
    pia->sim_dxP = sim_dxP;
    for (size_t v=0;v<3;v++)
      pia->s_xmax_img[v] = s_xmax_img[v];
    //calculate_pix_integration_pts(reinterpret_cast<void *>(pia));
    tp_addWork(&tp,calculate_pix_integration_pts,
                   reinterpret_cast<void *>(pia),
                   "image::add_integration_pts_to_pixels()");
#else // THREADS

    p = &(pix[i]);
    
    //
    // Set integration points dx = half the cell size.
    //
    double hh=0.5;
    p->int_pts.dx      = hh;
    p->int_pts.dx_phys = sim_dxP*hh;
    p->int_pts.npt     = static_cast<int>((s_xmax_img[ZZ]+0.5)/hh) +1;
    p->int_pts.p = mem.myalloc(p->int_pts.p, p->int_pts.npt);

    //
    // Set positions of each point along the line of sight, and 
    // assign neighbouring cells and their associated weights.
    //
    struct point_4cellavg *pt;
    double ppos_isim[3];
    //    cell *c = (*(p->cells.begin()));
    cell *c = p->inpixel;

    for (int ipt=0;ipt<p->int_pts.npt; ipt++) {
      pt = &(p->int_pts.p[ipt]);
      
      pt->pos[XX] = p->ix[XX] +0.5;
      pt->pos[YY] = p->ix[YY] +0.5;
      pt->pos[ZZ] = ipt*hh;
      //
      // Convert image position to a simulation position.
      //
      get_sim_Dpos(pt->pos, ppos_isim);
      //
      // pass in position, starting cell, and get out the four surrounding
      // cells (or some nulls if it's not surrounded), and the bilinear
      // interpolation weights associated with each cell.
      //
      find_surrounding_cells(ppos_isim, c, pt->ngb, pt->wt);
      //
      // Now for each point, we have set its position, neighbours, weights.
      //

//       if (i==33) {
// 	cout <<"\t-------------------\n";
// 	rep.printVec("AXES",sa,3);
// 	rep.printVec("IMG point",pt->pos,3);
// 	rep.printVec("SIM point",ppos_isim,3);
// 	if (pt->ngb[0]) rep.printVec("SIM cell0",pt->ngb[0]->pos,3);
// 	if (pt->ngb[0]) rep.printVec("IMG cell0",pt->ngb[0]->Ph,3);
// 	if (pt->ngb[1]) rep.printVec("SIM cell1",pt->ngb[1]->pos,3);
// 	if (pt->ngb[1]) rep.printVec("IMG cell1",pt->ngb[1]->Ph,3);
// 	if (pt->ngb[2]) rep.printVec("SIM cell2",pt->ngb[2]->pos,3);
// 	if (pt->ngb[2]) rep.printVec("IMG cell2",pt->ngb[2]->Ph,3);
// 	if (pt->ngb[3]) rep.printVec("SIM cell3",pt->ngb[3]->pos,3);
// 	if (pt->ngb[3]) rep.printVec("IMG cell3",pt->ngb[3]->Ph,3);
// 	rep.printVec("simxmax",sim_xmaxI,3);
// 	cout <<"\twt[] = ["<<pt->wt[0]<<", "<<pt->wt[1]<<", "<<pt->wt[2]<<", "<<pt->wt[3]<<"]";
// 	cout <<"\tsum = "<<pt->wt[0]+pt->wt[1]+pt->wt[2]+pt->wt[3]<<"\n";
// 	cout <<"\t-------------------\n";
//      }

    }

#endif // THREADS

    //cout <<"\t------------------- NEXT PIXEL --------------- \n";
    //
    // Now for each pixel, we have initialised its points, set npt, dx, and dx_phys
    //
  }
  return;
}



// ##################################################################
// ##################################################################



void image::find_surrounding_cells(const double x[3],
				   cell *c,
				   cell *ngb[4],
				   double wt[4]
				   )
{
  if (!point_in_Isim_domain(x)) {
    for (int v=0;v<4;v++) {
      ngb[v] = 0;
      wt[v]  = 0.0;
    }
  }
  else {
    //
    // move to cell nearest point
    // Note sa[YY] is the axis c should be already on, so no
    // need to move in this direction, only sa[XX] and sa[ZZ]
    //
    if (!pconst.equalD(x[sa[YY]],c->pos[sa[YY]])) {
      rep.error("WARNING: find_surrounding_cells() y-values not the same!",x[sa[YY]]-c->pos[sa[YY]]);
    }
    cell *seek=c;
    enum direction d1=NO, d2=NO;
    int sign=0;

    //
    // Move along Z-dir to get to the cell immediately past the position we are 
    // seeking, or the last cell in that direction, whichever is nearer.
    //
    if (x[sa[ZZ]]>seek->pos[sa[ZZ]]) {
      d1=get_posdir(sa[ZZ]); sign= 1;
    }
    else {
      d1=get_negdir(sa[ZZ]); sign=-1;
    }
    while ( gptr->NextPt(seek,d1)!=0 && sign*(x[sa[ZZ]]-seek->pos[sa[ZZ]])>0.0 )
      seek=gptr->NextPt(seek,d1);


    //
    // Move along X-dir now, and we should be at one of the nearest 4 cells.
    //
    if (x[sa[XX]]>seek->pos[sa[XX]]) {
      d1=get_posdir(sa[XX]); sign= 1;
    }
    else {
      d1=get_negdir(sa[XX]); sign=-1;
    }
    while ( gptr->NextPt(seek,d1)!=0 && sign*(x[sa[XX]]-seek->pos[sa[XX]])>0.0 )
      seek=gptr->NextPt(seek,d1);

    //
    // Make sure cell is less than dx*sqrt(2.0) from point in each direction.
    // (sanity check!)
    //
    double dist2 = (x[sa[XX]]-seek->pos[sa[XX]])*(x[sa[XX]]-seek->pos[sa[XX]]);
    dist2 += (x[sa[ZZ]]-seek->pos[sa[ZZ]])*(x[sa[ZZ]]-seek->pos[sa[ZZ]]);
    if (dist2 > 2*sim_dxI*sim_dxI) {
      rep.printVec("posn",x,3); rep.printVec("cell",seek->pos,3);
      rep.error("Nearest cell is more than root2*dx from point in xz plane",dist2);
    }

    // 
    // Figure out which cell we have, assign it, calculate offsets,
    // Then get the four weights and add the neighbouring cells as 
    // the other neighbours.
    //
    // if cx>px {ct=1}  else {ct=0}
    // if cz>pz {ct+=2} else {ct=ct}
    // ct: XNZN=0,XPZN=1,XNZP=2,XPZP=3
    // set ngb[ct] = cell;
    // 
    // d1,d2 are directions from seek-cell to the other corners.
    //
    int ct=0;
    if (x[sa[XX]]>seek->pos[sa[XX]]) {
      d1=get_posdir(sa[XX]); ct+=0;
    }
    else {
      d1=get_negdir(sa[XX]); ct+=1;
    }
    if (x[sa[ZZ]]>seek->pos[sa[ZZ]]) {
      d2=get_posdir(sa[ZZ]); ct+=0;
    }
    else {
      d2=get_negdir(sa[ZZ]); ct+=2;
    }
    ngb[ct] = seek;
    
    //
    // int pos[4][2];
    //   set pos[0] to [0,0]
    //   set pos[1] to [2,0]
    //   set pos[2] to [0,2]
    //   set pos[3] to [2,2]
    // int offset[2] = [cp[sa[XX]]-pos[ct][0], cp[sa[ZZ]]-pos[ct][1]]
    // Add offset to each position.  This sets pos[ct] to the cell's position,
    // and the other neighbours are dx away.
    //
    int pos[4][2] = {{0,0},{sim_dxI,0},{0,sim_dxI},{sim_dxI,sim_dxI}};
    //for (int i=0;i<2;i++)
    //  cout <<"i="<<i<<" and pos[i] = ["<<pos[i][0]<<","<<pos[i][1]<<"]\n";
    int offset[2] = {seek->pos[sa[XX]]-pos[ct][0],
		     seek->pos[sa[ZZ]]-pos[ct][1]};
    for (int v=0;v<4;v++) {
      pos[v][0] += offset[0];
      pos[v][1] += offset[1];
    }

    //
    // now weights can be got from bilinear interpolation formula:
    // FROM: http://en.wikipedia.org/wiki/Bilinear_interpolation
    //   w[0] = (pos3[x]-px)(pos3[z]-pz)/dx/dx  (so over 4)
    //   w[1] = (px-pos2[x])(pos2[z]-pz)/4
    //   w[2] = (pos1[x]-px)(pz-pos1[z])/4
    //   w[3] = (px-pos0[x])(pz-pos0[z])/4
    sign = sim_dxI*sim_dxI;
    wt[0] = (pos[3][0] - x[sa[XX]]) *(pos[3][1] - x[sa[ZZ]]) /sign;
    wt[1] = (x[sa[XX]] - pos[2][0]) *(pos[2][1] - x[sa[ZZ]]) /sign;
    wt[2] = (pos[1][0] - x[sa[XX]]) *(x[sa[ZZ]] - pos[1][1]) /sign;
    wt[3] = (x[sa[XX]] - pos[0][0]) *(x[sa[ZZ]] - pos[0][1]) /sign;

    //
    // Assign other neighbours.
    //
    if      (ct==3) {
      ngb[2] = gptr->NextPt(ngb[ct],d1);
      ngb[1] = gptr->NextPt(ngb[ct],d2);
      if (ngb[1])
	ngb[0] = gptr->NextPt(ngb[1],d1);
      else if (ngb[2])
	ngb[0] = gptr->NextPt(ngb[2],d2);	
      else
	ngb[0] = 0;
    }
    else if (ct==2) {
      ngb[3] = gptr->NextPt(ngb[ct],d1);
      ngb[0] = gptr->NextPt(ngb[ct],d2);
      if (ngb[0])
	ngb[1] = gptr->NextPt(ngb[0],d1);
      else if (ngb[3])
	ngb[1] = gptr->NextPt(ngb[3],d2);
      else
	ngb[1] = 0;
    }
    else if (ct==1) {
      ngb[0] = gptr->NextPt(ngb[ct],d1);
      ngb[3] = gptr->NextPt(ngb[ct],d2);
      if (ngb[0])
	ngb[2] = gptr->NextPt(ngb[0],d2);
      else if (ngb[3])
	ngb[2] = gptr->NextPt(ngb[3],d1);
      else
	ngb[2] = 0;
    }
    else {
      ngb[1] = gptr->NextPt(ngb[ct],d1);
      ngb[2] = gptr->NextPt(ngb[ct],d2);
      if (ngb[1])
	ngb[3] = gptr->NextPt(ngb[1],d2);
      else if (ngb[2])
	ngb[3] = gptr->NextPt(ngb[2],d1);
      else
	ngb[3] = 0;
    }

    //
    // Debug info:
    //
//     cout <<"\t-------------------\n";
//     rep.printVec("point",x,3);
//     if (ngb[0]) rep.printVec("cell0",ngb[0]->pos,3);
//     if (ngb[1]) rep.printVec("cell1",ngb[1]->pos,3);
//     if (ngb[2]) rep.printVec("cell2",ngb[2]->pos,3);
//     if (ngb[3]) rep.printVec("cell3",ngb[3]->pos,3);
//     rep.printVec("corner0",pos[0],2);
//     rep.printVec("corner1",pos[1],2);
//     rep.printVec("corner2",pos[2],2);
//     rep.printVec("corner3",pos[3],2);
//     cout <<"\twt[] = ["<<wt[0]<<", "<<wt[1]<<", "<<wt[2]<<", "<<wt[3]<<"]";
//     cout <<"\tsum = "<<wt[0]+wt[1]+wt[2]+wt[3];
//     cout <<"\t-------------------\n";

  } // else point is in domain...    
  //
  // Now we should have set ngb and weights.
  //
  return;
}



// ##################################################################
// ##################################################################



bool image::cell_is_in_pixel(double *cp, ///< Cell position (in image coordinates).
			     pixel  *p  ///< pixel in question.
			     )
{
  bool inside=true;

  //
  // Pixel indices are zero offset, and pixels are unit diameter,
  // so a pixel centre is (ix[v]+0.5), left edge is ix[v], right
  // edge is (ix[v]+1)
  //
  for (int v=0;v<2;v++) {
    if ( (cp[v] < p->ix[v]) || (cp[v]>(p->ix[v]+1)) )
      inside=false;
  }

  return inside;
}



// ##################################################################
// ##################################################################



void image::add_cells_to_pixels()
{
  if (!pix)
    rep.error("Can't add cells to uninitialised pixels!",pix);

  for (int v=0;v<im_npixels;v++) {
    pix[v].ncells=0;
    //pix[v].cells.clear();
    pix[v].inpixel = 0;
  }
  
  int ix, iy, ipix;
  //double TEMP_maxx=0.0;
  cell *c = gptr->FirstPt();
  do {
    //
    // Positions are guaranteed to be positive in image coordinates.
    //
    ix = static_cast<int>(c->Ph[XX]);
    //TEMP_maxx=max(c->Ph[XX],TEMP_maxx);
    iy = static_cast<int>(c->Ph[YY]);
    ipix = iy*im_npix[XX]+ix;
    if (ipix<0 || ipix>im_npixels)
      rep.error("cell with id following is outside image!",c->id);
    else {
      //      pix[ipix].cells.push_back(c);
      if (pix[ipix].ncells==0) {
	pix[ipix].inpixel = c;
	pix[ipix].ncells++;
      }
    }
    //  } while ( (c=gptr->NextPt(c)) !=0);
  } while ( (c=gptr->NextPt(c)) !=0);

  //
  // Checking to make sure that all pixels have cells in them, and
  // that we are not double counting cells, or missing some of them.
  //
  long int ct=0;
  for (int v=0;v<im_npixels;v++) {
    if (pix[v].ncells<=0) {
      rep.printVec("s_xmax_img",s_xmax_img,3);
      //cout <<" centre of right-most cell="<<TEMP_maxx<<endl;
      rep.error("Some pixels have no cells in them!",v);
    }
    ct += pix[v].ncells;
  }

  //if (ct != SimPM.Ncell)
  //  rep.error("Wrong number of cells in pixels",ct-SimPM.Ncell);

  return;
}



// ##################################################################
// ##################################################################



void image::set_cell_positions_in_image()
{
  //
  // First make sure we are using minimal cells.
  //
  if (!CI.query_minimal_cells()) {
    rep.error("image::set_cell_positions_in_image() needs minimal_cells to be set!",99);
  }
  
  //
  // Now for each cell, allocate Ph[3] and give it a coordinate
  // in the image coordinate system.
  //
  cell *c = gptr->FirstPt();
  do {
    c->Ph = mem.myalloc(c->Ph,3);
    get_image_Ipos(c->pos, ///< input position, sim coords, integer units.
		   c->Ph   ///< output: image coords and units.
		   );
    //cout <<c->id<<"\t\t"; rep.printVec("cell-pos",c->pos,3);
    //cout <<"\t\t"; rep.printVec("img-pos ",c->Ph,3);
  } while ((c=gptr->NextPt(c))!=0);

  //
  // Set flag indicating we have set positions.
  //
  cell_positions_set=true;
  return;
}



// ##################################################################
// ##################################################################



void image::delete_cell_positions()
{
  if (!cell_positions_set) {
    cout << "Don't call del_img_cl_pos() if positions aren't set!!!\n";
    return;
  }

  //
  // For each cell, delete Ph[3] and make it a null pointer.
  //
  cell *c = gptr->FirstPt();
  do {
    c->Ph = mem.myfree(c->Ph);
  } while ((c=gptr->NextPt(c))!=0);

  //
  // Unset cell positions flag.
  //
  cell_positions_set=false;
  return;
}



// ##################################################################
// ##################################################################



void image::initialise_pixels()
{
  pixel *p=0;
  for (int i=0;i<im_npixels;i++) {
    p = &(pix[i]);

    //
    // Set ix,iy,ipix
    //
    p->ipix = i;
    p->ix[XX] = i%im_npix[XX]; // remainder>0 when a,b>0 in a=x*b+r
    p->ix[YY] = i/im_npix[XX];
    p->int_pts.p = 0; // zero this pointer, so i know it's uninitialised.
  }

  return;
}



// ##################################################################
// ##################################################################



void image::delete_pixel_data(pixel *p)
{
  if (p->int_pts.p) {
    p->int_pts.p = mem.myfree(p->int_pts.p);
  }

  return;
}



// ##################################################################
// ##################################################################



double image::get_pt_density(struct point_4cellavg *pt)
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *pt->ngb[v]->P[RO];
  }
  return val;
}



// ##################################################################
// ##################################################################



double image::get_pt_neutral_numberdensity(
        struct point_4cellavg *pt, int ifrac
        )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *(pt->ngb[v]->P[RO]
                                      *(1.0-pt->ngb[v]->P[ifrac]));
  }
  val *= SimPM.EP.H_MassFrac/pconst.m_p();
  //cout <<"image::get_pt_neutral_numberdensity() n(H0) = "<<val<<"\n";
  return val;
}



// ##################################################################
// ##################################################################


//#define SQRT_DENSITY
//#define NO_DENSITY_WT
#define LINMAX_DENSITY
#define MAXDENS 25000.0

double image::get_pt_StokesQ(struct point_4cellavg *pt, const int ,
			     const int bx, const int by, const int bz,
			     const int signx, const int , const int signz,
			     const double sintht, const double costht)
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //

  //
  // New Q calculated from Bx,By according to
  // Q = sum_i wt[i]*|f(n_H)*(Bx^2-By^2)/sqrt(Bx^2+By^2)|_i
  //
  double val=0.0, bx2=0.0, by2=0.0, btot=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      //
      // Get Bx along line of sight:
      //
      bx2 = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      bx2 *= bx2;
      //
      // By is easier, and so then is btot.
      //
      by2 = pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      btot = sqrt(bx2+by2);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p())*(bx2-by2)/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *(bx2-by2)/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *(bx2-by2)/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  //
  // End of new calc
  //

  return val;
}



// ##################################################################
// ##################################################################



double image::get_pt_StokesU(struct point_4cellavg *pt, const int ,
			     const int bx, const int by, const int bz,
			     const int signx, const int signy, const int signz,
			     const double sintht, const double costht)
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //

  //
  // New Q calculated from Bx,By according to
  // Q = sum_i wt[i]*|f(n_H)*(2*Bx*By)/sqrt(Bx^2+By^2)|_i
  //
  double val=0.0, btot=0.0, bxy=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      //
      // Calculate Bx and add Bx^2 to btot
      //
      bxy = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      btot = bxy*bxy;
      //
      // Multipy Bx by By, add By^2 to btot
      //
      bxy *= signy*pt->ngb[v]->P[by];
      btot += pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      //
      // Btot = sqrt(Bx^2+By^2) now.  bxy=Bx*By.
      //
      btot = sqrt(btot);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p())*2.0*bxy/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *2.0*bxy/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *2.0*bxy/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  //
  // End of new calc
  //

  return val;
}



// ##################################################################
// ##################################################################



double image::get_pt_BXabs(struct point_4cellavg *pt, const int,
			   const int bx, const int by, const int bz,
			   const int signx, const int , const int signz,
			   const double sintht, const double costht)
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0, btot=0.0, bx2=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      bx2 = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      bx2 *= bx2;
      btot = sqrt(pt->ngb[v]->P[bx]*pt->ngb[v]->P[bx] +
		  pt->ngb[v]->P[by]*pt->ngb[v]->P[by] +
		  pt->ngb[v]->P[bz]*pt->ngb[v]->P[bz]);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p()) *bx2/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *bx2/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *bx2/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  return val;
}



// ##################################################################
// ##################################################################



double image::get_pt_BYabs(struct point_4cellavg *pt, const int ,
			   const int bx, const int by, const int bz,
			   const int , const int , const int ,
			   const double , const double )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0, btot=0.0, by2=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      by2 = pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      btot = sqrt(pt->ngb[v]->P[bx]*pt->ngb[v]->P[bx] +
		  by2 +
		  pt->ngb[v]->P[bz]*pt->ngb[v]->P[bz]);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p()) *by2/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *by2/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *by2/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  return val;
}



// ##################################################################
// ##################################################################



void image::calculate_pixel(struct pixel *px,                 ///< pointer to pixel
			    const struct vel_prof_stuff *vps, ///< struct with info for velocity binning.
			    const int what_to_integrate,      ///< flag for what to integrate.
			    double *im,                       ///< array of pixel data.
			    double *tot_mass                  ///< general purpose counter for stuff.
			    )
{
  //   cout <<"calculate_pixel: i="<<i<<" what="<<what_to_integrate;
  //   cout <<" mass="<<*tot_mass<<" vps->v_min="<<vps->v_min;
  //   cout <<" im="<<im<<" img="<<IMG<<" px="<<px<<endl;
  
  //
  // Simpson's Rule Integration along the pixel line of sight:
  // integral \simeq \frac{hh}{3} (f(0) +4f(1) +2f(2) +4f(3) +...+4f(n-1) +f(n))
  // for even n, we have npt=n+1 points.
  //
  int wt=2,
    npt=px->int_pts.npt;
  double hh=px->int_pts.dx_phys;
  double ans=0.0;
  
  if       (what_to_integrate==I_DENSITY) {
    ans += get_pt_density(&(px->int_pts.p[0]));
    for (int v=1; v<(npt-1); v++) {
      wt = 6-wt;
      ans += wt *get_pt_density(&(px->int_pts.p[v]));
    }
    ans += get_pt_density(&(px->int_pts.p[npt-1]));
    ans *= hh/3.0;

    *tot_mass += ans;
    im[px->ipix] = ans;
  }
  
  else if (what_to_integrate==I_NEUTRAL_NH) {
    ans += get_pt_neutral_numberdensity(&(px->int_pts.p[0]),SimPM.ftr);
    for (int v=1; v<(npt-1); v++) {
      wt = 6-wt;
      ans += wt *get_pt_neutral_numberdensity(&(px->int_pts.p[v]),SimPM.ftr);
    }
    ans += get_pt_neutral_numberdensity(&(px->int_pts.p[npt-1]),SimPM.ftr);
    ans *= hh/3.0;

    *tot_mass += ans;
    im[px->ipix] = ans;
  }
  
  else if (what_to_integrate==I_B_STOKESQ ||
	   what_to_integrate==I_B_STOKESU ||
	   what_to_integrate==I_BXabs     ||
	   what_to_integrate==I_BYabs) {
    int bx=0,by=0,bz=0;			       
    if      (sa[XX]==XX) bx=BX;
    else if (sa[XX]==YY) bx=BY;
    else if (sa[XX]==ZZ) bx=BZ;
    else rep.error("Bad axis from IMG[x]",sa[XX]);
    if      (sa[YY]==XX) by=BX;
    else if (sa[YY]==YY) by=BY;
    else if (sa[YY]==ZZ) by=BZ;
    else rep.error("Bad axis from IMG[y]",sa[YY]);
    if      (sa[ZZ]==XX) bz=BX;
    else if (sa[ZZ]==YY) bz=BY;
    else if (sa[ZZ]==ZZ) bz=BZ;
    else rep.error("Bad axis from IMG[z]",sa[ZZ]);
    int signx = ss[XX], signy = ss[YY], signz = ss[ZZ];
    double st = sin(static_cast<double>(theta));
    double ct = cos(static_cast<double>(theta));
    //cout <<"using element "<<bz<<" for LOS, and "<<bx<<" for x-component; theta="<<angle<<" deg.\n";

    if      (what_to_integrate==I_B_STOKESQ) {
      ans += get_pt_StokesQ(&(px->int_pts.p[0]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      for (int v=1; v<(npt-1); v++) {
	wt = 6-wt;
	ans += wt *get_pt_StokesQ(&(px->int_pts.p[v]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      }
      ans += get_pt_StokesQ(&(px->int_pts.p[npt-1]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      ans *= 1.0/3.0;
      *tot_mass += ans;
      im[px->ipix] = ans;
    }
    else if (what_to_integrate==I_B_STOKESU) {
      ans += get_pt_StokesU(&(px->int_pts.p[0]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      for (int v=1; v<(npt-1); v++) {
	wt = 6-wt;
	ans += wt *get_pt_StokesU(&(px->int_pts.p[v]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      }
      ans += get_pt_StokesU(&(px->int_pts.p[npt-1]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      ans *= 1.0/3.0;
      *tot_mass += ans;
      im[px->ipix] = ans;
    }
    else if (what_to_integrate==I_BXabs) {
      ans += get_pt_BXabs(&(px->int_pts.p[0]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      for (int v=1; v<(npt-1); v++) {
	wt = 6-wt;
	ans += wt *get_pt_BXabs(&(px->int_pts.p[v]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      }
      ans += get_pt_BXabs(&(px->int_pts.p[npt-1]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      ans *= 1.0/3.0;
      *tot_mass += ans;
      im[px->ipix] = ans;
    }
    else if (what_to_integrate==I_BYabs) {
      ans += get_pt_BYabs(&(px->int_pts.p[0]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      for (int v=1; v<(npt-1); v++) {
	wt = 6-wt;
	ans += wt *get_pt_BYabs(&(px->int_pts.p[v]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      }
      ans += get_pt_BYabs(&(px->int_pts.p[npt-1]),SimPM.ftr, bx,by,bz,signx,signy,signz,st,ct);
      ans *= 1.0/3.0;
      *tot_mass += ans;
      im[px->ipix] = ans;
    }
    else rep.error("Bad what to integrate -- Stokes Q/U or BX/BY",what_to_integrate);
  } // I_STOKES Q/U or BX/BY
  
  else if (what_to_integrate==I_VX) {
    //
    // This is quite simple -- VX is just a scalar quantity at each point, but we do want
    // a profile of it so it is more complicated than the density above.
    // Velocities are in km/s
    //
    int vx = VX;
    int vz = VZ;
    int Nbins = vps->npix[2];
    double bin_size = (vps->v_max-vps->v_min)/Nbins;
    double profile[Nbins];
    double temp_profile[Nbins];
    for (int v=0;v<Nbins;v++) profile[v] = temp_profile[v] = 0.0;
    class point_velocity VELX(vx, vz,
			      ss[XX], ss[ZZ],
			      theta,
			      vps->v_min, vps->v_max, Nbins);
    //
    // Add line broadening, if needed, second arg is FWHM in code's velocity units.
    // If type==1, then smooth everything with FWHM.
    // If type==2, use point temperature to do standard doppler broadening at each pt.
    //
    if      (vps->smooth==1) {
      VELX.set_broadening(1,vps->broadening);
      cout <<"Broadening cell velocities by "<<vps->broadening<<" (FWHM) in velocity units.\n";
    }
    else if (vps->smooth==2) {
      VELX.set_broadening(2,0.0);
      //cout <<"Broadening cell velocities by Doppler (thermal) broadening, on point-by-point basis.\n";
    }
    //
    // Get VX and density at each point.  Velocity determines which
    // bin this point contributes too; density determines the amplitude it adds to 
    // that bin.  The get_point_VX_profile() function does the binning for us.
    //  This is Simpson's Rule Integration.
    //
    // The function get_point_VX_profile() puts the mass in the right velocity bin
    // in the temp_profile[] array, and if smooth==2 it will smooth this with a 
    // Gaussian corresponding to the Doppler broadening from the point's temperature.
    //
    VELX.get_point_VX_profile(&(px->int_pts.p[0]), temp_profile, SimPM.ftr);
    for (int v=0;v<Nbins;v++) profile[v] +=  temp_profile[v];
    
    for (int ii=1; ii<(npt-1); ii++) {
      wt=6-wt;
      VELX.get_point_VX_profile(&(px->int_pts.p[ii]), temp_profile, SimPM.ftr);
      for (int v=0;v<Nbins;v++) profile[v] +=  wt*temp_profile[v];
      //cout <<"prof: ";
      //for (int v=0;v<Nbins;v++) cout <<" "<<temp_profile[v];
      //cout <<endl;
    }
    
    VELX.get_point_VX_profile(&(px->int_pts.p[npt-1]), temp_profile, SimPM.ftr);
    for (int v=0;v<Nbins;v++) profile[v] +=  temp_profile[v];
    
    for (int v=0;v<Nbins;v++) profile[v] *= hh/3.0;
    
    //
    // Smooth profile: This only does anything if smooth==1, so constant smoothing.
    // If smooth==2, then the line broadening is already done on point-by-point basis.
    // if smooth==0, then nothing happens either.
    //
    if      (vps->smooth==1) {
      VELX.smooth_profile_FFT(profile);
    }
    //
    // Now put each bin in the image, where velocity is the third dimension in the
    // data cube.  So we skip npixel elements for each profile element.
    //
    for (int v=0;v<Nbins;v++)
      im[px->ipix +v*im_npixels] = profile[v];
    
    *tot_mass = 0.0;
  }

  else if (what_to_integrate==I_VEL_LOS) {
    //
    // Line of Sight Velocity is more complicated to calculate, because I need 
    // to put the data in a range of bins.
    //
    // First we set the angles and normals, to choose how much of each velocity
    // component contributes to the line of sight.
    //
    //double ct = cos(static_cast<double>(angle)*M_PI/180.0);
    //double st = sin(static_cast<double>(angle)*M_PI/180.0);
    
    int vx=0,vz=0;
    
    if      (sa[ZZ]==XX) vz=VX;
    else if (sa[ZZ]==YY) vz=VY;
    else if (sa[ZZ]==ZZ) vz=VZ;
    else rep.error("Bad axis from IMG[z]",sa[ZZ]);
    if      (sa[XX]==XX) vx=VX;
    else if (sa[XX]==YY) vx=VY;
    else if (sa[XX]==ZZ) vx=VZ;
    else rep.error("Bad axis from IMG[x]",sa[XX]);
    //cout <<"using element "<<vz<<" for LOS, and "<<vx<<" for x-component; theta="<<angle<<" deg.\n";
    
    //
    // Should make this be dynamic:
    // Velocities are in km/s
    //
    int Nbins = vps->npix[2];
    double bin_size = (vps->v_max-vps->v_min)/Nbins;

    double profile[Nbins];
    double temp_profile[Nbins];
    for (int v=0;v<Nbins;v++) profile[v] = temp_profile[v] = 0.0;
    //double vel=0.0, norm=0.0;
    //int ibin;
    
    //
    // Set up class for calculating points.
    //
    class point_velocity VLOS(vx, vz,
			      ss[XX], ss[ZZ],
			      theta,
			      vps->v_min, vps->v_max, Nbins);
    //
    // Add line broadening, if needed, second arg is FWHM in code's velocity units.
    // If type==1, then smooth everything with FWHM.
    // If type==2, use point temperature to do standard doppler broadening at each pt.
    //
    if      (vps->smooth==1) {
      VLOS.set_broadening(1,vps->broadening);
      //cout <<"Broadening cell velocities by "<<vps->broadening<<" (FWHM) in velocity units.\n";
    }
    else if (vps->smooth==2) {
      VLOS.set_broadening(2,0.0);
      //cout <<"Broadening cell velocities by Doppler (thermal) broadening, on point-by-point basis.\n";
    }
    
    //
    // Get the LOS velocity and density at each point.  Velocity determines which
    // bin this point contributes too; density determines the amplitude it adds to 
    // that bin.  This is Simpson's Rule Integration.
    //
    // The function get_point_v_los_profile() puts the mass in the right velocity bin
    // in the temp_profile[] array, and if smooth==2 it will smooth this with a 
    // Gaussian corresponding to the Doppler broadening from the point's temperature.
    //
    VLOS.get_point_v_los_profile(&(px->int_pts.p[0]), temp_profile, SimPM.ftr);
    for (int v=0;v<Nbins;v++) profile[v] +=  temp_profile[v];
    
    for (int ii=1; ii<(npt-1); ii++) {
      wt=6-wt;
      VLOS.get_point_v_los_profile(&(px->int_pts.p[ii]), temp_profile, SimPM.ftr);
      for (int v=0;v<Nbins;v++) profile[v] +=  wt*temp_profile[v];
    }
    
    VLOS.get_point_v_los_profile(&(px->int_pts.p[npt-1]), temp_profile, SimPM.ftr);
    for (int v=0;v<Nbins;v++) profile[v] +=  temp_profile[v];
    
    for (int v=0;v<Nbins;v++) profile[v] *= hh/3.0;
    
    //
    // Smooth profile: This only does anything if smooth==1, so constant smoothing.
    // If smooth==2, then the line broadening is already done on point-by-point basis.
    // if smooth==0, then nothing happens either.
    //
    if      (vps->smooth==1) {
      VLOS.smooth_profile_FFT(profile);
    }
    
    //
    // Euler Integration:
    //
    //VLOS.get_point_v_los_profile(&(px->int_pts.p[0]), temp_profile);
    //for (int v=0;v<Nbins;v++) profile[v] +=  1.e19*0.5*temp_profile[v];
    //for (int ii=1; ii<(npt-1); ii++) {
    //  VLOS.get_point_v_los_profile(&(px->int_pts.p[ii]), temp_profile);
    //  for (int v=0;v<Nbins;v++) profile[v] +=  1.e19*temp_profile[v];
    //}
    //VLOS.get_point_v_los_profile(&(px->int_pts.p[npt-1]), temp_profile);
    //for (int v=0;v<Nbins;v++) profile[v] +=  1.e19*0.5*temp_profile[v];
    //for (int v=0;v<Nbins;v++) profile[v] *= hh;
    //VLOS.smooth_profile(profile);
    
    //
    // Now put each bin in the image, where velocity is the third dimension in the
    // data cube.  So we skip npixel elements for each profile element.
    //
    for (int v=0;v<Nbins;v++)
      im[px->ipix +v*im_npixels] = profile[v];
    
    *tot_mass = 0.0;
  }
  
  else if (what_to_integrate==I_EMISSION) {
    //
    // Set up class for calculating points.
    // None of the parameters mean anything for calculating emission;
    // they're just for velocity.
    //
    class point_velocity VLOS(1,1,1,1,1.0,1.0,2.0,1);
    double alpha=0.0, j=0.0, dtau=0.0;
    ans=0.0;
    //
    // Now do an Euler Integration with forward differencing.
    // This uses Rybicki & Lightman's basic solution for constant
    // source function, eq. 1.30, and ignores scattering.
    // 
    for (int v=0; v<npt; v++) {
       VLOS.get_point_Halpha_params(&(px->int_pts.p[v]), SimPM.ftr, &alpha, &j);
       dtau = alpha*hh;
       if (dtau < 1e-7) {
	 //
	 // We can approximate the exponential by Taylor Series to avoid roundoff errors.
	 // 
	 ans = ans*(1.0-dtau) +j*hh;
       }
       else {
	 //
	 // Don't approximate the exponentials because absorption is quite strong
	 // 
	 ans = ans*exp(-dtau) +j*(1.0-exp(-dtau))/alpha;
	 //cout <<" dtau = "<<dtau<<endl;
       }
    }

    *tot_mass += 0.0;  // this means nothing for emission integration...
    im[px->ipix] = ans; ///1.0e80;
    //cout <<"\t\tans = "<<ans<<endl;
  } // I_EMISSION
  
  else if (what_to_integrate==I_NII6584) {
    //
    // Set up class for calculating points.
    // None of the parameters mean anything for calculating [NII];
    // they're just for velocity.
    //
    class point_velocity VLOS(1,1,1,1,1.0,1.0,2.0,1);
    double alpha=0.0, j=0.0, dtau=0.0;
    ans=0.0;
    //
    // Now do an Euler Integration with forward differencing.
    // This uses Rybicki & Lightman's basic solution for constant
    // source function, eq. 1.30, and ignores scattering.
    // 
    for (int v=0; v<npt; v++) {
       VLOS.get_point_NII6584_params(&(px->int_pts.p[v]), SimPM.ftr, &alpha, &j);
       dtau = alpha*hh;
       if (dtau < 1e-7) {
	 //
	 // We can approximate the exponential by Taylor Series to avoid roundoff errors.
	 // 
	 ans = ans*(1.0-dtau) +j*hh;
       }
       else {
	 //
	 // Don't approximate the exponentials because absorption is quite strong
	 // 
	 ans = ans*exp(-dtau) +j*(1.0-exp(-dtau))/alpha;
	 //cout <<" dtau = "<<dtau<<endl;
       }
    }

    *tot_mass += 0.0;  // this means nothing for emission integration...
    im[px->ipix] = ans; ///1.0e80;
    //cout <<"\t\tans = "<<ans<<endl;
  } // I_NII6584
  
  else {
    rep.error("don't know what to integrate!",what_to_integrate);
  }
  //  cout <<"finished pixel "<<i<<" at address "<<px<<"\n";
  return;
}

// -------------------------------------------------------------
// *************************************************************
// ************************ IMAGE CLASS ************************
// *************************************************************
// -------------------------------------------------------------




// -------------------------------------------------------------
// *************************************************************
// ******************** POINT VELOCITY CLASS *******************
// *************************************************************
// -------------------------------------------------------------



// ##################################################################
// ##################################################################



point_velocity::point_velocity(const int vx_comp,  ///< velocity component perp. to LOS direction (contributing)
			       const int vz_comp,  ///< velocity componenet along LOS
			       const int signx,    ///< +1 if looking along +ve vx axis; -1 otherwise
			       const int signz,    ///< +1 if looking along +ve vz axis; -1 otherwise
			       const double theta, ///< Angle to LOS (radians)
			       const double velocity_min,      ///< minimum velocity in range
			       const double velocity_max,      ///< maximum velocity in range
			       const int velocity_Nbins        ///< Number of bins.
			       )
{
  vx = vx_comp;
  vz = vz_comp;
  sx = signx;
  sz = signz;
  ct = cos(theta);
  st = sin(theta);
  v_min      = velocity_min;
  v_max      = velocity_max;
  v_Nbins    = velocity_Nbins;
  v_binsize = (v_max-v_min)/v_Nbins;
  return;
}



// ##################################################################
// ##################################################################



void point_velocity::set_broadening(
        const int type,    ///< Type: 1=constant Gaussian broadening.
        const double fwhm ///< FWHM of Gaussian to smooth profile by.
        )
{
  broaden = type;
  //
  // Sigma is the smoothing length, in units of the velocity domain.
  //
  sigma   = 0.424661*fwhm/(v_max-v_min);

  FT_Ng = 1;
  while (FT_Ng<v_Nbins) FT_Ng *=2;
  //FT_Ng*=2;
  //cout <<"v_Nbins="<<v_Nbins<<" and FT_Ng="<<FT_Ng<<endl;

  return;
}




// ##################################################################
// ##################################################################



void point_velocity::four1(
        double *data,         ///< length 2*nn, for nn real and nn imaginary values.
        unsigned long int nn, ///< half the length of the array.
        int isign             ///< =1 for forward transform, =-1 for inverse transform.
        )
{
  //
  // This is the NR version of the FFT.  Instead of trying to convert it
  // all for zero-offset arrays, I copied Martin White, and access the
  // (n-1)th element of data[] whenever I need to.
  //
  unsigned long int n=0, m=0, j=0, i=0, istep=0, mmax=0;
  double theta, wtemp, wr, wpr, wpi, wi;
  double tempr, tempi;

  n = nn*2;
  j=1;
  for (i=1; i<n; i+=2) {
    if (j>i) {
      //cout <<"i="<<i<<" j="<<j<<" m="<<m<<endl;
      std::swap(data[j-1],data[i-1]);
      std::swap(data[  j],data[  i]);
    }
    m = n/2;
    while (m>=2 && j>m) {
      j -= m;
      m /= 2;
    }
    j += m;
  }


  mmax = 2;
  while (n>mmax) {
    istep = mmax*2;
    theta = isign*(M_PI*2.0/mmax);
    wtemp = sin(0.5*theta);
    wpr   = -2.0*wtemp*wtemp;
    wpi   = sin(theta);
    wr    = 1.0;
    wi    = 0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m; i<=n; i+=istep) {
	j = i+mmax;
	tempr = wr*data[j-1] -wi*data[j  ];
	tempi = wr*data[j  ] +wi*data[j-1];
	data[j-1] = data[i-1] -tempr;
	data[j  ] = data[i  ] -tempi;
	data[i-1] += tempr;
	data[i  ] += tempi;
      }
      wr = (wtemp=wr)*wpr -wi*wpi +wr;
      wi = wi*wpr +wtemp*wpi +wi;
    }
    mmax = istep;
  }

  return;
}



// ##################################################################
// ##################################################################



void point_velocity::smooth_profile_FFT(double *data ///< Array of velocity bins to smooth.
				    )
{
  //
  // This works with Fourier Transform methods, rather than doing the convolution
  // integral.  We remove the mean, zero pad up to a power of 2, FT, multiply by 
  // the FT of the Gaussian, IFT, and add back the mean.
  //
  if (broaden!=1) {
    //cout <<" no need to smooth profile; broaden="<<broaden<<endl;
    return;
  }

  //
  // Remove the mean
  //
  double mean=0.0;
  for (int v=0;v<v_Nbins;v++)
    mean += data[v];
  mean /= v_Nbins;
  for (int v=0;v<v_Nbins;v++)
    data[v] -= mean;

  //
  // Put data into Cx. array, with zero padding.
  //
  //cout <<"v_Nbins="<<v_Nbins<<" and FT_Ng="<<FT_Ng<<endl;
  double cxdata[2*FT_Ng];
  for (int v=0;v<2*FT_Ng;v++)
    cxdata[v]=0.0;
  for (int v=0;v<v_Nbins;v++)
    cxdata[2*v] = data[v];

  //
  // FT the data
  //
  four1(cxdata, static_cast<unsigned long>(FT_Ng), 1);
  //
  // Multiply by FT of gaussian.
  // The FT is in frequency units f=k/(2PI)
  // for each f, the FT of the Gaussian is exp(-2*PI*PI*sigma*sigma*f*f)
  // *BUT* remember the FFT returns frequencies in the range -N/2,N/2, so
  // I need to do the numbering carefully.
  // Martin does it in smoothmap.c, where if i>FT_Ng/2 f=FT_Ng-i; else f=i;
  // Because we are squaring, the sign doesn't matter.
  int f;
  double expfactor;
  for (int v=0; v<FT_Ng; v++) {
    if (v>FT_Ng/2) f = FT_Ng-v;
    else         f = v;
    expfactor = 2*M_PI*M_PI*sigma*sigma*f*f;
    //cout <<"expfactor="<<expfactor<<endl;
    expfactor = exp(-expfactor);
    //
    // We divide by FT_Ng because the IFT returns FT_Ng times the function.
    //
    cxdata[2*v+0] *= expfactor/FT_Ng;
    cxdata[2*v+1] *= expfactor/FT_Ng;
  }
  
  
  //
  // IFT the data
  //
  four1(cxdata, static_cast<unsigned long>(FT_Ng), -1);
 
  //
  // Move from Cx. array into original array, and add back mean.
  //
  for (int v=0;v<v_Nbins;v++)
    data[v] = cxdata[2*v] +mean;

  return;
}



// ##################################################################
// ##################################################################



void point_velocity::get_point_VX_profile(
        const struct point_4cellavg *pt, ///< point to add to profile.
        double *prof, ///< Array of velocity bins to put profile into.
        const int ifrac ///< index of i-fraction in P.V.
        )
{
  //
  // Gets the velocity profile, using VX only, ignoring projection
  // and orientation effects.  This is NOT a realistic observation,
  // but just used to give a picture of what VX is doing along the
  // length of a pillar.  It is identical to the more realistic
  // get_point_v_los_profile() function below except that it calls a
  // different function to get the point's velocity.
  //
  double vel, norm;
  int ibin;
  vel  = get_point_VX(pt);
  norm = get_point_neutralH_numberdensity(pt,ifrac);
  ibin = get_velocity_bin_number(vel);
  for (int v=0;v<v_Nbins;v++) prof[v] = 0.0;
  if (ibin>0 && ibin<v_Nbins)
    prof[ibin] = norm;
  //cout <<"vx="<<vel<<" norm="<<norm<<" prof["<<ibin<<"] = "<<prof[ibin]<<endl;
  if (broaden==2) {
    broaden_profile(pt,prof,ifrac,vel,norm);
  }
  //for (int v=0;v<v_Nbins;v++) cout <<" "<<prof[v];
  //cout <<endl;  
  return;
}



// ##################################################################
// ##################################################################



void point_velocity::get_point_v_los_profile(
        const struct point_4cellavg *pt, ///< point to add to profile.
        double *prof, ///< Array of velocity bins to put profile into.
        const int ifrac ///< index of i-fraction in P.V.
        )
{
  double vel, norm;
  int ibin;
  vel  = get_point_los_velocity(pt);
  norm = get_point_neutralH_numberdensity(pt,ifrac);
  ibin = get_velocity_bin_number(vel);
  for (int v=0;v<v_Nbins;v++) prof[v] = 0.0;
  if (ibin>0 && ibin<v_Nbins)
    prof[ibin] = norm;

  if (broaden==2) {
    broaden_profile(pt,prof,ifrac,vel,norm);
  }
  
  return;
}



// ##################################################################
// ##################################################################



void point_velocity::broaden_profile(const struct point_4cellavg *pt, ///< point to add to profile.
				     double *prof,    ///< velocity bins.
				     const int ifrac, ///< index of i-fraction in P.V.
				     double vel,      ///< LOS velocity of point.
				     double norm      ///< Normalisation of profile.
				     )
{
  //
  // This should be a delta-function input profile, and we will convolve this with
  // (i.e. multiply by) a Gaussian to get the doppler broadening.
  //
  if (fabs(vel)<TINYVALUE && norm<TINYVALUE)
    return;

  //cout <<"\n********** vel="<<vel<<" norm="<<norm<<endl;
  //print_array("profile, pre-broadening ",prof,v_Nbins);

  //
  // T  = gas temperature.
  // ma = mass of molecule/atom we are modelling.
  // sigma = the thermal broadening sigma in the gaussian.
  //
  double T = get_point_temperature(pt,ifrac);
  //double ma = 27.0*pconst.m_p(); // mass of 12CO molecule
  //double sigma2 = pconst.kB()*T/ma;
  //double sigma  = sqrt(sigma2);

  //
  // Here I've done the multiplications already: the prefactor is kB/2/27m_p
  // The factor of 2 is in the Gaussian function.
  // The factor of 27 is for Carbon Monoxide CO, being its atomic mass.
  //
  double sigma2 = 1.530e6*T;         // this includes the factor of 2 in the Gaussian exponential.
  double sigma  = sqrt(M_PI*sigma2); // this includes the sqrt(2pi) in the Gaussian normalisation.
  double imin, imid, imax;
  double temp;

  for (int v=0;v<v_Nbins;v++) {
    //
    // Get mean value in bin by using Simpson's 3-pt rule.
    //
    temp = 0.0;
    imin = v_min +v*v_binsize;
    imax = imin+v_binsize;
    imid = 0.5*(imin+imax);

    //temp +=     exp(-0.5*(imin-vel)*(imin-vel)/sigma2);
    //temp += 4.0*exp(-0.5*(imid-vel)*(imid-vel)/sigma2);
    //temp +=     exp(-0.5*(imax-vel)*(imax-vel)/sigma2);
    //temp /= 6.0;
    //temp /= sqrt(2.0*M_PI)*sigma;
    temp +=     exp(-(imin-vel)*(imin-vel)/sigma2);
    temp += 4.0*exp(-(imid-vel)*(imid-vel)/sigma2);
    temp +=     exp(-(imax-vel)*(imax-vel)/sigma2);
    temp /= 6.0;
    temp /= sigma;

    prof[v]  = norm*temp*v_binsize; // normalise correctly by multiplying by interval width.
  }
  //print_array("profile, post-broadening",prof,v_Nbins);

  //for (int v=0;v<v_Nbins;v++) {
    //
    // Add in temperature-dependent normalisation, which physically comes
    // from the emissivity of a certain line due to relative level-populations
    // at different temperatures
    //
    //prof[v] *= 1.0;
  //}
  return;
}



// ##################################################################
// ##################################################################



double point_velocity::get_point_temperature(
        const struct point_4cellavg *pt,
        const int ifrac
        )
{
  double val = 0.0;
  //
  // If microphysics is set up, then use MP->Temperature() to get the
  // temperature.  Otherwise assume a pure Hydrogen gas.
  //
  if (MP) {
    for (int v=0;v<4;v++) {
      if (pt->ngb[v]) {
        val += pt->wt[v] *MP->Temperature(pt->ngb[v]->P,SimPM.gamma);
      }
    }
  }
  else {
    rep.error("get_point_temperature(): no microphysics class",1);
    //  
    // First get the mean of p/rho/(1+x)
    //
    for (int v=0;v<4;v++) {
      if (pt->ngb[v]) {
        //cout <<"p,ro,if = "<<pt->ngb[v]->P[PG]<<", "<< pt->ngb[v]->P[RO] <<", "<<pt->ngb[v]->P[ifrac]<<endl;
        val += pt->wt[v] *(pt->ngb[v]->P[PG]/pt->ngb[v]->P[RO]/(1.0+pt->ngb[v]->P[ifrac]));
      }
    }
    //
    // multiply by m_p/k_B = 1.67e-24/1.38e-16 = 1.21e-8
    //
    val *= 1.21e-8;
    //  cout <<"Temperature="<<val<<endl;
  }
  return val;
}
    


// ##################################################################
// ##################################################################



double point_velocity::get_point_density(const struct point_4cellavg *pt)
{
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *pt->ngb[v]->P[RO];
  }
  return val;
}



// ##################################################################
// ##################################################################



double point_velocity::get_point_neutralH_numberdensity(
        const struct point_4cellavg *pt,
        const int ifrac
        )
{
  double val=0.0;
  //
  // First get neutral mass density.
  //
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      val += pt->wt[v] *(pt->ngb[v]->P[RO]*(1.0-pt->ngb[v]->P[ifrac]));
    }
  }
  //
  // Then scale by H mass fraction and divide by mass of H.
  //
  val *= SimPM.EP.H_MassFrac/pconst.m_p();
  return val;
}



// ##################################################################
// ##################################################################


///
/// Get the absorption and emission coefficients for H-alpha
/// recombination radiation, according to Hummer94 and Henney et al.
/// (2005)'s generic formulae.
/// 
/// I found emissivities according to Storey & Hummer
/// (1995,MNRAS,272,41), and they drop linearly with temperature.
/// Not sure why that is, but it is for Ha,Hb,Hg, etc.
///
/// UPDATE: Replaced with a fit to Ostebrock (1989)'s tables, which
/// seem to be more reliable.
///
void point_velocity::get_point_Halpha_params(
        const struct point_4cellavg *pt, ///< point in question.
        const int ifrac, ///< index of Prim.Vector with Ion. fraction.
        double *alpha,   ///< absorption coefficient (photons/cm)
        double *j        ///< emission coeff (phot/cm^3/s/ster)
        )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn;
  T  = get_point_temperature(pt,ifrac);
  //
  // get_point_neutralH_numberdensity() assume mean mass per particle
  // is m_p(), so we multiply by the H mass fraction to get H number
  // density.
  // get_point_density() is mass density, so we also divide by m_p().
  // Then the electron number density is the H+ number density times
  // (1+X(He)/4X(H)).
  //
  nn = get_point_neutralH_numberdensity(pt,ifrac);
  ni = get_point_density(pt)*SimPM.EP.H_MassFrac/pconst.m_p();
  ni = std::max(0.0, ni-nn);
  //
  // First absorption, from Henney et al. (2009) assuming the opacity
  // is from dust, so that neutrals and ions both count.
  //
  *alpha = (ni+nn) *5.0e-22; // cgs units hardcoded.
  //*alpha = 0.0; // zero absorption (cuts out simulation edges).
  //
  // Emissivity in H-alpha is
  //   j = 1.12e-22*n_e*n_p/pow(T,0.9) erg/cm3/s/sr,
  // from Osterbrock (book, edition 2006, table 4.4).
  // Converted to per square arcsec this is
  //   j = 2.63e-33*n_e*n_p/pow(T,0.9) erg/cm3/s/sq.arcsec.
  // Assume n_e=n_p (i.e. ignore electrons from Helium).
  //
  if (T<1.0) {
    // can get zero temperature if point is off-grid!!!
    *j = 0.0;
  }
  else {
    *j = 2.63e-33*ni*ni*exp(-0.9*log(T));
    // OLD VALUES *j = 2.056e-14*ni*ni/T;
    // Emissivity scales more steeply than recombination rate
  }
  return;
}



// ##################################################################
// ##################################################################



///
/// Get the absorption and emission coefficients for [N II] 6584AA
/// forbidden line emission, according to the formula in Dopita 
/// (1973,A&A,29,387)
/// 
void point_velocity::get_point_NII6584_params(
        const struct point_4cellavg *pt, ///< point in question.
        const int ifrac, ///< index of Prim.Vector with Ion. fraction.
        double *alpha,   ///< absorption coefficient (/cm)
        double *j        ///< emission coeff (erg/cm^3/s/sq.arcsec)
        )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn;
  T  = get_point_temperature(pt,ifrac);
  //
  // get_point_neutralH_numberdensity() assume mean mass per particle
  // is m_p(), so we multiply by the H mass fraction to get H number
  // density.
  // get_point_density() is mass density, so we also divide by m_p().
  // Then the electron number density is the H+ number density times
  // (1+X(He)/4X(H)).
  //
  nn = get_point_neutralH_numberdensity(pt,ifrac);
  ni = get_point_density(pt)*SimPM.EP.H_MassFrac/pconst.m_p();
  ni = std::max(0.0, ni-nn);
  //
  // First absorption, from Henney et al. (2009) assuming the opacity
  // is from dust, so that neutrals and ions both count.
  //
  *alpha = (ni+nn) *5.0e-22; // cgs units hardcoded.
  //*alpha = 0.0; // zero absorption (cuts out simulation edges).
  //
  // Emissivity in [N II] 6584AA is
  //  j([NII] ll 6584) =
  //   6.82e-18 n_e*n_p*f(N)*exp(-chi/kT)/(4*pi*sqrt(T))
  // in units of erg/cm3/s/sr, and use
  //  n_e*n_p = rho^2 y^2 (X_H/m_p)^2
  // and we convert to erg/cm2/s/sq.arcsec.
  //
  // Furthermore, we assume N has its solar ISM abundance of
  // A(N)=7.85 or f(N)=7.08e-5.
  //
  if (T<1.0) {
    // can get zero temperature if point is off-grid!!!
    *j = 0.0;
  }
  else {
    *j = 9.03e-34*ni*ni*exp(-2.1855e4/T)/sqrt(T);
  }
  return;
}




// ##################################################################
// ##################################################################



double point_velocity::get_point_los_velocity(
        const struct point_4cellavg *pt
        )
{
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *(sx*pt->ngb[v]->P[vx]*st
                                      +sz*pt->ngb[v]->P[vz]*ct);
  }
  //  cout <<"sx="<<sx<<" sz="<<sz;
  //  cout <<" ct="<<costheta<<" st="<<sintheta<<" val="<<val;
  //val /= 1.0e5; // convert to km/s
  //val *=15.0;
  //  cout <<" scaled val="<<val<<endl;
  //if (fabs(val)>20.0) 
  // cout <<"sx="<<sx<<" sz="<<sz<<" ct="<<costheta;
  // cout <<" st="<<sintheta<<" val="<<val<<" old_val="<<val/5.0<<"\n";
  return val;
}



// ##################################################################
// ##################################################################



double point_velocity::get_point_perp_velocity(const struct point_4cellavg *pt
					       )
{
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *(sx*pt->ngb[v]->P[vx]*ct
                                      -sz*pt->ngb[v]->P[vz]*st);
  }
  //cout <<"sx="<<sx<<" sz="<<sz<<" ct="<<ct<<" st="<<st<<" val="<<val;
  //val /= 1.0e5; // convert to km/s
  //val *=5.0;
  // cout <<" scaled val="<<val<<endl;
  return val;
}



// ##################################################################
// ##################################################################



double point_velocity::get_point_VX(const struct point_4cellavg *pt
				    )
{
  //
  // This does the bi-linear interpolation from the four nearest
  // cells in the plane, to the point in the centre of the pixel on
  // the line of sight.  The weights of each point are already
  // calculated in the setup of the image with its lines of sight.
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      val += pt->wt[v] *pt->ngb[v]->P[VX];
      //cout <<"vx="<<sx*pt->ngb[v]->P[VX]<<endl;
    }
  }
  return val;
}



// ##################################################################
// ##################################################################



int point_velocity::get_velocity_bin_number(
        const double vel ///< point's velocity
        )
{
  int ibin = static_cast<int>((vel - v_min)/v_binsize);
  if      (vel<v_min) {
    //cout <<"vel<v_min: "<<vel<<endl;
    ibin=-1;
  }
  else if (vel>v_max) {
    //cout <<"vel>v_max: "<<vel<<endl;
    ibin=v_Nbins;
  }
  else {
    if (vel < v_min +    ibin*v_binsize) {
      //cout <<"vel="<<vel<<" INcreasing bin.\n";
      //cout <<"ibin="<<ibin;
      //cout.width(20);
      //cout <<" v_min="<<v_min<<" v_max="<<v_max<<"\tvm+i*sze="<<v_min+ibin*v_binsize;
      //cout.width(20);
      //cout <<"\t(v-vm)/sze="<<(vel - v_min)/v_binsize<<endl;
      ibin--;
    }
    if (vel >= v_min +(ibin+1)*v_binsize) {
      //cout <<"vel="<<vel<<" DEcreasing bin.\n";
      //cout.width(20);
      //cout <<"v_min="<<v_min<<" vmax="<<v_max<<" binsize="<<v_binsize<<endl;
      ibin++;
    }
    if ((vel < v_min+ibin*v_binsize) ||
        (vel >= v_min+(ibin+1)*v_binsize)) {
      rep.error("Can't locate velocity in a velocity bin.",vel);
    }
  }
  return ibin;
}



// ##################################################################
// ##################################################################



// ------------------------------------------------------------
// ************************************************************
// ------------------------------------------------------------
//
// Dummy program, to test the coordinate transforms.
//
// int main(int argc, char **argv)
// {
//   //
//   // First initialise MPI, even though this is a single processor
//   // piece of code.
//   //
//   int err = COMM->init(&argc, &argv);

//   int angle=atoi(argv[1]);
//   //
//   // Normal vector to trace rays along (and measure angle from (if non-zero)).
//   //
//   string n=argv[2]; enum direction normal=NO;
//   if      (n=="XN") normal=XN;
//   else if (n=="XP") normal=XP;
//   else if (n=="YN") normal=YN;
//   else if (n=="YP") normal=YP;
//   else if (n=="ZN") normal=ZN;
//   else if (n=="ZP") normal=ZP;
//   else rep.error("Bad normal direction",n);
//   //
//   // Perpendicular direction, around which to rotate the view.
//   //
//   n=argv[3]; enum direction perpdir=NO;
//   if      (n=="XN") perpdir=XN;
//   else if (n=="XP") perpdir=XP;
//   else if (n=="YN") perpdir=YN;
//   else if (n=="YP") perpdir=YP;
//   else if (n=="ZN") perpdir=ZN;
//   else if (n=="ZP") perpdir=ZP;
//   else rep.error("Bad perp direction",n);

//   //*******************************************************************
//   //*******************************************************************
//   SimPM.ndim=3;
//   SimPM.Ncell=1;
//   for (int v=0;v<3;v++){
//     SimPM.Xmin[v] = 0.0;
//     SimPM.NG[v] = 4*v+4;
//     SimPM.Ncell *= SimPM.NG[v];
//     SimPM.Xmax[v] = static_cast<double>(SimPM.NG[v])/2.0; // dx=0.5
//     SimPM.Range[v] = SimPM.Xmax[v]-SimPM.Xmin[v];
//   }
//   SimPM.nvar=5;
//   SimPM.eqnNDim=3;
//   SimPM.solverType=2;
//   SimPM.coord_sys=COORD_CRT;
//   SimPM.eqntype=EQEUL;
//   SimPM.gridType=1;
//   //
//   // Set low-memory cells
//   //
//   CI.set_minimal_cell_data();

//   //
//   if (grid) rep.error("grid already setup, so bugging out",grid);
//   try {
//     grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
//   }
//   catch (std::bad_alloc) {
//     rep.error("(trunks::setup_grid) Couldn't assign data!", grid);
//   }
//   cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;

//   cout <<"-----------------------------------------------------------------\n";  
//   cout <<"*****************************************************************\n";
//   cout <<"-----------------------------------------------------------------\n";  
//   //class image *cc = new image (normal,angle,perpdir,grid);
//   class image cc(normal,angle,perpdir,grid);

//   int npix[2]; npix[0]=npix[1]=-1;
//   //cc->get_npix(npix);
//   cc.get_npix(npix);
//   cout <<"npix = ["<<npix[0]<<", "<<npix[1]<<"]\n"; 

//   //cc->set_cell_positions_in_image();
//   //cc->add_cells_to_pixels();
//   cc.set_cell_positions_in_image();
//   cc.add_cells_to_pixels();
//   cc.add_integration_pts_to_pixels();

//   cell *c = grid->FirstPt();
//   do {
//     cout <<c->id<<"\t\t"; rep.printVec("cell-pos",c->pos,3);
//     cout <<"\t\t"; rep.printVec("img-pos ",c->Ph,3);
//   } while ((c=grid->NextPt(c))!=0);

//   cout <<"-----------------------------------------------------------------\n";  
//   cout <<"*****************************************************************\n";
//   cout <<"-----------------------------------------------------------------\n";  
//   //  if (cc) {
//   //  delete cc; cc=0;
//   //}
//   cc.delete_cell_positions();

//   if(grid!=0) {
//     cout << "\t Deleting Grid Data..." << endl;
//     delete grid; grid=0;
//   }
//   COMM->finalise();
//   delete COMM; COMM=0;

//   return 0;
// }


// ##################################################################
// ##################################################################



