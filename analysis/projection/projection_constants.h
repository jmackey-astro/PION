
#ifndef PROJECTION_CONSTANTS_H
#define PROJECTION_CONSTANTS_H

#define PROJ_D     0  ///< projected density
#define PROJ_NtD   1  ///< projected neutral density
#define PROJ_InD   2  ///< projected ionised density
#define PROJ_EM    3  ///< projected emission measure.
#define PROJ_X00p1 4  ///< X-ray emission >0.1 keV
#define PROJ_X00p2 5  ///< X-ray emission >0.2 keV
#define PROJ_X00p3 6  ///< X-ray emission >0.3 keV
#define PROJ_X00p5 7  ///< X-ray emission >0.5 keV
#define PROJ_X01p0 8  ///< X-ray emission >1.0 keV
#define PROJ_X02p0 9  ///< X-ray emission >2.0 keV
#define PROJ_X05p0 10  ///< X-ray emission >5.0 keV
#define PROJ_X10p0 11 /// < X-ray emission >10.0 keV
#define PROJ_HA    12  ///< projected H-alpha emission.
#define PROJ_NII   13 ///< projected ionised metal-line.

#define OP_TEXT 0
#define OP_FITS 1
#define OP_SILO 2
#define OP_VTK  3

#define NIMG     14 ///< total number of images.
#define N_SCALAR 12  ///< number of scalar images.

//
// data read from file is put in a 2D array with these values for
// the first index.
//
#define DATA_R   0 ///< radius
#define DATA_D   1 ///< density
#define DATA_P   2 ///< pressure
#define DATA_V   3 ///< velocity
#define DATA_TR0 4 ///< tracer 0
#define DATA_TR1 5 ///< tracer 1
#define DATA_T   6 ///< Temperature


#endif // PROJECTION_CONSTANTS_H

