
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//
// This is the NR version of the FFT.  Instead of trying to convert it
// all for zero-offset arrays, I copied Martin White, and access the
// (n-1)th element of data[] whenever I need to.
//
void four1(double *data,
	   unsigned long int nn,       // half the length of the array.
	   int isign
	   )
{

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

void fourn(double data[], int nn[], int ndim, int isign)
/*
   Does the N-dimensional complex transform of a function.  The array
   layout for real and imaginary parts is e.g. (2D)
   data[2*(nn[1]*i+j)  ] = real part
   data[2*(nn[1]*i+j)+1] = imag part
   By MW, 2001.
*/
{
int		idim;
int		i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
int		ibit,k1,k2,n,nprev,nrem,ntot;
double		tempi,tempr;
double		theta,wi,wpi,wpr,wr,wtemp;

for (ntot=1, idim=1; idim<=ndim; idim++)
  ntot *= nn[idim-1];
nprev=1;
for (idim=ndim;idim>=1;idim--) {
  n=nn[idim-1];
  nrem=ntot/(n*nprev);
  ip1= 2*nprev;
  ip2= ip1*n;
  ip3= ip2*nrem;
  i2rev=1;
  for (i2=1; i2<=ip2; i2+=ip1) {
    if (i2 < i2rev) {
      for (i1=i2; i1<=i2+ip1-2; i1+=2) {
        for (i3=i1; i3<=ip3; i3+=ip2) {
          i3rev=i2rev+i3-i2;
          tempr= data[i3-1];
          tempi= data[i3  ];
          data[i3-1] = data[i3rev-1];
          data[i3  ] = data[i3rev  ];
          data[i3rev-1] = tempr;
          data[i3rev  ] = tempi;
        }
      }
    }
    ibit=ip2/2;
    while (ibit >= ip1 && i2rev > ibit) {
      i2rev -= ibit;
      ibit = ibit/2;
    }
    i2rev += ibit;
  }
  ifp1=ip1;
  while (ifp1 < ip2) {
    ifp2 = 2*ifp1;
    theta=isign*6.28318530717959/(ifp2/ip1);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (i3=1;i3<=ifp1;i3+=ip1) {
      for (i1=i3;i1<=i3+ip1-2;i1+=2) {
        for (i2=i1;i2<=ip3;i2+=ifp2) {
          k1=i2;
          k2=k1+ifp1;
          tempr=wr*data[k2-1]-wi*data[k2  ];
          tempi=wr*data[k2  ]+wi*data[k2-1];
          data[k2-1]  = data[k1-1]-tempr;
          data[k2  ]  = data[k1  ]-tempi;
          data[k1-1] += tempr;
          data[k1  ] += tempi;
        }
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    ifp1=ifp2;
  }
  nprev *= n;
  }
}

void print_array(string name, double *arr, int nel)
{
  cout <<name<<": [";
  for (int v=0;v<nel-1;v++) cout <<arr[v]<<", ";
  cout <<arr[nel-1]<<"]\n";
  return;
}


#define FWHM 0.424661

int main(int argc, char **argv)
{

  if  (argc != 3) {
    cout <<"usage: smooth1d <n-elements> <smooth-FWHM(pix)>\n";
    return 1;
  }

  int nel = atoi(argv[1]);
  int asz=1;
  while (asz<nel) asz *=2;
  //
  // now asz>nel, so multiply by 2 again for zero padding.
  //
  asz*=2;
  cout <<"nel="<<nel<<" and asz="<<asz<<endl;

  double sigma = atof(argv[2])*FWHM/asz;

  //
  // Put data in the array.
  //
  double *data = new double [nel];
  for (int v=0;v<nel;v++) {
    data[v] = 2.0;
    if (v>nel/5   && v<2*nel/5) data[v] *=2.5;
    if (v>4*nel/6 && v<5*nel/6) data[v] *=4.0;
  }
  print_array("original", data,nel);

  //
  // Remove the mean
  //
  double mean=0.0;
  for (int v=0;v<nel;v++)
    mean += data[v];
  mean /= nel;
  for (int v=0;v<nel;v++)
    data[v] -= mean;

  //
  // Put data into Cx. array, with zero padding.
  //

  double *cxdata = new double [2*asz];
  for (int v=0;v<2*asz;v++)
    cxdata[v]=0.0;
  //
  // The following bit of code works for the test problem.
  //
//   for (int v=0;v<asz;v++) {
//     cxdata[2*v+1]=0.0;
//     cxdata[2*v  ]=2.0-mean;
//   }
  
  for (int v=0;v<nel;v++)
    cxdata[2*v] = data[v];
  print_array("ORIG",cxdata,2*asz);

  //
  // FT the data
  //
  four1(cxdata, static_cast<unsigned long>(asz), 1);
  //fourn(cxdata, &asz, 1, 1);
  //print_array("  FT",cxdata,2*asz);


  //
  // Multiply by FT of gaussian.
  // The FT is in frequency units f=k/(2PI)
  // for each f, the FT of the Gaussian is exp(-2*PI*PI*sigma*sigma*f*f)
  // *BUT* remember the FFT returns frequencies in the range -N/2,N/2, so
  // I need to do the numbering carefully.
  // Martin does it in smoothmap.c, where if i>Ng/2 f=Ng-i; else f=i;
  // Because we are squaring, the sign doesn't matter.
  int f;
  double expfactor;
  for (int v=0; v<asz; v++) {
    if (v>asz/2) f = asz-v;
    else         f = v;
    expfactor = 2*M_PI*M_PI*sigma*sigma*f*f;
    //cout <<"expfactor="<<expfactor<<endl;
    expfactor = exp(-expfactor);
    //
    // We divide by asz because the IFT returns asz times the function.
    //
    cxdata[2*v+0] *= expfactor/asz;
    cxdata[2*v+1] *= expfactor/asz;
  }
  
  
  //
  // IFT the data
  //
  four1(cxdata, static_cast<unsigned long>(asz), -1);
  //fourn(cxdata, &asz, 1, -1);
  //  for (int v=0;v<2*asz;v++) cxdata[v] /= asz;
  //print_array(" IFT",cxdata,asz);

  //
  // Move from Cx. array into original array, and add back mean.
  //
  double *data2 = new double [nel];
  for (int v=0;v<nel;v++)
    data2[v] = cxdata[2*v] +mean;

  //
  // output the data.
  //
  print_array("original?",data2,nel);
  ofstream outf("temp.txt");
  for (int v=0;v<nel;v++)
    outf <<v<<"\t"<< data[v]+mean <<"\t"<< data2[v] <<endl;
  outf.close();

  delete [] data;
  delete [] data2;
  delete [] cxdata;
  return 0;
}
