#include "fourierCPP.h"
#include <cmath>
#include <iostream>
#include <string>

#ifdef USE_FFTW2

FOUR_FFTW2_ONE::FOUR_FFTW2_ONE(std::string prep, int nn[2]) : in_(NULL), out_(NULL) 
{
  plan_ = NULL;
  size_of_data_ = nn[0]*nn[1];
  in_ = (fftw_complex*) new double[size_of_data_*2];
  fftw_direction direction;
  if (prep == std::string("for2")) 
    {
      direction = FFTW_FORWARD;
    }
  else if (prep == std::string("back2")) 
    {
      direction = FFTW_BACKWARD;
    }
  else
    {
      std::cerr << " FOUR_FFTW2_ONE:: ERROR unknown type of fft preparation = " << prep << std::endl;
      exit(0);
    }
  plan_ = fftw2d_create_plan_specific(nn[1],nn[0],direction, FFTW_ESTIMATE | FFTW_IN_PLACE, in_,1,out_,1);
}
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FOUR_FFTW2_MANY::FOUR_FFTW2_MANY(std::string prep, int nn[2]) 
{
  plan_a_ = NULL;
  plan_b_ = NULL;
  size_of_data_ = nn[0]*nn[1];
  in_ = (fftw_complex*) new double[size_of_data_*2];
      	
  if (prep == std::string("for3"))
    {
      howmanya_ = nn[0]/2;
      howmanyb_ = nn[1];

      istridea_ = nn[0];
      istrideb_ = 1;

      idista_ = 1;
      idistb_ = nn[0];

      plan_a_ = fftw_create_plan_specific(nn[1],FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE,in_,nn[0],out_,1);
      plan_b_ = fftw_create_plan_specific(nn[0],FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE,in_,1,out_,1);
    }
  else if (prep ==  std::string("back3"))
    {
      howmanya_ = nn[1];
      howmanyb_ = nn[0]/2;

      istridea_ = 1;
      istrideb_ = nn[0];

      idista_ = nn[0];
      idistb_ = 1;

      plan_a_ = fftw_create_plan_specific(nn[0],FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE,in_,1,out_,1);
      plan_b_ = fftw_create_plan_specific(nn[1],FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE,in_,nn[0],out_,1);
    }
  else
    {
      std::cerr << " FOUR_FFTW2_MANY:: ERROR unknown type of fft preparation = " << prep << std::endl;
      exit(0);
    }
}
#endif

#ifdef USE_FFTW3

FOUR_FFTW3_ONE::FOUR_FFTW3_ONE(std::string prep, int nn[2]) 
{
  plan_ = NULL;
  
  size_of_data_ = nn[0]*nn[1];
  in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*size_of_data_);
  out_ = in_;
  
   
  if (prep == std::string("for2")) 
    {
      plan_ = fftw_plan_dft_2d(nn[1],nn[0],in_, out_, FFTW_FORWARD, FFTW_ESTIMATE);
    }
  else if (prep == std::string("back2")) 
    {
      plan_ = fftw_plan_dft_2d(nn[1],nn[0],in_, out_, FFTW_BACKWARD,FFTW_ESTIMATE);
    }
  else
    {
      std::cerr << " FOUR_FFTW3:: ERROR unknown type of fft preparation = " << prep << std::endl;
      exit(0);
    }  
}

FOUR_FFTW3_MANY::FOUR_FFTW3_MANY(std::string prep, int nn[2])
{
  plan_a_ = NULL;
  plan_a_ = NULL;
  
  size_of_data_ = nn[0]*nn[1];
  in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*size_of_data_);
  out_ = in_;

  if (prep == std::string("for3"))
    {
      plan_a_ = fftw_plan_many_dft( 1, &nn[1],nn[0]/2,in_, NULL,nn[0],1, out_, NULL,nn[0],1, FFTW_FORWARD,FFTW_ESTIMATE);
      plan_b_ = fftw_plan_many_dft( 1, &nn[0],nn[1],in_, NULL,1,nn[0], out_, NULL,1,nn[0], FFTW_FORWARD,FFTW_ESTIMATE);
    }
  else if (prep ==  std::string("back3"))
    {
      plan_a_ = fftw_plan_many_dft( 1, &nn[0],nn[1],in_, NULL,1, nn[0], out_, NULL,1,nn[0], FFTW_BACKWARD,FFTW_ESTIMATE);

      plan_b_ =  fftw_plan_many_dft( 1, &nn[1],nn[0]/2,in_, NULL,nn[0],1, out_, NULL,nn[0],1, FFTW_BACKWARD,FFTW_ESTIMATE);
    }
  else
    {
      std::cerr << " FOUR_FFTW3_MANY:: ERROR unknown type of fft preparation = " << prep << std::endl;
      exit(0);
    }
}


#endif


#ifdef USE_FFT_LOCAL


LOCAL_FOURIER::LOCAL_FOURIER(std::string prep,int nn[2]) 
    {
      int k;
      size_of_data_ = nn[0]*nn[1];
      in_ = new double[size_of_data_*2];
      for (k=0; k<size_of_data_*2; k++) in_[k] = 0.0;
      if (prep == std::string("for2") || prep == std::string("for3")) direction_ = 1;
      else 
	if (prep == std::string("back2") || prep == std::string("back3")) direction_ = -1;
	else 
	  {
	    std::cerr << " LOCAL_FOURIER:: ERROR unknown type of fft preparation = " << prep << std::endl;
	    exit(0);
	  }

      nn_[0] = nn[0];
      nn_[1] = nn[1];
    }


void LOCAL_FOURIER::fourtrans (double* data,int nn[],int isign)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4, i__5, i__6;
  double d__1;
  
  /* Builtin functions 
    double sin(); */
  
  /* Local variables */
  int idim, ibit, nrem, ntot, i2rev, i3rev, n;
  double theta, tempi, tempr;
  int nprev, i2;
  double wtemp;
  int i1, i3, k1, k2;
  double wi, wr;
  int ip1, ip2, ip3;
  double wpi, wpr;
  int ifp1, ifp2;
  int ndim=2;
  /* Parameter adjustments */
  --nn;
  --data;
  
  /* Function Body */
  ntot = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim)
    {
      ntot *= nn[idim];
    }
  nprev = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1;
    i__2 = ip2;
    i__3 = ip1;
    for (i2 = 1; i__3 < 0 ? i2 >= i__2 : i2 <= i__2; i2 += i__3) {
      if (i2 < i2rev) {
	i__4 = i2 + ip1 - 2;
	for (i1 = i2; i1 <= i__4; i1 += 2) {
	  i__5 = ip3;
	  i__6 = ip2;
	  for (i3 = i1; i__6 < 0 ? i3 >= i__5 : i3 <= i__5; i3 += 
	       i__6) {
	    i3rev = i2rev + i3 - i2;
	    tempr = data[i3];
	    tempi = data[i3 + 1];
			data[i3] = data[i3rev];
	    data[i3 + 1] = data[i3rev + 1];
	    data[i3rev] = tempr;
	    data[i3rev + 1] = tempi;
	  }
	}
      }
      ibit = ip2 / 2;
    L1:
      if (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit /= 2;
	goto L1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
  L2:
    if (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign * 6.28318530717959 / (ifp2 / ip1);
      /* Computing 2nd power */
      d__1 = sin(theta * .5);
      wpr = d__1 * d__1 * -2.;
      wpi = sin(theta);
      wr = 1.;
      wi = 0.;
      i__3 = ifp1;
      i__2 = ip1;
      for (i3 = 1; i__2 < 0 ? i3 >= i__3 : i3 <= i__3; i3 += i__2) {
	i__4 = i3 + ip1 - 2;
	for (i1 = i3; i1 <= i__4; i1 += 2) {
	  i__6 = ip3;
		    i__5 = ifp2;
	  for (i2 = i1; i__5 < 0 ? i2 >= i__6 : i2 <= i__6; i2 += 
	       i__5) {
	    k1 = i2;
	    k2 = k1 + ifp1;
	    tempr = wr * data[k2] - wi * data[k2 + 1];
	    tempi = wr * data[k2 + 1] + wi * data[k2];
	    data[k2] = data[k1] - tempr;
	    data[k2 + 1] = data[k1 + 1] - tempi;
	    data[k1] += tempr;
	    data[k1 + 1] += tempi;
	  }
	}
	wtemp = wr;
	wr = wr * wpr - wi * wpi + wr;
	wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
      goto L2;
    }
    nprev = n * nprev;
  }
  return;
}
#endif
