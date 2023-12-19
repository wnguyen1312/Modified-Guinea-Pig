#ifndef FIELD_CPP_SEEN
#define FIELD_CPP_SEEN
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "typeDefs.h"
#include "define.h"
#include "mathconst.h"
#include "fourierCPP.h"
#include "typeDefs.h"

class FIELD
{

  int nb_cells_x_, nb_cells_y_;
  int dist_size_;
  int integration_method_;
  PHI_FLOAT *dist_;

  PHI_FLOAT *phi1_;
  PHI_FLOAT *phi2_;

  ABSTRACT_FOURIER* fourier_transform_forward_;
  ABSTRACT_FOURIER* fourier_transform_backward_;


  inline void init_fourier_tools(FFT_SERVER* fourier)
  {
    int nn[2];
    //    int k;
    if (fourier == NULL) 
      {
	std::cerr << " FIELD::init_fourier_tools : you must give a fourier server with integration_method = 2 " << std::endl;
	exit(0); 
      }
    //int size = nb_cells_x_*nb_cells_y_*8;
    //  temp_ = new PHI_FLOAT[size];
    //  for (k=0; k<size; k++) temp_[k] = 0.0;
    nn[0]=2*nb_cells_y_;
    nn[1]=2*nb_cells_x_;
    fourier_transform_forward_ = fourier->new_fft(std::string("for3"),nn);
    fourier_transform_backward_ = fourier->new_fft(std::string("back3"),nn);
  }

    void set (int nx, int ny, int integration) 
    {
      //int k;
      int nbnodes = nx*ny;
      nb_cells_x_ = nx;
      nb_cells_y_ = ny;
      if (nbnodes >0) 
	{
	  phi1_=new PHI_FLOAT[nbnodes];
	  phi2_=new PHI_FLOAT[nbnodes];
	}
      else 
	{
	  phi1_ = NULL;
	  phi2_ = NULL;
	}
      integration_method_ = integration;
      if (integration_method_ == 2)
	{
	  dist_size_ = 8*nx*ny;
	}
      else
	{
	  if ( integration_method_ <= 0) dist_size_ =0;
	  else dist_size_ = nx*ny;
	}
      if (dist_size_ > 0) 
	{
	  dist_ = new PHI_FLOAT[dist_size_];
	}
      else dist_ = NULL;
    }

 public :


    FIELD() 
      { 
	set (0,0,0);
	fourier_transform_forward_ = NULL;
	fourier_transform_backward_ = NULL;
      }

    FIELD(int integration_method, int nx, int ny, FFT_SERVER* fourier=NULL)
      {
	set (nx,ny, integration_method);
	if (integration_method == 2) init_fourier_tools(fourier);
	else 
	  {
	    fourier_transform_forward_ = NULL;
	    fourier_transform_backward_ = NULL;
	  }
      }

    FIELD(const FIELD& f)
      {
	int k;
	set (f.nb_cells_x_, f.nb_cells_y_, f.integration_method_); 
	for (k = 0; k < dist_size_; k++) dist_[k] = (f.dist_)[k];
	if (integration_method_ == 2)
	  {
	    fourier_transform_forward_ = f.fourier_transform_forward_;
	    fourier_transform_backward_ = f.fourier_transform_backward_;
	  }
	else 
	  {
	    fourier_transform_forward_ = NULL;
	    fourier_transform_backward_ = NULL;
	  }
      }

    FIELD( FIELD& f)
      {
	int k;
	set (f.nb_cells_x_, f.nb_cells_y_, f.integration_method_); 
	for (k = 0; k < dist_size_; k++) dist_[k] = (f.dist_)[k];
	if (integration_method_ == 2)
	  {
	    fourier_transform_forward_ = f.fourier_transform_forward_;
	    fourier_transform_backward_ = f.fourier_transform_backward_;
	  }
	else 
	  {
	    fourier_transform_forward_ = NULL;
	    fourier_transform_backward_ = NULL;
	  }
      }

  ~FIELD() 
    {
      delete [] dist_;
      delete [] phi1_;
      delete [] phi2_;
    }

  inline FIELD& operator = (const FIELD& f)
    {
      if (this == &f) return *this; // protect against self-assignment
      int k;
      // delete old memory:
      delete [] dist_;
      delete [] phi1_;
      delete [] phi2_;
      // assign new memory:
      set (f.nb_cells_x_, f.nb_cells_y_, f.integration_method_); 
      if (integration_method_ == 2)
	{
	  fourier_transform_forward_ = f.fourier_transform_forward_;
	  fourier_transform_backward_ = f.fourier_transform_backward_;
	}
      for (k = 0; k < dist_size_; k++) dist_[k] = f.dist_[k];
      return *this;
    }

  inline const PHI_FLOAT* get_phi(int index) const 
  {
    if (index == 1) return phi1_;
    else
      if (index == 2) return phi2_;
      else 
	{
	  std::cerr << " FIELD::get_phi() : ERROR index = " << index << std::endl;
	  exit(0);
	}
  }

  void fold_fft(const PHI_FLOAT *rho1,const PHI_FLOAT *rho2,const PHI_FLOAT *dist,PHI_FLOAT *phi1,
		PHI_FLOAT *phi2,int n_x,int n_y, float charge_sign,   ABSTRACT_FOURIER* fourier_forward, ABSTRACT_FOURIER* fourier_backward);

  void foldfronts (const PHI_FLOAT *rho,const PHI_FLOAT *dist,PHI_FLOAT *phi,int n_x,int n_y, float charge_sign);
  
  void sor2 (const PHI_FLOAT *rho,PHI_FLOAT *phi,int n_x,int n_y,PHI_FLOAT *parameter, float charge_sign);
  

  /*! This routine is a subroutine for init_grid. 

    calculate the function F(x,y) of the thesis p. 20 (integral of the 2D Green function of the Poisson equation
  */
  // thesis p. 20 : F(x,y) function multiplied by 4.pi.eps0
  inline double f_potential(double x,double y)
  {
    return x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  }

  double f_potential_2(double x0,double y0,double dx,double dy);

  void foldfields (const PHI_FLOAT *rho, const PHI_FLOAT *dist,PHI_FLOAT *phi,int n_x,int n_y, float charge_sign);
  
  void dist_init(PHI_FLOAT factor, float deltax, float deltay,   FFT_SERVER* fourier);

  inline void foldfields (const PHI_FLOAT *rho1,const PHI_FLOAT *rho2,float charge_sign)
  {
    foldfields(rho1,dist_,phi1_, nb_cells_x_, nb_cells_y_, charge_sign);
    foldfields(rho2,dist_,phi2_, nb_cells_x_, nb_cells_y_, charge_sign);
  }

  inline void fold_fft(const PHI_FLOAT *rho1,const PHI_FLOAT *rho2,float charge_sign)
  {
    fold_fft(rho1, rho2, dist_, phi1_, phi2_, nb_cells_x_, nb_cells_y_, charge_sign, fourier_transform_forward_, fourier_transform_backward_);
  }

  inline  void foldfronts (const PHI_FLOAT *rho1,const PHI_FLOAT *rho2, PHI_FLOAT *parameter, float charge_sign)
  {
    foldfronts(rho1, dist_, phi1_, nb_cells_x_, nb_cells_y_, charge_sign);
    foldfronts(rho2, dist_, phi2_, nb_cells_x_, nb_cells_y_, charge_sign);
    sor2(rho1, phi1_, nb_cells_x_, nb_cells_y_,parameter, charge_sign);
    sor2(rho2, phi2_, nb_cells_x_, nb_cells_y_,parameter, charge_sign);
  }
};

#endif
