#ifndef FOURIER_SEEN
#define FOURIER_SEEN
//#include "typeDefs.h"
#include <cstdio>   //   exit function 
#include <cstdlib>  //   exit function
#include <iostream>
#include <string>
#include <vector>

#include "config.h"  // have fftw/sfftw/dfftw

class ABSTRACT_FOURIER
{

 protected:
  int size_of_data_;

 public:
 ABSTRACT_FOURIER()  : size_of_data_(0) {;}

  virtual  ~ABSTRACT_FOURIER() {;}


  virtual void make() = 0;

  virtual double* data_vector() = 0;

};

#ifdef USE_FFT_LOCAL
class LOCAL_FOURIER : public ABSTRACT_FOURIER
{

 protected : 

  double* in_;
  int direction_;
  int nn_[2];

 LOCAL_FOURIER() : in_(NULL) {;}

 private :
  /// assignment and copy constructor not implemented nor used
  LOCAL_FOURIER& operator=(const LOCAL_FOURIER&);
  LOCAL_FOURIER(LOCAL_FOURIER&);

  void fourtrans (double* data,int nn[],int isign);

 public :

  
  LOCAL_FOURIER(std::string prep,int nn[2]); 

  virtual  ~LOCAL_FOURIER() { delete [] in_;}

  virtual inline  double* data_vector() {   return in_;}

  virtual inline void make() {  fourtrans(in_,nn_,direction_);}

};
#endif

#ifdef USE_FFTW2
extern "C"
{
#ifdef HAVE_FFTW_H
#include "fftw.h"
#elif defined(HAVE_DFFTW_H)
#include "dfftw.h"
#elif defined(HAVE_SFFTW_H)
#include "sfftw.h"
#else // This is for old ./configure
#include "fftw.h"
#endif
}

class FOUR_FFTW2_ONE : public ABSTRACT_FOURIER
{

 protected : 

  fftwnd_plan plan_;

  fftw_complex* in_;
  fftw_complex* out_;

 FOUR_FFTW2_ONE() : in_(NULL), out_(NULL) {plan_ = NULL;}


 public :

  FOUR_FFTW2_ONE(std::string prep, int nn[2]);

  virtual  ~FOUR_FFTW2_ONE() 
    {
      if (in_ != NULL) delete [] in_;
      if (plan_ != NULL) { fftwnd_destroy_plan(plan_);}
    }

  virtual inline  double* data_vector() { return (double*)in_;}


  virtual inline void make() { fftwnd_one(plan_,in_,out_);}

};

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class FOUR_FFTW2_MANY : public ABSTRACT_FOURIER
{

 protected : 

  fftw_plan plan_a_,plan_b_;

  fftw_complex* in_;
  fftw_complex* out_;

  int howmanya_, howmanyb_;
  int istridea_, istrideb_;
  int idista_, idistb_;

 FOUR_FFTW2_MANY() : in_(NULL), out_(NULL) 
    {
      plan_a_ = NULL;
      plan_b_ = NULL;
    }


 public :

  FOUR_FFTW2_MANY(std::string prep, int nn[2]);




  virtual  ~FOUR_FFTW2_MANY() 
    {
      if (in_ != NULL) delete [] in_;
      if (plan_a_ != NULL) {fftw_destroy_plan(plan_a_);}
      if (plan_b_ != NULL) {fftw_destroy_plan(plan_b_);}
    }

  virtual inline  double* data_vector() { return (double*)in_;}


  virtual inline void make()
  {
    fftw(plan_a_,howmanya_,in_,istridea_,idista_,NULL,1,1);
    fftw(plan_b_,howmanyb_,in_,istrideb_,idistb_,NULL,1,1);
  }

};
#endif

#ifdef USE_FFTW3

extern "C"
{
#include "fftw3.h"
}



class FOUR_FFTW3_ONE : public ABSTRACT_FOURIER
{

 protected : 

  fftw_plan plan_;
  fftw_complex *in_;
  fftw_complex *out_;

 FOUR_FFTW3_ONE() : in_(NULL), out_(NULL) 
    {
      plan_ = NULL;
    }


 public :



  FOUR_FFTW3_ONE(std::string prep, int nn[2]);

  virtual  ~FOUR_FFTW3_ONE() 
    {
      fftw_free(in_);
      if (plan_ != NULL) { fftw_destroy_plan(plan_);}
    }

  virtual inline  double* data_vector() { return (double*)in_;}

  virtual inline void make() { fftw_execute(plan_);}
};
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class FOUR_FFTW3_MANY : public ABSTRACT_FOURIER
{

 protected : 

  fftw_plan plan_a_;
  fftw_plan plan_b_;
  fftw_complex *in_;
  fftw_complex *out_;

 FOUR_FFTW3_MANY() : in_(NULL), out_(NULL) 
    {
      plan_a_ = NULL;
      plan_b_ = NULL;
    }


 public :



  FOUR_FFTW3_MANY(std::string prep, int nn[2]);

  virtual  ~FOUR_FFTW3_MANY() 
    {
      fftw_free(in_);
      if (plan_a_ != NULL) { fftw_destroy_plan(plan_a_);}
      if (plan_b_ != NULL) {fftw_destroy_plan(plan_b_);}
    }

  virtual inline  double* data_vector() { return (double*)in_;}



  virtual inline void make()
  {
    fftw_execute(plan_a_);
    fftw_execute(plan_b_);
  }
};
#endif


class FFT_SERVER
{

 private : 

  std::vector<ABSTRACT_FOURIER*> pointer_fourier_;

 public : 

  FFT_SERVER() {;}

  ~FFT_SERVER()
    {
      unsigned int k;
      for (k=0; k < pointer_fourier_.size(); k++) delete pointer_fourier_[k];
    }

  inline ABSTRACT_FOURIER* new_fft(std::string prep, int nn[2])
  {
  
#ifdef USE_FFT_LOCAL
    pointer_fourier_.push_back(new LOCAL_FOURIER(prep,nn));
#else

#ifdef USE_FFTW2
    if (prep == std::string("for2") || prep == std::string("back2"))
      {
	pointer_fourier_.push_back(new FOUR_FFTW2_ONE(prep, nn));
      }
    else
      {
	pointer_fourier_.push_back(new FOUR_FFTW2_MANY(prep, nn));
      }
#else  

#ifdef USE_FFTW3
    if (prep == std::string("for2") || prep == std::string("back2"))
      {
	pointer_fourier_.push_back(new FOUR_FFTW3_ONE(prep, nn));
      }
    else
      {
	pointer_fourier_.push_back(new FOUR_FFTW3_MANY(prep, nn));
      }
#else
    std::cerr << " GUINEA:: error : a Fourier Transform must be assigned " << std::endl;
    exit(0);
#endif
#endif
#endif
    return pointer_fourier_.back();

  }


};

#endif
