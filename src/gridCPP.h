#ifndef GRID_SEEN
#define GRID_SEEN
#include <cmath>
#include <cstdio>
#include <list>
#include <vector>
#include "typeDefs.h"
#include "meshCPP.h"
#include "beamCPP.h"
#include "extraPhotonCPP.h"
#include "switchesCPP.h"
#include "resultsCPP.h"
#include "pairsCPP.h"
#include "lumiCPP.h"
#include "jetParameterCPP.h"
#include "backgroundCPP.h"
#include "rndmCPP.h"
#include "bhabhaSamplesCPP.h"
#include "abstractParticle.h"
#include "fourierCPP.h"
#include "fieldCPP.h"
#include "parametersCPP.h"
#include "mathematicalEntities.h"
#include "EDMwriter.h"

//! virtual class 
/*! 
 attributes of GENERAL_GRID are arrays for field computations.


 attributes of derived class GRID are arrays giving for each cell pointers to particles and photons which are affected to this cell. (as a result of 'distribute')

 inheriting class EXTRA_GRID has the same attributes as GENERAL_GRID
*/ 

class SLICE_ON_GRID
{

  float nb_part_per_macro_;
  float scal_step_;
  int nb_macroparticles_;
  BEAM* beam_;
  PHI_FLOAT *rho_;
  int dim_rho_;
  int nb_cells_y_;

  // not implemented
  SLICE_ON_GRID(SLICE_ON_GRID& s);

  inline void set_rho(int dim)
    {
      if (rho_ != NULL) delete [] rho_;
      if (dim > 0) 
	{
	  rho_ = new PHI_FLOAT[dim];
	}
      else
	{
	  rho_ = NULL;
	}
    }

  inline  void set_size(int nx, int ny)
    {
      nb_cells_y_ = ny;
      dim_rho_ = nx*ny;
      set_rho(dim_rho_);
    }
  
  
  inline void setZero()
    {
      nb_part_per_macro_ = 0.0;
      scal_step_ = 1.0;
      nb_macroparticles_ = 0;
      beam_ = NULL;
      rho_ = NULL;
      dim_rho_ = 0;
      nb_cells_y_ = 0;
    }
  
  public :
    
    SLICE_ON_GRID() 
    {
      setZero();
    }
  
  SLICE_ON_GRID(int nx, int ny)
    {
      setZero();
      set_size(nx,ny);
    } 
  
  SLICE_ON_GRID(const SLICE_ON_GRID& s) 
    {
      setZero();
      beam_ = s.beam_;
      nb_part_per_macro_ =  s.nb_part_per_macro_;
      scal_step_ = s.scal_step_ ;
      nb_macroparticles_ = s.nb_macroparticles_;
      nb_cells_y_ = s.nb_cells_y_;
      dim_rho_ = s.dim_rho_;
      set_rho( dim_rho_);
      /*       std::cout << " SLICE_ON_GRID : copy constructor " << std::endl; */
      /*       std::cout << " old address " << s.rho_ << " new address " << rho_ << std::endl; */
      /*       std::cout << " nb_macroparticles_ " <<  nb_macroparticles_ << std::endl; */
      /*       std::cout << " nb_part_per_macro_ " << nb_part_per_macro_ << " scal_step_ " << scal_step_ << " nb_cells_y_ " << nb_cells_y_ << " dim_rho_ " << dim_rho_ << std::endl; */
      /*       std::cout << " *********************************************** " << std::endl; */
    }
  
  ~SLICE_ON_GRID()
    {
      delete [] rho_;
    }
  
  inline SLICE_ON_GRID& operator = (const SLICE_ON_GRID& s)
    {
      if (this == &s) return *this; // protect against self-assignment
      int k;
      beam_ = s.beam_;
      nb_part_per_macro_ =  s.nb_part_per_macro_;
      scal_step_ = s.scal_step_ ;
      nb_macroparticles_ = s.nb_macroparticles_;
      nb_cells_y_ = s.nb_cells_y_;
      dim_rho_ = s.dim_rho_;
      // delete old memory:
      delete[] rho_;
      rho_ = NULL;
      // assign new memory:
      set_rho( dim_rho_);
      for (k=0; k < dim_rho_ ; k++) rho_[k] = s.rho_[k];
      return *this;
    }
  
  inline void connect_beam(BEAM* b) { beam_ = b;}
  
  
  inline void resize(int nx, int ny)
  {
    nb_cells_y_ = ny;
    dim_rho_ = nx*ny;
    set_rho( dim_rho_);
  }

  inline void set_macroparticles(int mac)
    {
      nb_macroparticles_ = mac;
    }
  
  inline int get_macroparticles() const 
    {
      return nb_macroparticles_;
    }
  
  inline const BEAM* get_beam() const 
    {
      return beam_;
    }
  
  inline BEAM* get_beam() 
    {
      return beam_;
    }
  
  inline void set_scal_step(float s)
    {
      scal_step_ = s;
    }
  
  inline float get_scal_step() const 
    {
      return scal_step_;
    }
  
  inline void update_nb_part_per_macro(float nb_macro_part)
    {
      if (nb_macroparticles_ > 0)
	{
	  nb_part_per_macro_ = nb_macro_part/(float)nb_macroparticles_;
	}
      else nb_part_per_macro_ = 0.0;
    }
  
  inline float get_nb_part_per_macro() const 
    {
      return nb_part_per_macro_;
    }
  
  inline const PHI_FLOAT* get_rho() const { return rho_;}
  
  inline void assignChargeToNGP(int i1,int i2, float chgSgn) 
    {  
      rho_[i1*nb_cells_y_+i2] +=  nb_part_per_macro_*chgSgn;
    }
  
  inline void assignChargeToCIC(int i1,int i2, float h_x, float h_y,  float chgSgn) 
    { 
      int j;
      PHI_FLOAT charge = nb_part_per_macro_*chgSgn;
      j=i1*nb_cells_y_+i2;
      rho_[j] += (1.0-h_x)*(1.0-h_y)*charge;
      j=(i1+1)*nb_cells_y_+i2;
      rho_[j] += h_x*(1.0-h_y)*charge;
      j=i1*nb_cells_y_+i2+1;
      rho_[j] += (1.0-h_x)*h_y*charge;
      j=(i1+1)*nb_cells_y_+i2+1;
      rho_[j] += h_x*h_y*charge;
    };
  
  inline void razRho()
    {
      int k;
      for (k=0; k < dim_rho_ ; k++)
	{
	  rho_[k]=0.0;	
	}
    }
  
  inline void symmetrizeCharges(int n_cell_x, int n_cell_y)
    {
      PHI_FLOAT tmp;
      int i1, i2;
      int j1, j2, j3, j4;
      for (i1=0;i1<n_cell_x/2;i1++){
	for (i2=0;i2<n_cell_y/2;i2++)
	  {
	    j1=i1*n_cell_y + i2;
	    j2=(n_cell_x - 1-i1)*n_cell_y +i2;
	    j3=i1*n_cell_y +(n_cell_y -1-i2);
	    j4=(n_cell_x-1-i1)*n_cell_y +(n_cell_y -1-i2);
	    tmp=0.25*(rho_[j1]+rho_[j2]+rho_[j3]+rho_[j4]);
	    rho_[j1]=tmp;
	    rho_[j2]=tmp;
	    rho_[j3]=tmp;
	    rho_[j4]=tmp;
	    
	    /* 	      tmp=0.25*(rho2_[j1]+rho2_[j2]+rho2_[j3]+rho2_[j4]); */
	    /* 	      rho2_[j1]=tmp; */
	    /* 	      rho2_[j2]=tmp; */
	    /* 	      rho2_[j3]=tmp; */
	    /* 	      rho2_[j4]=tmp; */
	  }
      }
    };
};

class GENERAL_GRID
{
  
 protected: 
  
  /*! number of cells in three coordinates */
  int n_cell_x_,n_cell_y_,n_cell_z_; 
  int timestep_;
  MESH mesh_;
  float min_x_,max_x_,min_y_,max_y_,min_z_,max_z_;
  float cut_x_,cut_y_,cut_z_;
  //  float nb_part_per_macro_[2];
  float step_;
  //    float scal_step[2];
  float rho_x_1_,rho_y_1_,rho_sum_1_,rho_x_2_,rho_y_2_,rho_sum_2_;

  float delta_x_inverse_,delta_y_inverse_;
  
  PHI_FLOAT rho_factor_;
  
  // int nb_macroparticles_[2];

  int integration_method_;

  //  PHI_FLOAT *dist_;

  FIELD champ_;

  SLICE_ON_GRID slice_of_beam_[2];

  //  BEAM* beam1_;
  // PHI_FLOAT *rho1_;

  //  BEAM* beam2_;
  //  PHI_FLOAT *rho2_;

  GENERAL_GRID();

  virtual void distribute_coherent_particles(int i_slice1,
					  int i_slice2,
				    int electron_distribution_rho) = 0;

  public :

  virtual ~GENERAL_GRID();

  inline void connect_beams(BEAM* b1, BEAM* b2)
    {
      slice_of_beam_[0].connect_beam(b1);
      slice_of_beam_[1].connect_beam(b2);
      /*   beam1_ = b1; */
      /*   beam2_ = b2; */
    } 
  
  inline int get_n_cell_x() const {return n_cell_x_;};
  inline int get_n_cell_y() const {return n_cell_y_;};
  inline int get_n_cell_z() const {return n_cell_z_;};
  
  // no_beam = 1,2
  inline int get_nb_macroparticles(int nobeam) const 
    {
      return  slice_of_beam_[nobeam-1].get_macroparticles();
      // return nb_macroparticles_[nobeam-1];
    }
  
  inline int get_timestep() const {return timestep_;};
  
  inline float get_min_x() const {return min_x_;};
  inline float get_max_x() const {return max_x_;};
  inline float get_min_y() const {return min_y_;};
  inline float get_max_y() const {return max_y_;};
  inline float get_min_z() const {return min_z_;};
  inline float get_max_z() const {return max_z_;};
  inline float get_delta_x() const {return mesh_.get_delta_x();};
  inline float get_delta_y() const {return mesh_.get_delta_y();};
  inline float get_delta_z() const {return mesh_.get_delta_z();};
  inline float get_scal_step(int index) const
    {
      return slice_of_beam_[index].get_scal_step();
      //return scal_step[index];
    }
  //  inline float get_cut_z() const {return cut_z_;};
  inline float get_step() const {return step_;};
  inline float get_offset_x() const {return mesh_.get_offset_x();};
  inline float get_offset_y() const {return mesh_.get_offset_y();};
  
  inline const MESH* get_mesh() const {return &mesh_;}
  
  inline float get_rho_x_1() const {return rho_x_1_;};
  inline float get_rho_y_1() const {return rho_y_1_;};
  inline float get_rho_sum_1() const {return rho_sum_1_;};
  inline float get_rho_sum_2() const {return rho_sum_2_;};
  inline float get_rho_x_2() const {return rho_x_2_;};
  inline float get_rho_y_2() const {return rho_y_2_;};
  
  
  inline PHI_FLOAT get_rho_factor() const {return rho_factor_;};
  
  inline const PHI_FLOAT* get_phi1() const {return champ_.get_phi(1);};
  inline const PHI_FLOAT* get_phi2() const {return champ_.get_phi(2);};

  //  void get_cuts(float& x,float& y,float&z) const {x = cut_x_; y = cut_y_; z = cut_z_;};

  // nobeam = 1,2
  inline void set_nb_macroparticles(int nobeam, int npart)
    { 
      slice_of_beam_[nobeam-1].set_macroparticles(npart);
      //nb_macroparticles_[nobeam-1] = npart;
    }
  
  /* inline void update_nb_part_per_macro(int nobeam, float n_particles) */
  /* { */
  /*   int index = nobeam -1; */
  /*   if (nb_macroparticles_[index] > 0) */
  /*     { */
  /*       nb_part_per_macro_[index] = n_particles/(float)nb_macroparticles_[index]; */
  /*     } */
  /* } */
  /* inline void razRhos() */
  /* { */
  /*   int i1, i2, j; */
  /*   for (i1=0;i1<n_cell_x_;i1++) */
  /*     { */
  /*       for (i2=0;i2<n_cell_y_;i2++) */
  /* 	{ */
  /* 	  j=i1*n_cell_y_+i2; */
  /* 	  rho1_[j]=0.0; */
  /* 	  rho2_[j]=0.0; */
  /* 	} */
  /*     } */
  /* }; */
  
  inline void absoluteCoordinates(float x, float y, PHI_FLOAT& h_x, PHI_FLOAT& h_y) const
    {
      h_x=x*delta_x_inverse_ + mesh_.get_offset_x();
      h_y=y*delta_y_inverse_ + mesh_.get_offset_y();
    }
  
  
  // given coordinates x,y : return the cell coordinates i1,i2 of the cell
  // containing (x,y) ; ai is the biggest integer that is lower than x
  // bi is the biggest integer lower than y-1/2; weight = y -1/2 -i2
  inline void cellLocalization(float x, float y, int& ai, int& bi, float& weight) const
    {
      ai=(int)floor(x);
      weight = y-0.5;
      bi=(int)floor(weight);
      weight -= bi;
    };
  
  inline void cellLocalization(PHI_FLOAT x, PHI_FLOAT y, int& ai, int& bi, PHI_FLOAT& weight) const
    {
      ai=(int)floor(x);
      weight = y-0.5;
      bi=(int)floor(weight);
      weight -= bi;
    };
  
  inline int slice_from_z(float z)
    {
      return (int)floor(z/get_delta_z()+0.5*(float)n_cell_z_);
    }
  
  inline void phiValuesX(int i1, int i2,register float h, const PHI_FLOAT * phi, PHI_FLOAT& phi1_x, PHI_FLOAT& phi2_x, PHI_FLOAT& phi3_x) const
    {
      int j;
      register float h_p;
      // 
      //
      h_p=1.0-h;
      j=i1*n_cell_y_+i2;
      phi1_x=h*phi[j+n_cell_y_+1]+h_p*phi[j+1];
      phi2_x=h*phi[j+n_cell_y_] + h_p*phi[j];
      phi3_x=h*phi[j+n_cell_y_-1] + h_p*phi[j-1];
    }

  inline void phiValuesY(int i1, int i2,register float h, const PHI_FLOAT * phi, PHI_FLOAT& phi1_y, PHI_FLOAT& phi2_y, PHI_FLOAT& phi3_y) const
    {
      int j;
      register float h_p;
      
      h_p=1.0-h;
      j=i1*n_cell_y_+i2;
      phi1_y=h*phi[j+n_cell_y_+1] + h_p*phi[j+n_cell_y_];
      phi2_y=h*phi[j+1] + h_p*phi[j];
      phi3_y=h*phi[j-n_cell_y_+1] + h_p*phi[j-n_cell_y_];
    }
  
  inline bool particleInGrid(float xpart, float ypart, float& h_x,float& h_y, int& i1, int& i2)
    { 
      h_x=(xpart*delta_x_inverse_ + mesh_.get_offset_x()-0.5);
      i1=(int)floor(h_x);
      h_x -= (float)i1;
      h_y=(ypart*delta_y_inverse_ +  mesh_.get_offset_y()-0.5);
      i2=(int)floor(h_y);
      h_y -= (float)i2;
      if ((i1>=0)&&(i1<n_cell_x_-1)&&(i2>=0)&&(i2<n_cell_y_-1)) return true;
      else return false;
    }
  
  inline bool particleInGrid(float xpart, float ypart, int& i1, int& i2)
    { 
      i1=(int)floor(xpart*delta_x_inverse_+ mesh_.get_offset_x());
      i2=(int)floor(ypart*delta_y_inverse_+ mesh_.get_offset_y());
      if ((i1>=0)&&(i1<n_cell_x_-1)&&(i2>=0)&&(i2<n_cell_y_-1)) return true;
      else return false;
    }
  
  inline void assignChargeToNGP(PHI_FLOAT* rho, int i1,int i2, PHI_FLOAT charge) {  rho[i1*n_cell_y_+i2] += charge;};
  
  inline void assignChargeToCIC(PHI_FLOAT* rho, int i1,int i2, float h_x, float h_y, PHI_FLOAT charge) 
    { 
      int j;
      j=i1*n_cell_y_+i2;
      rho[j] += (1.0-h_x)*(1.0-h_y)*charge;
      j=(i1+1)*n_cell_y_+i2;
      rho[j] += h_x*(1.0-h_y)*charge;
      j=i1*n_cell_y_+i2+1;
      rho[j] += (1.0-h_x)*h_y*charge;
      j=(i1+1)*n_cell_y_+i2+1;
      rho[j] += h_x*h_y*charge;
    };
  
  /* inline void symmetrizeCharges() */
  /* { */
  /*   PHI_FLOAT tmp; */
  /*   int i1, i2; */
  /*   int j1, j2, j3, j4; */
  /*       for (i1=0;i1<n_cell_x_/2;i1++){ */
  /* 	  for (i2=0;i2<n_cell_y_/2;i2++) */
  /* 	    { */
  /* 	      j1=i1*n_cell_y_+i2; */
  /* 	      j2=(n_cell_x_-1-i1)*n_cell_y_+i2; */
  /* 	      j3=i1*n_cell_y_+(n_cell_y_-1-i2); */
  /* 	      j4=(n_cell_x_-1-i1)*n_cell_y_+(n_cell_y_-1-i2); */
  /* 	      tmp=0.25*(rho1_[j1]+rho1_[j2]+rho1_[j3]+rho1_[j4]); */
  /* 	      rho1_[j1]=tmp; */
  /* 	      rho1_[j2]=tmp; */
  /* 	      rho1_[j3]=tmp; */
  /* 	      rho1_[j4]=tmp; */
  /* 	      tmp=0.25*(rho2_[j1]+rho2_[j2]+rho2_[j3]+rho2_[j4]); */
  /* 	      rho2_[j1]=tmp; */
  /* 	      rho2_[j2]=tmp; */
  /* 	      rho2_[j3]=tmp; */
  /* 	      rho2_[j4]=tmp; */
  /* 	  } */
  /*       } */
  /* }; */
  
  inline bool coordinatesInGridRange(float xpart, float ypart) const
    { 
      if ((xpart>min_x_)&&(xpart<max_x_)&&(ypart>min_y_)&&(ypart<max_y_)) return true;
      else return false;
    };
  
  void interpolePotential(float xpart,float ypart, PHI_FLOAT& h_x, PHI_FLOAT& h_y, PHI_FLOAT& phi1_x, PHI_FLOAT& phi2_x, PHI_FLOAT& phi3_x, PHI_FLOAT& phi1_y,PHI_FLOAT& phi2_y, PHI_FLOAT& phi3_y, const PHI_FLOAT *phi) const;
  
  TRIVECTOR ElectricFieldCIC(float xpart,float ypart,const PHI_FLOAT *phi) const;
  
  void computeFields(int integrationMethod, float charge_sign, PHI_FLOAT *sor_parameter); 
  
  virtual void distribute_particles(int i_slice1,
				    int i_slice2,
				    int electron_distribution_rho, 
				    int force_symmetric) = 0;
  
  
  void  init_grid_phys (float n_particles1, float n_particles2,float cut_x,float cut_y,float cut_z, float charge_sign,   FFT_SERVER* fourier);
  
};



class GRID : public GENERAL_GRID , public ABSTRACT_IO_CLASS
{

  ///////////// private struct ///////////////////////////////////////  
  struct DISTRIBUTE
  {
    long int in, out;
    float delta;
    long int tot;
    DISTRIBUTE():in(0),out(0),delta(0.0),tot(0){}
  };
  /////////////////////////////////////////////////////////////////////
  
  ///////////// private class ///////////////////////////////////////
  class PARTICLE_POINTER 
  {
    
  protected :
    float weight_;
    
  public :

    PARTICLE_POINTER() {;}
    PARTICLE_POINTER(float weight) {weight_ = weight;}


    virtual  ~PARTICLE_POINTER() {;}
    inline float weight() const {return weight_;}; 
    virtual   void velocities(float& vx, float& vy) const = 0;
    virtual   const TRIDVECTOR& getSpin() const = 0;
    virtual    float energy() const = 0;
    //   virtual  float q2() const = 0;
    //   virtual  float eorg() const = 0;
  };

  
  class BEAM_PARTICLE_POINTER : public PARTICLE_POINTER
  {
  
  private:
    
    const ABSTRACT_PARTICLE* partic_;
  
  public:
  
    BEAM_PARTICLE_POINTER() {;}
  
    BEAM_PARTICLE_POINTER(const ABSTRACT_PARTICLE* partic, float weight) : PARTICLE_POINTER(weight)
    {
      partic_ =partic;
      //	weight_ = weight;
    }
    
    ~BEAM_PARTICLE_POINTER() {;}
  
    virtual  inline  void velocities(float& vx, float& vy) const {partic_->velocities(vx, vy);}
  
    virtual inline const TRIDVECTOR& getSpin() const 
    { 
      return  partic_->getSpin();
    }
 
    virtual   inline   float energy() const { return partic_->energy();}
    
    //  virtual inline float q2() const { return 0.0;};
    //   virtual inline float eorg() const { return 0.0;};
    
    
    /*  inline float get_spin(int nc) const   */
    /*  {  */
    /*    //  std::cout << " BEAM_PARTICLE_POINTER:: get_spin() " << std::endl; */
    /*    return partic_->get_spin()(nc-1); */
    /*  } */
    
  };
  //////////////////
  
  class PHOTON_POINTER : public PARTICLE_POINTER
  {
    
  protected :
    
  public : 
    
    PHOTON_POINTER() {;}
    PHOTON_POINTER(float weight) : PARTICLE_POINTER(weight) {;}
    virtual ~PHOTON_POINTER() {;}
 
    virtual  float q2() const = 0;
    virtual  float eorg() const = 0;

  };

  class BEAM_PHOTON_POINTER : public PHOTON_POINTER
  {
    
  private:
    
    const PHOTON* partic_;
    
  public:
    
    BEAM_PHOTON_POINTER() {;}
    
    BEAM_PHOTON_POINTER(const PHOTON* partic, float weight) : PHOTON_POINTER(weight)
    {
      partic_ =partic;
    }
    
    virtual ~BEAM_PHOTON_POINTER() {;}
    
    virtual inline void velocities(float& vx, float& vy) const {partic_->velocities(vx, vy);}
    
    virtual inline const TRIDVECTOR& getSpin() const 
    { 
      std::cerr << " a photon has no spin " << std::endl;
      exit(0);
      /*       const TRIDVECTOR tv(0.0,0.0,0.0); */
      /*       return tv; */
    }
    
    virtual inline float energy() const { return partic_->energy();}
  
    virtual inline float q2() const { return 0.0;};
    virtual inline float eorg() const { return -1.0;};
    
    inline float spin(int nc) const  { return partic_->getSpin()(nc-1);}
    
  };
  //////////////////
  
  class EXTRA_PHOTON_POINTER : public PHOTON_POINTER
  {
    EXTRA* extra_phot_;
    
  public :
    
    EXTRA_PHOTON_POINTER() : extra_phot_(NULL) {;}
    ~EXTRA_PHOTON_POINTER()
      {
	//  std::cout << " destructor EXTRA_PHOTON_POINTER() pointer = "<< extra_phot_ << endl;
	delete extra_phot_;
      }
    
    EXTRA_PHOTON_POINTER& operator=(const EXTRA_PHOTON_POINTER& pt)
      {
	if (this == &pt) return *this; // Gracefully handle self assignment
	
	// call base class assignment operator
	PHOTON_POINTER::operator=(pt);
	
	float energy,vx,vy, q2, eorg;
	pt.extra_phot_->get_parameters(energy,vx,vy,q2,eorg);
	EXTRA* tmp = new EXTRA(energy, vx, vy, q2, eorg);
	delete extra_phot_;
	extra_phot_ = tmp;
	
	return *this;
      }

    EXTRA_PHOTON_POINTER(const EXTRA_PHOTON_POINTER& pt) : PHOTON_POINTER(pt)
    {
      float energy,vx,vy, q2, eorg;
      pt.extra_phot_->get_parameters(energy,vx,vy,q2,eorg);
      extra_phot_ = new EXTRA(energy, vx, vy,q2,eorg);
    }
    
    EXTRA_PHOTON_POINTER(float energy,float vx,float vy,float q2,float eorg,float weight) : PHOTON_POINTER(weight)
    {
      extra_phot_ = new EXTRA(energy, vx, vy,q2,eorg);
    }
    virtual inline float energy() const {return extra_phot_->energy();}
    virtual inline float q2() const {return extra_phot_->q2();}
    virtual inline float eorg() const {return extra_phot_->eorg();}
    virtual inline void velocities(float& vx, float& vy) const { extra_phot_->velocities(vx, vy);}
    virtual inline const TRIDVECTOR& getSpin() const 
    { 
      std::cerr << " no spin for extra photon " << std::endl;
      exit(0);
      /*       const TRIDVECTOR tv(0.0,0.0,0.0); */
      /*       return tv; */
    }
    
  };
  //////////////////////////////////////////////////////////////////////
  
  
  //! for each cell, a list of pointers to particles assigned to this celle (beam1)
  std::vector< std::list<BEAM_PARTICLE_POINTER> > part_pointer1_;
  //! for each cell, a list of pointers to particles assigned to this celle (beam2)
  std::vector< std::list<BEAM_PARTICLE_POINTER> > part_pointer2_;
  
  //! for each cell, a list of pointers to photons assigned to this celle (beam1)
  std::vector< std::list<BEAM_PHOTON_POINTER> > grid_photon1_;
  //! for each cell, a list of pointers to photons assigned to this celle (beam2)
  std::vector< std::list<BEAM_PHOTON_POINTER> > grid_photon2_;
  
  //! for each cell, a list of pointers to extra_photons assigned to this celle (beam2)
  std::vector< std::list<EXTRA_PHOTON_POINTER> > extra_photon_pointer1_;
  //! for each cell, a list of pointers to extra_photons assigned to this celle (beam2)
  std::vector< std::list<EXTRA_PHOTON_POINTER> > extra_photon_pointer2_;
  
  // A std::vector of the photons generated by incoherent pairs. They do not iteract with the beam.
  std::vector<TERTPHOTON> tertphot_;
  
  RNDM* rndm_generator_;
  
  
  RESULTS results_;
  COHERENT_RESULTS coherent_results_;
  // PHOTON_RESULTS photon_results_;
  TRIDENT_RESULTS trident_results_;
  
  
  COMPT compt_;
  
  BHABHA bhabhas_;
  
  
  FILE_IN_OUT* photon_file_;
  FILE_IN_OUT* hadron_file_;
  
  LUMI_HEAP_EE lumi_heap_ee_;
  LUMI_HEAP lumi_heap_eg_;
  LUMI_HEAP lumi_heap_ge_;
  LUMI_HEAP lumi_heap_gg_;
  
  
  DISTRIBUTE distribute1_;
  DISTRIBUTE distribute2_;
  
  ABSTRACT_MINIJETS* minijets_;
  
  GENERAL_CROSS* cross_;
  
  /* inline void connect_secondaries(PAIR_BEAM* secondaries) */
  /* { */
  /*   secondaries_ = secondaries; */
  /* } */
  
  void  lumi_init( const SWITCHES& switches);
  
  //void deltaVelocityFromFieldCIC(float xpart,float ypart, float energy,   PHI_FLOAT *phi, float distance,float& ax, float& ay);
  
  inline float spread_energy(float ener, int which_spread, float spread, RNDM& rndm_generator) const 
  {
    JET_FLOAT dx=3.4,xmin=-1.7;
    float energy = ener;
    switch(which_spread)
      {
      case 0:
	break;
      case 1:
	energy *= 1.0 + spread*(rndm_generator.rndm()*dx+xmin);
	break;
      case 2:
	if ( rndm_generator.rndm() < 0.5 )
	  {
	    energy *= 1.0 + spread;
	  }
	else
	  {
	    energy *= 1.0 - spread;
	  }
	break;
      case 3:
	energy *= 1.0 + spread*rndm_generator.gasdev();
	break;
      }
    return energy;
  }
  
void collide_ee(int cellx, int celly,float min_z, const BEAM_PARTICLE_POINTER& pointer1, const BEAM_PARTICLE_POINTER& pointer2, SWITCHES& switches,PAIR_BEAM& secondaries,int time_counter, int beamslice1, int beamslice2) ;
  
  void collide_ge(int cellx, int celly,float min_z, const BEAM_PHOTON_POINTER& photon_pointer, const BEAM_PARTICLE_POINTER& particle_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& rndm_generator);
  
  void collide_eg(int cellx, int celly,float min_z,const BEAM_PARTICLE_POINTER& particle_pointer, const BEAM_PHOTON_POINTER& photon_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& rndm_generator);

  void collide_gg(int cellx, int celly, float min_z, const BEAM_PHOTON_POINTER& photon_pointer1, const BEAM_PHOTON_POINTER& photon_pointer2, PAIR_BEAM& secondaries, PAIR_BEAM& muons, SWITCHES& switches, RNDM& rndm_generator, double *returnCrossSectionToAdd);

  // with dir = 1 and photon_q2=0., it is the old collide_ge_2
  // with dir = 2 and photon_q2=0., it is the old collide_eg_2
  // with dir = 1 and photon_q2 != 0., it is the old collide_ge_3
  // with dir = 2 and photon_q2!= 0., it is the old collide_eg_3

  // argument  index_of_process is old 'scdn1+scdn2'
  // 0 : Breit-Wheeler
  // 1 : Bethe-Heitler
  // 2 : Landau-Lifschitz

  inline void collide_compton(int index_of_process, int cellx, int celly,float min_z, const PHOTON_POINTER& extr_phot, const BEAM_PARTICLE_POINTER& particle,PAIR_BEAM& secondaries,SWITCHES& switches, RNDM& rndm_generator, int dir)
  {
    if (switches.get_do_compt())
      {
	float particleEnergy = particle.energy();
	float particleVx,particleVy;
	particle.velocities(particleVx, particleVy);

	float photonEnergy = extr_phot.energy();
	float photon_q2 =  extr_phot.q2(); 

	float weight = extr_phot.weight()*particle.weight();

	compt_.compt_do(mesh_, cellx, celly,min_z, secondaries,index_of_process,particleEnergy,photonEnergy,photon_q2,particleVx, particleVy,weight,dir,switches, rndm_generator);
      }
  }

  void collide_gg_XX(int index_of_process, int cellx, int celly,float min_z, const PHOTON_POINTER& photon1, const PHOTON_POINTER& photon2, SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator, double *returnCrossSectionToAdd);

  void make_hadrons_gg2(float min_z, float energy1,float q2_1,float energy2,float q2_2,float lumi, SWITCHES& switches,double *returnCrossSectionToAdd, RNDM& rndm_generator);

  void isr2(float e1,float e2,float *e1p,float *e2p, RNDM& rndm_generator);

  //void interpolePotential(float xpart,float ypart, PHI_FLOAT& h_x, PHI_FLOAT& h_y, PHI_FLOAT& phi1_x, PHI_FLOAT& phi2_x, PHI_FLOAT& phi3_x, PHI_FLOAT& phi1_y,PHI_FLOAT& phi2_y, PHI_FLOAT& phi3_y, PHI_FLOAT *phi);

  inline void store_coherent_particle(BEAM& beam,float energy,float vx,float vy, float x,float y,int slice)
  {
    beam.newCoherent(slice, x,y,vx,vy,energy);
    coherent_results_.addCountEnergy(energy);
  }

  inline void store_trident_particle(BEAM& beam,float energy,float vx,float vy, float x,float y,int slice)
  {
    beam.newTrident(slice, x,y,vx,vy,energy);
    //  std::cout << "GRID::create trident : vel "<< energy << " "  << x << " " << y << " " << vx << " " << vy << std::endl;	    
    trident_results_.addCountEnergy(energy);
  }

  int coherent_generate(BEAM& beam,int i_slice, float upsilon,PHOTON& phot, RNDM& rndm_generator);

  bool pick_coherent_energy(float ups,float& energy, RNDM& rndm_generator);

  void apply_electric_field_on_pair(float step_2, float ex, float ey, float& vx,float& vy,float e_inv) const
  {
    // apply electric field in one step : step_2
    // e_inv : inverse of energy

    vx += step_2*ex*e_inv;
    vy += step_2*ey*e_inv;
  }

  /*
    void scale_pair(float vold2, float& vx, float& vy, float& vz, float& vx0, float& vy0, float& vz0,float& eng, float& e_inv, float& e_inv2 ) const
    {
    #ifdef SCALE_ENERGY
    float scal;
    scal=sqrt((vold2*eng*eng + EMASS2)/((vx*vx+vy*vy+vz*vz)*eng*eng + EMASS2));
    vx*=scal;
    vy*=scal;
    vz*=scal;
    #ifdef PAIR_SYN
    vx0 *= scal;
    vy0 *= scal;
    vz0 *= scal;
    #endif
    eng/=scal;
    e_inv=1.0/eng;
    e_inv2=e_inv*e_inv;
    #endif
    }
  */

  void apply_magnetic_field_on_pair(float fac_theta,float step_q, float e_inv2, float bx, float by, float& vx,float& vy, float& vz,float& theta) const;

  int field_pair(const PAIR_PARTICLE& pair, const std::vector<GENERAL_GRID*>& grids,float& ex,float& ey, float& bx,float& by, int extra_grids, float charge_sign_0);

  //void step_pair_1X(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair,float step,int n, int extra_grids, float charge_sign_0, RNDM& rndm_generator);
  void step_pair_1(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair,double mass,float step,int n, int extra_grids, float charge_sign_0, RNDM& rndm_generator);
  void step_pair_1_tertphot(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair,double mass,float step,int n, int extra_grids, float charge_sign_0, RNDM& rndm_generator);
  // float  doElossParticle(PARTICLE& particle, int sokolov, TRIDVECTOR Efiekd, TRIDVECTOR Bfield, BEAM& beam, float initialEnergyForLoss, int i_slice,float radius_i, float dz, float emin, int do_prod, float charg_sign);

  //void computePhotonEnergies(std::vector<float>& photonEnergies, PARTICLE& particle, float  upsilon, int sokolov, TRIDVECTOR Efield, TRIDVECTOR Bfield,  float energy, float radius_i,  float dz, float charge_sign);

  void registerPhotons(const std::vector<float>& photonEnergies, PARTICLE& particle, int i_beam, int i_slice);
  // void referenceSpin(PARTICLE& part, TRIDVECTOR& e1, TRIDVECTOR& e2, TRIDVECTOR& e3, TRIDVECTOR Efield, TRIDVECTOR Bfield, float charge_sign);
  // void distributeElectronScatter1(int i_beam, int i_slice,float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);
  void distributeScatter1(const std::vector<PARTICLE*>& theParticles,float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);
  //void distributeElectronScatter2(int i_beam, int i_slice,float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);
  void distributeScatter2(const std::vector<PARTICLE*>& theParticles,float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);

  void electronScatter( std::vector< std::list<EXTRA_PHOTON_POINTER> >& extra_phot_ptr, int i_beam, int i_slice, float xmin, double s4,double lns4, float ratio, int ext_field, int geom, float ratio_i, float r_scal);

  //void distributePhotonInBeam(   std::vector< std::list<BEAM_PARTICLE_POINTER> >& grid_photon, BEAM* beam, int slice,float ratio, float ratio_i, RNDM& rndm_generator);
  void distributePhotonInBeam(   std::vector< std::list<BEAM_PHOTON_POINTER> >& grid_photon, int i_beam, int slice,float ratio, float ratio_i, RNDM& rndm_generator);

  //! Distributes the virtual photons 
  void distribute_virtual_photons(int i_slice1,
				  int i_slice2,
				  SWITCHES* switches,
				  double s4, double lns4);

  virtual void distribute_coherent_particles(int i_slice1,
					     int i_slice2,
					     int electron_distribution_rho);

  virtual void distribute_trident_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho);
  //! Distributes the particles for background calculation 
  //void distributeCoherentScatterInBeam1(int i_beam, int i_slice, float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);
  // void distributeCoherentScatterInBeam2(int i_beam, int i_slice, float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer);
  void distribute_particles_for_background(int i_slice1,
					   int i_slice2,  
					   float electron_ratio);

  //!  Distributes the particles from coherent pair creation for the background calculation

  void distribute_coherent_particles_for_background(int i_slice1,
						    int i_slice2, 
						    float electron_ratio);

  void distribute_trident_particles_for_background(int i_slice1,
						   int i_slice2, 
						   float electron_ratio);

  TRIVECTOR electric_field_out_of_main_grid(const std::vector<GENERAL_GRID*>& grids,int beam,PHI_FLOAT x,PHI_FLOAT y, int extra_grids) const;

  //int field_coherent(const std::vector<GENERAL_GRID*>& grids,int beam,PHI_FLOAT x,PHI_FLOAT y,float ener, float distance, float& deltaVx,float& deltaVy, int extra_grids);
  //void assignBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro, DISTRIBUTE& distribute);
  void assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  //void assignBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro, DISTRIBUTE& distribute);
  void assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  //void assignCoherentBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro, DISTRIBUTE& distribute);
  void assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  //void assignCoherentBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro, DISTRIBUTE& distribute);
  void assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  void assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  void assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute);

  void move_particles(const std::vector<GENERAL_GRID*>& grids, int i_beam, int i_slice, int interpolation,int do_beamstrahlung,int do_trident, int sokolov, float emin, int do_prod, int extra_grids, float charge_sign,  int bmt_rotate);

  void move_coherent_particles(const std::vector<GENERAL_GRID*>& grids, int i_beam, int i_slice, int interpolation, int do_beamstrahlung,float emin, int do_prod,int extra_grids, float charge_sign);

  void move_trident_particles(const std::vector<GENERAL_GRID*>& grids, int i_beam, int i_slice, int interpolation, int do_beamstrahlung,float emin, int do_prod,int extra_grids, float charge_sign);

  inline TRIVECTOR EBfieldOnParticle(const PARTICLE* particle, const std::vector<GENERAL_GRID*>& grids, const PHI_FLOAT *phi, int i_beam, int extra_grids) const
  {
    float xpart, ypart;
    TRIVECTOR EBfield;
    particle->XYposition(xpart, ypart);
    if ( coordinatesInGridRange(xpart, ypart) )
      {
	// electric field in GV/nm
	EBfield = ElectricFieldCIC(xpart, ypart, phi);
      }
    else
      {
	EBfield = electric_field_out_of_main_grid(grids,i_beam, (PHI_FLOAT)xpart, (PHI_FLOAT)ypart, extra_grids);
      }
    return EBfield;
  }

  //float moveSingleParticle(PARTICLE* particle, TRIVECTOR EBfield, float dzOnRadius,bool rotateSpin, const std::vector<GENERAL_GRID*>& grids, BEAM& beam, int i_beam, int i_slice, PHI_FLOAT *phi, float distance, int do_beamstrahlung, int sokolov, float emin,int do_prod, int extra_grids, float charge_sign);

  //void particleBeamstrahlung(PARTICLE* particle, TRIDVECTOR Efield, TRIDVECTOR Bfield, float dzOnRadius,  float emin, std::vector<float>& photonEnergies);
  //void particleBeamstrahlungSokolov(PARTICLE_WITH_SPIN* particle, TRIDVECTOR Efield, TRIDVECTOR Bfield, float dzOnRadius, float emin, float charge_sign, std::vector<float>& photonEnergies);

  void beamstrahlungSingleCoherentParticle(PARTICLE* particle, TRIVECTOR EBfield, float dzOnRadius, const std::vector<GENERAL_GRID*>& grids, int i_beam, int i_slice, const PHI_FLOAT *phi, float /*distance*/, float emin,int do_prod, int extra_grids, float charge_sign);

  //void advanceCoherentParticlesNGP();
  //void advanceCoherentParticlesCIC(PARTICLE& particle, const std::vector<GENERAL_GRID*>& grids, PHI_FLOAT *phi, int i_beam, float emin, int do_prod, int extra_grids, float scal_step, RNDM& rndm_generator);
  //float requiv(double lns4, float xmin,int iflag);
  //void mequiv (double s4, double lns4, float xmin,float e,int iflag,float *eph,float *q2,float *one_m_x, RNDM& rndm_generator);

  inline   void saveOnHadronFile(int num, float energy1, float energy2, float min_z, RNDM& rndm_generator, int DEF_DO_EDM4HEP )
  {
    int i;
    for (i=0;i<num;i++)
      {
	float temp_z = 1e-3*(min_z+rndm_generator.rndm())*mesh_.get_delta_z();

	hadron_file_->save_hadron(energy1,energy2,temp_z);
#ifdef USE_EDM4HEP
	if (do_edm4hep==1) EDM4HEPWriter::edmfile().save_hadron_EDM(energy1,energy2,temp_z);
#endif
      }
  }

  inline void clear_photon_pointer()
  {
    unsigned int j;
    for (j=0; j < grid_photon1_.size(); j++) grid_photon1_[j].clear();
    for (j=0; j < grid_photon2_.size(); j++) grid_photon2_[j].clear();
  }

 public :

  GRID();
  ~GRID();

  inline void connect_random_generator( RNDM* rndm_generator) 
  {
    rndm_generator_ = rndm_generator;
  }

  inline bool random_ok() const
  {
    return rndm_generator_ != NULL;
  }

  inline RESULTS& get_results() { return results_;}

  inline const COHERENT_RESULTS& get_coherent_results() const { return coherent_results_;}

  // inline const PHOTON_RESULTS& get_photon_results() const { return photon_results_;}

  inline const TRIDENT_RESULTS& get_trident_results() const { return trident_results_;}

  inline const COMPT& get_compton() const {return compt_;}

  /* inline  void get_distributes(DISTRIBUTE& distr1, DISTRIBUTE& distr2) const  */
  /* {  */
  /*   distr1 = distribute1_; */
  /*   distr2 = distribute2_; */
  /* } */

  inline void get_miss(float& miss1, long& out1, float& miss2, long& out2) const 
  {
    miss1 = distribute1_.delta;
    out1 = distribute1_.tot;
    miss2 = distribute2_.delta;
    out2 = distribute2_.tot;
  }

  inline const ABSTRACT_MINIJETS* get_minijets() const {return minijets_;}

  inline const BHABHA& get_bhabhas() const { return bhabhas_;}

  void read(const PARAMETERS& param, int automatic);

  inline void generalInit(const SWITCHES& switches, std::string photon_file,   std::string hadron_file, std::string c_phot_name, std::string bhabha_samples_file, std::string bhabha_photon_samples_file, GENERAL_CROSS* cross)
  {
    results_.connect_switches(&switches);
    if (switches.get_do_coherent()) coherent_results_.init();
    if (switches.get_do_trident()) trident_results_.init();
   
    lumi_init(switches);
    if ( !photon_file.empty() ) 
      {  
	photon_file_ = new FILE_IN_OUT();
	photon_file_->open_file(photon_file, "w");
      }
    else photon_file_ = NULL;
   
   
    if ( !hadron_file.empty() ) 
      {  
	hadron_file_ = new FILE_IN_OUT();
	hadron_file_->open_file(hadron_file, "w");
      }
    else hadron_file_ = NULL;
    cross_ = cross;
    if ( !c_phot_name.empty() ) 
      {
	compt_.connect_compt_phot_file(c_phot_name); 
      }
    if (switches.get_do_bhabhas()) bhabhas_.load_samples(switches.get_do_bhabhas(), bhabha_samples_file, bhabha_photon_samples_file);
  }

  void get_cuts(float& cx,float& cy,float& cz) const;
  //void init_grid_comp (int n_cell_x, int n_cell_y, int integration_method);
  void init_grid_comp (int n_cell_x, int n_cell_y, int integration_method,   FFT_SERVER* fourier = NULL);

  /*! Routine to calculate the parameters necessary for the iterative method of
    the potential calculation sor2 */
  void init_sor2 (PHI_FLOAT *parameter);

  inline void init_extra_photons()
  {
    int totalSize = n_cell_x_*n_cell_y_;
    //  extra_photons_cell1_ = std::vector< std::list<EXTRA> >(totalSize);
    //   extra_photons_cell2_ = std::vector< std::list<EXTRA> >(totalSize);
    extra_photon_pointer1_ = std::vector< std::list<EXTRA_PHOTON_POINTER> >(totalSize);
    extra_photon_pointer2_ = std::vector< std::list<EXTRA_PHOTON_POINTER> >(totalSize);
    clear_extra_photons();
    /*      extra_photons1_.init_extra_photons(totalSize); */
    /*      extra_photons2_.init_extra_photons(totalSize); */
    /*      extra_photons1_.clear_extra_photons(); */
    /*      extra_photons2_.clear_extra_photons(); */
  }

  inline void init_jet(float s,float ptmin,int iparam,int jet_pythia, int jet_select, std::string jetfileName)
  {
    if (!jet_pythia)
      {
	minijets_ = new MINIJETS(s, ptmin, iparam, jet_select,jetfileName);
      }
    else
      {
	minijets_ = new MINIJETS_PYTHIA(s, ptmin, iparam, jet_select,jetfileName);
      }
  };

  inline void all_distribute(int i1, int i2, SWITCHES& switches,double s4, double lns4)
  {
    // distribute the particle (electrons...) charges over the grid
    distribute_particles(i1, i2, switches.get_electron_distribution_rho(), switches.get_force_symmetric());
  
    // Distributes the particles for background calculation 
    distribute_particles_for_background(i1, i2, switches.get_electron_ratio());
  
    if (switches.get_do_pairs()||switches.get_do_hadrons()||switches.get_do_compt() || switches.get_do_muons())
      {
	distribute_virtual_photons(i1, i2, &switches,s4, lns4);
      }
  }

  virtual void distribute_particles(int i_slice1,
				    int i_slice2, 
				    int electron_distribution_rho, 
				    int force_symmetric);

  void distribute_photons(int slice_1,int slice_2, float photon_ratio, RNDM& rndm_generator);

  void check_distribute(int what);

  void save_lumi_on_files(SWITCHES& switches, std::string lumi_ee_out, std::string lumi_eg_out, std::string lumi_ge_out, std::string lumi_gg_out);

  void save_tertphot_on_file(std::string tertphot);

  inline void clear_extra_photons()
  {
    unsigned int j;
    for (j=0; j < extra_photon_pointer1_.size(); j++) extra_photon_pointer1_[j].clear();
    for (j=0; j < extra_photon_pointer2_.size(); j++) extra_photon_pointer2_[j].clear();
  }

  //! compute luminosity at each step
  /*!
    for each grid cell, call collide_ee(), computing all collisions between a particle from beam1 and a particle from beam2. The summed luminosity is added in 'results' (lumi_fine)
  */

  void step_lumi(float min_z,PAIR_BEAM& secondaries, int time_counter, SWITCHES& switches, int beamslice1, int beamslice2) ;

  inline void update_slice_charge(int slice1, int slice2)
  {
    /*  beam1_->meanPositionOfSlice(slice1, rho_x_1_, rho_y_1_); */
    /*  beam2_->meanPositionOfSlice(slice2, rho_x_2_, rho_y_2_); */
    slice_of_beam_[0].get_beam()->meanPositionOfSlice(slice1, rho_x_1_, rho_y_1_);
    slice_of_beam_[1].get_beam()->meanPositionOfSlice(slice2, rho_x_2_, rho_y_2_);
    float  rho_sum = slice_of_beam_[0].get_beam()->numberOfParticlesOfSlice(slice1);
    rho_sum_1_ = rho_sum*slice_of_beam_[0].get_nb_part_per_macro(); 
    rho_sum = slice_of_beam_[1].get_beam()->numberOfParticlesOfSlice(slice2); 
    rho_sum_2_ = rho_sum*slice_of_beam_[1].get_nb_part_per_macro();    
  }

  inline void moveAllParticles(const std::vector<GENERAL_GRID*>& grids,  int i_slice1, int i_slice2 , int interpolation,int do_beamstrahlung, int sokolov, float emin, int do_prod, int extra_grids, float charge_sign,int bmt_rotate,int do_trident)
  {
    move_particles(grids, 1, i_slice1,interpolation, do_beamstrahlung,do_trident,sokolov,  emin, do_prod, extra_grids, charge_sign,bmt_rotate);
    move_coherent_particles(grids, 1, i_slice1, interpolation, do_beamstrahlung, emin, do_prod, extra_grids, charge_sign);
  
    //  std::cout << " advance the second beam.... " << std::endl;
    move_particles(grids, 2, i_slice2,interpolation, do_beamstrahlung,do_trident,sokolov,  emin, do_prod, extra_grids, charge_sign,bmt_rotate);
    move_coherent_particles(grids, 2, i_slice2, interpolation, do_beamstrahlung, emin, do_prod, extra_grids, charge_sign);
    if(do_trident)
      {
	move_trident_particles(grids, 1, i_slice1, interpolation, do_beamstrahlung, emin, do_prod, extra_grids, charge_sign);
	move_trident_particles(grids, 2, i_slice2, interpolation, do_beamstrahlung, emin, do_prod, extra_grids, charge_sign);
      }
  };

  inline void move_photons(BEAM& beam,int i_beam,int i_slice)
  {
    float delta;
    delta=step_*    slice_of_beam_[i_beam-1].get_scal_step();
    beam.move_photon_beam(i_slice, delta);
  }

  void move_photons2(BEAM& beam,int ibeam,int i_slice,RNDM& rndm_generator);

  void move_pairs(const std::vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& rndm_generator);

  void move_pairs_tertphot(const std::vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& rndm_generator);

  void photon_lumi(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator);

  void photon_lumi_2(float min_z,SWITCHES& switches,PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator);

  void photon_lumi_3(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries,RNDM& rndm_generator);

  inline void set_bpm_signals() 
  {
    results_.bpm_signal(*slice_of_beam_[0].get_beam());
    results_.bpm_signal(*slice_of_beam_[1].get_beam());
    results_.bpm_signal_coherent(*slice_of_beam_[0].get_beam());
    results_.bpm_signal_coherent(*slice_of_beam_[1].get_beam());

  }

  virtual std::string output_flow() const;

#ifdef USE_EDM4HEP
  void EDM_output() const;
#endif
};

class EXTRA_GRID : public GENERAL_GRID
{

 private :

 EXTRA_GRID() : GENERAL_GRID() {;};
  
  virtual void distribute_coherent_particles(int i_slice1,
					     int i_slice2,
					     int electron_distribution_rho);
  
  virtual void distribute_trident_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho);
  
  //void assignCoherentBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro);
  void assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice);
  
  // void assignBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro);
  void assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice);
  
  // void assignBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro);
  void assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice);
  
  // void assignCoherentBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro);
  void assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice);
  
  void assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice);
  void assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice);
  
 public : 
 EXTRA_GRID(const GENERAL_GRID& grid) : GENERAL_GRID(grid) {;}
  ~EXTRA_GRID();
  
  // old name : void distribute_particles0_n(BEAM *beam1, int i_slice1,
  virtual void distribute_particles(int i_slice1,
				    int i_slice2,
				    int electron_distribution_rho, 
				    int force_symmetric);
};
#endif
