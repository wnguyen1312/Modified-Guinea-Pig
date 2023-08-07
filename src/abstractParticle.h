#ifndef ABSTRACTPARTICLE_SEEN
#define ABSTRACTPARTICLE_SEEN
#include <iostream>
#include <sstream>
#include <string>

#include "typeDefs.h"
#include "abstractIOclass.h"
#include "mathematicalEntities.h"

class ABSTRACT_BEAM
{

 protected:
 public: 
  ABSTRACT_BEAM() {;}
  virtual ~ABSTRACT_BEAM() {;}
  virtual int sizeOfPhotonSlice(int slice) const = 0;
  virtual int number_of_slices() const =0;
  virtual int numberOfParticlesOfSlice(int slice) const = 0;
  
};


class ABSTRACT_PARTICLE
{
  
 protected: 
  
  // unit of xpos_, ypos_ : nm 
  float xpos_, ypos_;
  
  float zpos_;
  
  // unit of velocity : radian (?)
  float vx_,vy_;

  // energy in GeV
  float energy_;

  inline void set(float energyi, float xi, float yi,float zi,float vxi,float vyi)
    {
      energy_ = energyi;
      xpos_= xi;
      ypos_ = yi;
      zpos_ = zi;
      vx_ = vxi;
      vy_ = vyi;
    }
  
  inline void set(const ABSTRACT_PARTICLE& part)
    {
      set (part.energy_, part.xpos_, part.ypos_, part.zpos_, part.vx_, part.vy_);
    }
  
  
  // internal units seems to be nm for x,y,z
  inline void transform_to_internal_units(float facPosition, float facVelocity)
    {
      xpos_ *= facPosition;
      ypos_ *= facPosition;
      zpos_ *= facPosition;
      
      vx_ *= facVelocity;
      vy_ *= facVelocity;
    }
  
  
  ABSTRACT_PARTICLE (float xi, float yi,float zi, float vxi,float vyi,float energyi)
    {
      set(energyi, xi, yi,zi, vxi, vyi);
    }

  
  public : 

	//added method to access the position of each particle --> should remove print function later

    float getXpos() const
    {
        //std::cout << "Xpos: " << xpos_ << std::endl;
        return xpos_;
    }

    float getYpos() const
    {
        //std::cout << "Ypos: " << ypos_ << std::endl;
        return ypos_;
    }

    float getZpos() const
    {
        //std::cout << "Zpos: " << zpos_ << std::endl;
        return zpos_;
    }



    ABSTRACT_PARTICLE () 
    {
      set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
    }
  
  virtual ~ABSTRACT_PARTICLE () {;}
  
  inline float z() const {return zpos_;}
  inline void setZ( float zi) { zpos_ = zi;};
  
  virtual const TRIDVECTOR& getSpin() const = 0;
  
  
  inline void razEnergy() {energy_ = 0.0;}
  
  
  inline void velocities(float& vxout, float& vyout) const 
    {
      vxout = vx_;
      vyout = vy_;
    };
  
  inline void velocities(double& vxout, double& vyout) const 
    {
      vxout = (double)vx_;
      vyout = (double)vy_;
    };
  
  inline void XYposition(PHI_FLOAT& xpos, PHI_FLOAT& ypos) const
    {
      xpos = (PHI_FLOAT)xpos_;
      ypos = (PHI_FLOAT)ypos_;
    };
  
  
  inline void XYposition(float& xpos, float& ypos) const
    {
      xpos = xpos_;
      ypos = ypos_;
    };
  
  inline float energy() const {return energy_;};
  
  inline void setEnergy(float en) {energy_ = en;}; 
  
  inline void advancePosition(float distance)
    {
      xpos_ += vx_*distance;
      ypos_ += vy_*distance;
    };
  
  inline void advanceVelocities(float Fx, float Fy, float step)
  {
    vx_ += step*Fx/energy_;
    vy_ += step*Fy/energy_;
  }
  
  inline void advanceVelocities(float Fx, float Fy, float step, float& deltaVx, float& deltaVy)
    {
      deltaVx = step*Fx/energy_;
      deltaVy = step*Fy/energy_;
      vx_ += deltaVx;
      vy_ += deltaVy;
    }

  
  inline void set_z_for_dump(int istep,float dz0, float max_z, int sign_label ) 
    {
      float z,dz;
      float sign = (float)sign_label;
      
      z = zpos_;
      dz = max_z-istep*dz0;
      //  std::cout << " set_z_for_dump istep= " << istep <<  " dz= " << dz << std::endl;
      z += dz;
      zpos_ = sign*(z + dz0);
    }
  
  inline void set_z_velocity_corrected_for_dump(int slice, int istep,int complement, float dz0, float max_z, int sign_label ) 
    {
      float dz;
      //      float vx, vy unused;
      set_z_for_dump(istep, dz0, max_z, sign_label);
      dz=(slice-istep + complement -1)*dz0;
      advancePosition(-dz);
    }
};



class ABSTRACT_BHABHA_PHOTON_SAMPLES
{
 public :

   ABSTRACT_BHABHA_PHOTON_SAMPLES() {;}
 virtual  ~ABSTRACT_BHABHA_PHOTON_SAMPLES() {;}
 virtual  int get_label() const = 0;
 virtual  unsigned int nb_samples() const = 0;
 
 virtual  void get_parameters_for_output(unsigned int number, int& number_bhabha, float& en,float& vx,float& vy, float&vz) const = 0;
 
 virtual  void add_bhabha_photon(int nbhabha, float px, float py, float pz, float en) = 0;
// virtual  void create_bhabha_photon(int nbhabha, float px, float py, float pz, float en) = 0;
 
 
};

class ABSTRACT_BHABHASAMPLES
{
  public : 
    ABSTRACT_BHABHASAMPLES() {;}
  virtual ~ABSTRACT_BHABHASAMPLES() {;}
  virtual unsigned int nb_samples() const = 0;
  
  virtual  void get_parameters_for_output(unsigned int number, unsigned int& evtIdx, float& eCM, float& mother1_en,float&e1,float&vx1,float& vy1, float&vz1, float& mother2_en, float& e2, float& vx2, float&vy2, float&vz2, int& nbphot) const = 0;
  
  virtual  void add_bhabha(unsigned int evtIdx, float px1, float py1, float pz1, float e1, float px2, float py2, float pz2, float e2, int nbphot) = 0;
  
};

class ABSTRACT_CROSS_DATA
{
  
  public :
    
    ABSTRACT_CROSS_DATA() {;}
  virtual ~ABSTRACT_CROSS_DATA() {;}
  virtual void resize(int n, int nval) = 0;
  virtual void add_data(float ener, const float* data) =0;
  
};


class ABSTRACT_LUMI_HEAP : public ABSTRACT_IO_CLASS
{

 public:
  ABSTRACT_LUMI_HEAP() {;}
  virtual ~ABSTRACT_LUMI_HEAP() {;}

  virtual  int nb_pairs() const = 0;
  virtual  void get_parameters_for_output(unsigned int number, float& e1,float& e2,float& x,float& y,float& z) const = 0;
  virtual  void get_parameters_for_output(unsigned int number, float& e1,float& e2,float& x,float& y,float& z, float& vx1,float& vy1,float& vx2,float& vy2, float& sx1, float& sy1, float& sz1, float& sx2, float& sy2, float& sz2,int& t)  const = 0;
  virtual std::string output_flow() const 
    {
      std::ostringstream out;
      out << " ABSTRACT_LUMI_HEAP:: no data for output file in abstract class " << std::endl;
      return out.str();
    }
};

class PARTICLE_INTERFACE : public ABSTRACT_PARTICLE
{
  friend class BEAM_FROM_FILE;
  protected : 
    
    
    TRIDVECTOR polar_;
  
  
  
 public:
  
  PARTICLE_INTERFACE() {;}
  virtual ~PARTICLE_INTERFACE() {;}
  
  inline bool good() const
    {
    bool test = true;
    // if the energy is 0 or negative the reading of the particle on the
    // file has probably failed. This test is made because the program
    // has to divide by energy.
    if ( energy_ < 1.9e-12 ) test = false;
    return test;
  }

 inline void init_from_input_file(float energy, float xpos, float ypos, float zpos, float vx, float vy, float polx, float poly, float polz)
   {
     set(energy, xpos, ypos,zpos,vx,vy);
     polar_.setComponents((double)polx, (double)poly, (double)polz);
   }



 inline void get_parameters(float& energy, float& xpos, float& ypos, float& zpos, float& vx, float& vy) const
 {
   xpos = xpos_;
   ypos = ypos_;
   zpos = zpos_;
   vx = vx_;
   vy = vy_;
   energy = energy_;
 }



 inline float get_helicity() const { return polar_(0);}
 
 virtual inline const TRIDVECTOR& getSpin() const 
   { 
     return polar_;
   }
 
};


#endif
