#ifndef PARTICLEBEAM_SEEN
#define PARTICLEBEAM_SEEN
#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <list>
#include <string>
#include <vector>
#include "particlesCPP.h"
#include "rndmCPP.h"
#include "physconst.h"
#include "mathconst.h"


#include "fileInputOutput.h"
#include "abstractParticle.h"
#include "EDMwriter.h"

#ifdef USE_EDM4HEP
class FILE_IN_OUT_EDM4HEP;
class EDM4HEPWriter;
#endif

class ABSTRACT_PARTICLE_BEAM
{
  
 protected :

   float gamma_;

 public : 
   
   ABSTRACT_PARTICLE_BEAM() 
   {
     gamma_ = 0.0;
   }
 virtual ~ABSTRACT_PARTICLE_BEAM() {;}
 
 inline float gamma() const { return gamma_;}
 
 virtual void beamXyRms(float& xmean, float& ymean, float& sigmaxRms, float& sigmayRms) const = 0;
 virtual void beamZRms(float& zmean, float& sigmazRms) const = 0;
 
 // emittances in mm.mrad
 virtual  void emittances(float& emittx, float& emitty) const =0;
 
};

class BEAM_FROM_FILE : public ABSTRACT_PARTICLE_BEAM
{
  
  std::vector<PARTICLE_INTERFACE> particles_;
  
  BEAM_FROM_FILE() {;}
  
  public : 

    virtual ~BEAM_FROM_FILE() {;}

  
  BEAM_FROM_FILE(std::string filename)
    {
      FILE_IN_OUT filin;
      PARTICLE_INTERFACE parttemp;
      filin.open_file(filename, "r");
      
      unsigned int counter = 0;
      
      while(filin.read_particle(parttemp))
	{
	  counter++;
	  parttemp.transform_to_internal_units(1e3, 1e-6);
	  if (parttemp.good() )
	    {
	      particles_.push_back(parttemp);
	    }
	  else 
	    {
	      std::cerr << " BEAM_FROM_FILE:: error reading file " << filename << " line " << counter << " : skipped " << std::endl;
	    }
	}
      filin.close();
      compute_gamma();
    }
  
  inline const PARTICLE_INTERFACE& pick_last()
  {
    return particles_.back();
  }
  inline void pick_last_particle(float& ener, float& x, float& y, float& z, float& vx, float& vy, TRIDVECTOR& polar)
  {
    particles_.back().get_parameters(ener, x, y, z, vx, vy);
    polar = particles_.back().getSpin();
  }

  inline void erase_last_particle()
    {
      particles_.pop_back();
    }
  
  inline bool not_empty() const
    {
      return (int) particles_.size();
    }
  
  virtual void emittances(float& emittx, float& emitty) const;
  
  virtual  void beamXyRms(float& x0, float& y0,float& sigmax, float& sigmay) const;
  virtual  void beamZRms(float& z0,float& sigmaz) const;
  


  inline void compute_gamma()
    {
      unsigned int k;
      gamma_ = 0.0;
      int n_particles = (int) particles_.size();
      if (n_particles <= 0) return;
      for (k=0; k < particles_.size(); k++)
	{
	  gamma_ += particles_[k].energy();
	}
      gamma_ /= EMASS*n_particles;
    }  
};


class PARTICLE_BEAM : public ABSTRACT_IO_CLASS, public ABSTRACT_PARTICLE_BEAM
{

  std::vector< std::vector<PARTICLE*> > particle_;
  std::vector< std::vector<PARTICLE*> > coherent_;
  std::vector< std::vector<PARTICLE*> > trident_;
  unsigned long int initial_number_of_particles_;
  unsigned long int number_of_particles_dispatched_in_slices_;
  RNDM* rndm_generator_;
  unsigned int typOfParticle_;
  float sigmax_, sigmay_, sigmaz_;
  TRIDVECTOR polarization_;

  inline void set(int bmt_rotate,TRIDVECTOR polar)
  {
    initial_number_of_particles_ = 0;
    number_of_particles_dispatched_in_slices_ = 0;
    gamma_ = 0.0;
    if (!ZPOS) 
      {
	typOfParticle_ = 0;
	return;
      }
    if (bmt_rotate)
      {
	typOfParticle_ = 3;
	polarization_ = polar;
      }
    else
      {
	/* 	if (EXT_FIELD)  */
	/* 	  { */
	typOfParticle_ = 2;
	polarization_ = 0.0;
	/*       } */
	/* 	else */
	/* 	  { */
	/* 	    typOfParticle_ = 1; */
	/* 	    polarization_ = 0.0; */
	/* 	  } */
      }
  }

  inline void compute_gamma()
    {
      unsigned long int j,k;
      gamma_ = 0.0;
      if (number_of_particles_dispatched_in_slices_ <= 0) return;
      for (j=0; j < particle_.size(); j++)
	{
	  for (k=0; k < particle_[j].size(); k++)
	    {
	      gamma_ += particle_[j][k]->energy();
	    }
	}
      gamma_ /= EMASS*number_of_particles_dispatched_in_slices_;
    }
  
  inline void backParticlesBefore(float max_z)
    {
      unsigned long int j,k;
      float dist=0.5*(max_z-SQRT3*sigmaz_);
      
      for (j = 0; j < particle_.size(); j++)
	{
	  for ( k = 0; k < particle_[j].size(); k++)
	    {
	      particle_[j][k]->advancePosition(-dist);
	    }
	}
    }
    
  inline void backSlicesWithZ(float max_z)
    {
      unsigned long int j,k;
      float z, dist;
      // int n_slice = particle_.size();
      for (k=0;k<particle_.size();k++)
	{
	  for (j=0;j< particle_[k].size();j++)
	    {
	      z = particle_[k][j]->z();
	      dist = 0.5*(max_z+z);
	      particle_[k][j]->advancePosition(-dist);
	    }
	}
    }
  inline void backSlicesWithZ2(float max_z)
    {
      unsigned long int j,k;
      float z, dist;
      //int n_slice = particle_.size();
      for (k=0;k<particle_.size();k++)
	{
	  for (j=0;j<particle_[k].size();j++)
	    {
	      z = particle_[k][j]->z();
	      dist = 0.5*(max_z - z);
	      particle_[k][j]->advancePosition(-dist);
	    }
	}
    }
  
  inline void backSlices(int beam, float step,  int timestep, float scalStep[2])
    {
      unsigned long int k, j;
      int n_slice = (int) particle_.size();
      float dist=step*timestep*0.5 *(scalStep[beam-1] - scalStep[2-beam])*(n_slice-1);
      dist+=0.5*(1.0-1.0/timestep)*step*timestep;
      for (k=0;k<particle_.size();k++)
	{
	  for (j=0;j<particle_[k].size();j++)
	    {
	      particle_[k][j]->advancePosition(-dist);
	    }
	  dist+=step*timestep*scalStep[2-beam];
	}
    }
    
  inline void backSlices2(int nbeam, float step,  int timestep, float scalStep[2])
    {
      unsigned long int k, j;
      int n_slice = (int) particle_.size();
      float dist=step*timestep*0.5*((n_slice+1)*scalStep[nbeam-1] + (n_slice-1)*scalStep[2-nbeam]);
      for (k=0;k<particle_.size();k++)
	{
	  for (j=0;j<particle_[k].size();j++)
	    {
	      particle_[k][j]->advancePosition(-dist);
	    }
	  dist -= step*timestep*scalStep[2-nbeam];
	}
    }
  
  void fill_beam(int dist_x, int dist_z, float delta_z, float sigma_x,float sigma_y, float sigma_z, float sigma_x_prime,float sigma_y_prime, float energy, int do_bds_spin_rotation);
  
  void fill_symmetric_beam(int dist_x, int dist_z, float delta_z, float sigma_x,float sigma_y, float sigma_z, float sigma_x_prime,float sigma_y_prime, float energy, int do_bds_spin_rotation);
  

  inline void get_rndm_xy_particle(float sigmax, float sigmay, float sigxp, float sigyp, float& x, float& y, float& vx, float& vy) const
    {
      x = rndm_generator_->gasdev()*sigmax;
      y = rndm_generator_->gasdev()*sigmay;
      vx = rndm_generator_->gasdev()*sigxp;
      vy = rndm_generator_->gasdev()*sigyp;
    }

  void assign_xyz_normal_distribution(int dist_z,float delta_z,float sigma_x,float sigma_y, float sigmaz, float sigma_x_prime,float sigma_y_prime, float energy, int do_bds_spin_rotation);
  
  void assign_symmetric_xyz_normal_distribution(int dist_z,float delta_z,float sigma_x,float sigma_y,float sigmaz, float sigma_x_prime,float sigma_y_prime, float energy, int do_bds_spin_rotation);
  
  //  inline void dispatch_random_particle_in_slices(float ztemp, float delta_z, float sigma_x,float sigma_y, float /*sigma_z*/, float sigma_x_prime,float sigma_y_prime, float energy, int nSlices)
  //    {
  //      float xtemp, ytemp, vxtemp, vytemp;
  //      int j=(int)floor(ztemp/delta_z+0.5* nSlices+1);
  //      TRIDVECTOR polar;
  //      if (typOfParticle_ == 3)
  //	{
  //	  polar = polarization_;
  //	}
  //      else polar = 0.;
  //      if ((j>0)&&(j<= nSlices))
  //	{
  //	  get_rndm_xy_particle(sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, xtemp, ytemp, vxtemp, vytemp);
  //	  set_new_particle_in_slice(j-1,xtemp, ytemp,ztemp, vxtemp, vytemp,energy, polar);
  //	}
  //    }

  //New method that approximately takes into account the spin rotation from the final focus system
  inline void dispatch_random_particle_in_slices(float ztemp, float delta_z, float sigma_x,float sigma_y, float /*sigma_z*/, float sigma_x_prime,float sigma_y_prime, float energy, int nSlices, int do_bds_spin_rotation)
    {
      float xtemp, ytemp, vxtemp, vytemp;
      int j=(int)floor(ztemp/delta_z+0.5* nSlices+1);
      TRIDVECTOR polar;
      if ((j>0)&&(j<= nSlices))
	{
	  get_rndm_xy_particle(sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, xtemp, ytemp, vxtemp, vytemp);
	  
	  if (typOfParticle_ == 3)
	    {
	      if(do_bds_spin_rotation)
		{
		  double polx, poly, polz, temp0, temp1, rotangle, theta,qi,qj,qr;
		  rotangle=2.26938288523942*energy; //\gamma*alpha/(2*pi) in rad
		  polarization_.getComponents(polx, poly, polz);

		  temp0=sqrt(vxtemp*vxtemp+vytemp*vytemp); // here temp0 is for vector normalization
		  theta=asin(temp0);
		  theta*=rotangle;

		  qi=-vytemp*sin(0.5*theta)/temp0;
		  qj=vxtemp*sin(0.5*theta)/temp0;
		  qr=cos(0.5*theta);

		  temp0=(1.0-2.0*qj*qj)*polx+2.0*qi*qj*poly+2.0*qj*qr*polz; //re-using temp0 for new polx
		  temp1=2.0*qi*qj*polx+(1.0-2.0*qi*qi)*poly-2.0*qi*qr*polz;
		  polz=-2.0*qj*qr*polx+2.0*qi*qr*poly+(1.0-2.0*(qi*qi+qj*qj))*polz;
		  polar.setComponents(temp0,temp1,polz);
		  //		  double polx, poly, polz, temp, rotangle;
		  //		  rotangle=2.26938288523942*energy; //\gamma*alpha/(2*pi) in rad
		  //		  polarization_.getComponents(polx, poly, polz);
		  //		  //Spin rotation in the small angle approximation
		  //		  temp=polz-rotangle*(vxtemp*polx+vytemp*poly);
		  //		  polx=polx+rotangle*vxtemp*polz;
		  //		  poly=poly+rotangle*vytemp*polz;
		  //		  polar.setComponents(polx,poly,temp);
		}
	      else polar = polarization_;
	    }
	  else polar = 0.;
	  set_new_particle_in_slice(j-1,xtemp, ytemp,ztemp, vxtemp, vytemp,energy, polar);
	}
    }
  

  //  inline void dispatch_symmetric_random_particle_in_slices(float ztemp, float delta_z, float sigma_x,float sigma_y, float /*sigma_z*/, float sigma_x_prime,float sigma_y_prime, float energy, int nSlices)
  //    {
  //      float xtemp, ytemp, vxtemp, vytemp;
  //      double polx, poly, polz;
  //      TRIDVECTOR polar;
  //      if (typOfParticle_ == 3)
  //	{
  //	  polar = polarization_;
  //	  polar.getComponents(polx, poly, polz);
  //	}
  //      else 
  //	{
  //	  polx = 0.0;
  //	  poly = 0.0;
  //	  polz = 0.0;
  //	}
  //      int j=(int)floor(ztemp/delta_z+0.5* nSlices+1);
  //      if ((j>0)&&(j<= nSlices))
  //	{
  //	  get_rndm_xy_particle(sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, xtemp, ytemp, vxtemp, vytemp);
  //	  j--;
  //	  polar.setComponents(polx, poly, polz);
  //	  set_new_particle_in_slice(j, xtemp, ytemp, ztemp, vxtemp, vytemp, energy, polar);
  //	  polar.setComponents(-polx, poly, polz);
  //	  set_new_particle_in_slice(j, -xtemp, ytemp, ztemp, -vxtemp, vytemp, energy, polar);
  //	  polar.setComponents(-polx, -poly, polz);
  //	  set_new_particle_in_slice(j, -xtemp, -ytemp, ztemp, -vxtemp, -vytemp, energy, polar);
  //	  polar.setComponents(polx, -poly, polz);
  //	  set_new_particle_in_slice(j, xtemp, -ytemp, ztemp, vxtemp, -vytemp, energy, polar);	
  //	  
  //	}
  //    }

  //New method that approximately takes into account the spin rotation from the final focus system
  inline void dispatch_symmetric_random_particle_in_slices(float ztemp, float delta_z, float sigma_x,float sigma_y, float /*sigma_z*/, float sigma_x_prime,float sigma_y_prime, float energy, int nSlices, int do_bds_spin_rotation)
    {
      float xtemp, ytemp, vxtemp, vytemp;
      double polx, poly, polz;
      TRIDVECTOR polar;

      int j=(int)floor(ztemp/delta_z+0.5* nSlices+1);
      if ((j>0)&&(j<= nSlices))
	{
	  get_rndm_xy_particle(sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, xtemp, ytemp, vxtemp, vytemp);
	  if (typOfParticle_ == 3)
	    {
	      if(do_bds_spin_rotation)
		{
		  double temp0 ,temp1, rotangle, theta,qi,qj,qr;
		  rotangle=2.26938288523942*energy; //\gamma*alpha/(2*pi) in rad
		  polarization_.getComponents(polx, poly, polz);
		  //Spin rotation in the small angle approximation
		  if(vytemp==0.0){
		    theta=vxtemp;
		  }else if(vxtemp==0.0){
		    theta=vytemp;
		  }else{
		    theta=asin(vxtemp/cos(atan(vxtemp/vytemp)));
		  }
		  theta*=rotangle;
		  temp0=sqrt(vxtemp*vxtemp+vytemp*vytemp);
		  qi=-vytemp*sin(0.5*theta)/temp0;
		  qj=vxtemp*sin(0.5*theta)/temp0;
		  qr=cos(0.5*theta);

		  temp0=(1.0-2.0*qj*qj)*polx+2.0*qi*qj*poly+2.0*qj*qr*polz; //re-using temp0
		  temp1=2.0*qi*qj*polx+(1.0-2.0*qi*qi)*poly-2.0*qi*qr*polz;
		  polz=-2.0*qj*qr*polx+2.0*qi*qr*poly+(1.0-2.0*(qi*qi+qj*qj))*polz;
		  polx=temp0;
		  poly=temp1;
		    //		  temp=polz-rotangle*(vxtemp*polx+vytemp*poly);
		    //		  polx=polx+rotangle*vxtemp*polz;
		    //		  poly=poly+rotangle*vytemp*polz;
		    //		  polz=temp;
		}
	      else
		{
		  polar = polarization_;
		  polar.getComponents(polx, poly, polz);
		}
	    }
	  else 
	    {
	      polx = 0.0;
	      poly = 0.0;
	      polz = 0.0;
	    }
	  j--;
	  polar.setComponents(polx, poly, polz);
	  set_new_particle_in_slice(j, xtemp, ytemp, ztemp, vxtemp, vytemp, energy, polar);
	  polar.setComponents(-polx, poly, polz);
	  set_new_particle_in_slice(j, -xtemp, ytemp, ztemp, -vxtemp, vytemp, energy, polar);
	  polar.setComponents(-polx, -poly, polz);
	  set_new_particle_in_slice(j, -xtemp, -ytemp, ztemp, -vxtemp, -vytemp, energy, polar);
	  polar.setComponents(polx, -poly, polz);
	  set_new_particle_in_slice(j, xtemp, -ytemp, ztemp, vxtemp, -vytemp, energy, polar);	
	}
    }
  
  void  set_particles_offset(int slice, float offset_x,float offset_y,float waist_x,float waist_y);

  inline void set_new_particle_in_slice(int slice, float x,float y,float z,float vx,float vy,float energy, TRIDVECTOR polar)
    {
      push_back_empty_particle(slice,1);
      int index = (int) particle_[slice].size() - 1;
      particle_[slice][index]->setParticle(x,y,z,vx,vy,energy, &polar);
      number_of_particles_dispatched_in_slices_++;
    }
  
  inline void set_new_particle_in_slice(int slice, const PARTICLE_INTERFACE& part)
    {
      push_back_empty_particle(slice,1);
      int index = (int) particle_[slice].size() - 1;
      particle_[slice][index]->setParticle(part);
      
      number_of_particles_dispatched_in_slices_++;
    }
  
  inline PARTICLE* newPart(const PARTICLE& part) const 
    {
      PARTICLE* newp;
      switch(typOfParticle_)
	{
	  /*        case 1: */
	  /* 	 newp = new PARTICLE(part); */
	  /* 	 break; */
       case 2:
	 newp=  new PARTICLE(part);  
	 break;
       case 3:
	 newp =  new PARTICLE_WITH_SPIN(dynamic_cast<const PARTICLE_WITH_SPIN&>(part) ); 
	 break; 
       default:
	 std::cerr << " PARTICLE_BEAM::newPart(const PARTICLE* part) : unknown type of particle " << std::endl;
	 newp = NULL;
	 exit(0);
	} 
     return newp;
    }

  inline PARTICLE* newPart() const 
    {
      PARTICLE* newp;
      switch(typOfParticle_)
	{
	  /*        case 1: */
	  /* 	 newp = new PARTICLE(); */
	  /* 	 break; */
	case 2:
	  newp=  new PARTICLE();  
	  break;
	case 3:
	  newp =  new PARTICLE_WITH_SPIN(); 
	  break; 
	default:
	  std::cerr << " PARTICLE_BEAM::newPart() : unknown type of particle " << std::endl;
	  newp = NULL;
	  exit(0);
	} 
      return newp;
    }
  
  inline void push_back_empty_particle(int slice, int number)
    {
      int k;
      switch(typOfParticle_)
	{
	  /*        case 1: */
	  /* 	 for (k=0; k<number; k++) particle_[slice].push_back(new PARTICLE()); */
	  /* 	 break; */
	case 2:
	  for (k=0; k<number; k++) particle_[slice].push_back(new PARTICLE());
	  break;
	case 3:
	  for (k=0; k<number; k++) particle_[slice].push_back(new PARTICLE_WITH_SPIN());
	  break;
	default:
	  std::cerr << " PARTICLE_BEAM::make : unknown type of particle " << std::endl;
	  exit(0);
	}
    }

  inline PARTICLE* generic_particle_for_dump(PARTICLE& part,int istep,float dz0, float max_z, int sign_label ) const
    {
      float z,dz;
      float sign = 1e-3*(float)sign_label;
      
      z = part.z();
      dz = max_z-istep*dz0;
      z += dz;
      z = sign*(z + dz0);
      part.setZ(z);
      return &part;  //BD 17/06/2010
    }
  
  inline  PARTICLE* generic_velocity_corrected_particle_for_dump(PARTICLE& part,int slice, int istep,int complement, float dz0, float max_z, int sign_label ) const
    {
      float dz;
      float vx, vy;
      PARTICLE* newp = generic_particle_for_dump(part, istep, dz0, max_z, sign_label);
      part.velocities(vx, vy);
      dz=(slice-istep + complement -1)*dz0;
      newp->apply_position_offset(-dz*vx, -dz*vy);
      return newp;
    }
  
  
  public : 


 PARTICLE_BEAM() : ABSTRACT_IO_CLASS(), ABSTRACT_PARTICLE_BEAM(),
    initial_number_of_particles_(0),number_of_particles_dispatched_in_slices_(0),
    typOfParticle_(0),sigmax_(0.0),sigmay_(0.0),sigmaz_(0.0)
    {
      rndm_generator_ = NULL;
    }
  
  ~PARTICLE_BEAM()
    {
      unsigned long int j,k;
      for (j=0; j < particle_.size(); j++)
	{
	  for (k=0; k < particle_[j].size(); k++)
	    {
	      delete particle_[j][k];
	    }
	}
      for (j=0; j < coherent_.size(); j++)
	{
	  for (k=0; k < coherent_[j].size(); k++)
	    {
	      delete coherent_[j][k];
	    }
	}
      for (j=0; j < trident_.size(); j++)
	{
	  for (k=0; k < trident_[j].size(); k++)
	    {
	      delete trident_[j][k];
	    }
	}
    }

 PARTICLE_BEAM(int nb_slices,int bmt_rotate,TRIDVECTOR polar, RNDM* rndm_generator) : ABSTRACT_IO_CLASS(), ABSTRACT_PARTICLE_BEAM(),
    initial_number_of_particles_(0),number_of_particles_dispatched_in_slices_(0),
    typOfParticle_(0),sigmax_(0.0),sigmay_(0.0),sigmaz_(0.0)

    {
      set(bmt_rotate, polar);
      particle_ = std::vector< std::vector<PARTICLE*> >(nb_slices);
      coherent_ = std::vector< std::vector<PARTICLE*> >(nb_slices);
      trident_ = std::vector< std::vector<PARTICLE*> >(nb_slices);
      rndm_generator_ = rndm_generator;
    }
  
  inline  std::vector<PARTICLE*>& getParticleVector(int slice)  { return particle_[slice];}
  
  inline const std::vector<PARTICLE*>& getParticleVector(int slice) const 
    {
      return particle_[slice];
    }

  inline  std::vector<PARTICLE*>& getCoherentVector(int slice) {return coherent_[slice];}

  inline  std::vector<PARTICLE*>& getTridentVector(int slice) {return trident_[slice];}
  
  inline const std::vector<PARTICLE*>& getCoherentVector(int slice) const {return coherent_[slice];}

  inline const std::vector<PARTICLE*>& getTridentVector(int slice) const {return trident_[slice];}
    
  inline void newCoherent(int slice,float x,float y, float vx,float vy, float energy)
    {
      switch(typOfParticle_)
	{
	  /*        case 1: */
	  /* 	 coherent_[slice].insert(coherent_[slice].begin(),new PARTICLE()); */
	  /* 	 break; */
	case 2:
	  coherent_[slice].insert(coherent_[slice].begin(), new PARTICLE());
	  break;
	case 3:
	  coherent_[slice].insert(coherent_[slice].begin(), new PARTICLE_WITH_SPIN());
	  break;
	default:
	  std::cerr << " PARTICLE_BEAM::newCoherent : unknown type of particle " << std::endl;
	  exit(0);
	}   
      
      coherent_[slice][0]->setParticle(x,y,0.0,vx,vy,energy);  
    }

  inline void newTrident(int slice,float x,float y, float vx,float vy, float energy)
    {
      switch(typOfParticle_)
	{
	  // case 1 missing
	case 2:
	  trident_[slice].insert(trident_[slice].begin(), new PARTICLE());
	  break;
	case 3:
	  trident_[slice].insert(trident_[slice].begin(), new PARTICLE_WITH_SPIN());
	  break;
	default:
	  std::cerr << " PARTICLE_BEAM::newTrident : unknown type of particle " << std::endl;
	  exit(0);
	}   
      trident_[slice][0]->setParticle(x,y,0.0,vx,vy,energy);  
    }

  void rotate_particles(float angle);
  void set_angle_particles(float x_angle, float y_angle,float delta_z);

  inline   void  set_particles_offset(float offset_x,float offset_y,float waist_x, float waist_y)
    {
      int j;
      for (j = 0; j < (int) particle_.size(); j++)
	{
	  set_particles_offset(j, offset_x, offset_y, waist_x, waist_y);   
	}
    }
  
  inline void adjustToMeanPosition()
    {
      int k;
      float xpart, ypart;
      for (k = 0; k < (int) particle_.size(); k++)
	{
	  meanPositionOfSlice(k,xpart, ypart);
	  set_particles_offset(k, -xpart, -ypart, 0.0, 0.0);
	}
    }
  
  void backstep (int beam,int trav_focus, float max_z, float step,  int timestep, float scal_step[2]);

  void backstep2 (int nbeam, float max_z, float step,  int timestep, float scal_step[2]);

  // after loading, the object pointed by bff is empty and 
  // the pointer bff is invalidated
  unsigned int load_particles(BEAM_FROM_FILE*& bff, float emin, float zmin, float deltaz, float sigx, float sigy, float sigz);

  void init_particles(unsigned long int nbpart, float sigma_x,float sigma_y, float sigma_z,int dist_x,int dist_z,float emx, float emy,float delta_z,float energy, int charge_symmetric, int do_bds_spin_rotation);
  
  virtual void beamXyRms(float& xmean, float& ymean, float& sigmaxRms, float& sigmayRms) const;
  virtual void beamZRms(float& zmean, float& sigmazRms) const;

  // emittances in mm.mrad
  virtual  void emittances(float& emittx, float& emitty) const;
  
  double meanLostEnergy(float ebeam) const;

  inline float  sigmaz() const { return sigmaz_;}
  inline void transverse_sigmas(float& sigx, float& sigy) const
    {
      sigx = sigmax_;
      sigy = sigmay_;
    }

  inline void index_slice(int number, int& index, int& slice) const
  {
    slice = -1;
    do
      {
	slice++;
	number -= (int) particle_[slice].size();
      }
    while (number >= 0);
    index = number + (int) particle_[slice].size();
  }

  void meanPositionOfSlice(int slice, float& x,float& y) const;
  
  inline int number_of_coherent_vectors() const {return (int) coherent_.size();}
  
  inline int number_of_slices() const { return (int) particle_.size();}

  unsigned long int numberOfParticles() const {return number_of_particles_dispatched_in_slices_;}
  
  inline unsigned long int numberOfParticlesOfSlice(int slice) const
    {
      return (int) particle_[slice].size();
    };
  
  void transverseRms(int slice,double& xmin,double& xmax,double& xmean,double&  ymin,double& ymax,double& ymean,double&  sigmaxRms,double& sigmayRms) const; 

  void ang_dis(unsigned int n_bin, std::vector< std::vector<float> >& bin ) const;
  
  virtual std::string output_flow() const; 

#ifdef USE_EDM4HEP
  void EDM_output(int beam_label) const;
#endif

  int store_beam(std::string name) const;

#ifdef USE_EDM4HEP
  void store_beam_EDM(BEAM_TYPE type) const;
  void store_coherent_beam_EDM(BEAM_TYPE type) const;
  void store_trident_beam_EDM(BEAM_TYPE type) const;
#endif

  void dump_beam(std::string name, int istep, int every_particle, int timestep, float step, float max_z, int  sign_label);
  void store_coherent_beam(std::string name) const;
  void store_trident_beam(std::string name) const;

};


class  PHOTON_BEAM 
{

 private:

  int n_slice_;
  std::vector< std::vector<PHOTON> > slice_photon_vector_;
  // AL: this seems useless std::list<PHOTON>::iterator end_;
  
  PHOTON_COUNTER photon_count_;

  inline void photon_info_slice(int i_slice,double& sum, int& n, int& n_tot) const
    {
      double s=0.0;
      long number=0;
      unsigned int k;
      for (k=0; k < slice_photon_vector_[i_slice].size(); k++)
	{
	  if (slice_photon_vector_[i_slice][k].energy() != 0.0)
	    {
	      s += fabs(slice_photon_vector_[i_slice][k].energy());
	      number++;
	    }
	  n_tot++;
	}
      sum = s;
      n = number;
    }
  
  inline void new_loaded_photon(int sliceDATA, float energy, float hel, float xt, float yt, float vx,float vy)
    {
      if( sliceDATA < 0 || sliceDATA >= n_slice_) return;
      // guineapig original, don't attribute a value to z for a photon
      photon_count_.addPhoton(energy);
      slice_photon_vector_[sliceDATA].insert(slice_photon_vector_[sliceDATA].begin(), PHOTON(xt, yt,0.0, vx, vy, energy, hel,photon_count_.getNumber()));
    }


  public:
    
 PHOTON_BEAM() : n_slice_(0) /* AL: this seems useless   , end_(NULL) */ {}
  
  PHOTON_BEAM(int n_sliceDATA);
  
  ~PHOTON_BEAM();
  
  
  inline std::vector<PHOTON>& getPhotonVector(int slice) { return slice_photon_vector_[slice];}
  
  inline const std::vector<PHOTON>& getPhotonVector(int slice) const { return slice_photon_vector_[slice];}
  
  inline const PHOTON_COUNTER& get_photon_counter() const { return photon_count_;}
  
  void load_photons(std::string filename, int type_of_beam, float delta_z, float max_z, int n_cell_z);
  
  inline int sizeOfSlice(int slice) const  {return (int) slice_photon_vector_[slice].size();}
  
  //int store_photon(std::string name) const;

  inline const PHOTON& new_photon(float energy, PARTICLE& particle,  int slice)
    {
      photon_count_.addPhoton(energy);
      slice_photon_vector_[slice].push_back(PHOTON(energy, particle, 0.0, photon_count_.getNumber()));
      return slice_photon_vector_[slice].back();
    }
  
  inline void photon_info(double& sum,int& number) const
    {
      int i,n,nt=0;
      double s;
      sum = 0.0;
      number = 0;
      for (i=0;i<n_slice_;i++)
	{
	  photon_info_slice(i,s,n,nt);
	  sum += s;
	  number += n;
	}
    }
  
  inline void move_photons(int i_slice, float delta)
    {
      unsigned int k;
      for (k=0; k < slice_photon_vector_[i_slice].size(); k++)
	{
	  slice_photon_vector_[i_slice][k].advancePosition(delta);
	}
    }
  
  void dump_photons(std::string name,int istep, int every_particle,int timestep, float step, float max_z, int  sign_label);
  
};

#endif
