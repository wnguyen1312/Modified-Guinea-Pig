#include "particlesCPP.h"

#include <iostream>

void PARTICLE_WITH_SPIN::rotateBMT(TRIDVECTOR Efield, TRIDVECTOR Bfield, float charge_sign, float dt) 
{
  // Efield : electric field (lab. frame) : GV/nm
  // Bfield : magnetic field X c (lab. frame) : GV/nm
  // dt : time step in nanometers (c.dt)

  //  std::cout << " PARTICLE_WITH_SPIN::rotateBMT VERIFIER a " << std::endl;
  int k;
  TRIDVECTOR Blong, Btrans, betaE, omega;
  
  // gamma = E[GeV] / M[Gev/c2]
  double gamma = (double)energy_ / (double)EMASS;

  //double beta2 = 1.0 - 1.0/gamma;  //unused!!!
  //   double betaz = sqrt(beta2 - (double)vx_*(double)vx_ - (double)vy_*(double)vy_);
  //   std::cout << " betaz = " << betaz << std::endl;

  double vitx = (double)vx_;
  double vity = (double)vy_;

  // the longitudinal betaz is assumed to be equal to 1
  // so, the module of the velocity beta is 1 too

  // projection of B on the momentum
  double BlongProjection = vitx*Bfield(0) + vity*Bfield(1) + Bfield(2);
  Blong.setComponents(BlongProjection*vitx, BlongProjection*vity, BlongProjection);
  Btrans.setComponents(Bfield(0) - Blong(0), Bfield(1) - Blong(1), Bfield(2) - Blong(2));
  // P X beta
  betaE(0) = vity*Efield(2) - Efield(1);
  betaE(1) = Efield(0) - vitx*Efield(2);
  betaE(2) = vitx*Efield(1) - vity*Efield(0);

  //  double anom = (float)anomalousMoment(Efield,Bfield);
  double anom = (float)PHYSTOOLS::BformFunction(ups_);
  double a1 = 1 + anom;
  double ga1 = 1 + gamma*anom;
  double ce1 = anom + 1./(1.0 + gamma);
  //  std::cout << " a= " << anom << " ce1= " << ce1 << std::endl;

  //  float ang = 0.0; // alternative 1
  for (k = 0; k< 3; k++)
    {
      omega(k) = -(double)charge_sign*(ga1*Btrans(k) + a1*Blong(k) - ce1*gamma*betaE(k));

      // alternative 1 
      // ang += omega(k)*omega(k);
    }

  // alternative 1
  //   std::cout << " alternative 1 " << std::endl;
  //   ang = sqrt(ang);
  //   omega *= -1.0/ang;
  //   ang *= dt/energy_;
  //   double r[3][3];
  //   double omegtemp[3];
  //   for (k= 0; k < 3; k++) omegtemp[k] = omega(k);
  
  //   TOOLS::rotmat(ang,omegtemp, r );
  //     double spintemp[3];
  //     double bidon;
  //     for (k= 0; k < 3; k++) spintemp[k] = spin_(k);
  //     for (k = 0; k < 3; k++)
  //       {
  // 	   bidon = 0.0;
  //         for (j=0; j < 3; j++)
  //  	 {
  // 	   bidon +=  r[k][j]*spintemp[j];
  // 	   // 	   spin_[k] += r[k][j]*spintemp[j];
  //  	 }
  // 	spin_.setComponent(k,bidon);
  //      }
  // end alternative
  
  // this vector has to be multiplied by the factor -e/(m0.gamma)
  // m0 : rest mass of the particles
  // e : charge of the particle 
  // for ELECTRONS (or POSITRONS) : factor = - e.c^2/(m0.c^2.gamma)
  // equal to c^2/E[GeV] * 10^(-9)
  // one has then to multiply by dt[s] = dt[m]/c = dt[nm]*10^(-9)/c
  // 
  // the final factor is : - c (1/E[GeV]) dt[nm] 10^(-18)
  // BUT : the Electric field was given in GV/nm (and B in the consistent unit), so that
  // the factor 10^(-18) has been already taken into account. 
  // the magnetic field is already multiplied by c, so : 

  // alternative 2
  //  std::cout << " alternative 2 " << std::endl;
  omega *= (double)dt/(double)energy_;
  TRIDVECTOR spintemp(spin_);
  
  spin_(0) += omega(1)*spintemp(2) - omega(2)*spintemp(1);
  spin_(1) += omega(2)*spintemp(0) - omega(0)*spintemp(2);
  spin_(2) += omega(0)*spintemp(1) - omega(1)*spintemp(0);
  // end alternative 
  
  
  // renormalisation 
  //  spin_.renormalize();
}

void PAIR_PARTICLE::apply_electric_field_on_pair(float step_2, float ex, float ey) 
{
  // apply electric field : step_2
  // e_inv : inverse of energy
  if (icharge_)
    {
      ex = -ex;
      ey = -ey;
    }
  advanceVelocities(ex, ey, step_2);
}

float PAIR_PARTICLE::apply_magnetic_field_on_pair(float fac_theta, float step_q, float bx, float by) 
{
  const float eps=1e-35;
  float theta;
  float b_norm, b_norm_i,a1,a2,a3,vb;
  if (icharge_)
    {
      bx = -bx;
      by = -by;
    }
  
  // normalising B
  b_norm=sqrt(bx*bx+by*by);
  b_norm_i=1.0/std::max(b_norm,eps);
  bx *= b_norm_i;
  by *= b_norm_i;
  
  // v x B
  vb = vx_*by-vy_*bx;
  
  theta=fac_theta*(velocityz_*velocityz_*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv_*e_inv_*step_q;
  a3=0.5*theta;

  // 1 - theta**2/2
  a1=1.0-a3;

  // theta
  theta=sqrt(theta);
  a2=theta*velocityz_;

  // a3 = (theta**/2).(v.B)  [B norm]
  a3*=vx_*bx+vy_*by;

  // apply magnetic field (??)
  vx_ = vx_*a1-a2*by+a3*bx;
  vy_ = vy_*a1+a2*bx+a3*by;
  velocityz_ = velocityz_*a1+theta*vb;
  return theta;
}

float  PAIR_PARTICLE::scale_pair(float vold2, double mass)
{
#ifdef SCALE_ENERGY
  float scal;
  scal=sqrt( (vold2*energy_*energy_ + mass*mass)/(velocity_q()*energy_*energy_ + mass*mass)   );
  vx_ *= scal;
  vy_ *= scal;
  velocityz_ *= scal;
  energy_ /= scal;
  e_inv_ *=  scal;
  return scal;
#endif
}

void PAIR_PARTICLE::synchrotron_radiation(float step, double mass, float step_fraction, float vx0, float vy0, float vz0, std::vector<float>* photon_e, RNDM& rndm_generator) 
{
#ifdef PAIR_SYN
  unsigned int j;
  float eng0, scal;
  std::vector<float> ph;
  float dz = step*1e-9;

  std::cout << "Value of dz: " << dz << std::endl;

  float dzOnRadiusFrac = sqrt((vx_-vx0)*(vx_-vx0)+(vy_-vy0)*(vy_-vy0)+(velocityz_-vz0)*(velocityz_-vz0))/step_fraction;
  float radius_i = dzOnRadiusFrac/dz;

  TRIDVECTOR vdummy;   
      float upsilon=LAMBDA_BAR*EMASS/(mass*mass*mass)*energy_*energy_*radius_i;
      dzOnRadiusFrac*=EMASS/mass;

      std::cout << "SYNCHROTRON RADIATION IS CALLED " << std::endl;
      PHYSTOOLS::synrad_no_spin_flip(upsilon,energy_, dzOnRadiusFrac,ph, rndm_generator);
      if ((int)ph.size() > 0) 
	{
	  eng0 = energy_;
	  for (j=0;j < ph.size() ;j++)
	    {
	      energy_ -= ph[j];
	      photon_e->push_back(ph[j]);
	    }
	  scal=sqrt( ( (energy_*energy_-mass*mass)*eng0*eng0)/((eng0*eng0-mass*mass)*energy_*energy_) );
	  vx_ *= scal;
	  vy_ *= scal;
	  velocityz_ *= scal;
	}
#endif
}
void PAIR_PARTICLE::synchrotron_radiation(float step, double mass, float step_fraction, float vx0, float vy0, float vz0, RNDM& rndm_generator) 
{
#ifdef PAIR_SYN
  unsigned int j;
  float eng0, scal;
  std::vector<float> ph;
  float dz = step*1e-9;
  float dzOnRadiusFrac = sqrt((vx_-vx0)*(vx_-vx0)+(vy_-vy0)*(vy_-vy0)+(velocityz_-vz0)*(velocityz_-vz0))/step_fraction;
  float radius_i = dzOnRadiusFrac/dz;

  TRIDVECTOR vdummy;   
      float upsilon=LAMBDA_BAR*EMASS/(mass*mass*mass)*energy_*energy_*radius_i;
      dzOnRadiusFrac*=EMASS/mass;

      
      PHYSTOOLS::synrad_no_spin_flip(upsilon,energy_, dzOnRadiusFrac,ph, rndm_generator);
      if ((int)ph.size() > 0) 
	{

std::cout << "SYNCHROTRON RADIATION IS CALLED " << std::endl;

	  eng0 = energy_;
	  for (j=0;j < ph.size() ;j++)
	    {
	      energy_ -= ph[j];
	    }
	  scal=sqrt( ( (energy_*energy_-mass*mass)*eng0*eng0)/((eng0*eng0-mass*mass)*energy_*energy_) );
	  vx_ *= scal;
	  vy_ *= scal;
	  velocityz_ *= scal;
	}
#endif
}

float PAIR_PARTICLE::apply_final_half_step_fields(float step, double mass, float ex,float ey, float bx, float by, float /*thetamx*/, RNDM& /*rndm_generator*/)
{
  float vold2;
  float theta;
  vold2 = velocity_q();
  apply_electric_field_on_pair(0.5*step, ex, ey); 
  scale_pair(vold2, mass);   
  theta = apply_magnetic_field_on_pair(0.25, step*step, bx, by); 
  update_energy(mass);
  return theta;
}

float PAIR_PARTICLE::apply_initial_half_step_fields(float step, double mass, float ex,float ey, float bx, float by, std::vector<float>* photon_e, RNDM& rndm_generator)
{
  float vx0,vy0,vz0, vold2;
  float step_2, step_q, scal;
  float theta;
  step_2=0.5*step;
  step_q=step*step;
#ifdef PAIR_SYN
  // store the speed in v0
  //  velocities(vx0, vy0, vz0);

  ABSTRACT_PARTICLE::velocities(vx0, vy0);
  vz0 = velocityz_;
#endif
  theta = apply_magnetic_field_on_pair(0.25, step_q, bx, by); 
  //  thetamax=2.0*theta;


  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0); 
  
  synchrotron_radiation(step, mass, 0.5, vx0, vy0, vz0, photon_e, rndm_generator);
  return theta;
}
float PAIR_PARTICLE::apply_initial_half_step_fields(float step, double mass, float ex,float ey, float bx, float by, RNDM& rndm_generator)
{
  float vx0,vy0,vz0, vold2;
  float step_2, step_q, scal;
  float theta;
  step_2=0.5*step;
  step_q=step*step;
#ifdef PAIR_SYN
  // store the speed in v0
  //  velocities(vx0, vy0, vz0);

  ABSTRACT_PARTICLE::velocities(vx0, vy0);
  vz0 = velocityz_;
#endif
  theta = apply_magnetic_field_on_pair(0.25, step_q, bx, by); 
  //  thetamax=2.0*theta;


  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0); 
  
  synchrotron_radiation(step, mass, 0.5, vx0, vy0, vz0, rndm_generator);
  return theta;
}
float PAIR_PARTICLE::apply_full_step_fields(float step, double mass, float ex,float ey, float bx, float by, std::vector<float>* photon_e, RNDM& rndm_generator) 
{
  float vx0,vy0,vz0, vold2;
  float step_2, step_q, scal;
  float theta;
  step_2=0.5*step;
  step_q=step*step;

#ifdef PAIR_SYN
  // store the speed in v0
  //  velocities(vx0, vy0, vz0);
  ABSTRACT_PARTICLE::velocities(vx0, vy0);
  vz0 = velocityz_;
#endif
  /* scd new */
  // save the norm of the speed
  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0);   
  
  theta = apply_magnetic_field_on_pair(1.0, step_q, bx, by); 
  
  /* scd new */
  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0);  
  
  synchrotron_radiation(step, mass, 1.0, vx0, vy0, vz0, photon_e, rndm_generator);
  return theta;
}

float PAIR_PARTICLE::apply_full_step_fields(float step, double mass, float ex,float ey, float bx, float by, RNDM& rndm_generator) 
{
  float vx0,vy0,vz0, vold2;
  float step_2, step_q, scal;
  float theta;
  step_2=0.5*step;
  step_q=step*step;

#ifdef PAIR_SYN
  // store the speed in v0
  //  velocities(vx0, vy0, vz0);
  ABSTRACT_PARTICLE::velocities(vx0, vy0);
  vz0 = velocityz_;
#endif
  /* scd new */
  // save the norm of the speed
  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0);   
  
  theta = apply_magnetic_field_on_pair(1.0, step_q, bx, by); 
  
  /* scd new */
  vold2 = velocity_q();
  apply_electric_field_on_pair(step_2, ex, ey); 
  scal = scale_pair(vold2, mass);   
  scale_velocities_sync(scal, vx0, vy0, vz0);  
  
  synchrotron_radiation(step, mass, 1.0, vx0, vy0, vz0, rndm_generator);
  return theta;
}


void PAIR_PARTICLE::update_energy(double mass)
{
  float scal;
  last_rescaling_ok_ = true;
#ifdef SCALE_ENERGY
  scal=sqrt( (energy_*energy_-mass*mass)/(  (energy_*energy_)* velocity_q()   ));
  if (fabs(scal-1.0)>0.01)
    {
      last_rescaling_ok_ = false;
      std::cerr << " PAIR_PARTICLE::scale_energy : scale should be 1 : " << std::endl;
      std::cerr << " unsigned_energy_ " << energy_  <<  " module of velocity " << sqrt( velocity_q() ) << " scale " << scal << std::endl;
    }
  vx_ *= scal;
  vy_ *= scal;
  velocityz_ *= scal;
#endif
}
