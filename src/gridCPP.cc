#include <iostream>
#include "gridCPP.h"
#include "mathconst.h"
#include "physconst.h"
#include "mathematicalTools.h"
#include "physicalTools.h"
#include <fstream>
GENERAL_GRID::GENERAL_GRID()
{
  n_cell_x_=32;
  n_cell_y_=32;
  n_cell_z_=32;
  timestep_ =10;

  min_x_ = 0.0;
  max_x_ = 0.0; 
  min_y_ = 0.0;
  max_y_ = 0.0;
  min_z_ = 0.0;
  max_z_ = 0.0;

  cut_x_ = 0.0;
  cut_y_ = 0.0;
  cut_z_ = 0.0;

  integration_method_ = 0;
  step_ = 0.0;


  rho_x_1_ = 0.0;
  rho_y_1_ = 0.0;
  rho_sum_1_ = 0.0;
  rho_x_2_ = 0.0;
  rho_y_2_ = 0.0;
  rho_sum_2_ = 0.0;

  delta_x_inverse_ = 0.0;
  delta_y_inverse_ = 0.0;

  rho_factor_ = 0.0;

}

void GENERAL_GRID::interpolePotential(float xpart,float ypart, PHI_FLOAT& h_x, PHI_FLOAT& h_y, PHI_FLOAT& phi1_x, PHI_FLOAT& phi2_x, PHI_FLOAT& phi3_x, PHI_FLOAT& phi1_y,PHI_FLOAT& phi2_y, PHI_FLOAT& phi3_y, const PHI_FLOAT *phi) const
{
  int i1, i2;
  float h;
  absoluteCoordinates(xpart, ypart,h_x, h_y);

  // locate in Y
  cellLocalization(h_x, h_y, i1, i2, h);

  phiValuesY(i1, i2,h, phi, phi1_y, phi2_y, phi3_y);

  // locate in X
  cellLocalization(h_y, h_x,i2,i1,h); 


  phiValuesX(i1, i2, h, phi, phi1_x, phi2_x,phi3_x );


  h_x -= floor(h_x);
  h_y -= floor(h_y);
}

// electric field in GV/nm
TRIVECTOR GENERAL_GRID::ElectricFieldCIC(float xpart,float ypart, const PHI_FLOAT *phi) const
{
  //  int k;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  TRIVECTOR Efield;
  interpolePotential(xpart, ypart, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);

  Efield(0) = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_;
  Efield(1) = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_;

  Efield(2) = 0.0;
  // the potential seems to have been multiplied by 2 (for taking into 
  // account the magnetic field)
  Efield *= 0.5;
  return Efield;
}

void GENERAL_GRID::computeFields(int integrationMethod, float charge_sign, PHI_FLOAT *sor_parameter) 
    {
      // (see commentaries in init_grid_phys() )

      //int nn[2];
      switch (integrationMethod)
	{
	case 1:	  
	  champ_.foldfields(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho(),charge_sign);
	  break;
	case 2: 
	  
	  champ_.fold_fft(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho(),charge_sign);
	  break; 
	       
	case 3:
	  champ_.foldfronts(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho() , sor_parameter, charge_sign);
	}
    }


GENERAL_GRID::~GENERAL_GRID()
{
  //  delete [] rho1_;
  //  delete [] rho2_;
}



void  GENERAL_GRID::init_grid_phys (float n_particles1, float n_particles2, float cut_x,float cut_y,float cut_z, float charge_sign,   FFT_SERVER* fourier)
{
  //int i1,i2,j,j0,i;
  // PHI_FLOAT factor,phi0;
   PHI_FLOAT factor;
   //double x0,y0;
   //int nn[2];

  float offsetx, offsety;
  float deltax, deltay, deltaz;
  
  slice_of_beam_[0].update_nb_part_per_macro(n_particles1);
  slice_of_beam_[1].update_nb_part_per_macro(n_particles2);
    /////
  cut_x_ = cut_x;
  cut_y_ = cut_y;
  cut_z_ = cut_z;
  min_x_=-((float)n_cell_x_-2)/((float)n_cell_x_)*cut_x;
  max_x_=((float)n_cell_x_-2)/((float)n_cell_x_)*cut_x;
  min_y_=-((float)n_cell_y_-2)/((float)n_cell_y_)*cut_y;
  max_y_=((float)n_cell_y_-2)/((float)n_cell_y_)*cut_y;
  offsetx=0.5*(float)n_cell_x_;
  offsety=0.5*(float)n_cell_y_;
  deltax=2.0*cut_x/((float)n_cell_x_);
  deltay=2.0*cut_y/((float)n_cell_y_);
  deltaz=2.0*cut_z/((float)n_cell_z_);
  mesh_ = MESH(deltax, deltay, deltaz, offsetx, offsety);
  delta_x_inverse_ = 1.0/deltax;
  delta_y_inverse_ = 1.0/deltay;


  min_z_ = -cut_z;
  max_z_ = cut_z;
  step_=deltaz/(2.0*(float)timestep_);
  //

  factor=-charge_sign*2.0*RE/(deltaz)*EMASSeV;
  rho_factor_ = 2.0*factor;
  factor/=deltax*deltay;
  champ_.dist_init(factor, deltax, deltay,fourier);
}


GRID::GRID() : GENERAL_GRID()
{
  photon_file_ = NULL;
  hadron_file_ = NULL;
  minijets_ = NULL;
  cross_ = NULL;
  rndm_generator_ = NULL;
}

GRID::~GRID() 
{
  if (photon_file_ != NULL) 
    {
      photon_file_->close();
      delete photon_file_;
    }
  if (hadron_file_ != NULL) 
    {
      hadron_file_->close();
      delete hadron_file_;
    }
  delete minijets_;
}




void  GRID::lumi_init(const SWITCHES& switches)
{
  int nmax;
  if(switches.get_do_lumi()){
    if (switches.get_do_lumi()&1) 
      {
	nmax = switches.get_num_lumi();
	if (nmax < 100) nmax = 100;
	lumi_heap_ee_ = LUMI_HEAP_EE(nmax, switches.get_lumi_p(), rndm_generator_);
      }
    if (switches.get_do_lumi()&2) 
      {
	nmax = switches.get_num_lumi_eg();
	if (nmax < 100) nmax = 100;
	lumi_heap_eg_ = LUMI_HEAP(nmax,switches.get_lumi_p_eg(), rndm_generator_);
	lumi_heap_ge_ = LUMI_HEAP(nmax, switches.get_lumi_p_eg(), rndm_generator_);
      }
    if (switches.get_do_lumi()&4) 
      {
	nmax = switches.get_num_lumi_gg();
	if (nmax < 100) nmax = 100;
	lumi_heap_gg_ = LUMI_HEAP( nmax,switches.get_lumi_p_gg(), rndm_generator_);
      }
  }
}

void GRID::save_lumi_on_files(SWITCHES& switches, std::string lumi_ee_out, std::string lumi_eg_out, std::string lumi_ge_out, std::string lumi_gg_out)
{
  if(switches.get_do_lumi())
    {
      if (switches.get_do_lumi()&1) 
	{
	  lumi_heap_ee_.saveLumi(lumi_ee_out);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      lumi_heap_ee_.saveLumi_EDM(EDM4HEPWriter::edmfile());
	    }
#endif
	}
      if (switches.get_do_lumi()&2) 
	{	
	  lumi_heap_eg_.saveLumi(lumi_eg_out);
	  lumi_heap_ge_.saveLumi(lumi_ge_out);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      lumi_heap_eg_.saveLumi_eg_EDM(EDM4HEPWriter::edmfile());
	      lumi_heap_ge_.saveLumi_ge_EDM(EDM4HEPWriter::edmfile());
	    }
#endif
	}
      if (switches.get_do_lumi()&4) 
	{
	  lumi_heap_gg_.saveLumi(lumi_gg_out);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      lumi_heap_gg_.saveLumi_gg_EDM(EDM4HEPWriter::edmfile());
	    }	
#endif
	}
    }
}

void GRID::save_tertphot_on_file(std::string tertphotfile)
{
  FILE_IN_OUT filout;
  filout.open_file(tertphotfile,"w");
  for(unsigned int i =0; i<tertphot_.size();i++)
    {
      filout.save_object_on_persistent_file(&tertphot_[i]);
    }
  filout.close();
}

void GRID::read(const PARAMETERS& param, int automatic)
{
  int n_n;
  float n_f;
  if (automatic != 1) 
    {
      n_n = param.readIValue("n_x");
      if (n_n > 0) n_cell_x_= n_n;
      n_n = param.readIValue("n_y");
      if (n_n > 0) n_cell_y_ = n_n;
      n_n = param.readIValue("n_z");
      if (n_n > 0) n_cell_z_ = n_n;
    }
  n_n = param.readIValue("n_t");
  if (n_n > 0) timestep_ = n_n;
  n_n = param.readIValue("n_m.1");
  //  if (n_n > 0) nb_macroparticles_[0] = n_n;
  if (n_n > 0) slice_of_beam_[0].set_macroparticles(n_n);
  n_n = param.readIValue("n_m.2");
  //  if (n_n > 0) nb_macroparticles_[1] = n_n;
  if (n_n > 0) slice_of_beam_[1].set_macroparticles(n_n);

  // the indexation is unusual, but it is like that in the original guinea-pig C
  n_f = param.readFValue("scale_step.1");
  //  if (n_f > 0.0) scal_step[1] = n_f;
  if (n_f > 0.0) slice_of_beam_[1].set_scal_step(n_f);  
  n_f = param.readFValue("scale_step.2");
  //  if (n_f > 0.0) scal_step[0] = n_f;
  if (n_f > 0.0) slice_of_beam_[0].set_scal_step(n_f);
}


std::string GRID::output_flow() const 
{
  std::ostringstream out;
  out << title(std::string("grid parameters"));
  out << "n_x = " << n_cell_x_ << " n_y = " << n_cell_y_ << std::endl;

  out << "n_z = " << n_cell_z_ << " n_t = " << timestep_ << std::endl;

  out << "n_m.1 = " << slice_of_beam_[0].get_macroparticles()  << " n_m.2 = " << slice_of_beam_[1].get_macroparticles() << std::endl;

  out << "cut_x = " << cut_x_ << " nm ; cut_y = " << cut_y_ << " nm ; cut_z = " << cut_z_*1e-3 << " micrometers " << std::endl;

  out << std::endl;
  out << " ...................................................... " << std::endl;
  out << " relative amount of interacting particles that were outside the grid during one time step : " << std::endl;
  out << "beam 1 : miss = " << distribute1_.delta <<  "  beam 2 : miss = " << distribute2_.delta << std::endl;
  out << " number of interacting part. that were outside the grid during one time step : " << std::endl;
  out << "out.1=" << distribute1_.tot <<  ";out.2=" << distribute2_.tot << ";" << std::endl; 


 return out.str();
}

#ifdef USE_EDM4HEP
void GRID::EDM_output() const
{
  EDM4HEPWriter::edmfile().write_Metadata("n_x",n_cell_x_);
  EDM4HEPWriter::edmfile().write_Metadata("n_y",n_cell_y_);
  EDM4HEPWriter::edmfile().write_Metadata("n_z",n_cell_z_);
  EDM4HEPWriter::edmfile().write_Metadata("n_t",timestep_);
  EDM4HEPWriter::edmfile().write_Metadata("n_m.1",slice_of_beam_[0].get_macroparticles());
  EDM4HEPWriter::edmfile().write_Metadata("n_m.2",slice_of_beam_[1].get_macroparticles());
  EDM4HEPWriter::edmfile().write_Metadata("cut_x",cut_x_);
  EDM4HEPWriter::edmfile().write_Metadata("cut_y",cut_y_);
  EDM4HEPWriter::edmfile().write_Metadata("cut_z",cut_z_*1e-3);
  EDM4HEPWriter::edmfile().write_Metadata("beam_1_miss",distribute1_.delta);
  EDM4HEPWriter::edmfile().write_Metadata("beam_2_miss",distribute2_.delta);
  EDM4HEPWriter::edmfile().write_Metadata("out.1",distribute1_.tot);
  EDM4HEPWriter::edmfile().write_Metadata("out.2",distribute2_.tot);
}
#endif

void GRID::init_grid_comp(int n_cell_x, int n_cell_y, int integration_method,   FFT_SERVER* fourier)
{
  n_cell_x_ = n_cell_x;
  n_cell_y_ = n_cell_y;
  integration_method_ = integration_method;
   champ_ = FIELD(integration_method, n_cell_x_, n_cell_y_, fourier);
   slice_of_beam_[0].resize(n_cell_x_, n_cell_y_);
   slice_of_beam_[1].resize(n_cell_x_, n_cell_y_);

  part_pointer1_ =  std::vector< std::list<BEAM_PARTICLE_POINTER> >(n_cell_x_*n_cell_y_);
  part_pointer2_ =  std::vector< std::list<BEAM_PARTICLE_POINTER> >(n_cell_x_*n_cell_y_);


  grid_photon1_ = std::vector< std::list<BEAM_PHOTON_POINTER> >(n_cell_x_*n_cell_y_);
  grid_photon2_ = std::vector< std::list<BEAM_PHOTON_POINTER> >(n_cell_x_*n_cell_y_);
}


/*! Routine to calculate the parameters necessary for the iterative method of
   the potential calculation sor2 */

void GRID::init_sor2 (PHI_FLOAT *parameter)
{
  PHI_FLOAT factor;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();
  //  float deltaz = mesh_.get_delta_z();

  factor=1.0/(8.0*PI*RE/mesh_.get_delta_z()*1e9/2000.0);
  parameter[0]=factor*(deltay/deltax);
  parameter[1]=parameter[0];
  parameter[2]=factor*(deltax/deltay);
  parameter[3]=parameter[2];
  parameter[4]=-2.0*(parameter[0]+parameter[2]);
  parameter[5]=(deltay*deltay*cos(PI/(float)n_cell_x_)
		+deltax*deltax*cos(PI/(float)n_cell_y_))
                        /(deltax*deltax+deltay*deltay);
}



void GRID::check_distribute(int what)
{
    float tmp;
    if (what==0){
	distribute1_.delta=0.0;
	distribute2_.delta=0.0;
	distribute1_.tot=0;
	distribute2_.tot=0;
    }
    if (what<=1){
	distribute1_.in=0; distribute1_.out=0;
	distribute2_.in=0; distribute2_.out=0;
    }
    if (what==2){
	tmp=(float)distribute1_.out
	    /((float)(distribute1_.in+distribute1_.out));
	if (tmp>distribute1_.delta) distribute1_.delta=tmp;
	tmp=(float)distribute2_.out
	    /((float)(distribute2_.in+distribute2_.out));
	if (tmp>distribute2_.delta) distribute2_.delta=tmp;
	if (distribute1_.out>distribute1_.tot){
	  distribute1_.tot=distribute1_.out;
	}
	if (distribute2_.out>distribute2_.tot){
	  distribute2_.tot=distribute2_.out;
	}
    }
    if (what==3){
      	printf("miss_1=%f;miss_2=%f;\n",distribute1_.delta,distribute2_.delta);
    }
}




/*
  Distributes the beam particles for field calculation
*/

void GRID::distribute_particles(int i_slice1,
				 int i_slice2, 
				 int electron_distribution_rho, 
				 int force_symmetric)
{
  
  slice_of_beam_[0].razRho();
  slice_of_beam_[1].razRho();
  //  razRhos(); 

  switch (electron_distribution_rho)
    {
    case 1:
      assignBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
      assignBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
      break;
    case 2:

      assignBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
      assignBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
      break;

    }

  distribute_coherent_particles(i_slice1, i_slice2, electron_distribution_rho);
  distribute_trident_particles(i_slice1, i_slice2, electron_distribution_rho);  
  if (force_symmetric) 
    {
      slice_of_beam_[0].symmetrizeCharges(n_cell_x_,n_cell_y_);
      slice_of_beam_[1].symmetrizeCharges(n_cell_x_,n_cell_y_);
    }
}




/*
  Distribute the coherent pair particles for the calculation of the fields
 */

void GRID::distribute_coherent_particles(int i_slice1,
					  int i_slice2,
					  int electron_distribution_rho)
{
  switch (electron_distribution_rho){
    case 1:


      assignCoherentBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
      assignCoherentBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
      break;
  case 2:

    assignCoherentBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
    assignCoherentBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  }
  
}

/*
  Distribute the trident pair particles for the calculation of the fields
 */
void GRID::distribute_trident_particles(int i_slice1,int i_slice2,int electron_distribution_rho)
{
  switch (electron_distribution_rho){
  case 1:
    assignTridentBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
    assignTridentBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  case 2:
    assignTridentBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
    assignTridentBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  }
}

//! Distributes the particles for background calculation 

void GRID::distribute_particles_for_background(int i_slice1,
				 int i_slice2,  
				 float electron_ratio)
{
  int j,i1,i2;
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;
  
    for (i1=0;i1<n_cell_x_;i1++)
      {
	for (i2=0;i2<n_cell_y_;i2++)
	  {
	    j=i1*n_cell_y_+i2;
	    part_pointer1_[j].clear();
	    part_pointer2_[j].clear();
	  }
      }

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  
    distributeScatter1( slice_of_beam_[0].get_beam()->getParticleVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
   distributeScatter1(slice_of_beam_[1].get_beam()->getParticleVector(i_slice2), ratio, ratio_i_2, part_pointer2_);
  distribute_coherent_particles_for_background(i_slice1, i_slice2, electron_ratio);
  distribute_trident_particles_for_background(i_slice1, i_slice2, electron_ratio);
}

//
void GRID::distributeScatter1(const std::vector<PARTICLE*>& theParticles, float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer)
{

  float xpart, ypart;
  int i1, i2;
  //  int k,j;
  int j;
  unsigned int k;
  //    const std::vector<PARTICLE*>& theParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < theParticles.size(); k++)
    {
      if (rndm_generator_->rndm()<ratio)
	{
	theParticles[k]->XYposition(xpart, ypart);
	  if (particleInGrid(xpart, ypart, i1, i2))
	    {
	      j=i1*n_cell_y_+i2;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(theParticles[k], ratio_i));
	    }
	}
    }
}

// This method is not used, uses time to calculate.
// Don't delete for the moment

void GRID::distributeScatter2(const std::vector<PARTICLE*>& theParticles,float ratio, float ratio_i, std::vector< std::list<BEAM_PARTICLE_POINTER> >& part_pointer)
{

  float xpart, ypart,poids,h_x,h_y;
  int i1, i2;
  int j;
  unsigned int k;

  //  const std::vector<PARTICLE*>& theParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < theParticles.size(); k++)
      {
	theParticles[k]->XYposition(xpart, ypart);
	if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	  {
	    if (rndm_generator_->rndm()<ratio)
	      {
		j=i1*n_cell_y_+i2;
		poids = (1.0-h_x)*(1.0-h_y)*ratio_i;

	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(theParticles[k], poids));
		j=(i1+1)*n_cell_y_+i2;
		poids =  h_x*(1.0-h_y)*ratio_i;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(theParticles[k], poids));
		j=i1*n_cell_y_+i2+1;
		poids =  (1.0-h_x)*h_y*ratio_i;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(theParticles[k], poids));
		j=(i1+1)*n_cell_y_+i2+1;
		poids = h_x*h_y*ratio_i;
			      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(theParticles[k], poids));
	      }
	  }
      }
}


//  Distributes the particles from coherent pair creation for the background
//  calculation
 

void GRID::distribute_coherent_particles_for_background(int i_slice1,
					  int i_slice2,
					  float electron_ratio)

{

  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
    distributeScatter1(slice_of_beam_[0].get_beam()->getCoherentVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
        distributeScatter1(slice_of_beam_[1].get_beam()->getCoherentVector(i_slice2),ratio, ratio_i_2, part_pointer2_);
}

//  Distributes the particles from trident pair creation for the background
//  calculation
void GRID::distribute_trident_particles_for_background(int i_slice1,int i_slice2,float electron_ratio)
{
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  distributeScatter1(slice_of_beam_[0].get_beam()->getTridentVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
  distributeScatter1(slice_of_beam_[1].get_beam()->getTridentVector(i_slice2), ratio, ratio_i_2, part_pointer2_);
}

//! Distributes the virtual photons 

void GRID::distribute_virtual_photons(int i_slice1,
				 int i_slice2,
				 SWITCHES* switches,
				 double s4, double lns4)
{
   const float   emass2=EMASS*EMASS;
  float xmin,r_scal;
  float ratio,ratio_i_1,ratio_i_2;
    int geom;

  r_scal=switches->get_r_scal();
  geom=switches->get_geom();
  ratio=switches->get_electron_ratio();
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();
  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  clear_extra_photons();
  //   extra_photons1_.clear_extra_photons();
  //   extra_photons2_.clear_extra_photons();
  //  switch (switches->get_electron_distribution_scatter())
  //   {
  //    case 1:
  xmin=switches->get_compt_x_min()*emass2/s4;

  electronScatter(extra_photon_pointer1_, 1, i_slice1, xmin, s4, lns4, ratio, switches->get_ext_field(), geom, ratio_i_1, r_scal);
  
  electronScatter( extra_photon_pointer2_,  2, i_slice2, xmin, s4, lns4, ratio, switches->get_ext_field(), geom, ratio_i_2, r_scal);
  

	  
  //     break;
  // case 2:
  
  //   distributeElectronScatter2(beam1_, i_slice1, ratio, ratio_i_1, part_pointer1_, rndm_generator);
  //   distributeElectronScatter2(beam2_, i_slice2, ratio, ratio_i_2, part_pointer2_, rndm_generator);
  //   break;
  //   }
}

void GRID::electronScatter(std::vector< std::list<EXTRA_PHOTON_POINTER> >& extra_phot_ptr, int i_beam, int i_slice, float xmin, double s4,double lns4, float ratio, int ext_field, int geom, float ratio_i, float r_scal)
{
  unsigned int k;
  int i1, i2,j;
  int i_equiv=6;
  int n_phot, i_phot;
  float energy, r_phot;
  float e_phot,q2,one_m_x,x=0.0,y=0.0, radius, theta;
  float xVelocity, yVelocity;
  //PHYSTOOLS phys;

    const std::vector<PARTICLE*>& theParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < theParticles.size(); k++)
      //  for (i=beam->firstParticleOfSlice(i_slice);i<beam->firstParticleOfSlice(i_slice+1);i++)
    {
      //      energy=beam->energyOfParticle(i);
      energy=theParticles[k]->energy();
      //r_phot=phys.requiv(lns4,xmin,i_equiv)*ratio;
      r_phot=PHYSTOOLS::requiv(lns4,xmin,i_equiv)*ratio;
      n_phot=(int)floor(r_phot);
      r_phot -= n_phot;
      if(rndm_generator_->rndm()<r_phot) n_phot++;
      for (i_phot=0;i_phot<n_phot;i_phot++)
	{

	  //phys.mequiv(s4, lns4,xmin,energy,i_equiv,&e_phot,&q2,&one_m_x, *rndm_generator_);
	  PHYSTOOLS::mequiv(s4, lns4,xmin,energy,i_equiv,&e_phot,&q2,&one_m_x, *rndm_generator_);
	  //#ifdef EXT_FIELD
	  if (ext_field && (pow(q2/(EMASS*EMASS),1.5)*energy*energy < e_phot*e_phot*theParticles[k]->getUps()) )
	    {
	      e_phot=-1.0;
	    }
	  //#endif
	  if (e_phot>0.0)
	    {
	      switch(geom)
		{
		case 0:
		  theParticles[k]->XYposition(x, y);
		  break;
		case 1:
		  radius=HBAR*Cvelocity/sqrt(q2*one_m_x)*r_scal*1e9;
		  radius=std::min(radius,float(1.0e5));

		  theParticles[k]->XYposition(x, y);
		  x += rndm_generator_->rndm_sincos(&theta)*radius;
		  y += theta*radius;
		  break;
		case 2:
		  radius=HBAR*Cvelocity/sqrt(q2*one_m_x)*r_scal*1e9;
		  radius=std::min(radius,float(1e5));
		  theParticles[k]->XYposition(x, y);
		  x += rndm_generator_->gasdev()*radius;
		  y += rndm_generator_->gasdev()*radius;
		  break;
		}
	      if (particleInGrid(x, y, i1, i2))
		{
		  j=i1*n_cell_y_+i2;
		  theParticles[k]->velocities(xVelocity, yVelocity);
		  extra_phot_ptr[j].push_back(EXTRA_PHOTON_POINTER(e_phot, xVelocity, yVelocity,q2,energy,ratio_i));
		  //		  photon->store_vir_photon(e_phot,xVelocity, yVelocity,q2,energy,ratio_i,j);
		}
	      
	    }
	}
    }
}





void GRID::assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{

  unsigned int k;
  int i1, i2;
  float xpart, ypart;
  const std::vector<PARTICLE*>& theParticles = sog.get_beam()->getParticleVector(i_slice);    
  for (k = 0; k < theParticles.size(); k++)
    {
      theParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, i1, i2))
	{
	  sog.assignChargeToNGP( i1, i2, 1.0);
	  distribute.in++;
	}
      else   distribute.out++;		
    }
}


void GRID::assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  // n_macro : nb of particles per macroparticle
  unsigned int k;
  int i1, i2;
  float h_x,h_y;
  float xpart, ypart;
const std::vector<PARTICLE*>& theParticles = sog.get_beam()->getParticleVector(i_slice);
 for (k = 0; k < theParticles.size(); k++)
	{
	theParticles[k]->XYposition(xpart, ypart);
	  if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	    {
	      sog.assignChargeToCIC( i1, i2, h_x, h_y, 1.0);
	      distribute.in++;	    
	    }
	  else
	    {
	      //	      std::cout << " particle no "  << k+1 << " slice " << i_slice << " out of grid : x= " << xpart << " y= " << ypart << std::endl;
	      distribute.out++;
	    }
	  
	}
}


void GRID::assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const std::vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      
      if (particleInGrid(xpart, ypart,i1, i2))
	{
	  
	  sog.assignChargeToNGP( i1, i2,ch);
	  distribute.in++;
	}
      else
	{
	  distribute.out++;
	}
    }
}

void GRID::assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const std::vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;      
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	{
	  sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);
	  distribute.in++;
	}
      else{
	distribute.out++;
      }
    }
}

void GRID::assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const std::vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      
      if (particleInGrid(xpart, ypart,i1, i2))
	{
	  
	  sog.assignChargeToNGP( i1, i2,ch);
	  distribute.in++;
	}
      else
	{
	  distribute.out++;
	}
    }
}

void GRID::assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const std::vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;      
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	{
	  sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);
	  distribute.in++;
	}
      else{
	distribute.out++;
      }
    }
}

/*! This routine calculates the luminosity from the collision of the two
slices i1 and i2 of grid. The dimensions of the grids have to be the same. */

void GRID::step_lumi(float min_z, PAIR_BEAM& secondaries,int time_counter, SWITCHES& switches, int beamslice1=-1, int beamslice2=-1)
{
  float sum=0.0;
  int i1,i2,j;
  std::list<BEAM_PARTICLE_POINTER>::const_iterator pointer1,pointer2;

  const PHI_FLOAT* rho1 = slice_of_beam_[0].get_rho();
  const PHI_FLOAT* rho2 = slice_of_beam_[1].get_rho();
  for (i1=0;i1 < n_cell_x_;i1++)
    {
      for (i2=0;i2 < n_cell_y_;i2++)
	{
	  j=i1*n_cell_y_+i2;
	  sum += rho1[j]*rho2[j];
	  for (pointer1=part_pointer1_[j].begin(); pointer1 != part_pointer1_[j].end(); pointer1++)
	    {
	      for( pointer2 = part_pointer2_[j].begin(); pointer2 != part_pointer2_[j].end(); pointer2++)
		{
		  collide_ee(i1, i2, min_z, *pointer1, *pointer2,switches, secondaries,time_counter, beamslice1, beamslice2 ) ;
		}
	    }
	}
    }

  results_.add_lumi_fine( sum*1e18/(mesh_.get_delta_x()*mesh_.get_delta_y()*timestep_) );	    
}

void GRID::collide_ee(int cellx, int celly,float min_z, const BEAM_PARTICLE_POINTER& pointer1, const BEAM_PARTICLE_POINTER& pointer2, SWITCHES& switches,PAIR_BEAM& secondaries, int time_counter, int beamslice1=-1, int beamslice2=-1)
{
  //  JET_FLOAT bhabhan;
  float help,ecm,e1,e2;
  //float ecmratio;
  float weight = pointer1.weight() * pointer2.weight();
  int j1=0,j2=1;
  int bmt_precession = switches.get_bmt_precession();

  e1 =  fabs(pointer1.energy());
  e2 =  fabs(pointer2.energy());
  if (switches.get_do_espread())
    {
      e1 = spread_energy(e1, switches.get_which_espread1(), switches.get_espread1(), *rndm_generator_);
      e2 = spread_energy(e2, switches.get_which_espread2(), switches.get_espread2(), *rndm_generator_);
    }

  if (switches.get_do_isr()) isr2(e1,e2,&e1,&e2, *rndm_generator_);
  
  // Modification of the ecm by Strahinja to be safer with CLIC
  float vx1, vy1, vz1, vx2, vy2, vz2;
  float cosinus;
  pointer1.velocities(vx1, vy1);
  pointer2.velocities(vx2, vy2);
  vz1 = sqrt(1. - vx1*vx1 - vy1*vy1);
  vz2 = sqrt(1. - vx2*vx2 - vy2*vy2);
  cosinus = (vx1*vx2 + vy1*vy2 - vz1*vz2); // The minus sign because the particles have opposite vz
  ecm=sqrt(2.0*e1*e2*(1. - cosinus));

  float energy1 = pointer1.energy();
  float energy2 = pointer2.energy();
  if (switches.get_do_lumi()&1)
    {
      float p1Vx, p1Vy, p2Vx, p2Vy;
      pointer1.velocities(p1Vx, p1Vy);
      pointer2.velocities(p2Vx, p2Vy);
      if ( bmt_precession ) 
      //      if (SPIN)
		{
		  lumi_heap_ee_.lumi_store_ee(mesh_, cellx, celly,min_z, energy1, p1Vx, p1Vy, energy2, p2Vx, p2Vy, weight, pointer1.getSpin(),pointer2.getSpin() , time_counter);
		}
      else lumi_heap_ee_.lumi_store_ee(mesh_, cellx, celly,min_z, energy1, p1Vx, p1Vy, energy2, p2Vx, p2Vy, weight,time_counter);
    }

  if (energy1 < 0.0) j1=1;
  if (energy2 < 0.0) j2=0;
  results_.add_lumi(j1, j2, weight);

  if (energy1*energy2 > 0.0)
    {
      if (switches.get_do_cross())
		{
		  cross_->cross_add(e1,e2,weight);
		}
      results_.add_lumi_ee(weight);
      results_.add_step_lumi_ee(weight);

      if ( bmt_precession ) 
	//     if (SPIN)
		{
		  float spin1= pointer1.getSpin()(2);
		  float spin2= pointer2.getSpin()(2);
		  help=0.5*(1.0+spin1*spin2);
		  results_.add_lumi_pp(weight*help);
		}

      if (ecm > switches.get_ecm_min()) {
	results_.add_lumi_ee_high(weight);
        results_.add_step_lumi_ee_high(weight);
      }
      help= weight*ecm;
      results_.add_lumi_ecm(help);
      results_.add_lumi_ecm2(help*ecm);
    
    }
  
  // here the bhabhas

  if (switches.get_do_bhabhas())
    {
      float part1Vx, part1Vy, part2Vx, part2Vy;
      pointer1.velocities(part1Vx, part1Vy);
      pointer2.velocities(part2Vx, part2Vy);
      // The following line changed by SL to use signed energies (needed in boost_bhabha())
	bhabhas_.make_bhabha(secondaries, part1Vx, part1Vy, part2Vx, part2Vy, e1, e2, ecm, weight, mesh_, cellx, celly,min_z, switches, *rndm_generator_, beamslice1, beamslice2 );

// The following alternative tested whether it matters which beam moves in which direction
//      bhabhas_.make_bhabha(secondaries, part1Vx, part1Vy, part2Vx, part2Vy, -e1, e2, ecm, weight, mesh_, cellx, celly,min_z, switches, *rndm_generator_, beamslice1, beamslice2);
    }

  if (switches.get_do_jets())
    {
      minijets_->mkjll_(secondaries.get_pair_parameters(),pointer1.energy(),pointer2.energy(), weight, switches, *rndm_generator_);

    }
  if (switches.get_do_prod()==1)
    {

      std::cerr << " GRID::collide_ee : do_prod not implemented! " << std::endl;
      exit(0);
    }  
}


void GRID::isr2(float e1,float e2,float *e1p,float *e2p, RNDM& rndm_generator)
{
    double c=2.0*ALPHA_EM/PI,emass2=EMASS*EMASS;
    double x,s,beta,corr,tmp;

    s=4.0*e1*e2;
    beta=c*(log(s/emass2)-1.0);
    do{
	x=pow(1.0-rndm_generator.rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(rndm_generator.rndm()>corr);
    *e1p=e1*(1.0-x);
    do{
	x=pow(1.0-rndm_generator.rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(rndm_generator.rndm()>corr);
    *e2p=e2*(1.0-x);
}


/*! This routine moves the particles one timestep. */

void GRID::move_particles(const std::vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice,int interpolation,int do_beamstrahlung,int do_trident, int sokolov, float emin,int /*do_prod*/, int extra_grids, float charge_sign, int bmt_rotate)
{
  float vx,vy,xpos,ypos;

  // ultrarelivistic particles 
  //      Bfield *= 1.0/Cvelocity;
  // electric field: GV/nm
  // for the magnetic field, the product c.B
  // also in GV/nm
  //      Efield = EBfield;

  unsigned int k,i;
  const PHI_FLOAT *phi;
  //  BEAM* beam;
  if (i_beam==1) 
    {
      //     beam = slice_of_beam_[0].get_beam();
      //      phi=phi2_;
      phi= champ_.get_phi(2);
    }
  else 
    {
      //     beam = slice_of_beam_[1].get_beam();
      //     phi=phi1_;
      phi= champ_.get_phi(1);
    }
  const std::vector<PARTICLE*>& theParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);

  int charge_label = slice_of_beam_[i_beam-1].get_beam()->sign_label(); //get charge of the beam


  switch (interpolation)
    {
    case 1:
      std::cerr << " GRID::move_particles interpolation = 1, a traiter " << std::endl;
      exit(0);
      break;
    case 2:
      {
	TRIVECTOR EBfield;
	TRIDVECTOR Efield, Bfield;
	float dzOnRadius,oldener;
	////////////////////////////////////////////////
	// in case of do_beamstrahlung, see the implementation in the case of (do_prod>1)&&(i_beam==1) 
	///////////////////////////////////////////////
	PARTICLE* particleIterator;
	for (k = 0; k < theParticles.size(); k++)
	  {
	    particleIterator = theParticles[k];
	    oldener=particleIterator->getEnergy();

 //Storing the position of the beam particle into a txt 

std::ofstream outfile("beam_particle_position.txt", std::ios_base::app); // Open the file in append mode
if (outfile.is_open())
{
    //for (k = 0; k < theParticles.size(); k++)
    //{
        //particleIterator = theParticles[k];
        //oldener = particleIterator->getEnergy();

        // Access x, y, z positions using the particleIterator
        float x_position = particleIterator->getXpos();
        float y_position = particleIterator->getYpos();
        float z_position = particleIterator->getZpos();
	
        // Store x, y, z positions in the file
        outfile << x_position << "\t" << y_position << "\t" << z_position << "\t" << charge_label << std::endl;

    outfile.close(); // Close the file after writing
}

else
{
    std::cerr << "Unable to create/open file 'beam_particle_position.txt'" << std::endl;
    exit(1);
}



	    EBfield = EBfieldOnParticle( particleIterator, grids, phi, i_beam, extra_grids);


	    dzOnRadius = particleIterator->advanceDueToEBfield(EBfield, step_,slice_of_beam_[i_beam-1].get_scal_step());

	    Efield = EBfield;
	    Bfield = TRIDVECTOR(EBfield(1), -EBfield(0), 0.0);  
    //Storing field values


    // Create or open 'field.txt' for writing
    //std::cout << "Charge Label: " << charge_label << std::endl;

std::ofstream fieldFile("field.txt", std::ios::app);  // Use a different name

if (!fieldFile) {
    std::cerr << "Error opening field.txt for writing!" << std::endl;
    return;  // Just return without a value
}

// Write the 6 components into 6 columns using (0) notation
//Only write 4 components. z components are zero
fieldFile << Efield(0) << "\t"
          << Efield(1) << "\t"
          //<< Efield(2) << "\t"
          << Bfield(0) << "\t"
          //<< Bfield(1) << "\t"
          << Bfield(1) << std::endl;

// Close the file
fieldFile.close();


	    if (bmt_rotate)
	      {
		particleIterator->rotateBMT(Efield, Bfield, charge_sign, step_);
		if (do_beamstrahlung)
		  { 
		    std::vector<float> photonEnergies;
		    results_.updateUpsmax(particleIterator->getUps());
            std::cout << "Case 1" << std::endl; // Print statement
		    if (sokolov) 
		      {
			particleIterator->beamstrahlungSokolov(Efield, Bfield, dzOnRadius, emin, charge_sign,photonEnergies, *rndm_generator_);
		      }
		    else 
		      {
			particleIterator->beamstrahlung(Efield, Bfield, dzOnRadius,emin,photonEnergies, *rndm_generator_);
		      }
		    registerPhotons(photonEnergies, *(particleIterator), i_beam, i_slice);    
		  }
	      }
	    else
	      {
		if(do_beamstrahlung)
		  {
		    std::vector<float> photonEnergies;
		    results_.updateUpsmax(particleIterator->getUps());
            //std::cout << "Case 2" << std::endl; 
		    // 		    particleBeamstrahlung(particleIterator,Efield, Bfield, dzOnRadius, emin, photonEnergies);
		    particleIterator->beamstrahlung(Efield, Bfield, dzOnRadius, emin, photonEnergies, *rndm_generator_);
		    registerPhotons(photonEnergies, *(particleIterator), i_beam, i_slice);    
		  }
	      }
	    if(do_trident)
	      {
		std::vector<float> electrons,positrons,virt;
		if(do_beamstrahlung)
		  {
		    particleIterator->setUps(particleIterator->getUps()*particleIterator->getEnergy()/oldener);
		  }
		particleIterator->createTridents(step_*slice_of_beam_[i_beam-1].get_scal_step(),&electrons,&positrons,&virt,*rndm_generator_);
		for(i=0;i<electrons.size();i++)
		  {
		    particleIterator->XYposition(xpos,ypos);
		    particleIterator->velocities(vx,vy);
		    store_trident_particle(*slice_of_beam_[i_beam-1].get_beam(),electrons[i],vx,vy,xpos,ypos,i_slice);
		    store_trident_particle(*slice_of_beam_[i_beam-1].get_beam(),positrons[i],vx,vy,xpos,ypos,i_slice);
		  }
	      }
	  }
	break;
      }
    default:
      std::cerr << " GRID::move_particles : unknown value of interpolation type :  " << interpolation << std::endl;
      exit(0);
    }
}


void GRID::move_pairs(const std::vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& rndm_generator)
{
  float stepLocal,d;
  int n_pair_steps;
  double mass=pairs.get_pair_parameters().get_mass();

  std::vector<PAIR_PARTICLE>& the_pairs = pairs.get_pairs(i_slice);
  //  std::list<PAIR_PARTICLE>::iterator itr;
  unsigned int k;
  //  for (itr = the_pairs.begin(); itr !=  the_pairs.end(); itr++)
  for (k = 0; k <  the_pairs.size(); k++)
    {      
      //      d=sqrt((rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/fabs(itr->energy()));
      d=sqrt( (rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/the_pairs[k].unsigned_energy() );
      n_pair_steps=(int)d+1;
      pairs.add_pair_steps(n_pair_steps);
      stepLocal=step_/(float)n_pair_steps;
      step_pair_1(grids,the_pairs[k],mass,stepLocal,n_pair_steps, extra_grids, charge_sign_0, rndm_generator);
    }
}

void GRID::move_pairs_tertphot(const std::vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& rndm_generator)
{
  float stepLocal,d;
  int n_pair_steps;
  double mass=pairs.get_pair_parameters().get_mass();

  std::vector<PAIR_PARTICLE>& the_pairs = pairs.get_pairs(i_slice);
  //  std::list<PAIR_PARTICLE>::iterator itr;
  unsigned int k;
  //  for (itr = the_pairs.begin(); itr !=  the_pairs.end(); itr++)
  for (k = 0; k <  the_pairs.size(); k++)
    {      
      //      d=sqrt((rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/fabs(itr->energy()));
      d=sqrt( (rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/the_pairs[k].unsigned_energy() );
      n_pair_steps=(int)d+1;
      pairs.add_pair_steps(n_pair_steps);
      stepLocal=step_/(float)n_pair_steps;
      step_pair_1_tertphot(grids,the_pairs[k],mass,stepLocal,n_pair_steps, extra_grids, charge_sign_0, rndm_generator);
    }
}

/*
void GRID::deltaVelocityFromFieldCIC(float xpart,float ypart, float energy, PHI_FLOAT *phi, float distance, float& deltavx, float& deltavy)
{
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;

  interpolePotential(xpart, ypart, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);

  deltavx = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_/energy;
  deltavy = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_/energy;
  //  std::cout << " :deltaVelocityFromFieldCIC Ex/g= " << deltavx << std::endl;
      deltavx = -deltavx*distance;
      deltavy = -deltavy*distance;
}
*/

TRIVECTOR GRID::electric_field_out_of_main_grid(const std::vector<GENERAL_GRID*>& grids,int beam,PHI_FLOAT x,PHI_FLOAT y, int extra_grids) const
{
  int k;
  PHI_FLOAT tmp,dx,dy;
  const PHI_FLOAT *phi;
  // float temp;
  float ax, ay;
  TRIVECTOR EBfield;
  ax = 0.0;
  ay = 0.0;
  for(k = 0;k <= extra_grids;k++)
    {
      if (grids[k]->coordinatesInGridRange(x,y))
	{
	  if (beam == 1 ) phi = grids[k]->get_phi2();
	  else phi = grids[k]->get_phi1();

	  // electric field in V/m
	  EBfield = grids[k]->ElectricFieldCIC(x, y, phi);
	  return EBfield;
	}
    }
  if (extra_grids<2) 
    {
      EBfield.setComponents(0.0, 0.0, 0.0);
      return EBfield;
    }
  // particle is not in largest grid 
  if (beam==1)
    {
      dx=(x-grids[0]->get_rho_x_2());
      dy=(y-grids[0]->get_rho_y_2());
      tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_2()/(dx*dx+dy*dy);
      ax=dx*tmp;
      ay=dy*tmp;
    }
  else
    {
      dx=(x-grids[0]->get_rho_x_1());
      dy=(y-grids[0]->get_rho_y_1());
      tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_1()/(dx*dx+dy*dy);
      ax=dx*tmp;
      ay=dy*tmp;
    }

      EBfield.setComponents(0.5*ax, 0.5*ay, 0.0);
      return EBfield;
}


void GRID::registerPhotons(const std::vector<float>& photonEnergies, PARTICLE& particle, int i_beam, int i_slice)
{
  unsigned int k;
  BEAM* beam = slice_of_beam_[i_beam-1].get_beam();
  //  int  i_beam = beam.label();
  //  std::cout << " GRID::createPhotons : nb photons recuperes " << photonEnergies.size() << std::endl;
  for (k=0; k < photonEnergies.size() ;k++)
    {
      
      //      const PHOTON& foton = 
      beam->new_photon(photonEnergies[k],particle, i_slice);

      //     energy -= photonEnergies[k];
      //      const PHOTON& foton = beam->new_photon(photonEnergies[k],particle, i_slice);
//       if (photon_file_ != NULL) 
// 	{
// 	  if (i_beam != 1) 
// 	    {
// 	      PHOTON temp(foton);
// 	      float energyTemp = temp.energy();
// 	      temp.setEnergy(-energyTemp);
// 	      photon_file_->save_object_on_persistent_file(&temp);
// 	    }
// 	  else photon_file_->save_object_on_persistent_file(&foton);	  
// 	}
    }
}


void GRID::beamstrahlungSingleCoherentParticle(PARTICLE* particle, TRIVECTOR /*EBfield*/, float dzOnRadius, const std::vector<GENERAL_GRID*>& /*grids*/, int i_beam,  int i_slice, const PHI_FLOAT */*phi*/, float /*distance*/, float emin,int do_prod, int /*extra_grids*/, float /*charge_sign*/)

{	

  float upsilon = particle->getUps();


  float energy = 0.0;
  float initialEnergy = fabs(particle->energy());
  results_.updateUpsmax(upsilon);
  std::cout << "Upsilon value: " << upsilon << std::endl;
  std::vector<float> photonEnergies;
  TRIDVECTOR ev1, ev2, ev3; 
  energy = PHYSTOOLS::synrad_no_spin_flip(upsilon, initialEnergy, dzOnRadius,photonEnergies, *rndm_generator_);
  registerPhotons(photonEnergies, *(particle), i_beam, i_slice);
  std::cout << "Calling beamstrahlungSingleCoherentParticle()" << std::endl; 
  if (energy < emin)
    {
      std::cout << "PARTICLE:: e_low2: " << energy << std::endl;
      energy = emin;
    }
  if ((do_prod>1)&&( i_beam ==1))
    {
      std::cerr << " GRID:: do_prod not implemented "  << std::endl;
      exit(0);
    }
  if (particle->energy() < 0.0)
    {
      particle->setEnergy( -energy);
    }
  else
    {
      particle->setEnergy(energy);
    }	  
    
}

void GRID::move_coherent_particles(const std::vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice, int interpolation, int do_beamstrahlung, float emin, int do_prod, int extra_grids,float charge_sign)
{
  unsigned int k; 
  //  float initialEnergyForLoss;
  const PHI_FLOAT *phi;

  //  BEAM* beam;

  if (i_beam==1) 
    {
      //     beam = slice_of_beam_[0].get_beam();
      //      phi=phi2_;  
      phi = champ_.get_phi(2);
    }
  else
    { 
      //     beam = slice_of_beam_[1].get_beam();
      //     phi=phi1_;
      phi = champ_.get_phi(1);
    }
  const std::vector<PARTICLE*>& coherent = slice_of_beam_[i_beam-1].get_beam()->getCoherentVector(i_slice);

  switch (interpolation)
    {
    case 1:
      std::cerr << " GRID::move_coherent_particles interpolation = 1, a traiter " << std::endl;
      exit(0);
      break;
    case 2:
      {
	for (k = 0; k < coherent.size(); k++)
	  {

	    TRIVECTOR EBfield = EBfieldOnParticle(coherent[k], grids, phi, i_beam, extra_grids);
	    float dzOnRadius = coherent[k]->advanceDueToEBfield(EBfield, step_,  slice_of_beam_[i_beam-1].get_scal_step());
	    // (no spin rotation for coherent particles, nor sokolov ternov spin flip)
	    if(do_beamstrahlung)
	      {
		beamstrahlungSingleCoherentParticle(coherent[k], EBfield, dzOnRadius,grids, i_beam, i_slice, phi, step_, emin, do_prod, extra_grids, charge_sign);	
	      }  
	  }      
	break;
      }
    default:
      std::cerr << " GRID::move_coherent_particles : unknown value of interpolation type :  " << interpolation << std::endl;
      exit(0);
    }
}

void GRID::move_trident_particles(const std::vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice, int interpolation, int do_beamstrahlung, float emin, int do_prod, int extra_grids,float charge_sign)
{
  unsigned int k; 
  const PHI_FLOAT *phi;
  if (i_beam==1) 
    {
      phi = champ_.get_phi(2);
    }
  else
    { 
      phi = champ_.get_phi(1);
    }
  const std::vector<PARTICLE*>& trident = slice_of_beam_[i_beam-1].get_beam()->getTridentVector(i_slice);

  switch (interpolation)
    {
    case 1:
      std::cerr << " GRID::move_trident_particles interpolation = 1, a traiter " << std::endl;
      exit(0);
      break;
    case 2:
      {
	for (k = 0; k < trident.size(); k++)
	  {
	    TRIVECTOR EBfield = EBfieldOnParticle(trident[k], grids, phi, i_beam, extra_grids);
	    float dzOnRadius = trident[k]->advanceDueToEBfield(EBfield, step_, slice_of_beam_[i_beam-1].get_scal_step());
	    /* This routine includer transverse momentum rotation due to non-transversality if the fields */
	    //	    float dzOnRadius = trident[k]->advanceTridentDueToEBfield(EBfield, step_, slice_of_beam_[i_beam-1].get_scal_step());
	    // (no spin rotation for coherent particles, nor sokolov ternov spin flip)
	    if(do_beamstrahlung)
	      {
		beamstrahlungSingleCoherentParticle(trident[k], EBfield, dzOnRadius,grids, i_beam, i_slice, phi, step_, emin, do_prod, extra_grids, charge_sign);	
	      }  
	  }      
	break;
      }
    default:
      std::cerr << " GRID::move_trident_particles : unknown value of interpolation type :  " << interpolation << std::endl;
      exit(0);
    }
}


/* Distributes the photons for background calculation */

void GRID::distribute_photons(int slice_1,int slice_2, float photon_ratio, RNDM& rndm_generator)
{
  //int k;
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  //  std::cout << " GRID::distribute_photons " << std::endl;

  ratio = photon_ratio;
  clear_photon_pointer();
 
  if (ratio<eps) return;

  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);

  distributePhotonInBeam(grid_photon1_, 1, slice_1, ratio, ratio_i_1, rndm_generator);
  distributePhotonInBeam(grid_photon2_, 2, slice_2, ratio, ratio_i_2, rndm_generator);

}

//void GRID::distributePhotonInBeam(   std::vector< std::list<BEAM_PARTICLE_POINTER> >& grid_photon, BEAM* beam, int slice,float ratio, float ratio_i, RNDM& rndm_generator)
void GRID::distributePhotonInBeam(   std::vector< std::list<BEAM_PHOTON_POINTER> >& grid_photon,int i_beam, int slice,float ratio, float ratio_i, RNDM& rndm_generator)
{
  float xphot, yphot;
  int i1, i2, j;
   

  float temp;

  const std::vector<PHOTON>& thePhotons =  slice_of_beam_[i_beam-1].get_beam()->getPhotonVector(slice);
  unsigned int k;
  for(k=0; k < thePhotons.size(); k++)
    {
      temp = rndm_generator.rndm();
      if (temp<ratio)
	{
	  thePhotons[k].XYposition(xphot, yphot);

	  if (particleInGrid(xphot, yphot, i1, i2))
	    {
	      j = i1*n_cell_y_ + i2;
	      grid_photon[j].push_back(BEAM_PHOTON_POINTER(&thePhotons[k],ratio_i ));

		
	    }
	}
    }
}

/*! This routine calculates the electron-photon, positron-photon and
   photon-photon luminosities */

void GRID::photon_lumi(float min_z, SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator)
{

  int i1,i2,j,n_x,n_y;
  //  std::list<BEAM_PARTICLE_POINTER>::iterator photon_pointer, photon_pointer1,photon_pointer2;
  std::list<BEAM_PHOTON_POINTER>::const_iterator photon_pointer; 
  std::list<BEAM_PHOTON_POINTER>::const_iterator photon_pointer1,photon_pointer2;
  std::list<BEAM_PARTICLE_POINTER>::const_iterator particle_pointer;
  n_x = n_cell_x_;
  n_y = n_cell_y_;

  const PAIR_PARAMETER& pair_parameter = secondaries.get_pair_parameters();

  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;
	  for (photon_pointer=grid_photon1_[j].begin();photon_pointer != grid_photon1_[j].end(); photon_pointer++)
	    {
	      for (particle_pointer = part_pointer2_[j].begin(); particle_pointer != part_pointer2_[j].end();  particle_pointer++)
		{
		  collide_ge(i1, i2, min_z,*photon_pointer,*particle_pointer, secondaries, switches, pair_parameter, rndm_generator);
		}
	    }
	  
	  for (photon_pointer=grid_photon2_[j].begin(); photon_pointer!= grid_photon2_[j].end(); photon_pointer++)
	    {
	      for (particle_pointer= part_pointer1_[j].begin(); particle_pointer != part_pointer1_[j].end(); particle_pointer++)
		{
		  collide_eg( i1, i2, min_z,*particle_pointer, *photon_pointer, secondaries, switches,pair_parameter,rndm_generator);		  
		}
	    }
	  
	  
	  double returnCrossSectionToAdd[3];
	  double tempor[3];
	  int k;
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k] = 0.0;
	  }
	  
	  for ( photon_pointer1 = grid_photon1_[j].begin(); photon_pointer1 != grid_photon1_[j].end(); photon_pointer1++)
	    {
	      for( photon_pointer2 = grid_photon2_[j].begin(); photon_pointer2 != grid_photon2_[j].end(); photon_pointer2++)
		{
		  collide_gg(i1, i2, min_z,*photon_pointer1, *photon_pointer2, secondaries, muons, switches, rndm_generator, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	      
	    }	  
	  results_. cumulate_hadrons_gg(tempor);
	}
    }
}


void GRID::collide_ge(int cellx, int celly,float min_z, const BEAM_PHOTON_POINTER& photon_pointer, const BEAM_PARTICLE_POINTER& particle_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& rndm_generator)
{
   float photonEnergy = photon_pointer.energy();
  if (photonEnergy <= 0.0) return; 

  float particleEnergy = particle_pointer.energy(); 
  int j2=1;
  if (particleEnergy <0.0 ) j2=0;

   float weight = photon_pointer.weight()*particle_pointer.weight(); 

  results_.add_lumi(2, j2, weight);

  if (switches.get_do_lumi()&2)
    {
      lumi_heap_ge_.lumi_store(mesh_, cellx, celly,min_z, photonEnergy,particleEnergy, weight);
    }

  results_.add_lumi_ge(weight);

  if (switches.get_do_jets())
    minijets_->mkjbh1_(pair_parameter, photonEnergy,particleEnergy,weight, switches, rndm_generator);
  // previously collide_ge_2
  collide_compton(0, cellx, celly,min_z, photon_pointer, particle_pointer, secondaries, switches,rndm_generator, 1);		    
}

void GRID::collide_eg(int cellx, int celly,float min_z, const BEAM_PARTICLE_POINTER& particle_pointer, const BEAM_PHOTON_POINTER& photon_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& rndm_generator)
{
  float photonEnergy = photon_pointer.energy();
  if (photonEnergy <= 0.0) return;

  float particleEnergy = particle_pointer.energy();  
  int j1=0;
  if (particleEnergy < 0.0) j1=1;
 
  float weight = photon_pointer.weight()*particle_pointer.weight();

  results_.add_lumi(j1, 2, weight);

  if (switches.get_do_lumi()&2)
    {
      lumi_heap_eg_.lumi_store(mesh_, cellx, celly,min_z,particleEnergy,photonEnergy, weight);
    }

  results_.add_lumi_eg(weight);

  if (switches.get_do_jets())
    minijets_->mkjbh2_(pair_parameter, particleEnergy,photonEnergy,weight, switches, rndm_generator);
  // previously collide_eg_2 
  collide_compton(0,cellx, celly, min_z, photon_pointer, particle_pointer, secondaries, switches, rndm_generator, 2);		    
}

void GRID::collide_gg(int cellx, int celly, float min_z, const BEAM_PHOTON_POINTER& photon_pointer1, const BEAM_PHOTON_POINTER& photon_pointer2, PAIR_BEAM& secondaries ,PAIR_BEAM& muons, SWITCHES& switches, RNDM& rndm_generator, double *returnCrossSectionToAdd)
{
  JET_FLOAT ecm;
  //double beta_x,beta_y;
  // float particleVx1,particleVy1,particleVx2,particleVy2 ;
  float photonEnergy1 = photon_pointer1.energy();
  float photonEnergy2 = photon_pointer2.energy();
  float weight  = photon_pointer1.weight()*photon_pointer2.weight();

  ecm=sqrt(4.0*photonEnergy1*photonEnergy2);

  results_.add_lumi(2, 2, weight);
  if (ecm<=0.0) {
    return;
  }

  if (switches.get_do_lumi()&4){

    lumi_heap_gg_.lumi_store(mesh_, cellx, celly,min_z,photonEnergy1, photonEnergy2, weight);
  }
  results_.add_lumi_gg(weight);
  if (ecm*ecm > switches.get_gg_smin())
    {
      results_.add_lumi_gg_high(weight);
    }
  if (switches.get_do_jets())
    minijets_->mkjbw_(photonEnergy1,photonEnergy2,weight, switches, rndm_generator );

  collide_gg_XX(0,cellx, celly, min_z,photon_pointer1, photon_pointer2, switches, secondaries, muons, rndm_generator, returnCrossSectionToAdd);
}

//void GRID::collide_gg_XX(int index_of_process, int cellx, int celly,float min_z, const PARTICLE_POINTER& photon1, const PARTICLE_POINTER& photon2, SWITCHES& switches, PAIR_BEAM& secondaries, RNDM& rndm_generator, double returnCrossSectionToAdd[3])
void GRID::collide_gg_XX(int index_of_process, int cellx, int celly,float min_z, const PHOTON_POINTER& photon1, const PHOTON_POINTER& photon2, SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator, double *returnCrossSectionToAdd)
{
  double beta_x,beta_y;
  float particle1Vx,particle1Vy,particle2Vx,particle2Vy ;
  float phot1_q2, phot1_eorg;
  float phot2_q2, phot2_eorg;
  float weight;
  float energy1 = photon1.energy();
  float energy2 = photon2.energy();
  if (energy1*energy2 <= 0.0) return;

  photon1.velocities(particle1Vx,particle1Vy);
  photon2.velocities(particle2Vx,particle2Vy);
  beta_x=(particle1Vx*energy1 + particle2Vx*energy2)/(energy1+energy2);
  beta_y=(particle1Vy*energy1 + particle2Vy*energy2)/(energy1+energy2);

  phot1_q2 =  photon1.q2();
  phot2_q2 =  photon2.q2();
  phot1_eorg = photon1.eorg();
  phot2_eorg = photon2.eorg();

  weight  = photon1.weight()*photon2.weight();

  if (switches.get_do_pairs())
    {
      secondaries. make_pair_bw(mesh_, cellx, celly,min_z,index_of_process, energy1,phot1_q2,phot1_eorg,energy2,phot2_q2,phot2_eorg, weight,beta_x,beta_y, switches,rndm_generator);
    }
  if (switches.get_do_muons())
    {
      muons.make_muon(mesh_, cellx, celly,min_z, index_of_process, energy1,phot1_q2,energy2,phot2_q2,weight,beta_x,beta_y, switches, rndm_generator);
    }
  if (switches.get_do_hadrons()) make_hadrons_gg2(min_z, energy1,phot1_q2,energy2,phot2_q2,weight, switches, returnCrossSectionToAdd, rndm_generator);
}

/**********************************************************************/
/* Routines for the creation and storage of hadrons                   */
/**********************************************************************/

void GRID::make_hadrons_gg2(float min_z, float energy1,float q2_1,float energy2,float q2_2,float lumi, SWITCHES& switches, double *returnCrossSectionToAdd, RNDM& rndm_generator)
{
  double s,h;
  double cross_section=0.0;
  int num,k;
  double eps=0.0808,mu=0.4525;

  for (k=0; k<3; k++) returnCrossSectionToAdd[k] = 0.0;

  s=energy1*energy2;
  h=std::max(1.0,pow(s/100.0,0.43));
  if ((q2_1>h)||(q2_2>h)) {
    return;
  }
  s*=4.0;
  if (s<4.0) return;
  switch (switches.get_do_hadrons())
    {
    case 1:
      cross_section=(211.0*pow(s,eps)+297.0*pow(s,-mu))*1e-37*lumi;
      break;
    case 2:
      cross_section=lumi*200e-37*(1.0+6.3e-3*pow(log(s),2.1)+1.96*pow(s,-0.37));
      break;
    case 3:
      cross_section=(211.0*pow(s,eps)+215.0*pow(s,-mu))*1e-37*lumi;
      break;
    case 4:
      cross_section=(-0.99244+0.0989203*log(s+22865.9))*1e-34*lumi;
      break;
    }
  if (hadron_file_ != NULL)
    {
      h=cross_section*switches.get_hadron_ratio();
      num=(int)floor(h);
      h-=num;
      if(h>rndm_generator.rndm()) num++;
      saveOnHadronFile(num, energy1, energy2,min_z, rndm_generator, switches.get_do_edm4hep());
  }
  returnCrossSectionToAdd[0] = cross_section;
  if (s<25.0) return;
  returnCrossSectionToAdd[1] = cross_section;
  if (s<100.0) return;
  returnCrossSectionToAdd[2] = cross_section;
}

/*! This routine calculates the background using the virtual photons with the
   apropriate impact parameter. In the moment it just does pair creation. */

void GRID::photon_lumi_2(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& rndm_generator)
{
  int i1,i2,j,n_x,n_y;
  std::list<BEAM_PHOTON_POINTER>::iterator photon_pointer;
  //double beta_x,beta_y,e1,e2;
  n_x = n_cell_x_;
  n_y = n_cell_y_;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;

	  double returnCrossSectionToAdd[3];
	  double tempor[3];
	  int k;
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }


	  std::list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot;
	  //float extr_phot_e, extr_phot_weight, extr_phot_q2, extr_phot_eorg;
	  for (photon_pointer=grid_photon1_[j].begin(); photon_pointer != grid_photon1_[j].end(); photon_pointer++)
	    {
	      const std::list<EXTRA_PHOTON_POINTER>& extr_phot_list = extra_photon_pointer2_[j];
	      
	      for (extr_phot = extr_phot_list.begin() ; extr_phot!=extr_phot_list.end(); extr_phot++)
		{

		  collide_gg_XX(1, i1,i2,min_z, *photon_pointer, *extr_phot,switches, secondaries, muons, rndm_generator, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }

	  results_. cumulate_hadrons_ge(tempor);
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }
	  for( photon_pointer=grid_photon2_[j].begin(); photon_pointer!= grid_photon2_[j].end(); photon_pointer++)
	    {
	      const std::list<EXTRA_PHOTON_POINTER>& extr_phot_list = extra_photon_pointer1_[j];
	      for ( extr_phot = extr_phot_list.begin(); extr_phot != extr_phot_list.end(); extr_phot++)
		{
		  collide_gg_XX(1,i1, i2,min_z,  *extr_phot, *photon_pointer,switches,secondaries, muons, rndm_generator, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }
	  results_. cumulate_hadrons_eg(tempor);
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }

	  const std::list<EXTRA_PHOTON_POINTER>& extr_phot1_list = extra_photon_pointer1_[j];
	  std::list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot1;
	  std::list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot2;
	  for( extr_phot1 = extr_phot1_list.begin(); extr_phot1!=extr_phot1_list.end(); extr_phot1++)
	    {

	      const std::list<EXTRA_PHOTON_POINTER>& extr_phot2_list = extra_photon_pointer2_[j];

	      for ( extr_phot2 = extr_phot2_list.begin(); extr_phot2!=extr_phot2_list.end(); extr_phot2++)
		{

		  collide_gg_XX(2, i1, i2,min_z, *extr_phot1, *extr_phot2, switches,secondaries, muons, rndm_generator, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }
	  results_. cumulate_hadrons_ee(tempor);
	}
    }
}


void GRID::photon_lumi_3(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries,  RNDM& rndm_generator)
{
  int i1,i2,j,n_x,n_y;

  std::list<BEAM_PARTICLE_POINTER>::iterator particle_pointer;

  std::list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot1, extr_phot2;

  n_x = n_cell_x_;
  n_y = n_cell_y_;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;
	  const std::list<EXTRA_PHOTON_POINTER>& extr_phot1_list = extra_photon_pointer1_[j];

	  for( extr_phot1 = extr_phot1_list.begin(); extr_phot1!=extr_phot1_list.end(); extr_phot1++)
	    {
	      for ( particle_pointer = part_pointer2_[j].begin(); particle_pointer != part_pointer2_[j].end(); particle_pointer++)
		{
		  collide_compton(1,i1, i2,min_z, *extr_phot1, *particle_pointer, secondaries, switches, rndm_generator, 1);
		}
	    }
	  
	  const std::list<EXTRA_PHOTON_POINTER>& extr_phot2_list = extra_photon_pointer2_[j];

	  for (extr_phot2 = extr_phot2_list.begin(); extr_phot2!=extr_phot2_list.end(); extr_phot2++)
	    {
	      for (particle_pointer = part_pointer1_[j].begin(); particle_pointer != part_pointer1_[j].end(); particle_pointer++)
		{

		  collide_compton(1,i1, i2,min_z, *extr_phot2, *particle_pointer, secondaries, switches, rndm_generator,2);
		}
	    }
	}
    }
}


void GRID::move_photons2(BEAM& beam,int ibeam,int i_slice, RNDM& rndm_generator)
{
  float ax,ay,wgt;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  float radius_i;
  const PHI_FLOAT *phi;
  PHI_FLOAT upsilon,upsilon0,tmp;
  PHI_FLOAT x,y;
  int i, nd;


  float scal_step_local;
  //  scal_step_local = scal_step[ibeam-1];
  scal_step_local = slice_of_beam_[ibeam-1].get_scal_step();
  if (ibeam==1){
    //    phi=phi2_;
    phi = champ_.get_phi(2);
    wgt= slice_of_beam_[0].get_nb_part_per_macro();
  }
  else{
    //   phi=phi1_;
    phi = champ_.get_phi(1);
    wgt= slice_of_beam_[1].get_nb_part_per_macro();
  }
  std::vector<PHOTON>& thePhotons =  beam.getPhotonVector(i_slice);
  unsigned int k;
  for(k=0; k < thePhotons.size(); k++)
    {
      thePhotons[k].XYposition(x, y);
      if ( coordinatesInGridRange(x, y) )
	{ 
	  interpolePotential(x, y, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);
	  ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_;
	  ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_;
	}
      else
	{
	  ax = 0.0;
	  ay = 0.0;
	}
      float factor = step_*scal_step_local;
      //      photItr->move(factor);
      thePhotons[k].advancePosition(factor);

      radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step_local;
      
      upsilon0=LAMBDA_BAR/EMASS*radius_i;
      //     upsilon=upsilon0*photItr->get_energy()/EMASS;
      upsilon=upsilon0*thePhotons[k].energy()/EMASS;
      
      if (upsilon>1e-10)
	{
	  tmp=ALPHA_EM*ALPHA_EM/RE*step_*1e-9*scal_step_local*upsilon0*0.23*exp(-8.0/(3.0*upsilon));
	  tmp*=pow(1.0+0.22*upsilon,-1.0/3.0);
	  coherent_results_.updateProbmax(tmp);
	  tmp*=1.36;
	  if (rndm_generator.rndm()<tmp) {
	    if (tmp<0.1) {
	      if (coherent_generate(beam,i_slice,upsilon,thePhotons[k], rndm_generator)) {
		coherent_results_.updateSumeng(wgt*thePhotons[k].energy());
		thePhotons[k].razEnergy();//photon->energy=0.0;
		tmp=1.0;
		//nyes++;
	      }
	      else {
		tmp=0.0;
		//nno++;
	      }
	    }
	    else {
	      nd=(int)(tmp/0.1+1.0);
	      float ds=tmp/nd;
	      tmp=0.0;
	      for (i=0;i<nd;++i) {
		if (rndm_generator.rndm()<ds) {
		  if (coherent_generate(beam,i_slice,upsilon,thePhotons[k], rndm_generator)) {
		    coherent_results_.updateSumeng(wgt*thePhotons[k].energy());
		    thePhotons[k].razEnergy();//photon->energy=0.0;
		    tmp=1.0;
		    //nyes++;
		    i=nd;
		  }
		}
	      }
	    }
	  }
	  else {
	    tmp=0.0;
	  }
	  /*
	    if (tmp>1.0)
	    {
	    //	      coherent_results.updateSumeng(wgt*photItr->get_energy());
	    coherent_results_.updateSumeng(wgt*thePhotons[k].energy());
	    coherent_generate(beam,i_slice,upsilon,thePhotons[k], rndm_generator);
	    //	      photItr->razEnergy();
	    thePhotons[k].razEnergy();
	    tmp=1.0;
	    }
	    else
	    {
	    if (rndm_generator.rndm()<tmp*1.36) 
	    {
	    coherent_generate(beam,i_slice,upsilon,thePhotons[k], rndm_generator);
	    //		  photItr->razEnergy();
	    thePhotons[k].razEnergy();
	    }
	    }
	    // end test
	  */ 
	  tmp*=wgt;
	  //	  coherent_results.updateSums(tmp, photItr->get_energy(), upsilon0);
	  coherent_results_.updateSums(tmp, thePhotons[k].energy(), upsilon0);
	}
      //      photItr++;
    }
}




//void GRID::coherent_generate(BEAM& beam,int i_slice, float upsilon,PHOTON& phot, RNDM& rndm_generator)
int GRID::coherent_generate(BEAM& beam,int i_slice, float upsilon,PHOTON& phot, RNDM& rndm_generator)
{
  float x,y;
  float vx,vy;
  float energy;
  energy =  phot.energy();
  // number adjustment to be done for total cross section
  if (pick_coherent_energy(upsilon,energy, rndm_generator))
    {
      phot.XYposition(x, y);
      phot.velocities(vx,vy);
      store_coherent_particle(beam,energy,vx,vy,x,y,i_slice);
      store_coherent_particle(beam,-(phot.energy()-energy),vx,vy,x,y,
		     i_slice);
      return 1;
    }
  return 0;
}

bool GRID::pick_coherent_energy(float ups,float& energy, RNDM& rndm_generator)
{
  float a=0.13513,eta,dxdy,x,h,hp,tmp,tmp2,y,coh;
  y=2.0*rndm_generator.rndm()-1.0;
  eta=8.0/(3.0*ups);
  tmp=fabs(y/(1.0+(1.0-y*y)/(2.0*sqrt(eta))));
  tmp2=1.0-tmp*tmp;
  h=-log(tmp2);
  hp=tmp/tmp2*(1.0+(1.0+y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)));
  if (y>0.0){
    x=0.5*(1.0+sqrt(h/(eta+h)));
  }
  else{
    x=0.5*(1.0-sqrt(h/(eta+h)));
  }
  dxdy=hp/sqrt(h*(h+eta))*(1.0-h/(h+eta));
  coh=PHYSTOOLS::u(ups);
  tmp=PHYSTOOLS::fcp(ups,x)*a/coh*dxdy;
  if (rndm_generator.rndm()<tmp){
     energy *= x;
    return true;
  } else {
    return false;
  }
}


// void GRID::step_pair_1X(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& rndm_generator)
// {
//   const float eps=1e-35,emass2=EMASS*EMASS;
//   float x,y,z,vx,vy,vz,e_inv2,e_inv,step_2,step_q,vold2,scal,thetamax;
//   float energy;
//   float ex,ey,bx,by,b_norm,b_norm_i,theta,a1,a2,a3,vb,eng;
//   float ph[1000],vx0,vy0,vz0,eng0;
//   int nph;
//   int i,icharge,j;
  
//   pair.get_parameters(x,y,z,vx,vy,vz,energy);
//   if (energy >0.0)
//     {
//       eng = energy;
//       icharge=0;
//     }
//   else
//     {
//       eng=-energy;
//       icharge=1;
//     }
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
//   step_2=0.5*step;
//   step_q=step*step;
  
//   // initial half step 
  
//   field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//   if (icharge)
//     {
//       ex=-ex;
//       ey=-ey;
//       bx=-bx;
//       by=-by;
//     }
//   b_norm=sqrt(bx*bx+by*by);
//   b_norm_i=1.0/std::max(b_norm,eps);
//   bx*=b_norm_i;
//   by*=b_norm_i;
//   vb=vx*by-vy*bx;
// #ifdef PAIR_SYN
//   vx0=vx;
//   vy0=vy;
//   vz0=vz;
// #endif
//   theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//   a3=0.5*theta;
//   a1=1.0-a3;
//   theta=sqrt(theta);
//   thetamax=2.0*theta;
//   a2=theta*vz;
//   a3*=vx*bx+vy*by;
//   vz=vz*a1+theta*vb;
//   vx=vx*a1-a2*by+a3*bx;
//   vy=vy*a1+a2*bx+a3*by;
//   vold2=vx*vx+vy*vy+vz*vz;
//   vx+=step_2*ex*e_inv;
//   vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//   scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//   vx*=scal;
//   vy*=scal;
//   vz*=scal;
// #ifdef PAIR_SYN
//   vx0*=scal;
//   vy0*=scal;
//   vz0*=scal;
// #endif
//   eng/=scal;
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
// #endif
// #ifdef PAIR_SYN
//   synrad(eng,
// 	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))
// 	 /(0.5*step*1e-9),
// 	 step*1e-9,ph,&nph, rndm_generator);
//   if (nph>0) {
//     eng0=eng;
//     for (j=0;j<nph;j++){
//       eng-=ph[j];
//     }
//     scal=sqrt(((eng*eng-emass2)*eng0*eng0)
// 	      /((eng0*eng0-emass2)*eng*eng));
//     vx*=scal;
//     vy*=scal;
//     vz*=scal;
//   }
// #endif
//   /* loop over steps */
//   for(i=1;i<nbSteps;i++)
//     {
//       x+=vx*step;
//       y+=vy*step;
//       z+=vz*step;
//       field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//       if (icharge){
// 	ex=-ex;
// 	ey=-ey;
// 	bx=-bx;
// 	by=-by;
//       }
      
// #ifdef PAIR_SYN
//       vx0=vx;
//       vy0=vy;
//       vz0=vz;
// #endif
//       /* scd new */
//       vold2=vx*vx+vy*vy+vz*vz;
//       vx+=step_2*ex*e_inv;
//       vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//       scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//       vx*=scal;
//       vy*=scal;
//       vz*=scal;
// #ifdef PAIR_SYN
//       vx0*=scal;
//       vy0*=scal;
//       vz0*=scal;
// #endif
//       eng/=scal;
//       e_inv=1.0/eng;
//       e_inv2=e_inv*e_inv;
// #endif
//       b_norm=sqrt(bx*bx+by*by);
//       b_norm_i=1.0/std::max(b_norm,eps);
//       bx*=b_norm_i;
//       by*=b_norm_i;
//       vb=vx*by-vy*bx;
//       /*vb=0.0;*/
//       theta=(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//       a3=0.5*theta;
//       a1=1.0-a3;
//       theta=sqrt(theta);
//       thetamax=std::max(thetamax,theta);
//       a2=theta*vz;
//       a3*=vx*bx+vy*by;
//       /*a3=0.0;
// 	a1=1.0;*/
//       vz=vz*a1+theta*vb;
//       vx=vx*a1-a2*by+a3*bx;
//       vy=vy*a1+a2*bx+a3*by;
      
//       /* scd new */
//       vold2=vx*vx+vy*vy+vz*vz;
//       vx+=step_2*ex*e_inv;
//       vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//       scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//       vx*=scal;
//       vy*=scal;
//       vz*=scal;
// #ifdef PAIR_SYN
//       vx0*=scal;
//       vy0*=scal;
//       vz0*=scal;
// #endif
//       eng/=scal;
//       e_inv=1.0/eng;
//       e_inv2=e_inv*e_inv;
// #endif
// #ifdef PAIR_SYN
//       /* new for test */
//       synrad(eng,
// 	     sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))/step*1e9, step*1e-9,ph,&nph, rndm_generator);
//       if (nph>0) {
// 	eng0=eng;
// 	for (j=0;j<nph;j++){
// 	  eng-=ph[j];
// 	}
// 	scal=sqrt(((eng*eng-emass2)*eng0*eng0)
// 		  /((eng0*eng0-emass2)*eng*eng));
// 	vx*=scal;
// 	vy*=scal;
// 	vz*=scal;
//       }
// #endif
//     }
//   /* last half step */
//   x+=vx*step;
//   y+=vy*step;
//   z+=vz*step;
//   field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//   if (icharge){
//     ex=-ex;
//     ey=-ey;
//     bx=-bx;
//     by=-by;
//   }
  
//   /* scd new */
//   vold2=vx*vx+vy*vy+vz*vz;
//   vx+=step_2*ex*e_inv;
//   vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//   scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//   vx*=scal;
//   vy*=scal;
//   vz*=scal;
//   eng/=scal;
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
// #endif
//   b_norm=sqrt(bx*bx+by*by);
//   b_norm_i=1.0/std::max(b_norm,eps);
//   bx*=b_norm_i;
//   by*=b_norm_i;
//   vb=vx*by-vy*bx;
//   /*vb=0.0;*/
//   theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//   a3=0.5*theta;
//   a1=1.0-a3;
//   theta=sqrt(theta);
//   thetamax=std::max(thetamax,float(2.0)*theta);
//   a2=theta*vz;
//   a3*=vx*bx+vy*by;
//   /*a3=0.0;
//     a1=1.0;*/
//   vz=vz*a1+theta*vb;
//   vx=vx*a1-a2*by+a3*bx;
//   vy=vy*a1+a2*bx+a3*by;
// #ifdef SCALE_ENERGY
//   scal=sqrt((eng*eng-emass2)/((eng*eng)*(vx*vx+vy*vy+vz*vz)));
//   if (fabs(scal-1.0)>0.01)
//     printf("> %g %g %g %g %g\n",eng,eng-fabs(energy),thetamax,
// 	   sqrt(vx*vx+vy*vy+vz*vz),scal);
// #else
//   scal=1.0;
// #endif
  
//   if (icharge) energy = -eng;
//   else energy =eng; 
//   pair.set(x,y,z,vx*scal,vy*scal,vz*scal,energy);
// }
void GRID::step_pair_1(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair, double mass,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& rndm_generator)
{
  float thetamax;
  float ex,ey,bx,by;
  int i;
  // initial half step 
 
  // retrieve E and B fields
  field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = 2.0*pair.apply_initial_half_step_fields(step, mass, ex,ey, bx, by, rndm_generator);
  /* loop over steps */
  for(i=1;i<nbSteps;i++)
    {
      pair.advancePosition(step);

      // retrieve E and B fields
      field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);    
      thetamax = std::max (thetamax, pair.apply_full_step_fields(step, mass, ex,ey, bx, by, rndm_generator)); 
    }
  /* last half step */
  pair.advancePosition(step);
  field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = std::max( thetamax, (float)2.0*pair.apply_final_half_step_fields(step, mass, ex,ey, bx, by, thetamax,rndm_generator));
  if (	!pair.last_rescaling_ok() )
    {
      std::cerr << " GRID::step_pair_1() : " << " thetamax " << thetamax << std::endl;
    }	  
}

void GRID::step_pair_1_tertphot(const std::vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& pair, double mass,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& rndm_generator)
{
  float thetamax;
  float ex,ey,bx,by;
  int i;
  unsigned int j;
  std::vector<float> photon_e;
  // initial half step 
 
  // retrieve E and B fields
  field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = 2.0*pair.apply_initial_half_step_fields(step, mass, ex,ey, bx, by, &photon_e, rndm_generator);
  /* loop over steps */
  for(i=1;i<nbSteps;i++)
    {
      pair.advancePosition(step);

      // retrieve E and B fields
      field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);    
      thetamax = std::max (thetamax, pair.apply_full_step_fields(step, mass, ex,ey, bx, by, &photon_e, rndm_generator)); 
    }
  /* last half step */
  pair.advancePosition(step);
  field_pair(pair, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = std::max( thetamax, (float)2.0*pair.apply_final_half_step_fields(step, mass, ex,ey, bx, by, thetamax,rndm_generator));
  for(j=0;j<photon_e.size();j++)
    {
      tertphot_.push_back(TERTPHOTON(photon_e[j],pair));
    }
  if (	!pair.last_rescaling_ok() )
    {
      std::cerr << " GRID::step_pair_1_tertphot() : " << " thetamax " << thetamax << std::endl;
    }	  
}

void GRID::apply_magnetic_field_on_pair(float fac_theta, float step_q, float e_inv2, float bx, float by, float& vx,float& vy, float& vz, float& theta) const
{
  const float eps=1e-35;
  float b_norm, b_norm_i,a1,a2,a3,vb;

  // on normalise B
  b_norm=sqrt(bx*bx+by*by);
  b_norm_i=1.0/std::max(b_norm,eps);
  bx *= b_norm_i;
  by *= b_norm_i;

  // v x B
  vb = vx*by-vy*bx;


  theta=fac_theta*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
  a3=0.5*theta;

  // 1 - theta**2/2
  a1=1.0-a3;

  // theta
  theta=sqrt(theta);
  a2=theta*vz;

  // a3 = (theta**/2).(v.B)  [B norme]
  a3*=vx*bx+vy*by;

  // application du champ magnetique (??)
  vx=vx*a1-a2*by+a3*bx;
  vy=vy*a1+a2*bx+a3*by;
  vz=vz*a1+theta*vb;
}


int GRID::field_pair(const PAIR_PARTICLE& pair, const std::vector<GENERAL_GRID*>& grids,float& ex,float& ey, float& bx,float& by, int extra_grids, float charge_sign_0)
{
  PHI_FLOAT ax_1,ay_1,ax_2,ay_2;
  PHI_FLOAT tmp,dx,dy;
  float Ex1, Ey1,Ex2, Ey2;
  float x,y;
  int i;
  TRIVECTOR Efield1, Efield2;

  pair.XYposition(x, y);
  for(i=0;i<=extra_grids;i++)
    {
      if (grids[i]->coordinatesInGridRange(x,y))
	{
	  Efield1 = grids[i]->ElectricFieldCIC(x, y, grids[i]->get_phi1());
	  Efield2 = grids[i]->ElectricFieldCIC(x, y, grids[i]->get_phi2());

	  Efield1 *= -charge_sign_0;

	  Ex1 = Efield1(0);
	  Ey1 = Efield1(1);
	  Ex2 = Efield2(0);
	  Ey2 = Efield2(1);
	  ex = -(Ex2-Ex1);
	  ey = -(Ey2-Ey1);
	  by =  (Ex2 + Ex1);
	  bx = -(Ey2 + Ey1);

	  return i;
	}
    }
/* particle is not in largest grid */
  dx=(x-grids[0]->get_rho_x_1());
  dy=(y-grids[0]->get_rho_y_1());
  tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_1()/(dx*dx+dy*dy);
  tmp*=-charge_sign_0;

/*  tmp=grids[0].rho_factor*grids[0].rho_sum_1/(dx*dx+dy*dy);*/
  ax_1=dx*tmp;
  ay_1=dy*tmp;
  dx=(x-grids[0]->get_rho_x_2());
  dy=(y-grids[0]->get_rho_y_2());
  tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_2()/(dx*dx+dy*dy);
  ax_2=dx*tmp;
  ay_2=dy*tmp;

  ex = -0.5*(ax_2-ax_1);
  ey = -0.5*(ay_2-ay_1);
  by = 0.5*(ax_2+ax_1);
  bx = -0.5*(ay_2+ay_1);
  return -1;
}



EXTRA_GRID::~EXTRA_GRID()
{
}

//same as distribute_particles0 but without counting particles as being on or off the grid 

void EXTRA_GRID::distribute_particles(int i_slice1,
				   int i_slice2,
				   int electron_distribution_rho, 
				   int force_symmetric)
{

  slice_of_beam_[0].razRho();
  slice_of_beam_[1].razRho();

  //  razRhos(); 
  switch (electron_distribution_rho)
    {
    case 1:

//       assignBeamSliceNGP(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignBeamSliceNGP(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignBeamSliceNGP(slice_of_beam_[1], i_slice2);

      break;
    case 2:

//        assignBeamSliceCIC(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignBeamSliceCIC(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
  distribute_coherent_particles(i_slice1, i_slice2, electron_distribution_rho);
  distribute_trident_particles(i_slice1, i_slice2, electron_distribution_rho);
  if (force_symmetric)
    {
      slice_of_beam_[0].symmetrizeCharges(n_cell_x_,n_cell_y_);
      slice_of_beam_[1].symmetrizeCharges(n_cell_x_,n_cell_y_);
    }

    // symmetrizeCharges();
}


//  Distribute the coherent pair particles for the calculation of the fields
//  but does not count if they are on the grid or not (this routine is used for
//  the larger grids)
 

void EXTRA_GRID::distribute_coherent_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho)
{
  switch (electron_distribution_rho)
    {
    case 1:


//       assignCoherentBeamSliceNGP(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignCoherentBeamSliceNGP(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignCoherentBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignCoherentBeamSliceNGP(slice_of_beam_[1], i_slice2);
      break;
    case 2:


      assignCoherentBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignCoherentBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
}

void EXTRA_GRID::distribute_trident_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho)
{
  switch (electron_distribution_rho)
    {
    case 1:
      assignTridentBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignTridentBeamSliceNGP(slice_of_beam_[1], i_slice2);
      break;
    case 2:
      assignTridentBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignTridentBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
}

// void EXTRA_GRID::assignCoherentBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int i1, i2;
//   float xpart, ypart;  
//   float ch;
//   unsigned int particle;
//   const std::vector<PARTICLE*>& coherent = beam->getCoherentVector(i_slice);
//   for (particle = 0; particle < coherent.size(); particle++)
//     {
//       coherent[particle]->XYposition(xpart, ypart);
//       if (coherent[particle]->energy() < 0.0) ch=-1.0;
//       else ch=1.0;
//       if (particleInGrid(xpart, ypart,i1, i2)) assignChargeToNGP(rho, i1, i2, n_macro*ch);	    
//     }
// }
void EXTRA_GRID::assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const std::vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart,i1, i2)) sog.assignChargeToNGP(i1, i2,ch);	    
    }
}
void EXTRA_GRID::assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const std::vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart,i1, i2)) sog.assignChargeToNGP(i1, i2,ch);	    
    }
}
// void EXTRA_GRID::assignBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int k, i1, i2;
//   float xpart, ypart;
//     std::vector<PARTICLE*>& theParticles = beam->getParticleVector(i_slice);
//     for (k = 0; k < theParticles.size(); k++)
//    {
// 	theParticles[k]->XYposition(xpart, ypart);
//       if (particleInGrid(xpart, ypart, i1, i2)) assignChargeToNGP(rho, i1, i2, n_macro);
//     }
// }
void EXTRA_GRID::assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  unsigned int k;
  int i1, i2;
  float xpart, ypart;
    const std::vector<PARTICLE*>& theParticles = sog.get_beam()->getParticleVector(i_slice);
    for (k = 0; k < theParticles.size(); k++)
   {
	theParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, i1, i2)) sog.assignChargeToNGP(i1, i2,1.0);
    }
}
// void EXTRA_GRID::assignBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int k, i1, i2;
//   float h_x,h_y;
//   float xpart, ypart;
//     std::vector<PARTICLE*>& theParticles = beam->getParticleVector(i_slice);
//     for (k = 0; k < theParticles.size(); k++)
// 	{
// 	theParticles[k]->XYposition(xpart, ypart);
// 	  if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) assignChargeToCIC(rho, i1, i2, h_x, h_y, n_macro);	    	  
// 	}
// }
void EXTRA_GRID::assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  unsigned int k;
  int i1, i2;
  float h_x,h_y;
  float xpart, ypart;
  const std::vector<PARTICLE*>& theParticles = sog.get_beam()->getParticleVector(i_slice);
  for (k = 0; k < theParticles.size(); k++)
    {
      theParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, 1.0);	    	  
    }
}
// void EXTRA_GRID::assignCoherentBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int i1, i2;
//   float xpart, ypart;  
//   float ch;
//   float h_x,h_y;
//   unsigned int particle;
//   const std::vector<PARTICLE*>& coherent = beam->getCoherentVector(i_slice);
//   for (particle = 0; particle < coherent.size(); particle++)
//     {
//       coherent[particle]->XYposition(xpart, ypart);
//       if (coherent[particle]->energy() < 0.0) ch=-1.0;
//       else ch=1.0;
//       if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) assignChargeToCIC(rho, i1, i2, h_x, h_y, n_macro*ch);	
//     }
// }
void EXTRA_GRID::assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const std::vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);	
    }
}

void EXTRA_GRID::assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const std::vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);	
    }
}
