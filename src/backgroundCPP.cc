#include "backgroundCPP.h"
#include <iostream>


double COMPT::x_compt_ = 0.0;

/* px and py are beta_x and beta_y at input */
void COMPT::compt_do(const MESH& mesh, int cellx, int celly,float min_z, PAIR_BEAM& secondaries, int index_of_process, float epart,float ephot,float q2,float vx,float vy,float wgt, int dir, SWITCHES& switches, RNDM& rndm_generator)
{
  double tmp,y,scal,theta_g,theta_e,x,e,pz,pt,phi_e,s,px,py;
  int n,i;
  double eps=1e-5;

  if (q2>EMASS*EMASS) return;
  if ( index_of_process > 1) 
    {
      std::cerr << " ERROR in COMPT::compt_do : only BW or BH processes are allowed, index_of_process = " << index_of_process << std::endl;
      exit(0);
    }
  lsum_ += wgt; 
  ncall_++;
  s=4.0*ephot*epart;
  if (q2>s) return;
  if (s < switches.get_compt_x_min()*EMASS*EMASS*4.0) return;
  tmp=compt_tot(s)*wgt*switches.get_compt_scale();
  sum_ += tmp;
  n=(int)floor(tmp)+1;
  scal=tmp/n;
  x=4.0*ephot*epart/(EMASS*EMASS);
  for (i=0;i<n;i++) 
    {
      y=compt_select(s, rndm_generator);
      if (scal>eps)
	{
	  theta_g=EMASS/epart*sqrt((x-(x+1.0)*y)/y);
	  theta_e=theta_g*y/(1.0-y);
	  if (((1.0-y)*epart<switches.get_compt_emax())&&
	      ((1.0-y)*epart>switches.get_pair_ecut())) 
	    {
	      sum2_ += scal;
	      e=(1.0-y)*epart;
	      px=vx*e;
	      py=vy*e;
	      pt=theta_e*e;
	      phi_e=2.0*PI*rndm_generator.rndm();
	      px += pt*sin(phi_e);
	      py += pt*cos(phi_e);
	      pz=sqrt(e*e-px*px-py*py-EMASS*EMASS);
	      if (dir!=1) 
		{
		  e = -e;
		  pz = -pz;
		}
	      if(compt_results_.store_compt(index_of_process, e,px,py,pz,scal,rndm_generator)) 
		{
// 		 secondaries.new_pair(mesh, cellx, celly,min_z,index_of_process, e,px,py,pz, switches.get_pair_ratio(), switches.get_track_secondaries(), switches.get_store_secondaries(), rndm_generator);
		 secondaries.new_pair(mesh, cellx, celly,min_z,index_of_process, e,px,py,pz, switches.get_pair_ratio(), switches.get_track_pairs(), switches.get_store_pairs(), rndm_generator);
		  if (theta_e>0.1) 
		    {
		      pt=0.0;
		      std::cout << " COMPT::compt_do, warning : theta_e= " << theta_e << std::endl;
		    }
		  if(compton_phot_file_ != NULL)
		    {
		      compton_phot_file_->save_compton_photon(y*epart,vx*epart-px,vy*epart-py);
		    }
#ifdef USE_EDM4HEP
		  if (switches.get_do_edm4hep()) EDM4HEPWriter::edmfile().save_comptPhoton_EDM(px,py,pz);
#endif
		}
	    }
	  if (theta_g>1e-3) 
	    {
	      sum3_ += scal;
	      sume3_ += y*epart*scal;
	    }
	  if (theta_e>1e-3) 
	    {
	      sum4_ += scal;
	      sume4_ += (1.0-y)*epart*scal;
	    }
	}
    }
}

