#include "pairsCPP.h"
#include <cmath>
#include <list>
#include <sstream>
#include <string>

#include "physconst.h"


void PAIR_BEAM::distribute_pairs(float delta_z,unsigned int n)
{
  float startLocal,delta_i;
  unsigned int k;
  startLocal = 0.5*delta_z*n;
  delta_i=1.0/delta_z;
  if ( n > active_pairs_.size())
    {
      std::cout << " PAIR_BEAM::distribute_pairs vector too small " << std::endl;
      exit(0);
    }
  for(k=0;k<n;k++) 
    {
      active_pairs_[k].clear();
    }
  std::list<PAIR_PARTICLE>::iterator itr = reserve_.begin();
  int m;
  while(itr != reserve_.end())
    {
      m=(int)floor((startLocal-itr->z())*delta_i);
      if ((m>=0)&&(m< (int)n))
	{
	  active_pairs_[m].push_back(*itr);
	  reserve_.erase(itr++);
	}
      else
	{
	  itr++;
	}
    }
}

void PAIR_BEAM::move_unactive_pairs(float step)
{
  std::list<PAIR_PARTICLE>::iterator itr = reserve_.begin();
  while(itr != reserve_.end())
    {
      itr->advancePosition(step);
      itr++;
    }
}


void PAIR_BEAM::new_pair(const MESH& mesh, int cellx, int celly,float min_z, int index_of_process, float energy,float px,float py,float pz, float ratio, int tracking, int saving, RNDM& rndm_generator )
{
  // test if particle energy is above the required minumum 
  if (fabs(energy) < pair_parameter_.get_ecut() )   
    {
      return;
    }    

  // reduce the number of stored particles if requested (to speed up tracking) 
  if (rndm_generator.rndm_pairs() > ratio)
    {
      return;
    }

  // store particles for tracking?
  if (!tracking  && !saving  ) return;

  float xx,yy,zz,vxx,vyy,vzz;

  mesh.pair_guess_position_in_cell(cellx, celly, min_z, xx, yy, zz,rndm_generator);

  vxx=px/fabs(energy);
  vyy=py/fabs(energy);
  vzz=pz/fabs(energy);

  count_pairs_++;
  PAIR_PARTICLE pair_temp = PAIR_PARTICLE(count_pairs_, index_of_process,xx,yy,zz,vxx,vyy,vzz, energy);
  // if (store_pairs>1) store particles at production time for comparison
  if (saving > 1)
    {
      pairs0_.push_back(pair_temp);
    }

  if (tracking)
    {
      reserve_.push_back(pair_temp);
    }
}

void PAIR_BEAM::new_pair(const unsigned int evtIndex, const MESH& mesh, int cellx, int celly,float min_z, int index_of_process, float energy,float px,float py,float pz, float ratio, int tracking, int saving, RNDM& rndm_generator, int beamslice1=-1, int beamslice2=-1 )
{

  // test if particle energy is above the required minumum
  if (fabs(energy) < pair_parameter_.get_ecut() )
    {
      return;
    }

  // reduce the number of stored particles if requested (to speed up tracking)
  if (rndm_generator.rndm_pairs() > ratio) 
    {
      return;
    }    
  
  // store particles for tracking? 
  if (!tracking  && !saving  ) return;
  
  float xx,yy,zz,vxx,vyy,vzz;
  
  mesh.pair_guess_position_in_cell(cellx, celly, min_z, xx, yy, zz,rndm_generator);
  

  vxx=px/fabs(energy);
  vyy=py/fabs(energy);
  vzz=pz/fabs(energy);
  
  count_pairs_++;
  PAIR_PARTICLE pair_temp = PAIR_PARTICLE(evtIndex, index_of_process,xx,yy,zz,vxx,vyy,vzz, energy);

  // for retrieval of timing
  pair_temp.set_slices12( beamslice1 , beamslice2 );

  // if (store_pairs>1) store particles at production time for comparison 
  if (saving > 1)
    {
      pairs0_.push_back(pair_temp);
    }
 
  if (tracking) 
    {
      reserve_.push_back(pair_temp);
    }
//std::cout << " manu :  end of in pairsCPP:new_pair with evtIndex " << std::endl;

}

// incoherent pair creation : Breit-Wheeler processes
void PAIR_BEAM::make_pair_bw(const MESH& mesh, int cellx, int celly,float min_z,int index_of_process, float eph1,float q2_1,float eorg1, float eph2,float q2_2,float eorg2, float flum,float beta_x,float beta_y, SWITCHES& switches,RNDM& rndm_generator)
{
  const double emass_2=EMASS*EMASS,one=1.0;
  double beta,phi;
  double ptot, gam2i, c;
  double pz1,pz2,e1,e2,px1,px2,py1,py2;
  int i;
  double sigma;
  double sigmap = 0.0;
  int niter;
  double pt;
  double ecm2;

  ecm2 = eph1 * eph2;
  if (ecm2 < emass_2) return;
  gam2i = emass_2 / ecm2;
  beta = sqrt(one - gam2i);
  // formula Thesis, p. 39
  sigma = gam2i * 1.23088e-29 * flum * ((gam2i * 2. + 2. - gam2i * gam2i) *
	  log((one + beta) / (one - beta)) - (gam2i * 2. * gam2i / 
          (one - beta * beta) + 2.) * beta);

  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i)
    {
      PHYSTOOLS::mkit(gam2i, c, rndm_generator);
      ptot = sqrt(ecm2 - emass_2);
      pt = sqrt(one - c * c) * ptot;
      pz1 = c * ptot;
      pz2=-pz1;
      e1 = sqrt(ecm2);
      e2=e1;
      PHYSTOOLS::lorent_pair(eph1,eph2,e1,pz1);
      switch(switches.get_pair_q2())
	{
	case 0:
	  if ((q2_1>emass_2)||(q2_2>emass_2))
	    {
	      sigmap=0.0;
	    }
	  else
	    {
	      sigmap=sigma;
	    }
	  break;
	case 1:
	  if ((q2_1>(pt*pt+emass_2))||(q2_2>(pt*pt+emass_2)))
	    {
	      sigmap=0.0;
	    }
	  else
	    {
	      sigmap=sigma;
	    }
	  break;
	case 2:
	  if ((q2_1>ecm2)||(q2_2>ecm2))
	    {
	      sigmap=0.0;
	    }
	  else
	    {
	      sigmap=sigma;
	    }
	  break;
	}
      phi=2.0*PI*rndm_generator.rndm_pairs();
      px1=sin(phi)*pt;
      py1=cos(phi)*pt;
      px2=-px1;
      py2=-py1;
   
      PHYSTOOLS::lorent(e1,px1,beta_x);
   
      PHYSTOOLS::lorent(e1,py1,beta_y);
      
      PHYSTOOLS::lorent_pair(eph1,eph2, e2, pz2);
      PHYSTOOLS::lorent(e2,px2,beta_x);
      PHYSTOOLS::lorent(e2,py2,beta_y);
      e2= -e2;
    
      book_keeping(mesh, index_of_process, e1,px1,py1,pz1, e2, px2,py2,pz2, sigmap,cellx, celly,min_z, switches, rndm_generator );
      
      if (switches.get_beam_pair())
	{
	  if (eorg1>0.0)
	    {
	      book_keeping_p(mesh, index_of_process, eorg1-eph1,sigmap,cellx, celly,min_z,switches, rndm_generator );
	    }
	  if (eorg2>0.0)
	    {
	      book_keeping_p(mesh,index_of_process, -(eorg2-eph2),sigmap,cellx, celly,min_z,switches,rndm_generator );
	    }
	}
    }
} /* pair_bw */


void  PAIR_BEAM::make_muon(const MESH& mesh, int cellx, int celly,float min_z, int index_of_process, float eph1,float q2_1,float eph2,float q2_2, float flum,float beta_x,float beta_y, SWITCHES& switches,RNDM& rndm_generator)
{
  const double mumass_2=MUMASS*MUMASS,one=1.0;
  double beta,phi;
  double ptot,gam2i,c,e1,e2;
  int i;
  double sigma;
  double sigmap = 0.0;
  int niter;
  double pt,px1,py1,pz1,px2,py2,pz2;
  double ecm2;
  const double emass2=EMASS*EMASS;
  const double fact=1.23088e-29*EMASS*EMASS/MUMASS/MUMASS;
  ecm2 = eph1 * eph2;
  if (ecm2 < mumass_2) return;
  gam2i = mumass_2 / ecm2;
  beta = sqrt(one - gam2i);
  sigma = gam2i * fact * flum * switches.get_muon_scale() * ((gam2i * 2. + 2. - gam2i * gam2i) *
				 log((one + beta) / (one - beta)) 
				 - (gam2i * 2. * gam2i 
				    / (one - beta * beta) + 2.) * beta);
  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i)
    {
      PHYSTOOLS::mkit(gam2i, c, rndm_generator);
      ptot = sqrt(ecm2 - mumass_2);
      pt = sqrt(one - c * c) * ptot;
      pz1 = c * ptot;
      pz2=-pz1;
      e1 = sqrt(ecm2);
      e2=e1;
      PHYSTOOLS::lorent_pair(eph1,eph2,e1,pz1);
      
      switch(switches.get_pair_q2()){
      case 0:
	if ((q2_1>emass2)||(q2_2>emass2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      case 1:
	if ((q2_1>(pt*pt+mumass_2))||(q2_2>(pt*pt+mumass_2))){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	    }
	break;
      case 2:
	if ((q2_1>ecm2)||(q2_2>ecm2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      }
      phi=2.0*PI*rndm_generator.rndm_pairs();
      px1=sin(phi)*pt;
      py1=cos(phi)*pt;
      px2=-px1;
      py2=-py1;
      PHYSTOOLS::lorent(e1,px1,beta_x);
      PHYSTOOLS::lorent(e1,py1,beta_y);

      PHYSTOOLS::lorent_pair(eph1,eph2, e2, pz2);
      PHYSTOOLS::lorent(e2,px2,beta_x);
      PHYSTOOLS::lorent(e2,py2,beta_y);
      e2= -e2;

      book_keeping_muon(mesh,index_of_process, e1,px1,py1,pz1, e2, px2,py2,pz2, sigmap,cellx, celly,min_z, switches, rndm_generator );
      
    }
} /* make_muon */



void PAIR_BEAM::load_events(int time_counter,float ratio, int tracking, RNDM& rndm_generator)
{
  float e,x,y,z,vx,vy,vz;
  std::string line;
  int  t;
  if (file_of_events_ == NULL) return;
       
  if (!event_to_store_.empty())
    {
      std::istringstream tt(event_to_store_);
      tt >> t;
      if (t<=time_counter)
	{
	  if (tt >> e >> x >> y >> z >> vx >> vy >> vz)
	    {
	      new_event(e, x, y, z, vx, vy, vz,ratio, tracking, rndm_generator);
	    }
	  else 
	    {
	      std::cerr << " PAIR_BEAM::load_events: error reading load event on file " << std::endl;
	      std::cerr <<  " e= " << e << " x = " << x << " y= " << y << " z= " << z << " vx = " << vx << " vy= " << vy << " vz= " << vz << std::endl; 
	    }
	}
      else return;
    }
  while (file_of_events_->read_line(line)) 
    {
      std::istringstream ss(line);
      ss >> t;
      if (t > time_counter)
	{
	  event_to_store_ = line;
	  break;
	}
      else 
	{
	  if (ss >> e >> x >> y >> z >> vx >> vy >> vz)
	      new_event(e, x, y, z, vx, vy, vz, ratio, tracking, rndm_generator);
	  else
	    {
	      std::cerr << " PAIR_BEAM::load_events: error reading load event on file " << std::endl;
	      std::cerr <<  " e= " << e << " x = " << x << " y= " << y << " z= " << z << " vx = " << vx << " vy= " << vy << " vz= " << vz << std::endl; 
	    }
	}
    }
}

std::string PAIR_BEAM::output_flow() const 
{
  // unused method : concerns only particular pairs (see 
  // PAIR_BEAM::compute_pairs_calls)
  double e1,e2;
  long int n1,n2;
  std::ostringstream out;
  out << title(std::string("pair beams"));
  compute_pairs_calls(n1, e1, n2, e2); 
  out << "pairs_ncal.1 = " << n1 << " energy (GeV) : pairs_cal.1 = " << e1 << std::endl;
  out << "pairs_ncal.2 = " << n2 << " energy (GeV) : pairs_cal.2 = " << e2 << std::endl;
  out << "mstp: " << count_pairs_ << " " << pair_track_.call() << " " <<  (float)pair_track_.step()/ std::max( (float)pair_track_.call(),float(1.0)) << std::endl;
  return out.str();
}

void PAIR_BEAM::compute_pairs_calls(long int& n1, double& e1, long int& n2, double& e2) const
{
  std::list<PAIR_PARTICLE>::const_iterator point;
  float vx, vy,vz, energy, abs_energy;
  float z,r,b;
  float pt,pz,r0,phi;
 
  z=2300.0;
  r=18.0;
  b=3.0;
  e1=0.0;
  e2=0.0;
  n1=0;
  n2=0;
  point= reserve_.begin();
  while(point != reserve_.end())
    {
      //      point->velocities(vx,vy,vz);
      point->velocities(vx,vy);
      vz = point->Zvelocity();
      abs_energy = point->unsigned_energy();
      energy = point->signed_energy();
	pt=sqrt(vx*vx+vy*vy)*abs_energy;
	pz=vz*energy;
	r0=3.3333e3*pt/b;
	phi=z/r0*pt/std::abs(pz);
	if (sqrt(2.0*(1.0-cos(phi)))*r0>r){
	  if (vz>0.0) {
	    e1+=abs_energy;
	    n1++;
	  }
	  else{
	    e2+=abs_energy;
	    n2++;
	  }
	}
	point++;
    }
}

void  PAIR_PARAMETER::init(BEAM& beam1, BEAM& beam2, int massflag, float pair_ecut, float pair_step, float step, int timestep)
{

  float sigma_x1, sigma_y1;
  float sigma_x2, sigma_y2;


  beam1.transverse_sigmas(sigma_x1, sigma_y1);
  beam2.transverse_sigmas(sigma_x2, sigma_y2);

  ecut=pair_ecut;

  s4 =  beam1.get_ebeam()  *  beam2.get_ebeam();
  lns4=log(s4/(EMASS*EMASS));
  d_eps_1_ = 2.0*RE*1e9* step *EMASS* pair_step / (( sigma_x1+ sigma_y1 )* sigma_y1 * timestep);
  d_eps_2_ = 2.0*RE*1e9*step*EMASS*pair_step / ((sigma_x2+sigma_y2)*sigma_y2*timestep);
  if(massflag==0){
    mass_=EMASS;
  }
  else if(massflag==1){
    mass_=MUMASS;
  }
}

void PAIR_PARAMETER::jet_equiv (float xmin,float e,int iflag,float& eph,float& q2,float& wgt, RNDM& rndm_generator) const
{
  const float emass2=EMASS*EMASS,eps=1e-30;
  float help,q2max,q2min,lnx,x,qxmin;
  wgt=1.0;
  switch (iflag)
    {
    case 1:
      q2max=1.0;
      eph=(1.0-xmin)*rndm_generator.rndm_equiv()+xmin;
      help=1.0-eph;
      if (help<=eps){
	  eph=0.0;
	  wgt=0.0;
	  q2=0.0;
	  return;
      }
      qxmin=eph*eph/help;
      q2min=emass2*qxmin;
      wgt=-0.5*(1.0+(help*help)) / (eph*log(xmin))
	   *log(q2max/q2min)
	   /lns4*(1.0-xmin);
      wgt=std::max(float(0.0),wgt);
      q2= q2min*pow(q2max/q2min,rndm_generator.rndm_equiv());
      eph *= e;
      return;
    case 2:
      eph=(1.0-xmin)*rndm_generator.rndm_equiv()+xmin;
      help=1.0 - eph;
      q2min=emass2;
      wgt=-0.5*(1.0+(help*help)) / (eph*log(xmin))*(1.0-xmin);
      wgt=std::max(float(0.0),wgt);
      q2= q2min*pow(e*e/q2min,rndm_generator.rndm_equiv());
      eph *= e;
      return;
    case 3:
      q2max=1.0;
      eph = (1.0-xmin)*rndm_generator.rndm_equiv()+xmin;
      help=1.0 - eph;
      if (help<=eps){
	  eph=0.0;
	  wgt=0.0;
	  q2=0.0;
	  return;
      }
      qxmin=eph*eph/help;
      q2min=emass2*qxmin;
      wgt=-0.5*(1.0+(help*help)) / (eph*log(xmin))
	   *log(q2max/q2min)
	   /lns4*(1.0-xmin);
      wgt=std::max(float(0.0),wgt);
      q2= q2min*pow(q2max/q2min,rndm_generator.rndm_equiv());
      eph *= e;
      return;
    case 4:
      eph=(1.0-xmin)*rndm_generator.rndm_equiv()+xmin;
      help=1.0 - eph;
      if (help<=0.0){
	  eph=0.0;
	  wgt=0.0;
	  q2=0.0;
	  return;
      }
      qxmin=eph*eph/help;
      q2min=emass2*qxmin;
      wgt=-0.5*(1.0+(help*help)) / (eph*log(xmin))
	   *(lns4-log(qxmin))
	   /lns4*(1.0-xmin);
      wgt=std::max(float(0.0),wgt);
      q2= q2min*pow(e*e/q2min,rndm_generator.rndm_equiv());
      eph*=e;
      return;
    case 5:
	if(rndm_generator.rndm_equiv()<0.5){
	    lnx=-sqrt(rndm_generator.rndm_equiv())
		*lns4;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_generator.rndm_equiv()*lns4;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=s4;
	}
	if((1.0+(1.0-x)*(1.0-x))*0.5<rndm_generator.rndm_equiv()){
	    eph=0.0;
	    q2=0.0;
	}
	else{
	    eph=e*x;
	    q2=q2min*pow(q2max/q2min,rndm_generator.rndm_equiv());
	}
	if (q2*(1.0-x)<x*x*emass2) eph=0.0;
	return;
    }
} 

