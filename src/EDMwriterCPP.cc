#ifdef USE_EDM4HEP

#include "particlesCPP.h"
#include "switchesCPP.h"
#include "beamParametersCPP.h"
#include "gridCPP.h"

#include <iostream>
#include <string>
#include <fstream> 
#include <iomanip>

#include "EDMwriter.h"

FILE_IN_OUT_EDM4HEP* EDM4HEPWriter::edmfile_ = nullptr;

FILE_IN_OUT_EDM4HEP& EDM4HEPWriter::edmfile() 
{
  if(not edmfile_) {
    throw std::runtime_error("EDM4HEP File was not initialised!");
  }
  return *edmfile_;
}

void EDM4HEPWriter::initialize(std::string const& name, SWITCHES& switches)
{
  edmfile_ = new FILE_IN_OUT_EDM4HEP(name,switches);
}


//create pair collection members
  void FILE_IN_OUT_EDM4HEP::save_pair_EDM(const PAIR_PARTICLE& obj, PAIR_PARTICLE_TYPE type) {
  float  vx, vy, xpos, ypos;
  // int process, slice1, slice2;

  obj.XYposition(xpos, ypos);
  obj.velocities(vx, vy);
  
  auto p =
    (type == electronPair ) ? m_pairsCollection->create():
    (type == electronPair0 ) ? m_pairs0Collection->create():
    (type == muonPair ) ? m_muonsCollection->create():
    (type == muonPair0) ? m_muons0Collection->create():
    (type == bhabhaPair ) ? m_bhabhasCollection->create():
    (type == bhabhaPair0) ? m_bhabhas0Collection->create():
    throw std::logic_error("Bad Particle Type used");

  if ( type == muonPair || type == muonPair0 ) {
    (obj.signed_energy() < 0) ? p.setPDG(-13): p.setPDG(13);
    p.setMass(0.1056583755);
  } else {
    (obj.signed_energy() < 0) ? p.setPDG(-11) : p.setPDG(11);
    p.setMass(0.00051099895000);
  };
  p.setVertex({xpos , ypos , obj.z()});
  p.setMomentum({vx*obj.energy(), vy*obj.energy(), obj.Zvelocity()*obj.energy()});

}

void FILE_IN_OUT_EDM4HEP::save_hadron_EDM(float& energy1, float& energy2, float& z)
{
  m_hadron_energy1Collection->push_back(energy1);
  m_hadron_energy2Collection->push_back(energy2);
  m_hadron_zCollection->push_back(z);
}

void FILE_IN_OUT_EDM4HEP::save_jet_EDM(float eph1, float eph2, float pz1, float pz2, float pt, int event)
{
  m_jets_eph1Collection->push_back(eph1);
  m_jets_eph2Collection->push_back(eph2);
  m_jets_pz1Collection->push_back(pz1);
  m_jets_pz2Collection->push_back(pz2);
  m_jets_ptCollection->push_back(pt);
  m_jets_eventCollection->push_back(event);
}

void FILE_IN_OUT_EDM4HEP::save_lumi_ee_EDM(float& energy1,float& energy2, float& xpos, float& ypos, float& zpos, float& vx1, float& vy1, float& vx2, float& vy2, float& sx1, float& sy1, float& sz1, float& sx2, float& sy2, float& sz2, int& t)
{
  m_lumiee_energy1Collection->push_back(energy1);
  m_lumiee_energy2Collection->push_back(energy2);
  m_lumiee_xposCollection->push_back(xpos);
  m_lumiee_yposCollection->push_back(ypos);
  m_lumiee_zposCollection->push_back(zpos);
  m_lumiee_vx1Collection->push_back(vx1);
  m_lumiee_vy1Collection->push_back(vy1);
  m_lumiee_vx2Collection->push_back(vx2);
  m_lumiee_vy2Collection->push_back(vy2);
  m_lumiee_sx1Collection->push_back(sx1);
  m_lumiee_sy1Collection->push_back(sy1);
  m_lumiee_sz1Collection->push_back(sz1);
  m_lumiee_sx2Collection->push_back(sx2);
  m_lumiee_sy2Collection->push_back(sy2);
  m_lumiee_sz2Collection->push_back(sz2);
  m_lumiee_tCollection->push_back(t);
}

void FILE_IN_OUT_EDM4HEP::save_lumi_EDM(float& energy1, float& energy2, float& xpos, float& ypos, float& zpos, LUMI_TYPE type)
{

  if(type == lumi_eg) {
    m_lumieg_energy1Collection->push_back(energy1);
    m_lumieg_energy2Collection->push_back(energy2);
    m_lumieg_xposCollection->push_back(xpos);
    m_lumieg_yposCollection->push_back(ypos);
    m_lumieg_zposCollection->push_back(zpos);
  } else if ( type == lumi_ge) {
    m_lumige_energy1Collection->push_back(energy1);
    m_lumige_energy2Collection->push_back(energy2);
    m_lumige_xposCollection->push_back(xpos);
    m_lumige_yposCollection->push_back(ypos);
    m_lumige_zposCollection->push_back(zpos);
  }else if(type == lumi_gg) {
    m_lumigg_energy1Collection->push_back(energy1);
    m_lumigg_energy2Collection->push_back(energy2);
    m_lumigg_xposCollection->push_back(xpos);
    m_lumigg_yposCollection->push_back(ypos);
    m_lumigg_zposCollection->push_back(zpos);
  }
}

void FILE_IN_OUT_EDM4HEP::save_beam_particle_EDM(const PARTICLE& obj, BEAM_TYPE type) {

  float energy = abs(obj.getEnergy());
  float  vx, vy, xpos, ypos;
  obj.XYposition(xpos, ypos);
  obj.velocities(vx, vy);
  float vz = sqrt(1-vx*vx-vy*vy);

  auto b=
    (type == beam1) ?  m_beam1Collection->create():
    (type == beam2) ?  m_beam2Collection->create():
    (type == coh1) ?  m_coh1Collection->create():
    (type == coh2) ?  m_coh2Collection->create():
    (type == tri1) ?  m_tri1Collection->create():
    (type == tri2) ?  m_tri2Collection->create():
    throw std::logic_error("Bad BEAM TYPE used");

  // if beams are the same charge, we assume they are both electrons, no-one would collide positron vs positron
  bool isElectron = (type == beam1 or type == beam2) and m_beams_are_same_charge;

  if(obj.getEnergy() > 0 or isElectron) {
     b.setPDG(11);//particles are electrons
  } else {
    b.setPDG(-11);//particles are positrons
  }

  b.setVertex({ (xpos*1e-3) , (ypos*1e-3) , (obj.z()*1e-3) });
  b.setMomentum({static_cast<float>(vx*1e6*energy), static_cast<float>(vy*1e6*energy), static_cast<float>(vz*1e6*energy)});
  b.setMass(0.00051099895000);
}

void FILE_IN_OUT_EDM4HEP::save_photons1_EDM(const PHOTON& obj)
{
  float energy = obj.energy();
  float  vx, vy;
  obj.velocities(vx, vy);
  float vz=sqrt(1-vx*vx-vy*vy);

  auto ph = m_photon1Collection->create();
  ph.setPDG(22);
  ph.setMomentum( { (vx*energy) , (vy*energy) , (vz*energy)}  ) ;
}

void FILE_IN_OUT_EDM4HEP::save_photons2_EDM(const PHOTON& obj)
{
  float energy = obj.energy();
  float  vx, vy;
  obj.velocities(vx, vy);
  float vz=sqrt(1-vx*vx-vy*vy);

  auto ph = m_photon2Collection->create();
  ph.setPDG(22);
  ph.setMomentum( { (vx*energy) , (vy*energy) , (vz*energy)}  ) ;

}

void FILE_IN_OUT_EDM4HEP::save_bhabhasamples_EDM(const ABSTRACT_BHABHASAMPLES* const bhabhas)
{
  unsigned int k, evtIndex;
  float eCM, mother1_en, e1, vx1, vy1, vz1, mother2_en, e2, vx2, vy2, vz2;
  int nbphot;
  for (k=0; k < bhabhas->nb_samples(); k++)
    {
      bhabhas->get_parameters_for_output(k, evtIndex, eCM, mother1_en, e1, vx1, vy1, vz1, mother2_en, e2, vx2, vy2, vz2, nbphot);

      m_bhabhas_prod_eventIndexCollection->push_back(evtIndex);
      m_bhabhas_prod_eCMCollection->push_back(eCM);
      m_bhabhas_prod_Energy1Collection->push_back(e1);
      m_bhabhas_prod_Energy2Collection->push_back(e2);
      m_bhabhas_prod_MotherEnergy1Collection->push_back(mother1_en);
      m_bhabhas_prod_MotherEnergy2Collection->push_back(mother2_en);
      m_bhabhas_prod_vx1Collection->push_back(vx1);
      m_bhabhas_prod_vy1Collection->push_back(vy1);
      m_bhabhas_prod_vz1Collection->push_back(vz1);
      m_bhabhas_prod_vx2Collection->push_back(vx2);
      m_bhabhas_prod_vy2Collection->push_back(vy2);
      m_bhabhas_prod_vz2Collection->push_back(vz2);
      m_bhabhas_prod_nbphotCollection->push_back(nbphot);
    }
}


void FILE_IN_OUT_EDM4HEP::save_bhabhaPhotonSamples_EDM(const ABSTRACT_BHABHA_PHOTON_SAMPLES*  bhabhaPhot)
{
  unsigned int k;
  float en, vx, vy, vz;
  int evtIndex;

  for (k=0; k < bhabhaPhot->nb_samples(); k++)
    {
      bhabhaPhot->get_parameters_for_output(k,evtIndex, en, vx, vy, vz);

      auto bhph =  m_bhPhotonProdCollection->create();
      bhph.setPDG(22);
      bhph.setMomentum({(vx*en),(vy*en),(vz*en)});
    }
}

void FILE_IN_OUT_EDM4HEP::save_bhabhaPhotons_EDM(const ABSTRACT_BHABHA_PHOTON_SAMPLES*  bhabhaPhot)
{
  unsigned int k;
  float en, vx, vy, vz;
  int evtIndex;

  for (k=0; k < bhabhaPhot->nb_samples(); k++)
    {
      bhabhaPhot->get_parameters_for_output(k,evtIndex, en, vx, vy, vz);

      auto bhph =  m_bhPhotonsCollection->create();
      bhph.setPDG(22);
      bhph.setMomentum({(vx*en),(vy*en),(vz*en)});
    }
}

void FILE_IN_OUT_EDM4HEP::save_comptPhoton_EDM(float px, float py, float pz)
{
  auto cph = m_comptPhotonCollection->create();
  cph.setPDG(22);
  cph.setMomentum({(px),(py),(pz)});
}

#endif
