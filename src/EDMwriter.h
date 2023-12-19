#ifndef EDMWRITER_SEEN 
#define EDMWRITER_SEEN

//to avoid unused parameter warnings
#define DEF_DO_EDM4HEP

#ifdef USE_EDM4HEP

//to avoid unused parameter warnings
#undef DEF_DO_EDM4HEP
#define DEF_DO_EDM4HEP do_edm4hep

#include "edm4hep/MCParticleCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"
#include "podio/UserDataCollection.h"
#include "switchesCPP.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <iomanip>



class PAIR_PARTICLE;
class PARTICLE;
class PHOTON;
class ABSTRACT_BHABHASAMPLES;
class ABSTRACT_BHABHA_PHOTON_SAMPLES;
class BEAM_PARAMETERS;
class GRID;
class RESULTS;

enum LUMI_TYPE {
  lumi_eg,
  lumi_ge,
  lumi_gg
};

enum BEAM_TYPE {
  beam1,
  beam2,
  coh1,
  coh2,
  tri1,
  tri2
};

enum PAIR_PARTICLE_TYPE {
  electronPair = 0,
  electronPair0,
  muonPair,
  muonPair0,
  bhabhaPair,
  bhabhaPair0
};

class FILE_IN_OUT_EDM4HEP
{
  podio::EventStore m_store;

  podio::ROOTWriter* m_writer;

  edm4hep::MCParticleCollection* m_pairsCollection;
  edm4hep::MCParticleCollection* m_pairs0Collection;

  edm4hep::MCParticleCollection* m_bhabhasCollection;
  edm4hep::MCParticleCollection* m_bhabhas0Collection;

  edm4hep::MCParticleCollection* m_muonsCollection;
  edm4hep::MCParticleCollection* m_muons0Collection;

  podio::UserDataCollection<float>* m_hadron_energy1Collection;
  podio::UserDataCollection<float>* m_hadron_energy2Collection;
  podio::UserDataCollection<float>* m_hadron_zCollection;

  podio::UserDataCollection<float>* m_jets_eph1Collection;
  podio::UserDataCollection<float>* m_jets_eph2Collection;
  podio::UserDataCollection<float>* m_jets_pz1Collection;
  podio::UserDataCollection<float>* m_jets_pz2Collection;
  podio::UserDataCollection<float>* m_jets_ptCollection;
  podio::UserDataCollection<int>* m_jets_eventCollection;

  edm4hep::MCParticleCollection* m_bhPhotonProdCollection;
  edm4hep::MCParticleCollection* m_bhPhotonsCollection;

  podio::UserDataCollection<float>* m_lumiee_energy1Collection;
  podio::UserDataCollection<float>* m_lumiee_energy2Collection;
  podio::UserDataCollection<float>* m_lumiee_xposCollection;
  podio::UserDataCollection<float>* m_lumiee_yposCollection;
  podio::UserDataCollection<float>* m_lumiee_zposCollection;
  podio::UserDataCollection<float>* m_lumiee_vx1Collection;
  podio::UserDataCollection<float>* m_lumiee_vy1Collection;
  podio::UserDataCollection<float>* m_lumiee_vx2Collection;
  podio::UserDataCollection<float>* m_lumiee_vy2Collection;
  podio::UserDataCollection<float>* m_lumiee_sx1Collection;
  podio::UserDataCollection<float>* m_lumiee_sy1Collection;
  podio::UserDataCollection<float>* m_lumiee_sz1Collection;
  podio::UserDataCollection<float>* m_lumiee_sx2Collection;
  podio::UserDataCollection<float>* m_lumiee_sy2Collection;
  podio::UserDataCollection<float>* m_lumiee_sz2Collection;
  podio::UserDataCollection<int>* m_lumiee_tCollection;
 
  podio::UserDataCollection<float>* m_lumieg_energy1Collection;
  podio::UserDataCollection<float>* m_lumieg_energy2Collection;
  podio::UserDataCollection<float>* m_lumieg_xposCollection;
  podio::UserDataCollection<float>* m_lumieg_yposCollection;
  podio::UserDataCollection<float>* m_lumieg_zposCollection;

  podio::UserDataCollection<float>* m_lumige_energy1Collection;
  podio::UserDataCollection<float>* m_lumige_energy2Collection;
  podio::UserDataCollection<float>* m_lumige_xposCollection;
  podio::UserDataCollection<float>* m_lumige_yposCollection;
  podio::UserDataCollection<float>* m_lumige_zposCollection;

  podio::UserDataCollection<float>* m_lumigg_energy1Collection;
  podio::UserDataCollection<float>* m_lumigg_energy2Collection;
  podio::UserDataCollection<float>* m_lumigg_xposCollection;
  podio::UserDataCollection<float>* m_lumigg_yposCollection;
  podio::UserDataCollection<float>* m_lumigg_zposCollection;

  edm4hep::MCParticleCollection* m_beam1Collection;
  edm4hep::MCParticleCollection* m_beam2Collection;
  edm4hep::MCParticleCollection* m_coh1Collection;
  edm4hep::MCParticleCollection* m_coh2Collection;
  edm4hep::MCParticleCollection* m_tri1Collection;
  edm4hep::MCParticleCollection* m_tri2Collection;
  
  edm4hep::MCParticleCollection* m_photon1Collection;
  edm4hep::MCParticleCollection* m_photon2Collection;

  edm4hep::MCParticleCollection* m_comptPhotonCollection;

  podio::UserDataCollection<int>* m_bhabhas_prod_eventIndexCollection;
  podio::UserDataCollection<float>* m_bhabhas_prod_eCMCollection;
  podio::UserDataCollection<float>* m_bhabhas_prod_Energy1Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_Energy2Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_MotherEnergy1Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_MotherEnergy2Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vx1Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vy1Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vz1Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vx2Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vy2Collection;
  podio::UserDataCollection<float>* m_bhabhas_prod_vz2Collection;
  podio::UserDataCollection<int>* m_bhabhas_prod_nbphotCollection;

  bool m_beams_are_same_charge;

 public :

 FILE_IN_OUT_EDM4HEP(std::string const& name, SWITCHES& switches): m_store() ,

    m_writer(new podio::ROOTWriter(name, &m_store)),

    m_pairsCollection(&m_store.create<edm4hep::MCParticleCollection>("Pairs")),
    m_pairs0Collection(&m_store.create<edm4hep::MCParticleCollection>("Pairs0")),

    m_bhabhasCollection(&m_store.create<edm4hep::MCParticleCollection>("Bhabhas")),
    m_bhabhas0Collection(&m_store.create<edm4hep::MCParticleCollection>("Bhabhas0")),

    m_muonsCollection(&m_store.create<edm4hep::MCParticleCollection>("Muons")),
    m_muons0Collection(&m_store.create<edm4hep::MCParticleCollection>("Muons0")),

    m_hadron_energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Hadrons_energy1")),
    m_hadron_energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Hadrons_energy2")),
    m_hadron_zCollection(&m_store.create<podio::UserDataCollection<float>>("Hadrons_z")),

    m_jets_eph1Collection(&m_store.create<podio::UserDataCollection<float>>("Jets_phEnergy1")),
    m_jets_eph2Collection(&m_store.create<podio::UserDataCollection<float>>("Jets_phEnergy2")),
    m_jets_pz1Collection(&m_store.create<podio::UserDataCollection<float>>("Jets_pz1")),
    m_jets_pz2Collection(&m_store.create<podio::UserDataCollection<float>>("Jets_pz2")),
    m_jets_ptCollection(&m_store.create<podio::UserDataCollection<float>>("Jets_pt")),
    m_jets_eventCollection(&m_store.create<podio::UserDataCollection<int>>("Jets_event")),

    m_bhPhotonProdCollection(&m_store.create<edm4hep::MCParticleCollection>("Bhabha_photons")),
    m_bhPhotonsCollection(&m_store.create<edm4hep::MCParticleCollection>("Bhabha_photons0")),

    m_lumiee_energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_energy1")),
    m_lumiee_energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_energy2")),
    m_lumiee_xposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_xpos")),
    m_lumiee_yposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_ypos")),
    m_lumiee_zposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_zpos")),
    m_lumiee_vx1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_vx1")),
    m_lumiee_vy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_vy1")),
    m_lumiee_vx2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_vx2")),
    m_lumiee_vy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_vy2")),
    m_lumiee_sx1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sx1")),
    m_lumiee_sy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sy1")),
    m_lumiee_sz1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sz1")),
    m_lumiee_sx2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sx2")),
    m_lumiee_sy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sy2")),
    m_lumiee_sz2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ee_sz2")),
    m_lumiee_tCollection(&m_store.create<podio::UserDataCollection<int>>("Lumi_ee_t")),

    m_lumieg_energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_eg_energy1")),
    m_lumieg_energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_eg_energy2")),
    m_lumieg_xposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_eg_xpos")),
    m_lumieg_yposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_eg_ypos")),
    m_lumieg_zposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_eg_zpos")),

    m_lumige_energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ge_energy1")),
    m_lumige_energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ge_energy2")),
    m_lumige_xposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ge_xpos")),
    m_lumige_yposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ge_ypos")),
    m_lumige_zposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_ge_zpos")),

    m_lumigg_energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_gg_energy1")),
    m_lumigg_energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Lumi_gg_energy2")),
    m_lumigg_xposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_gg_xpos")),
    m_lumigg_yposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_gg_ypos")),
    m_lumigg_zposCollection(&m_store.create<podio::UserDataCollection<float>>("Lumi_gg_zpos")),
    
    m_beam1Collection(&m_store.create<edm4hep::MCParticleCollection>("Beam1")),
    m_beam2Collection(&m_store.create<edm4hep::MCParticleCollection>("Beam2")),

    m_coh1Collection(&m_store.create<edm4hep::MCParticleCollection>("CoherentBeam1")),
    m_coh2Collection(&m_store.create<edm4hep::MCParticleCollection>("CoherentBeam2")),

    m_tri1Collection(&m_store.create<edm4hep::MCParticleCollection>("TridentBeam1")),
    m_tri2Collection(&m_store.create<edm4hep::MCParticleCollection>("TridentBeam2")),

    m_photon1Collection(&m_store.create<edm4hep::MCParticleCollection>("Photons1")),
    m_photon2Collection(&m_store.create<edm4hep::MCParticleCollection>("Photons2")),

    m_comptPhotonCollection(&m_store.create<edm4hep::MCParticleCollection>("ComptonPhotons")),

    m_bhabhas_prod_eventIndexCollection(&m_store.create<podio::UserDataCollection<int>>("Bhabhas_prod_EvtIndex")),
    m_bhabhas_prod_eCMCollection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_eCM")),
    m_bhabhas_prod_Energy1Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_Energy1")),
    m_bhabhas_prod_Energy2Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_Energy2")),
    m_bhabhas_prod_MotherEnergy1Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_MotherEnergy1")),
    m_bhabhas_prod_MotherEnergy2Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_MotherEnergy2")),
    m_bhabhas_prod_vx1Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vx1")),
    m_bhabhas_prod_vy1Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vy1")),
    m_bhabhas_prod_vz1Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vz1")),
    m_bhabhas_prod_vx2Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vx2")),
    m_bhabhas_prod_vy2Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vy2")),
    m_bhabhas_prod_vz2Collection(&m_store.create<podio::UserDataCollection<float>>("Bhabhas_prod_vz2")),
    m_bhabhas_prod_nbphotCollection(&m_store.create<podio::UserDataCollection<int>>("Bhabhas_prod_nbphot")),
    m_beams_are_same_charge(switches.get_charge_sign() == 1)
      {
	if (switches.get_do_bhabhas())
	  {
	    m_writer->registerForWrite("Bhabhas_prod_EvtIndex");
	    m_writer->registerForWrite("Bhabhas_prod_eCM");
	    m_writer->registerForWrite("Bhabhas_prod_Energy1");
	    m_writer->registerForWrite("Bhabhas_prod_Energy2");
	    m_writer->registerForWrite("Bhabhas_prod_MotherEnergy1");
	    m_writer->registerForWrite("Bhabhas_prod_MotherEnergy2");
	    m_writer->registerForWrite("Bhabhas_prod_vx1");
	    m_writer->registerForWrite("Bhabhas_prod_vy1");
	    m_writer->registerForWrite("Bhabhas_prod_vz1");
	    m_writer->registerForWrite("Bhabhas_prod_vx2");
	    m_writer->registerForWrite("Bhabhas_prod_vy2");
	    m_writer->registerForWrite("Bhabhas_prod_vz2");
	    m_writer->registerForWrite("Bhabhas_prod_nbphot");

	    m_writer->registerForWrite("Bhabha_photons");
	    m_writer->registerForWrite("Bhabha_photons0");

	    if (switches.get_track_pairs() || switches.get_store_pairs()) m_writer->registerForWrite("Bhabhas");
	    if (switches.get_store_pairs() > 1)	m_writer->registerForWrite("Bhabhas0");
	  }
	else
	  {
	    if (switches.get_track_pairs() || switches.get_store_pairs()) m_writer->registerForWrite("Pairs");
	    if (switches.get_store_pairs() > 1) m_writer->registerForWrite("Pairs0");
	  }

	if (switches.get_track_muons() || switches.get_store_muons()) m_writer->registerForWrite("Muons");
	if (switches.get_store_muons() > 1) m_writer->registerForWrite("Muons0");

     
	if (switches.get_do_hadrons())
	  {
	    m_writer->registerForWrite("Hadrons_energy1");
	    m_writer->registerForWrite("Hadrons_energy2");
	    m_writer->registerForWrite("Hadrons_z");
	  }

	if (switches.get_store_beam()) 
	  {
	    m_writer->registerForWrite("Beam1");
	    m_writer->registerForWrite("Beam2");
	    if ( switches.get_do_coherent() )
	      {
		m_writer->registerForWrite("CoherentBeam1");
		m_writer->registerForWrite("CoherentBeam2");
	      }
	    if ( switches.get_do_trident() )
	      {
		m_writer->registerForWrite("TridentBeam1");
		m_writer->registerForWrite("TridentBeam2");
	      }
	    if ( switches.get_write_photons() )
	      {
		m_writer->registerForWrite("Photons1");
		m_writer->registerForWrite("Photons2");
	      }
	  }

	m_writer->registerForWrite("ComptonPhotons");

	if (switches.get_do_jets())
	  {
	    m_writer->registerForWrite("Jets_phEnergy1");
	    m_writer->registerForWrite("Jets_phEnergy2");
	    m_writer->registerForWrite("Jets_pz1");
	    m_writer->registerForWrite("Jets_pz2");
	    m_writer->registerForWrite("Jets_pt");
	    m_writer->registerForWrite("Jets_event");
	  }

	if(switches.get_do_lumi())
	  {
	    if (switches.get_do_lumi()&1) 
	      {
		m_writer->registerForWrite("Lumi_ee_energy1");
		m_writer->registerForWrite("Lumi_ee_energy2");
		m_writer->registerForWrite("Lumi_ee_xpos");
		m_writer->registerForWrite("Lumi_ee_ypos");
		m_writer->registerForWrite("Lumi_ee_zpos");
		m_writer->registerForWrite("Lumi_ee_vx1");
		m_writer->registerForWrite("Lumi_ee_vy1");
		m_writer->registerForWrite("Lumi_ee_vx2");
		m_writer->registerForWrite("Lumi_ee_vy2");
		m_writer->registerForWrite("Lumi_ee_sx1");
		m_writer->registerForWrite("Lumi_ee_sy1");
		m_writer->registerForWrite("Lumi_ee_sz1");
		m_writer->registerForWrite("Lumi_ee_sx2");
		m_writer->registerForWrite("Lumi_ee_sy2");
		m_writer->registerForWrite("Lumi_ee_sz2");
		m_writer->registerForWrite("Lumi_ee_t");
	      }
	    if (switches.get_do_lumi()&2) 
	      {	
		m_writer->registerForWrite("Lumi_eg_energy1");
		m_writer->registerForWrite("Lumi_eg_energy2");
		m_writer->registerForWrite("Lumi_eg_xpos");
		m_writer->registerForWrite("Lumi_eg_ypos");
		m_writer->registerForWrite("Lumi_eg_zpos");
	     
		m_writer->registerForWrite("Lumi_ge_energy1");
		m_writer->registerForWrite("Lumi_ge_energy2");
		m_writer->registerForWrite("Lumi_ge_xpos");
		m_writer->registerForWrite("Lumi_ge_ypos");
		m_writer->registerForWrite("Lumi_ge_zpos");
	      }
	    if (switches.get_do_lumi()&4) 	
	      {
		m_writer->registerForWrite("Lumi_gg_energy1");
		m_writer->registerForWrite("Lumi_gg_energy2");
		m_writer->registerForWrite("Lumi_gg_xpos");
		m_writer->registerForWrite("Lumi_gg_ypos");
		m_writer->registerForWrite("Lumi_gg_zpos");
	      }
	 
	  }

      }
  virtual ~FILE_IN_OUT_EDM4HEP() {;}

  virtual void save_pair_EDM(const PAIR_PARTICLE& obj, PAIR_PARTICLE_TYPE type);

  virtual void save_hadron_EDM(float& energy1, float& energy2, float& z);
  virtual void save_jet_EDM(float eph1, float eph2, float pz1, float pz2, float pt, int event);
  virtual void save_lumi_ee_EDM(float& energy1,float& energy2, float& xpos, float& ypos, float& zpos, float& vx1, float& vy1, float& vx2, float& vy2, float& sx1, float& sy1, float& sz1, float& sx2, float& sy2, float& sz2, int& t);
  virtual void save_lumi_EDM(float& energy1, float& energy2, float& xpos, float& ypos, float& zpos, LUMI_TYPE type);

  virtual void save_beam_particle_EDM(const PARTICLE& obj, BEAM_TYPE type);

  virtual void save_photons1_EDM(const PHOTON& obj);
  virtual void save_photons2_EDM(const PHOTON& obj);
  virtual void save_bhabhasamples_EDM(const ABSTRACT_BHABHASAMPLES* const bhabhas);
  virtual void save_bhabhaPhotonSamples_EDM(const ABSTRACT_BHABHA_PHOTON_SAMPLES*  bhabhaPhot);
  virtual void save_bhabhaPhotons_EDM(const ABSTRACT_BHABHA_PHOTON_SAMPLES*  bhabhaPhot);
  virtual void save_comptPhoton_EDM(float px, float py, float pz);


  template<typename T>
  void write_Metadata(std::string const& name, T value) {
    auto& evtMD = m_store.getEventMetaData();
    evtMD.setValue(name, (T)(value));
  }

 /* virtual void write_grid_MD(const PAIRS_RESULTS obj) */
 /* { */
 /*  auto& evtMD1 =m_store.getEventMetaData() ; */
 /*   evtMD1.setValue( "n_pairs" , obj.number()) ; */
 /* } */


 virtual void close()
 {
   m_writer->writeEvent();
   m_writer->finish();
 }



};

template<>
inline void FILE_IN_OUT_EDM4HEP::write_Metadata<double>(std::string const& name, double value) {
  auto& evtMD = m_store.getEventMetaData();
  evtMD.setValue(name, static_cast<float>(value));
}

template<>
inline void FILE_IN_OUT_EDM4HEP::write_Metadata<long unsigned int>(std::string const& name, long unsigned int value) {
  auto& evtMD = m_store.getEventMetaData();
  evtMD.setValue(name, static_cast<int>(value));
}

template<>
inline void FILE_IN_OUT_EDM4HEP::write_Metadata<long int>(std::string const& name, long int value) {
  auto& evtMD = m_store.getEventMetaData();
  evtMD.setValue(name, static_cast<int>(value));
}

class EDM4HEPWriter 
{
 private:
  static FILE_IN_OUT_EDM4HEP* edmfile_;

 public:

  static void initialize(std::string const& name, SWITCHES& switches);
  static FILE_IN_OUT_EDM4HEP& edmfile();
  
 };


#endif

#endif
