#include <iostream>
#include <cstdlib>
#include <sstream>
#include "switchesCPP.h"
#include "physconst.h"
#include "EDMwriter.h"

SWITCHES::SWITCHES() 
{
  electron_distribution_rho=2;
  //    electron_distribution_scatter=1;
  electron_ratio=1.0;
  do_lumi=0;
  num_lumi=100000; num_lumi_eg=100000; num_lumi_gg=100000;
  do_cross=0;
  lumi_p=1e-29; lumi_p_eg=1e-29; lumi_p_gg=1e-29;
  do_photons_[0]=0;
  do_photons_[1]=0;
  write_photons=0; //=store_photons
  photon_distribution=1;
  photon_ratio=1.0;
  do_hadrons=0;
  store_hadrons=0;
  hadron_ratio=1e5;
  do_jets=0;
  jet_store=0; //=store_jets
  jet_pstar=2.0; //=jet_ptmin
  jet_ratio=1e5;
  do_pairs=0;
  load_event=0;
  //track_secondaries=0;
  track_pairs=0;
  track_muons=0;
  store_pairs=0;
  store_muons=0;
  do_tertphot=0;
  pair_ratio=1.0;
  muon_ratio=1.0;
  muon_scale=1.0;
  bhabha_ratio=1.0;
  pair_ecut=EMASS;
  muon_ecut=MUMASS;
  integration_method=2;
  extra_grids=0;
  time_order=2;
  interpolation=2;
  adjust=0;
  geom=1;
  r_scal=1.0;
  jet_pythia=0;
  jet_select=1;
  pair_q2=1;
  load_photon=0;
  load_beam=0;
  cuts_from_loaded_beam = 0;
  bmt_precession_ = 0;
  ST_spin_flip_ = 0;
  automatic_grid_sizing = 0;
  emin=1.0;
  charge_sign=-1.0; charge_sign_0=-1.0;
  do_beamstrahlung_=1;
  store_beam=0;
  do_cross_gg=0;
  force_symmetric=0;
  do_isr=0;
  do_espread=0;
  espread1=0.0; espread2=0.0; which_espread1=0; which_espread2=0;
  do_coherent=0;
  do_trident=0;
  gg_smin=4.0*150.0*150.0;
  twobeam=0;

  do_edm4hep=0;

  // default values from readData.cc
  beam_pair=0;
  do_bhabhas=0;
  bhabha_ecmload=500;
  bhabha_scal=1.e-29;

  charge_symmetric=0; // to be checked
  do_bds_spin_rotation=0;

  do_compt=0;
  do_compt_phot=0;
  compt_emax=100;
  compt_scale=1.0;
  compt_x_min=1.0;
  
  do_dump=0;
  do_muons=0;
  // do_prod is not implemented
  do_prod=0; prod_e=0.0; prod_scal=1e-29;
  do_size_log=0;

  dump_particle=1;
  dump_step=1;
  ecm_min=0.0;
  ext_field=0;
  pair_step=1.0;
  // CLIC defaults:
  f_rep = 50; n_b=312;

  rndm_load=1; rndm_save=1; 
  rndm_seed=1;
}


// void SWITCHES::readFirstBeamParameters(const PARAMETERS& param)
// {
//    n_b= param.readDValue("n_b");
//    f_rep=param.readDValue("f_rep");
//    charge_sign=param.readFValue("charge_sign");

//   if(charge_sign>0.0){
//     charge_sign_0=1.0;
//   }
//   if(charge_sign<0.0){
//     charge_sign_0=-1.0;
//   }

// }

void SWITCHES::read(const PARAMETERS& param)
{

   n_b= param.readDValue("n_b");
   f_rep=param.readDValue("f_rep");
   charge_sign=param.readFValue("charge_sign");

  if(charge_sign>0.0){
    charge_sign_0=1.0;
  }
  if(charge_sign<0.0){
    charge_sign_0=-1.0;
  }


  integration_method = param.readIValue("integration_method");
  //  silent= param.readIValue("silent");


  extra_grids = param.readIValue("grids")-1;


  if (extra_grids<0) extra_grids=0;

  load_photon = param.readIValue("load_photons");
  write_photons = param.readIValue("store_photons");
  do_photons_[0] = param.readIValue("do_photons.1");

  do_photons_[1] = param.readIValue("do_photons.2");
  
  do_beamstrahlung_ = param.readIValue("do_eloss");

  ecm_min = param.readFValue("ecm_min"); 

  float temp = param.readFValue("ecm_min_gg");
  gg_smin=4*temp*temp;

  do_hadrons = param.readIValue("do_hadrons");

  store_hadrons = param.readIValue("store_hadrons");

  hadron_ratio = param.readFValue("hadron_ratio");

  do_jets = param.readIValue("do_jets");

  do_pairs = param.readIValue("do_pairs");

  beam_pair = param.readIValue("beam_pair");

  //store_secondaries = param.readIValue("store_secondaries");
  store_pairs = param.readIValue("store_pairs");

  do_muons = param.readIValue("do_muons");

  store_muons = param.readIValue("store_muons");

  do_coherent = param.readIValue("do_coherent");

  do_trident = param.readIValue("do_trident");

  emin = param.readFValue("emin");

  //track_secondaries = param.readIValue("track_secondaries");
  track_pairs = param.readIValue("track_pairs");

  track_muons = param.readIValue("track_muons");

  do_tertphot = param.readIValue("do_tertphot");

  pair_ecut = param.readFValue("pair_ecut");

  muon_ecut = param.readFValue("muon_ecut");

  pair_step = param.readFValue("pair_step");

  electron_ratio = param.readFValue("electron_ratio");

  do_lumi = param.readIValue("do_lumi");

  bhabha_scal = param.readFValue("bhabha_scal");

  bhabha_ecmload = param.readFValue("bhabha_ecmload");

  do_bhabhas = param.readIValue("do_bhabhas");

  rndm_save = param.readIValue("rndm_save");

  rndm_load = param.readIValue("rndm_load");

  rndm_seed = param.readIValue("rndm_seed");

  // do_lumi_ee_2 = param.readIValue("do_lumi_ee_2");

  do_size_log = param.readIValue("do_size_log");

  // lumi_ee_2_n = param.readIValue("lumi_ee_2_n");

  // lumi_ee_2_xmin=param.readFValue("lumi_ee_2_min"); 

  // lumi_ee_2_xmax=param.readFValue("lumi_ee_2_max"); 

  // do_lumi_eg_2=param.readIValue("do_lumi_eg_2");

  // lumi_eg_2_n=param.readIValue("lumi_eg_2_n");

  // lumi_eg_2_xmin=param.readFValue("lumi_eg_2_min"); 

  // lumi_eg_2_xmax=param.readFValue("lumi_eg_2_max"); 

  // do_lumi_ge_2=param.readIValue("do_lumi_ge_2");

  // lumi_ge_2_n=param.readIValue("lumi_ge_2_n");

  // lumi_ge_2_xmin=param.readFValue("lumi_ge_2_min"); 

  // lumi_ge_2_xmax=param.readFValue("lumi_ge_2_max"); 

  // do_lumi_gg_2=param.readIValue("do_lumi_gg_2");

  // lumi_gg_2_n=param.readIValue("lumi_gg_2_n");

  // lumi_gg_2_xmin=param.readFValue("lumi_gg_2_min"); 

  // lumi_gg_2_xmax=param.readFValue("lumi_gg_2_max"); 

  do_cross=param.readIValue("do_cross");

  do_prod=param.readIValue("do_prod");

  load_event=param.readIValue("load_events");

  prod_e=param.readFValue("prod_e");

  prod_scal=param.readFValue("prod_scal");

  do_compt=param.readIValue("do_compt");

  do_compt_phot=param.readIValue("do_compt_phot");

  compt_x_min=param.readFValue("compt_x_min");

  compt_scale=param.readFValue("compt_scale");

  compt_emax=param.readFValue("compt_emax");

  do_isr=param.readIValue("do_isr");

  do_espread=param.readIValue("do_espread");

  espread1=param.readFValue("espread.1");

  which_espread1=param.readIValue("which_espread.1");

  espread2=param.readFValue("espread.2");

  which_espread2=param.readIValue("which_espread.2");

  num_lumi=param.readIValue("num_lumi");

  num_lumi_eg=param.readIValue("num_lumi_eg");

  num_lumi_gg=param.readIValue("num_lumi_gg");

  lumi_p=param.readFValue("lumi_p");

  lumi_p_eg=param.readFValue("lumi_p_eg");

  lumi_p_gg=param.readFValue("lumi_p_gg");

  photon_ratio=param.readFValue("photon_ratio");

  pair_ratio=param.readFValue("pair_ratio");

  muon_ratio=param.readFValue("muon_ratio");

  muon_scale=param.readFValue("muon_scale");

  bhabha_ratio=param.readFValue("bhabha_ratio");

  jet_ratio=param.readFValue("jet_ratio");

  jet_pstar=param.readFValue("jet_ptmin");

  jet_store=param.readIValue("store_jets");

  jet_select=param.readIValue("jet_log");

  geom=param.readIValue("beam_size");

  r_scal=param.readFValue("beam_size_scale"); 

  ext_field=param.readIValue("ext_field");

//   if (!EXT_FIELD) 
//     {
//       if (ext_field){
// 	//	fprintf(stderr,"EXT_FIELD = false (LesDifines.h) \n");
// 	fprintf(stderr,"Cannot use flag ext_field\n");
// 	//	fprintf(stderr,"Please recompile with EXT_FIELD = true \n");
// 	exit(1);
//       }
      //    }
  pair_q2=param.readIValue("pair_q2");

  store_beam=param.readIValue("store_beam");

  load_beam=param.readIValue("load_beam");

  cuts_from_loaded_beam = param.readIValue("cuts_from_loaded_beam");

  bmt_precession_   = param.readIValue("bmt_precession");
  ST_spin_flip_   = param.readIValue("ST_spin_flip");

  automatic_grid_sizing = param.readIValue("automatic_grid_sizing");

  load_photon=param.readIValue("load_photons");

  jet_pythia=param.readIValue("jet_pythia");

  force_symmetric=param.readIValue("force_symmetric");

  charge_symmetric=param.readIValue("charge_symmetric");

  do_bds_spin_rotation=param.readIValue("do_bds_spin_rotation");

  // beam_vx_min=param.readFValue("beam_vx_min");

  // beam_vx_max=param.readFValue("beam_vx_max");

  // beam_vx_interval=param.readIValue("beam_vx_interval");

  // beam_vy_min=param.readFValue("beam_vy_min");

  // beam_vy_max=param.readFValue("beam_vy_max");

  // beam_vy_interval=param.readIValue("beam_vy_interval");

  do_dump=param.readIValue("do_dump");

  dump_step=param.readIValue("dump_step");

  dump_particle=param.readIValue("dump_particle");
  if ( dump_step <= 0 ) dump_step = 1;
  if ( dump_particle <= 0 ) dump_particle = 1;
  check_consistency();

  do_edm4hep=param.readFValue("do_edm4hep");
    }

void SWITCHES::check_consistency() const
{
  if ( do_bhabhas && do_pairs)
    {
      std::cout << " do_pairs= " << do_pairs << " do_bhabhas = " << do_bhabhas << std::endl;
      std::cout << " ERROR : it is not allowed to have do_bhabhas= 1 together with do_pairs = 1 " << std::endl;
      exit(0);
    } 
  if ( do_bhabhas && do_compt)
    {
      std::cout << " do_pairs= " << do_pairs << " do_compt = " << do_compt << std::endl;
      std::cout << " ERROR : it is not allowed to have do_bhabhas= 1 together with do_compt = 1 " << std::endl;
      exit(0);
    } 
  //  if (store_secondaries && !track_secondaries)
  if (store_pairs==1 && !track_pairs)
    {
      //std::cerr << " store_secondaries = " << store_secondaries << " track_secondaries = " << track_secondaries << std::endl;
      //std::cerr << " WARNING : it is not very consistent to store secondaries without tracking them! " << std::endl;
      std::cout << " store_pairs = " << store_pairs << " track_pairs = " << track_pairs << std::endl;
      std::cout << " WARNING : it is not very consistent to store pairs without tracking them! " << std::endl;
      exit(0);
    }
  if (do_muons==0 && (track_muons || store_muons))
    {
      std::cout << " do_muons = " << do_muons << " store_muons = " << store_muons << " track_muons = " << track_muons << std::endl;
      std::cout << " WARNING : it is not very consistent to track or store muons without generating them! " << std::endl;
      exit(0);
    }
  if (store_muons==1 && !track_muons)
    {
      std::cout << " store_muons = " << store_muons << " track_muons = " << track_muons << std::endl;
      std::cout << " WARNING : it is not very consistent to store muons without tracking them! " << std::endl;
      exit(0);
    }
  if (track_muons==1 && muon_ecut<10*MUMASS)
    {
      std::cout << " track_muons = " << track_muons << " muon_ecut = "<< muon_ecut <<std::endl;
      std::cout << " WARNING : Tracking for low energy muons not verified functional " << std::endl;

    }

  if (cuts_from_loaded_beam > 0  && load_beam == 0) 
    {
      std::cout << " WARNING : the switch cuts_from_loaded_beam is without effect with load_beam = 0 " << std::endl;
    }
  if(do_tertphot && !(track_pairs || track_muons))
    {
      std::cout << " do_tertphot = " << do_tertphot << " track_pairs = "<< track_pairs << " track_muons = " << track_muons << std::endl;
      std::cout << " WARNING : Tertphots are produced by incoherent particles. They should be tracked." << std::endl;
    }
}

void SWITCHES::readTWOBEAM(const PARAMETERS& param)
{
  //  VALUE value;
  //  load_variable("twobeam",&value);
  //  twobeam=CONTENTS(value);
  twobeam=param.readIValue("twobeam");
}

// void SWITCHES::readCharge_sign_2(const PARAMETERS& param)
// {
//   //  VALUE value;
//   //  load_variable("charge_sign_2",&value);
//   //  charge_sign_2=CONTENTS(value);
//   charge_sign_2=param.readFValue("charge_sign_2");
// }

std::string SWITCHES::output_flow() const 
{
  std::ostringstream out;
  out << title(std::string("SWITCHES : "));
  out << "charge_sign = " << charge_sign << std::endl;
  out << "bmt_precession = " << bmt_precession_ << std::endl;
  out << "ST_spin_flip = " << ST_spin_flip_ << std::endl;
  out << "do_bds_spin_rotation = " << do_bds_spin_rotation << std::endl;
  out << "automatic_grid_sizing = " << automatic_grid_sizing << std::endl;
  out << "integration_method = " << integration_method << " force_symmetric = " << force_symmetric << std::endl;
  out <<  "rndm_load = " << rndm_load << "rndm_save = " << rndm_save << " rndm_seed = " << rndm_seed << std::endl;
  out << "do_photons.1 = " << do_photons_[0] << " do_photons.2 = " << do_photons_[1] << std::endl;
  out << "write_photons = " << write_photons << std::endl;
  out << "do_comp = " << do_compt << " do_prod = " << do_prod << std::endl;
  out << "electron_ratio = " << electron_ratio << std::endl;
  out << "compt_x_min = " << compt_x_min << " compt_emax = " << compt_emax << " GeV ; compt_scale = " << compt_scale << std::endl;
  out << "do_lumi = " << do_lumi << "num_lumi = " << num_lumi << " lumi_p = " << lumi_p << std::endl;
  out << "do_cross = " << do_cross << " do_isr = "<< do_isr << " do_espread = " << do_espread << std::endl;
  out << "photon_ratio = " << photon_ratio << std::endl;
  out << "do_hadrons = " << do_hadrons << " store_hadrons = " << store_hadrons << " hadron_ratio = " << hadron_ratio << std::endl;
  out << "do_jets = " << do_jets << " store_jets = " << jet_store << std::endl;
  out << "do_pairs = " << do_pairs << " load_events = " << load_event << std::endl;
  //out << "track_secondaries = " << track_secondaries << " pair_step = " << pair_step << std::endl;
  //out << "store_secondaries = " << store_secondaries << std::endl;
  out << "track_pairs = " << track_pairs << " pair_step = " << pair_step << std::endl;
  out << "store_pairs = " << store_pairs << std::endl;
  out << "bhabha_scal = " << bhabha_scal << " bhabha_ecmload = " << bhabha_ecmload << " GeV " <<  std::endl;
  out << "do_muons = " << do_muons << " track_muons = " << track_muons << " store_muons = " << store_muons << std::endl;
  out << "muon_ratio = " << muon_ratio << " muon_scale = " << muon_scale << " muon_ecut = " << muon_ecut <<std::endl;
  out << "do_coherent = " << do_coherent << std::endl;
  out << "do_trident = " << do_trident << std::endl;
  out << "emin = " << emin << std::endl;
  out << "grids = " << extra_grids+1 << std::endl;
  out << "pair_ecut = " << pair_ecut << " GeV " << std::endl;
  out << "pair_ratio = " << pair_ratio << std::endl;
  out << "bhabha_ratio = " << bhabha_ratio << std::endl;
  out << "pair_q2 = " << pair_q2 << std::endl;
  out << "beam_pair = " << beam_pair << std::endl;
  out << "jet_ratio = " << jet_ratio << std::endl;
  out << "jet_ptmin = " << jet_pstar << std::endl;
  out << "jet_pythia = " << jet_pythia << std::endl;
  out << "jet_log = " << jet_select << std::endl;
  out << "beam_size = " << geom << " beam_size_scale = " << r_scal << " ext_field = " << ext_field << std::endl;
  out << "espread.1 = " << espread1 << " which_espread.1 = " << which_espread1 << " espread.2 = " << espread2 << " which_espread.2 = " << which_espread2 << std::endl;
  out << "f_rep = " << f_rep << " n_b = " << n_b << std::endl;
  out << "do_edm4hep = " << do_edm4hep << std::endl;

  return out.str();
}

#ifdef USE_EDM4HEP
void SWITCHES::EDM_output() const
{
  EDM4HEPWriter::edmfile().write_Metadata("charge_sign",charge_sign);
  EDM4HEPWriter::edmfile().write_Metadata("bmt_precession",bmt_precession_);
  EDM4HEPWriter::edmfile().write_Metadata("ST_spin_flip",ST_spin_flip_);
  EDM4HEPWriter::edmfile().write_Metadata("do_bds_spin_rotation",do_bds_spin_rotation);
  EDM4HEPWriter::edmfile().write_Metadata("automatic_grid_sizing",automatic_grid_sizing);
  EDM4HEPWriter::edmfile().write_Metadata("integration_method",integration_method);
  EDM4HEPWriter::edmfile().write_Metadata("force_symmetric",force_symmetric);
  EDM4HEPWriter::edmfile().write_Metadata("rndm_load",rndm_load);
  EDM4HEPWriter::edmfile().write_Metadata("rndm_save",rndm_save);
  EDM4HEPWriter::edmfile().write_Metadata("rndm_seed",rndm_seed);
  EDM4HEPWriter::edmfile().write_Metadata("do_photons.1",do_photons_[0]);
  EDM4HEPWriter::edmfile().write_Metadata("do_photons.2",do_photons_[1]);
  EDM4HEPWriter::edmfile().write_Metadata("write_photons",write_photons);
  EDM4HEPWriter::edmfile().write_Metadata("do_comp",do_compt);
  EDM4HEPWriter::edmfile().write_Metadata("do_prod",do_prod);
  EDM4HEPWriter::edmfile().write_Metadata("electron_ratio",electron_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("compt_x_min",compt_x_min);
  EDM4HEPWriter::edmfile().write_Metadata("compt_emax",compt_emax);
  EDM4HEPWriter::edmfile().write_Metadata("compt_scale",compt_scale);
  EDM4HEPWriter::edmfile().write_Metadata("do_lumi",do_lumi);
  EDM4HEPWriter::edmfile().write_Metadata("num_lumi",num_lumi);
  EDM4HEPWriter::edmfile().write_Metadata("lumi_p",lumi_p);
  EDM4HEPWriter::edmfile().write_Metadata("do_cross",do_cross);
  EDM4HEPWriter::edmfile().write_Metadata("do_isr",do_isr);
  EDM4HEPWriter::edmfile().write_Metadata("do_espread",do_espread);
  EDM4HEPWriter::edmfile().write_Metadata("photon_ratio",photon_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("do_hadrons",do_hadrons);
  EDM4HEPWriter::edmfile().write_Metadata("store_hadrons",store_hadrons);
  EDM4HEPWriter::edmfile().write_Metadata("hadron_ratio",hadron_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("do_jets",do_jets);
  EDM4HEPWriter::edmfile().write_Metadata("store_jets",jet_store);
  EDM4HEPWriter::edmfile().write_Metadata("do_pairs",do_pairs);
  EDM4HEPWriter::edmfile().write_Metadata("load_events",load_event);
  EDM4HEPWriter::edmfile().write_Metadata("track_pairs",track_pairs);
  EDM4HEPWriter::edmfile().write_Metadata("pair_step",pair_step);
  EDM4HEPWriter::edmfile().write_Metadata("store_pairs",store_pairs);
  EDM4HEPWriter::edmfile().write_Metadata("bhabha_scal",bhabha_scal);
  EDM4HEPWriter::edmfile().write_Metadata("bhabha_ecmload",bhabha_ecmload);
  EDM4HEPWriter::edmfile().write_Metadata("do_muons",do_muons);
  EDM4HEPWriter::edmfile().write_Metadata("track_muons",track_muons);
  EDM4HEPWriter::edmfile().write_Metadata("store_muons",store_muons);
  EDM4HEPWriter::edmfile().write_Metadata("muon_ratio",muon_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("muon_scale",muon_scale);
  EDM4HEPWriter::edmfile().write_Metadata("muon_ecut",muon_ecut);
  EDM4HEPWriter::edmfile().write_Metadata("do_coherent",do_coherent);
  EDM4HEPWriter::edmfile().write_Metadata("do_trident",do_trident);
  EDM4HEPWriter::edmfile().write_Metadata("emin",emin);
  EDM4HEPWriter::edmfile().write_Metadata("grids",extra_grids+1);
  EDM4HEPWriter::edmfile().write_Metadata("pair_ecut",pair_ecut);
  EDM4HEPWriter::edmfile().write_Metadata("pair_ratio",pair_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("bhabha_ratio",bhabha_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("pair_q2",pair_q2);
  EDM4HEPWriter::edmfile().write_Metadata("beam_pair",beam_pair);
  EDM4HEPWriter::edmfile().write_Metadata("jet_ratio",jet_ratio);
  EDM4HEPWriter::edmfile().write_Metadata("jet_ptmin",jet_pstar);
  EDM4HEPWriter::edmfile().write_Metadata("jet_pythia",jet_pythia);
  EDM4HEPWriter::edmfile().write_Metadata("jet_log",jet_select);
  EDM4HEPWriter::edmfile().write_Metadata("beam_size",geom);
  EDM4HEPWriter::edmfile().write_Metadata("beam_size_scale",r_scal);
  EDM4HEPWriter::edmfile().write_Metadata("ext_field",ext_field);
  EDM4HEPWriter::edmfile().write_Metadata("espread.1",espread1);
  EDM4HEPWriter::edmfile().write_Metadata("which_espread.1",which_espread1);
  EDM4HEPWriter::edmfile().write_Metadata("espread.2",espread2);
  EDM4HEPWriter::edmfile().write_Metadata("which_espread.2",which_espread2);
  EDM4HEPWriter::edmfile().write_Metadata("f_rep",f_rep);
  EDM4HEPWriter::edmfile().write_Metadata("n_b",n_b);
  EDM4HEPWriter::edmfile().write_Metadata("do_edm4hep",do_edm4hep);
}
#endif
