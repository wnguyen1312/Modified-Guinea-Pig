#ifndef SWITCHES_SEEN
#define SWITCHES_SEEN

#include <string>

#include "define.h"
#include "abstractIOclass.h"
#include "parametersCPP.h"

/** 
 * class that holds global parameters
 * most parameters are explained in the manual
 */

class SWITCHES : public ABSTRACT_IO_CLASS
{
  //  int silent;
  int electron_distribution_rho;
  //   int electron_distribution_scatter;
  int do_lumi,num_lumi,num_lumi_eg,num_lumi_gg;
  float lumi_p,lumi_p_eg,lumi_p_gg;
  float electron_ratio;

  // 2 dimensions was enough, but the original C had 3 dimensions
  // to be seen..
  int do_photons_[2],write_photons;
  int photon_distribution,do_beamstrahlung_;
  float photon_ratio;
  int do_hadrons,store_hadrons;
  float hadron_ratio;
  int do_jets;
  float jet_pstar;
  float jet_ratio;
  //  int do_pairs,track_secondaries,do_muons,store_secondaries;
  int do_pairs,track_pairs,store_pairs;
  int do_muons,track_muons,store_muons;
  int do_tertphot;
  float muon_ratio,muon_scale,muon_ecut;
  float pair_ratio;
  float bhabha_ratio;
  /// pair_ecut: minimal energy in GeV the particles from pair creation need to have to be stored
  /// pair_step: scaling factor for the step size of the pairs; if the value is increased the step size is decreased
  float pair_ecut,pair_step;
  int integration_method; /*!< integration method in transverse plane 1=direct, 2= fast Fourier (default), 3 = iterative */
  int extra_grids;
  int time_order;
  int interpolation;
  int adjust;
  /// ext_field, if not 0 the program takes into account the effect of the magnetic field for the equivalent photon approximation
  int geom,ext_field,beam_pair;
  int twobeam;
  float r_scal;
  int jet_pythia,jet_select,jet_store;
  int pair_q2;
  int load_photon,load_beam,load_event;
  int cuts_from_loaded_beam;
  int bmt_precession_;
  int ST_spin_flip_;
  int automatic_grid_sizing;
  /// charge_sign is the relative sign of the charge of the two beams -1 is for e+e- and +1 for e-e-. If set to 0 no beam-beam force is assumed
  float emin,charge_sign,charge_sign_0;
  // float charge_sign_2;
  int store_beam;
  int do_cross,do_cross_gg,do_isr,do_espread,do_prod;
  /// RMS value of energy spread of beam 1 and 2, distribution set with which_espread
  float espread1,espread2;
  /// energy spreading method, used in gridCPP.h
  /// 0: no energy spread 1: flat distribution 2: two peaks 3: Gaussian
  int which_espread1,which_espread2;
  /// if force_symmetric is not equal to 0, the beams are assumed to be up-down and left-right symmetric
  int force_symmetric,charge_symmetric;
  int do_bds_spin_rotation;
  int do_coherent,do_compt,do_compt_phot,do_trident;
  /// load and save random number generator state
  int rndm_load,rndm_save;
  /// initial random seed
  unsigned long rndm_seed;
  float gg_smin,compt_x_min,compt_emax;
  /* int do_lumi_ee_2,lumi_ee_2_n; */
  /* float lumi_ee_2_xmin,lumi_ee_2_xmax; */
  /* int do_lumi_eg_2,lumi_eg_2_n; */
  /* float lumi_eg_2_xmin,lumi_eg_2_xmax; */
  /* int do_lumi_ge_2,lumi_ge_2_n; */
  /* float lumi_ge_2_xmin,lumi_ge_2_xmax; */
  /* int do_lumi_gg_2,lumi_gg_2_n; */
  /* float lumi_gg_2_xmin,lumi_gg_2_xmax; */
  int do_bhabhas;
  float bhabha_scal, bhabha_ecmload,ecm_min;
  /// values related to do_prod, not implemented
  float prod_e,prod_scal;
  float compt_scale;
  //  int hist_ee_bins,hist_espec_bins;
  //  float hist_ee_min,hist_ee_max,hist_espec_min,hist_espec_max;
  int do_size_log;
  //  float beam_vx_min,beam_vx_max,beam_vy_min,beam_vy_max;
  /// rep rate and number of bunches (only used in some output results)
  double f_rep,n_b;
  int do_dump,dump_step,dump_particle;
  //  int beam_vx_interval,beam_vy_interval;
  int do_edm4hep;

  void check_consistency() const;

 public:

  SWITCHES();

  //  void readFirstBeamParameters(const PARAMETERS& param);
  // void readBeamParametersContinued(const PARAMETERS& param);
  void read(const PARAMETERS& param);
  void readTWOBEAM(const PARAMETERS& param);
  //  void readCharge_sign_2(const PARAMETERS& param);

  virtual  std::string output_flow() const;

#ifdef USE_EDM4HEP
  void EDM_output() const;
#endif

  inline int getTWOBEAM() const {return twobeam;};
  inline int get_rndm_load() const {return rndm_load;};
  inline int get_rndm_save() const {return rndm_save;};
  inline unsigned long   get_rndm_seed() const {return rndm_seed;};
  inline int get_write_photons() const {return write_photons;};
  inline int get_do_lumi() const {return do_lumi;};
  /* inline int get_do_lumi_ee_2() const {return do_lumi_ee_2;}; */
  /* inline int get_do_lumi_ge_2() const {return do_lumi_ge_2;}; */

  inline int get_do_bhabhas() const {return do_bhabhas;};

inline int get_do_edm4hep() const {return do_edm4hep;};


  /* inline int get_do_lumi_eg_2() const {return do_lumi_eg_2;}; */
  inline int get_num_lumi() const {return num_lumi;};
  inline int get_num_lumi_eg() const {return num_lumi_eg;};
  inline int get_num_lumi_gg() const {return num_lumi_gg;};
  inline int get_do_compt_phot() const {return do_compt_phot;};
  inline int get_do_size_log() const {return do_size_log;};
  inline int get_do_pairs() const {return do_pairs;};
  inline int get_do_hadrons() const {return do_hadrons;};
  inline int get_do_compt() const {return do_compt;};
  inline int get_do_muons() const {return do_muons;};
  inline int get_integration_method() const {return integration_method;};
  inline int get_extra_grids() const { return extra_grids;};
  inline int get_load_photon() const {return load_photon;};
  inline int get_load_beam() const {return load_beam;};
  inline int get_automatic_grid_sizing() const {return automatic_grid_sizing;}
  inline int get_cuts_from_loaded_beam() const {return cuts_from_loaded_beam;}
  inline int get_bmt_precession() const {return bmt_precession_;}
  inline int get_ST_spin_flip() const {return ST_spin_flip_;}
  inline int get_charge_symmetric() const {return charge_symmetric;};
  inline int get_do_bds_spin_rotation() const {return do_bds_spin_rotation;};
  //  inline int get_electron_distribution_scatter() const { return electron_distribution_scatter;};
  inline int get_adjust() const {return adjust;};
  inline int get_time_order() const {return time_order;};
  //  inline int get_track_secondaries() const {return track_secondaries;};
  inline int get_track_pairs() const {return track_pairs;};
  inline int get_track_muons() const {return track_muons;};
  inline int get_do_tertphot() const {return do_tertphot;};
  inline int get_load_event() const {return load_event;};
  inline int get_electron_distribution_rho() const {return electron_distribution_rho;};
  inline int get_force_symmetric() const {return force_symmetric;};
  inline int get_geom() const {return geom;};
  inline int get_ext_field() const {return ext_field;};
  inline int get_do_espread() const { return do_espread;};
  inline int get_which_espread1() const { return which_espread1;};
  inline int get_which_espread2() const { return which_espread2;};
  inline int get_do_isr() const { return do_isr;};
  inline int get_do_cross() const {return do_cross;};
  inline int get_do_coherent() const {return do_coherent;};
  inline int get_do_trident() const {return do_trident;};
  inline int get_do_jets() const {return do_jets;};
  inline int get_do_prod() const {return do_prod;};
  inline int get_interpolation() const {return interpolation;};
  inline int get_do_beamstrahlung() const {return do_beamstrahlung_;};
  inline int get_do_photons(int nbeam) const {return do_photons_[nbeam-1];};
  inline int get_jet_pythia() const {return jet_pythia;};
  inline int get_jet_select() const {return jet_select;};
  inline int get_jet_store() const {return jet_store;};
  //  inline int get_store_secondaries() const {return store_secondaries;}
  inline int get_store_pairs() const {return store_pairs;}
  inline int get_store_muons() const {return store_muons;}
  /* inline int get_do_lumi_gg_2() const {return do_lumi_gg_2;} */
  inline int get_store_hadrons() const {return store_hadrons;}
  inline int get_pair_q2() const {return pair_q2;}
  inline int get_beam_pair() const {return beam_pair;}
  inline int get_store_beam() const { return store_beam;}
  inline int get_do_dump() const { return do_dump;}
  inline int get_dump_step() const { return dump_step;}
  inline int get_dump_particle() const { return dump_particle;}
  inline float get_lumi_p() const {return lumi_p;};
  inline float get_lumi_p_eg() const {return lumi_p_eg;};
  inline float get_lumi_p_gg() const {return lumi_p_gg;};
  inline float get_charge_sign() const {return charge_sign;};
  inline float get_charge_sign_0() const {return charge_sign_0;}
  inline float get_pair_ecut() const {return pair_ecut;};
  inline float get_muon_ecut() const {return muon_ecut;};
  inline float get_pair_step() const {return pair_step;};
  inline float get_electron_ratio() const {return electron_ratio;};
  inline float get_r_scal() const {return r_scal;};
  inline float get_compt_x_min() const {return compt_x_min;};
  inline float get_espread1() const { return espread1;};
  inline float get_espread2() const { return espread2;};
  inline float get_ecm_min() const {return ecm_min;};
  inline float get_emin() const {return emin;};
  inline float get_photon_ratio() const {return photon_ratio;};
  inline float get_jet_pstar() const {return jet_pstar;};
  inline float get_jet_ratio() const {return jet_ratio;};
  inline float get_pair_ratio() const {return pair_ratio;}
  inline float get_muon_ratio() const {return muon_ratio;}
  inline float get_muon_scale() const {return muon_scale;}
  inline float get_bhabha_ratio() const {return bhabha_ratio;}
  inline float get_compt_scale() const {return compt_scale;}
  inline float get_compt_emax() const {return compt_emax;}
  inline float get_gg_smin() const {return gg_smin;}
  inline float get_hadron_ratio() const {return hadron_ratio;}
  inline float get_bhabha_scal() const {return bhabha_scal;};
  inline float get_bhabha_ecmload() const {return bhabha_ecmload;};
  inline double get_f_rep() const {return f_rep;}
  inline double get_n_b() const {return n_b;}

};

#endif
