#include "guineapigCPP.h"
#include "option_args.h"
#include "EDMwriter.h"

GUINEA::GUINEA(char *name):cross_(NULL),number_of_stored_particles1_(0),number_of_stored_particles2_(0)
{
  time_counter_ = 0;
  std::cout << header();
  init_input_file_names();
  beam_parameters1_.setLabel('1');
  beam_parameters2_.setLabel('2');
  secondaries_.set_name(std::string("pairs"));
  muons_.set_name(std::string("muons"));
  params_.read_accelerator(name); 
}

void GUINEA::run( char *par,char *prot)
{
  //  the order of the following sequence must not be changed : 
  // reading ACCELERATOR parameters
  //  params_.read_accelerator(name); 
  // setting, from previously read data, of beam characteristics (emittance,
  // sigmas a.s;
  beam_parameters1_.read(params_);
  beam_parameters2_.read(params_);

  // reading 'PARAMETERS' : switches and grid parameters
  params_.read_parameters(par );
  switches.read(params_);
#ifdef USE_EDM4HEP
  if (switches.get_do_edm4hep()) EDM4HEPWriter::initialize("output.root", switches);
#endif
  // the next method needs the switches
  init_output_file_names();

#ifdef TWOBEAM
  std::cerr << " THE TWOBEAM IS NOT PROGRAMMED " << std::endl;
  exit(0);
#endif

  grid_.read(params_, switches.get_automatic_grid_sizing());
  //      std::cout << " guineapig : grid parameters "  << " nmacros " << grid_.get_nb_macroparticles(1) << std::endl;

  beam1_.connect_parameters(&beam_parameters1_);
  beam2_.connect_parameters(&beam_parameters2_);
  grid_.connect_beams(&beam1_, &beam2_);
  generator_ = RNDM(switches.get_rndm_seed());
  std::cout << "  the random number seed is : " << switches.get_rndm_seed() << std::endl;
  grid_.connect_random_generator(&generator_);
  if ( !check_parameters() ) exit(1);
  simulate();
  grid_.set_bpm_signals();
  outputs (std::string(prot) );
  if (switches.get_rndm_save()) generator_.rndm_save();
#ifdef USE_EDM4HEP
  if (switches.get_do_edm4hep()) EDM4HEPWriter::edmfile().close();
#endif
}

void GUINEA::save_results_on_files()
{
  grid_.save_lumi_on_files(switches, lumi_ee_file_, lumi_eg_file_, lumi_ge_file_, lumi_gg_file_);
  if(switches.get_do_tertphot())  grid_.save_tertphot_on_file(tertphot_file_);
#ifdef USE_EDM4HEP
  if(switches.get_do_jets()) std::cout << "Saving jets to event data model" <<std::endl;
#endif
  if (switches.get_store_beam()) 
    {
      float gridMaxZ = grid_.get_max_z();
      float gridStep = grid_.get_step();
      int gridTimestep = grid_.get_timestep();
      float gridScalStep[2];
      gridScalStep[0] = grid_.get_scal_step(0);
      gridScalStep[1] = grid_.get_scal_step(1);
      beam1_.backstep2(1, gridMaxZ, gridStep,  gridTimestep, gridScalStep);
      beam2_.backstep2(2, gridMaxZ, gridStep,  gridTimestep, gridScalStep);
      number_of_stored_particles1_ = beam1_.store_beam(beam1_file_);
      number_of_stored_particles2_ = beam2_.store_beam(beam2_file_);
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep())
	{
	  beam1_.store_beam_EDM(beam1);
	  beam2_.store_beam_EDM(beam2);
	}
#endif
      if ( switches.get_do_coherent() )
	{
	  beam1_.store_coherent_beam(coh1_file_);
	  beam2_.store_coherent_beam(coh2_file_);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      beam1_.store_beam_EDM(coh1);
	      beam2_.store_beam_EDM(coh2);
	    }
#endif
	}
      if ( switches.get_do_trident() )
	{
	  beam1_.store_trident_beam(tri1_file_);
	  beam2_.store_trident_beam(tri2_file_);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      beam1_.store_beam_EDM(tri1);
	      beam2_.store_beam_EDM(tri2);
	    }
#endif
	}
      if ( switches.get_write_photons() )
        {
	  FILE_IN_OUT filout;
	  filout.open_file(photon_file_, "w");
	  for (int h=0; h < grid_.get_n_cell_z(); h++)
	    {	    
	      std::vector<PHOTON>& thePhotons1 =  beam1_.getPhotonVector(h);
	      for(unsigned int k=0;k<thePhotons1.size();k++)
		{
		  if (thePhotons1[k].energy()>0.0){
		    filout.save_object_on_persistent_file( &thePhotons1[k]);
#ifdef USE_EDM4HEP
		    if (switches.get_do_edm4hep()) 
		      {
			EDM4HEPWriter::edmfile().save_photons1_EDM(thePhotons1[k]);
		      }
#endif
		  }
		}
	      std::vector<PHOTON>& thePhotons2 =  beam2_.getPhotonVector(h);
	      for(unsigned int k=0;k<thePhotons2.size();k++){
		float ene = thePhotons2[k].energy();
		if (thePhotons2[k].energy()>0.0)
		  {
		    thePhotons2[k].setEnergy(-ene);
		    filout.save_object_on_persistent_file(&thePhotons2[k]);
#ifdef USE_EDM4HEP
		    if (switches.get_do_edm4hep())
		      {
			EDM4HEPWriter::edmfile().save_photons2_EDM(thePhotons2[k]);
		      }
#endif
		  }
	      }
	    // photon_vector_.store_photon(photon_file_);
	    }
	  filout.close();
        }
    }
  //   if (switches.get_track_secondaries() || switches.get_store_secondaries()) secondaries_.save_pairs_on_file(secondaries_file_);
  //   if (switches.get_store_secondaries() > 1) secondaries_.save_pairs0_on_file(secondaries0_file_);
  if (switches.get_do_bhabhas())
    {
      grid_.get_bhabhas().save_on_files(switches.get_do_bhabhas(), switches.get_do_edm4hep(), bhabha_prod_, bhphoton_prod_,bhphotons_);
      if (switches.get_track_pairs() || switches.get_store_pairs()) 
	{
	  secondaries_.save_bhabhas_on_file(secondaries_file_);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep()) 
	    {
	      secondaries_.save_pairs_on_EDMcollection(bhabhaPair);
	    }
#endif
	}
      if (switches.get_store_pairs() > 1){
	secondaries_.save_bhabhas0_on_file(secondaries0_file_);
#ifdef USE_EDM4HEP
	if (switches.get_do_edm4hep())
	  {
	    secondaries_.save_pairs0_on_EDMcollection(bhabhaPair0);
	  }
#endif
      }
    }
  else
    {
      if (switches.get_track_pairs() || switches.get_store_pairs()) 
	{
	  secondaries_.save_pairs_on_file(secondaries_file_);
#ifdef USE_EDM4HEP
	  if (switches.get_do_edm4hep())
	    {
	      secondaries_.save_pairs_on_EDMcollection(electronPair);
	    }
#endif
	}
      if (switches.get_store_pairs() > 1) {
	secondaries_.save_pairs0_on_file(secondaries0_file_);
#ifdef USE_EDM4HEP
	if (switches.get_do_edm4hep())
	  {
	    secondaries_.save_pairs0_on_EDMcollection(electronPair0);
	  }
#endif
      }
    }
  if (switches.get_track_muons() || switches.get_store_muons()) 
    {
      muons_.save_pairs_on_file(muons_file_);
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep())
	{
	  muons_.save_pairs_on_EDMcollection(muonPair);
	}
#endif
    }
  if (switches.get_store_muons() > 1) 
    {
      muons_.save_pairs0_on_file(muons0_file_);
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep())
	{
	  muons_.save_pairs0_on_EDMcollection(muonPair0);
	}
#endif
    }
}

void GUINEA::outputs(std::string nameOfProtokoll)
{
  float miss1, miss2;
  long out1, out2;
  grid_.get_miss(miss1, out1, miss2, out2);
  if (miss1 > 0.1 ||   miss2 > 0.1)
    {
      std::cerr << " ---------------------------------------------------------- " << std::endl;
      std::cerr << " WARNING : there is a quite important quantity of particles travelling out of grid " << std::endl;
      std::cerr  << "miss.1 = " << miss1 << " out.1 = " << out1 << " miss.2 = " << miss2 << " out.2 = " << out2 << std::endl; 
      std::cerr << " maybe the grid dimensions should be rematched " << std::endl;
      std::cerr << " ---------------------------------------------------------- " << std::endl;
    }
  save_results_on_files();
  print_program_outputs(nameOfProtokoll);
  //for(int i1=1;i1<=6;i1++) time_.print_timer(i1);
  //time_.print_timer_all(6);
}

bool GUINEA::check_parameters() const
{
  int k;
  bool check = true;
  //std::cerr << " checking parameters... " << std::endl;
#ifdef USE_FFT_LOCAL
  int vnx = TOOLS::check_power_of_2(grid_.get_n_cell_x() );
  int vny = TOOLS::check_power_of_2(grid_.get_n_cell_y() );
  if (switches.get_integration_method() == 2)
    {
      if (!vnx || !vny  )
	{
	  std::cerr << " ERROR : with integration_method = 2 (FFT_LOCAL) the cell numbers must be power of 2 " << std::endl;
	  std::cerr << " Consider compiling with FFTW2 or FFTW3 (also for cpu-speed)" << std::endl;
	  std::cerr << " n_x = " << grid_.get_n_cell_x() << " n_y = " << grid_.get_n_cell_y()  << std::endl;
	  check = false;
	}
    }
#endif
  if ( switches.get_do_prod() )
    {
      std::cerr << " ERROR : key_word do_prod not completely implemented " << std::endl;
      // see grid::move_particles for the do_beamstrahlung case
      check = false;
    }
  if (beam_parameters1_.dist_x() != 0 || beam_parameters2_.dist_x() != 0)
    {
      std::cerr << " dist_x other than 0 is not available " << std::endl;
      std::cerr << " dist_x.1 = " << beam_parameters1_.dist_x() << " dist_x.2 = "
	   <<  beam_parameters2_.dist_x() << std::endl;
      check = false;
    }

  bool do_photons = switches.get_do_photons(1)>0 || switches.get_do_photons(2)>0;
  bool load_photon = switches.get_load_photon()>0;
  if (switches.get_do_coherent() && !do_photons && !load_photon)
    {
      std::cerr << " GUINEA::check_parameters:: it is not useful to have do_coherent=1 with no do_photons nor load_photons " << std::endl;
    }
  //if ( switches.get_track_secondaries() > 0 )
  if ( switches.get_track_pairs() > 0 )
    {
      bool do_pair_compt_bha = switches.get_do_pairs()>0 || switches.get_do_compt()>0 || switches.get_do_bhabhas() >0;
      if (!do_pair_compt_bha) 
	{
	  //std::cerr << " with track_secondaries > 0 there must be  either do_pairs = 1 or do_compt = 1 or do_bhabhas = 1 " << std::endl;
	  std::cerr << " with track_pairs > 0 there must be  either do_pairs = 1 or do_compt = 1 or do_bhabhas = 1 " << std::endl;
	  check = false;

	}
    }

  const BEAM_PARAMETERS* bpPtr[2]  = {&beam_parameters1_, &beam_parameters2_};  
  float cutz = params_.readFValue("cut_z");
  if ( !switches.get_load_beam() )
    {
      // program generated beam
      //       if ( cutz <= 0 ) 
      // 	{
      // 	  std::cerr << "  GUINEA::check_parameters: the data of cut_z is mandatory. program stop " << std::endl;
      // 	  check = false;      
      // 	}
      if ( switches.get_cuts_from_loaded_beam() ) 
	{
	  std::cerr << " GUINEA::check_parameters: WARNING  with beams generated by program the keyword cuts_from_loaded_beam is without effect "<< std::endl;
	}

      for (k=0; k < 2; k++)
	{
	  if ( bpPtr[k]->sigma_z() <= 0.0 )
	    {
	      std::cerr << "  GUINEA::check_parameters: the data of sigma_z is mandatory. program stop " << std::endl;
	      check = false;
	    }
	}

      for (k=0; k < 2; k++)
	{
	  int number=0;
	  if (bpPtr[k]->beta_x() > 0.0 ) number++;
	  if (bpPtr[k]->em_x() > 0.0 ) number++;
	  if (bpPtr[k]->sigma_x() > 0.0 ) number++;
	  if (number < 2 )
	    {
	      std::cerr << " GUINEA::check_parameters : ERROR among beta_x, emitt_x, sigma_x only one parameter is given for beam " << k+1 << " . Two are needed " << std::endl;
	      check = false;
	    }
	  number = 0;
	  if (bpPtr[k]->beta_y() > 0.0 ) number++;
	  if (bpPtr[k]->em_y() > 0.0 ) number++;
	  if (bpPtr[k]->sigma_y() > 0.0 ) number++;
	  if (number < 2 )
	    {
	      std::cerr << " GUINEA::check_parameters : ERROR among beta_y, emitt_y, sigma_y only one parameter is given. Two are needed  " << std::endl;
	      check = false;
	    }
	}
    }
  else
    {
      // loaded beam 
      if (switches.get_automatic_grid_sizing() )
	{
	  if ( switches.get_cuts_from_loaded_beam() ) 
	    {
	      std::cerr << " GUINEA::check_parameters: WARNING  with automatic_grid_sizing,  the keyword cuts_from_loaded_beam is without effect "<< std::endl;
	    }
	}
      else
	{
	  if ( params_.readFValue("cut_x") <= 0.0 || params_.readFValue("cut_y") <= 0.0 )
	    {
	      bool test = false;
	      for (k=0; k < 2; k++)
		{
		  int number=0;
		  if (bpPtr[k]->beta_x() > 0.0 ) number++;
		  if (bpPtr[k]->em_x() > 0.0 ) number++;
		  if (bpPtr[k]->sigma_x() > 0.0 ) number++;
		  if (number < 2 )
		    {
		      test = true;
		      check = false;
		    }
		  number = 0;
		  if (bpPtr[k]->beta_y() > 0.0 ) number++;
		  if (bpPtr[k]->em_y() > 0.0 ) number++;
		  if (bpPtr[k]->sigma_y() > 0.0 ) number++;
		  if (number < 2 )
		    {
		      test = true;
		      check = false;
		    }
		}
	      if (test ) 
		{
		  std::cerr << " GUINEA::check_parameters : ERROR cut_x or cut_y is missing. In that case the beam transverse characteristics (beta, sigma, emitt) must be defined. Else use the option automatic_grid_sizing = 1  " << std::endl;
		}
	    }
	  if ( cutz <= 0.0) 
	    {
	      bool test = false;
	      for (k=0; k < 2; k++)
		{
		  if ( bpPtr[k]->sigma_z() <= 0.0 )
		    {
		      test = true;
		      check = false;
		    }
		}
	      if (test) 
		{
		  std::cerr << "  GUINEA::check_parameters: one of sigma_z and cut_z must be given. program stop " << std::endl;
		}
	    }
	}
    }
  if ( !grid_.random_ok() )
    {
      std::cerr << " GUINEA::check_parameters: GRID seems not to be connected to the random generator " << std::endl;
      check = false;
    }
  if (switches.get_ST_spin_flip()  != 0) 
    {
      if (switches.get_bmt_precession() == 0) 
	{
	  std::cerr << " with ST_spin_flip, bmt_precession = 1 is MANDATORY " << std::endl;
	  check = false;
	}
    }

#ifndef USE_EDM4HEP
  if (switches.get_do_edm4hep()) {
    std::cerr << "GUINEA::check_parameters: EDM4hep was not compiled but do_edm4hep is enabled"  << std::endl;
    check = false;
  }
#endif

  //std::cerr << " end of checking ... " << std::endl;
  return check;
}

// void GUINEA::close()
// {
//     if (switches.get_rndm_save()) generator_.rndm_save();
// }

void GUINEA::set_simulation()
{
  grid_.check_distribute(0);
  int n_slice=grid_.get_n_cell_z();
  beam1_.make_beam(n_slice, switches.get_bmt_precession(), &generator_);
  beam2_.make_beam(n_slice, switches.get_bmt_precession(), &generator_);

  if (switches.get_do_size_log())
    {
      beam1_.write_size_init(std::string("beamsize"));
      beam2_.write_size_init(std::string("beamsize"));
    }

  if (switches.get_rndm_load()) 
    {
      generator_.rndm_load();
    }
  if (switches.get_load_event())
    {
      secondaries_.set_load_file(event_input_file_);
    }


  if (switches.get_do_jets())
    {
      grid_.init_jet(4.0*beam_parameters1_.ebeam()*beam_parameters2_.ebeam(), switches.get_jet_pstar(), switches.get_do_jets(),switches.get_jet_pythia(), switches.get_jet_select(), jet_file_);
    }
  if (switches.get_do_cross())
    {
      if (MCROSS) 
	{
	  if (AVERCROSS) cross_ = new maverCROSS(cross_input_file_);
	  else cross_ = new mCROSS(cross_input_file_);
	}
      else 
	{
	  if (AVERCROSS) cross_ = new averCROSS(cross_input_file_);
	  else cross_ = new CROSS(cross_input_file_);
	}
    }
  else cross_ = NULL;
  grid_.generalInit(switches, photon_file_, hadronfile_, compton_phot_file_,bhabha_input_file_, bhabha_photon_input_file_, cross_);
  //  if(switches.get_do_pairs()||switches.get_do_hadrons()||switches.get_do_compt()
  //       ||switches.get_do_muons())
  //    {
  //     grid_.init_extra_photons();
  //    }
}

void GUINEA::set_beams_and_grids()
{
  std::vector<float> size_x, size_y;
  float size_z;
  int updated_n_cell_x = grid_.get_n_cell_x();
  int updated_n_cell_y = grid_.get_n_cell_y();
  int load_beam = switches.get_load_beam();
  if (load_beam) 
    {
      switch (load_beam) 
	{
	case 1:
	  std::cerr << " GUINEA::set_beams_and_grids(): load_beam = 1 not yet implemented " << std::endl;
	  exit(0);
	case 2:
	  std::cerr << " GUINEA::set_beams_and_grids(): load_beam = 2 not yet implemented " << std::endl;
	  exit(0);
	case 3:
	  {
	    // these pointers will be deleted a the end of loading particles, 
	    // by the beam objects.
	    BEAM_FROM_FILE* bffe = new BEAM_FROM_FILE(get_elfilename());
	    BEAM_FROM_FILE* bffp = new BEAM_FROM_FILE(get_posfilename());
	    xycuts_for_grids(*bffe, *bffp, switches.get_extra_grids()+1, size_x, size_y, updated_n_cell_x, updated_n_cell_y);
	    size_z = zcut_for_grids(*bffe, *bffp);

	    transfer_loaded_particles_in_beam( beam1_, beam_parameters1_, bffe, size_z); 
	    transfer_loaded_particles_in_beam( beam2_, beam_parameters2_, bffp, size_z); 
	    init_grid_phys(beam_parameters1_.n_particles(), beam_parameters2_.n_particles(),size_x, size_y, size_z, updated_n_cell_x, updated_n_cell_y);
	       beam_displacements_from_data(grid_.get_delta_z()); 
	  }
	  break;

	default:
	  std::cerr << "load_beam set to incorrect value: " << std::endl;
	  std::cerr << " load_beam= " << load_beam << std::endl;
	  std::cerr << "must be within [0..2] " << std::endl;
	  exit(1);

	} 
    }
  else 
    {


      //      int  n_cell_z= params_.readIValue("n_z");
      int  n_cell_z = grid_.get_n_cell_z();
      if (n_cell_z <= 0) 
	{
	  std::cerr << " GUINEA::set_beams_and_grids() : ERROR : n_cell_z = " << n_cell_z << std::endl;
	  exit(0);
	}
      size_z = cut_z_from_data(switches.get_automatic_grid_sizing());
      float deltaz=2.0*size_z/((float)n_cell_z);
      beam1_.init_particles(grid_.get_nb_macroparticles(1),deltaz, switches.get_charge_symmetric(), switches.get_do_bds_spin_rotation());
      beam2_.init_particles(grid_.get_nb_macroparticles(2),deltaz, switches.get_charge_symmetric(), switches.get_do_bds_spin_rotation());
      beam_displacements_from_data(deltaz); 


      xycuts_for_grids( beam1_.particle_beam(), beam2_.particle_beam(),switches.get_extra_grids()+1, size_x, size_y, updated_n_cell_x, updated_n_cell_y);
      init_grid_phys(beam_parameters1_.n_particles(), beam_parameters2_.n_particles(), size_x,size_y,size_z,updated_n_cell_x, updated_n_cell_y);
    }

  if(switches.get_load_photon()) 
    {
      beam1_.load_photons(photon_input_file_, 1, grid_.get_delta_z(), grid_.get_max_z(), grid_.get_n_cell_z());
      beam2_.load_photons(photon_input_file_,2,  grid_.get_delta_z(), grid_.get_max_z(), grid_.get_n_cell_z());
    }

}

void GUINEA::printInitialBeam(const BEAM& beam)
{
  float bidonx, bidony, bidz, bidsigz;
  float bidsigx, bidsigy;
  float bidbetax, bidbetay;
  std::cout << " ******* initial beam characteristics : " << beam.label() << " ************** " << std::endl;
  beam.transverse_sigmas(bidsigx, bidsigy);
  beam.emittances(bidonx, bidony);
  beam.beamZRms(bidz, bidsigz);
  std::cout << " sigmax = " << bidsigx << " sigmay = " << bidsigy << std::endl;
  std::cout << " emittx =" << bidonx << " emitty = " << bidony << std::endl;
  std::cout << " Z0 = " << bidz << " sigZ = " << bidsigz << std::endl;
  bidbetax =  bidsigx* bidsigx*1e-9/(bidonx*EMASS/beam1_.get_ebeam());
  bidbetay =  bidsigy* bidsigy*1e-9/(bidony*EMASS/beam1_.get_ebeam());
  std::cout << " betax = " << bidbetax << " betay = " << bidbetay << std::endl;	
  std::cout << " ------------------------------------------------- " << std::endl;
}


void GUINEA::set_output_data_and_files()
{
  if (switches.get_do_pairs() || switches.get_do_compt() || switches.get_do_bhabhas()  )
    {
      secondaries_.resize(grid_.get_n_cell_z());
    }
  if (switches.get_do_muons() )
    {
      muons_.resize(grid_.get_n_cell_z());
    }
}

void GUINEA::simulate()
{
  //int i1;
  int gridTimestep;
  float gridScalStep[2], gridMaxZ, gridStep;  
  PHI_FLOAT sor_parameter[6];
  set_simulation();
  set_beams_and_grids();
  set_output_data_and_files();
  grid_.init_sor2(sor_parameter);
  gridMaxZ = grid_.get_max_z();
  gridStep = grid_.get_step();
  gridTimestep = grid_.get_timestep();
  gridScalStep[0] = grid_.get_scal_step(0);
  gridScalStep[1] = grid_.get_scal_step(1);

  beam1_.backstep(1, gridMaxZ, gridStep, gridTimestep, gridScalStep);
  beam2_.backstep(2, gridMaxZ, gridStep, gridTimestep, gridScalStep);

  /* Third parameter: 0 for electrons, 1 for muons*/
  secondaries_.set_pair_parameters(beam1_, beam2_, 0, switches.get_pair_ecut(), switches.get_pair_step(), grid_.get_step(), grid_.get_timestep());
  muons_.set_pair_parameters(beam1_, beam2_, 1, switches.get_muon_ecut(), switches.get_pair_step(), grid_.get_step(), grid_.get_timestep());

  switch (switches.get_time_order())
    {
    case 1:
      std::cerr << " SIMULATE time_order = 1 to be done " << std::endl;
      exit(0);
      //      break;
    case 2:
      if (switches.get_track_pairs() && switches.get_track_muons())
		{
		  beam_interaction_with_trackpair_muon(secondaries_, muons_, sor_parameter );
		}
      else if (switches.get_track_pairs()) 
		{
		  beam_interaction_with_trackpair(secondaries_, sor_parameter );
		}
      else if (switches.get_track_muons()) 
		{
		  beam_interaction_with_trackpair(muons_, sor_parameter );
		}
	  else  beam_interaction( sor_parameter );
      break;
    }
}

std::string GUINEA::output_flow() const 
{
  std::ostringstream out;
  //double esum1,esum2;
  out << title(std::string("other information "));
  //   esum1 = beam1_.meanEnergy();
  //   esum2 = beam2_.meanEnergy();
  //   out << " de1 = " <<  esum1 << " de2 = " <<  esum2 << std::endl;
  //if (switches.get_store_secondaries()) 
  if (switches.get_store_pairs()) 
    {
      out << " initial pairs are saved on file : " << secondaries0_file_ << std::endl;
    }
  //if (switches.get_track_secondaries()) 
  if (switches.get_track_pairs()) 
   {
     out << " pairs are saved on file : " << secondaries_file_ << std::endl;
   }
  if (switches.get_store_muons()) 
    {
      out << " initial muons are saved on file : " << muons0_file_ << std::endl;
    }
  //if (switches.get_track_secondaries()) 
  if (switches.get_track_muons()) 
   {
     out << " muons are saved on file : " << muons_file_ << std::endl;
   }
  if (switches.get_do_bhabhas())
    {
      int nbhabha_ini, nbhabha_photon_ini;
      grid_.get_bhabhas().numbers_of_loaded(nbhabha_ini, nbhabha_photon_ini); 
      out << " " << nbhabha_ini << " bhabhas have been loaded " << nbhabha_photon_ini << " bhabha_photons have been loaded " <<std::endl; 
      out << " read bhabha samples are saved on file : " << bhabha_prod_ << std::endl;
    }
  if (switches.get_do_jets())
    {
      out << " jets are saved on file : " << jet_file_ << std::endl;
    }
    if (switches.get_store_beam()) 
    {
      out << number_of_stored_particles1_  << " particles are saved on file : " << beam1_file_ << std::endl;
      out << number_of_stored_particles2_  << " particles are saved on file : " << beam2_file_ << std::endl;
      out << " coherent particles of beam1 are saved on file : " << coh1_file_ << std::endl;
      out << " coherent particles of beam2 are saved on file : " << coh2_file_ << std::endl;
      out << " trident particles of beam1 are saved on file : " << tri1_file_ << std::endl;
      out << " trident particles of beam2 are saved on file : " << tri2_file_ << std::endl;
    }
   return out.str();
}

void GUINEA::make_step(int i1,int i2,PHI_FLOAT *sor_parameter)
{
  int i_grid;
  int i_offset;

  float min_z = 0.5*(i2-i1-1);

  const PAIR_PARAMETER& ppar = secondaries_.get_pair_parameters();
  const PAIR_PARAMETER& mpar = muons_.get_pair_parameters();
  if(ppar.get_s4()!=mpar.get_s4()){
    std::cout << "s4 must be the same for muons and pairs" << std::endl;
  }
  grid_.all_distribute(i1, i2, switches,ppar.get_s4(), ppar.get_lns4());
				
  for (i_grid=1;i_grid<=switches.get_extra_grids();i_grid++)
    gridsPtr_[i_grid]->distribute_particles(i1, i2, switches.get_electron_distribution_rho(), switches.get_force_symmetric());
  
  time_.add_timer(1);
  grid_.step_lumi(min_z, secondaries_, time_counter_, switches, i1, i2 ) ;
  time_counter_++;

  time_.add_timer(2); 

  grid_.update_slice_charge(i1, i2);

  grid_.computeFields(switches.get_integration_method(), switches.get_charge_sign(),sor_parameter);
  for (i_grid=1;i_grid<=switches.get_extra_grids();i_grid++)
    {
      gridsPtr_[i_grid]->computeFields(switches.get_integration_method(), switches.get_charge_sign(),sor_parameter);
    }
  
  time_.add_timer(3);
  //int nbeam = 1;

  //   grid_.moveAllParticles( gridsPtr_, beam1_, nbeam, i1, switches.get_interpolation(), switches.get_do_beamstrahlung(),switches.get_ST_spin_flip(), switches.get_emin(), switches.get_do_prod(),switches.get_extra_grids(), switches.get_charge_sign(),switches.get_bmt_precession());

  grid_.moveAllParticles( gridsPtr_, i1,i2, switches.get_interpolation(), switches.get_do_beamstrahlung(),switches.get_ST_spin_flip(), switches.get_emin(), switches.get_do_prod(),switches.get_extra_grids(), switches.get_charge_sign(),switches.get_bmt_precession(),switches.get_do_trident());
  
  // nbeam = 2;
  //  std::cout << " advance the second beam.... " << std::endl;
  //  grid_.moveAllParticles( gridsPtr_, beam2_, nbeam, i2, switches.get_interpolation(), switches.get_do_beamstrahlung(),switches.get_ST_spin_flip(),switches.get_emin(), switches.get_do_prod(),switches.get_extra_grids(), switches.get_charge_sign(), switches.get_bmt_precession());
  
  time_.add_timer(4);
 
  if ( switches.get_do_photons(1) || switches.get_do_photons(2) || switches.get_load_photon())
    {
      grid_.distribute_photons(i1,i2, switches.get_photon_ratio(),generator_ );    
      grid_.photon_lumi(min_z,switches,secondaries_, muons_, generator_);

      if(switches.get_do_pairs()||switches.get_do_hadrons()||switches.get_do_compt() || switches.get_do_muons())
		{
		  // here the pairs are generated
		  grid_.photon_lumi_2(min_z,switches, secondaries_, muons_, generator_);
		}
      if (switches.get_do_compt()) 
		{
		  grid_.photon_lumi_3(min_z,switches, secondaries_, generator_);
		}
      if (switches.get_do_coherent())
		{
		  grid_.move_photons2(beam1_,1,i1, generator_);
		  grid_.move_photons2(beam2_,2,i2,  generator_);
		}
      else
		{
		  grid_.move_photons(beam1_,1,i1);
		  grid_.move_photons(beam2_,2,i2);
		}
      
    } // end if (store photons)
  time_.add_timer(5);
  i_offset=i1+i2-grid_.get_n_cell_z();
  //if(switches.get_track_secondaries())
  if(switches.get_track_pairs())
    {
      double d_eps_1, d_eps_2;

      secondaries_.get_pair_parameters().get_d_eps(d_eps_1, d_eps_2);
      if(switches.get_do_tertphot())
	{
	  if (i_offset<0)
	    {
	      grid_.move_pairs_tertphot(gridsPtr_, secondaries_, i1, d_eps_1, d_eps_2, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	  else
	    {
	      grid_.move_pairs_tertphot(gridsPtr_,secondaries_, i1-i_offset-1, d_eps_1, d_eps_2, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	}
      else
	{
	  if (i_offset<0)
	    {
	      grid_.move_pairs(gridsPtr_, secondaries_, i1, d_eps_1, d_eps_2, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	  else
	    {
	      grid_.move_pairs(gridsPtr_,secondaries_, i1-i_offset-1, d_eps_1, d_eps_2, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	}
    }
  if(switches.get_track_muons())
    {
      double d_eps_1_mu, d_eps_2_mu;

      muons_.get_pair_parameters().get_d_eps(d_eps_1_mu, d_eps_2_mu);
      if(switches.get_do_tertphot())
	{
	  if (i_offset<0)
	    {
	      grid_.move_pairs_tertphot(gridsPtr_, muons_, i1, d_eps_1_mu, d_eps_2_mu, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	  else
	    {
	      grid_.move_pairs_tertphot(gridsPtr_, muons_, i1-i_offset-1, d_eps_1_mu, d_eps_2_mu, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	}
      else
	{
	  if (i_offset<0)
	    {
	      grid_.move_pairs(gridsPtr_, muons_, i1, d_eps_1_mu, d_eps_2_mu, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }
	  else
	    {
	      grid_.move_pairs(gridsPtr_, muons_, i1-i_offset-1, d_eps_1_mu, d_eps_2_mu, switches.get_extra_grids(), switches.get_charge_sign_0(), generator_);
	    }	  
	}
    }
  // std::cout << " for the extra photons " << std::endl;
  if (switches.get_do_pairs() || 
      switches.get_do_hadrons() ||
      switches.get_do_compt() ||
      switches.get_do_muons())
    {
      grid_.clear_extra_photons();
    }
  time_.add_timer(6);
}



void GUINEA::iteration_on_overlaping_slices(int firstSliceOfBeam1, int lastSliceOfBeam2,PHI_FLOAT* sor_parameter)
{
  int i0;
  grid_.get_results().set_step_lumi_ee (0.0);
  grid_.get_results().set_step_lumi_ee_high(0.0);
  for (i0=0;i0<grid_.get_timestep();i0++)
    {
      make_time_step_on_slices(firstSliceOfBeam1, lastSliceOfBeam2, sor_parameter);
    }
   grid_.get_results().add_lumi_ee_step(grid_.get_results().get_step_lumi_ee());
   grid_.get_results().add_lumi_ee_high_step(grid_.get_results().get_step_lumi_ee_high());	
  if (switches.get_do_dump())
    {
      int istep;
      istep = firstSliceOfBeam1 + lastSliceOfBeam2;
      if (istep%switches.get_dump_step() ==0)
		{
		  dump_beams(istep,switches.get_dump_step(), switches.get_dump_particle());
		}
    } 
}

void GUINEA::iteration_on_overlaping_slices_with_trackpair(PAIR_BEAM& pair_beam_ref, int firstSliceOfBeam1, int lastSliceOfBeam2, PHI_FLOAT* sor_parameter)
{
  int i0;
  unsigned int numberToDistribute = lastSliceOfBeam2-firstSliceOfBeam1+1;
  grid_.get_results().set_step_lumi_ee (0.0);
  grid_.get_results().set_step_lumi_ee_high(0.0);
  if (switches.get_load_event()) 
    {
      //pair_beam_ref.load_events(time_counter_, switches.get_pair_ratio(), switches.get_track_secondaries(), generator_);
      pair_beam_ref.load_events(time_counter_, switches.get_pair_ratio(), switches.get_track_pairs(), generator_);
    }
  for (i0=0;i0<grid_.get_timestep();i0++)
    {
      pair_beam_ref.distribute_pairs(grid_.get_delta_z(),numberToDistribute);
      pair_beam_ref.move_unactive_pairs(grid_.get_step());
      make_time_step_on_slices(firstSliceOfBeam1, lastSliceOfBeam2, sor_parameter);
      // once the current pairs have been moved the have to be redistributed
      // and before that, desactived
      pair_beam_ref.desactive_current_pairs(numberToDistribute);
    }
  grid_.get_results().add_lumi_ee_step(grid_.get_results().get_step_lumi_ee());
  grid_.get_results().add_lumi_ee_high_step(grid_.get_results().get_step_lumi_ee_high());
  if (switches.get_do_dump())
    {
      int istep;
      istep = firstSliceOfBeam1 + lastSliceOfBeam2;
      if (istep%switches.get_dump_step() ==0)
	dump_beams(istep,switches.get_dump_step(), switches.get_dump_particle());
    } 
}

void GUINEA::iteration_on_overlaping_slices_with_trackpair_muon(PAIR_BEAM& pair_beam_ref, PAIR_BEAM& muon_beam_ref, int firstSliceOfBeam1, int lastSliceOfBeam2, PHI_FLOAT* sor_parameter)
{
  int i0;
  unsigned int numberToDistribute = lastSliceOfBeam2-firstSliceOfBeam1+1;
  grid_.get_results().set_step_lumi_ee (0.0);
  grid_.get_results().set_step_lumi_ee_high(0.0);
  if (switches.get_load_event()) 
    {
      //pair_beam_ref.load_events(time_counter_, switches.get_pair_ratio(), switches.get_track_secondaries(), generator_);
      pair_beam_ref.load_events(time_counter_, switches.get_pair_ratio(), switches.get_track_pairs(), generator_);
    }
  for (i0=0;i0<grid_.get_timestep();i0++)
    {
      pair_beam_ref.distribute_pairs(grid_.get_delta_z(),numberToDistribute);
      pair_beam_ref.move_unactive_pairs(grid_.get_step());
      muon_beam_ref.distribute_pairs(grid_.get_delta_z(),numberToDistribute);
      muon_beam_ref.move_unactive_pairs(grid_.get_step());

      make_time_step_on_slices(firstSliceOfBeam1, lastSliceOfBeam2, sor_parameter);
      // once the current pairs have been moved the have to be redistributed
      // and before that, desactived
      pair_beam_ref.desactive_current_pairs(numberToDistribute);
      muon_beam_ref.desactive_current_pairs(numberToDistribute);
    }
   grid_.get_results().add_lumi_ee_step(grid_.get_results().get_step_lumi_ee());
   grid_.get_results().add_lumi_ee_high_step(grid_.get_results().get_step_lumi_ee_high());
  if (switches.get_do_dump())
    {
      int istep;
      istep = firstSliceOfBeam1 + lastSliceOfBeam2;
      if (istep%switches.get_dump_step() ==0)
	dump_beams(istep,switches.get_dump_step(), switches.get_dump_particle());
    } 
}

void GUINEA::print_program_outputs(std::string nameOfProtokoll)
{

  FILE_IN_OUT filout;

  filout.open_file(nameOfProtokoll, "w");
  filout.set_header(header());
  filout.save_object_on_output_listing(&beam_parameters1_);
  filout.save_object_on_output_listing(&beam_parameters2_);
  filout.save_object_on_output_listing(&switches);
  filout.save_object_on_output_listing(&grid_);
  filout.save_object_on_output_listing(&beam1_);
  filout.save_object_on_output_listing(&beam2_);
#ifdef USE_EDM4HEP
  if (switches.get_do_edm4hep())
    {
      beam_parameters1_.EDM_output_1();
      beam_parameters2_.EDM_output_2();
      switches.EDM_output();
      grid_.EDM_output();
      beam1_.EDM_output();
      beam2_.EDM_output();
    }
#endif
  if (switches.get_do_coherent())
    {
      filout.save_object_on_output_listing(&grid_.get_coherent_results());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) grid_.get_coherent_results().EDM_output();
#endif
    }
  if (switches.get_do_trident())
    {
      filout.save_object_on_output_listing(&grid_.get_trident_results());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) grid_.get_trident_results().EDM_output();
#endif
    }
  if (switches.get_do_pairs() || switches.get_do_bhabhas() || switches.get_do_compt())
    {
      //     filout.save_object_on_output_file(&secondaries_);
      filout.save_object_on_output_listing(secondaries_.get_results());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) secondaries_.get_results()->EDM_output();
#endif
    }
  if (switches.get_do_muons())
    {
      filout.save_object_on_output_listing(muons_.get_results());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) muons_.get_results()->EDM_output();
#endif
    }
  if (switches.get_do_jets())
    {
      filout.save_object_on_output_listing(grid_.get_minijets());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) grid_.get_minijets()->EDM_output();
#endif
    }
  if (switches.get_do_cross())
    {
      filout.save_object_on_output_listing(cross_);
    }
  if (switches.get_do_compt())
    {
      filout.save_object_on_output_listing(&grid_.get_compton());
#ifdef USE_EDM4HEP
      if (switches.get_do_edm4hep()) grid_.get_compton().EDM_output();
#endif
    }
  //  filout.save_object_on_output_file(& muon_results);
  filout.save_object_on_output_listing(&grid_.get_results());
#ifdef USE_EDM4HEP
  if (switches.get_do_edm4hep()) grid_.get_results().EDM_output();
#endif
  filout.save_object_on_output_listing(this);
  filout.close();
}

void GUINEA::xycuts_for_grids(const ABSTRACT_PARTICLE_BEAM& bff1,const ABSTRACT_PARTICLE_BEAM& bff2 , int nbgrids, std::vector<float>& size_x, std::vector<float>& size_y, int& updated_n_cell_x, int& updated_n_cell_y ) const
{
  int counter;
  float tmp;
   //float cutx, cuty;
  if (nbgrids > 7) 
    {
      std::cerr << " GUINEA::cuts_from_sigmas : ERROR the number of grids is limited to 6 , nbgrids = " << nbgrids << std::endl;
      exit(0);
    }
  size_x.clear();
  size_y.clear();
  size_x.resize(nbgrids);
  size_y.resize(nbgrids);
  counter = 0;
  if (switches.get_automatic_grid_sizing() ) 
    {
      // automatic generation of cuts
      main_grid_automatic_xycuts(bff1,bff2, size_x[0], size_y[0], updated_n_cell_x, updated_n_cell_y);
    }
  else 
    {
      // cuts from user's data : 
      main_grid_xycuts_from_user(bff1,bff2, size_x[0], size_y[0]);
    }
  counter++;
  if (counter >= nbgrids) return;
  size_x[counter]=2.0*size_x[0];
  size_y[counter]=2.0*size_y[0];
  counter++;
  if (counter >= nbgrids) return;
  size_x[counter]=size_x[1];
  tmp=pow((double)(size_x[0]/size_y[0]),(double)0.333);
  size_y[counter]=size_y[1]*tmp;
  counter++;
  if (counter >= nbgrids) return;

  size_x[counter]=size_x[1];
  size_y[counter]=size_y[1]*tmp*tmp;
  counter++;
  if (counter >= nbgrids) return;

  size_x[counter]=size_x[1];
  size_y[counter]=size_x[1];
  counter++;
  if (counter >= nbgrids) return;

  size_x[counter]=size_x[1]*2.0;
  size_y[counter]=size_x[1]*2.0;
  counter++;
  if (counter >= nbgrids) return;

  size_x[counter]=size_x[1]*6.0;
  size_y[counter]=size_x[1]*6.0;
}

float GUINEA::zcut_for_grids(const ABSTRACT_PARTICLE_BEAM& bff1,const ABSTRACT_PARTICLE_BEAM& bff2 ) const
{
  //float  cutz, 
  float size_z;
  if (switches.get_automatic_grid_sizing() ) 
    {
      // automatic generation of cuts
      size_z = default_zcut_from_beams(bff1,bff2);
    }
  else 
    {
      // cuts from user's data : 
      size_z = main_grid_zcut_from_user(bff1,bff2);
    }
  return size_z;
}


void GUINEA::main_grid_xycuts_from_loaded_beam(const BEAM_FROM_FILE& bff1,const BEAM_FROM_FILE& bff2, float& size_x, float& size_y ) const
{
  float cutx, cuty;
      float dumx, dumy;
      float sigx, sigy, sigx2, sigy2; 
      cutx = params_.readFValue("cut_x_factor");
      cuty = params_.readFValue("cut_y_factor");
      bff1.beamXyRms(dumx, dumy, sigx, sigy);
      bff2.beamXyRms(dumx, dumy, sigx2, sigy2);
      if (sigx < sigx2) sigx = sigx2;
      if (sigy < sigy2) sigy = sigy2;
      if(cutx>0.0) size_x = cutx*sigx;
      else size_x = 3.0*sigx;
      if(cuty>0.0) size_y = cuty*sigy;
      else size_y = 3.0*sigy;  
}
float GUINEA::main_grid_zcut_from_loaded_beam(const BEAM_FROM_FILE& bff1,const BEAM_FROM_FILE& bff2) const
{
  float  cutz;
  float dumz, size_z;
      float sigz, sigz2; 
      cutz = params_.readFValue("cut_z_factor");
      bff1.beamZRms(dumz, sigz);
      bff2.beamZRms( dumz,sigz2);
      if (sigz < sigz2) sigz = sigz2;
      if(cutz>0.0) size_z = cutz*sigz;
      else size_z = default_zcut(sigz);
      return size_z;
}



void GUINEA::main_grid_automatic_xycuts(const ABSTRACT_PARTICLE_BEAM& bff1,const ABSTRACT_PARTICLE_BEAM& bff2 ,float& size_x, float& size_y, int& new_n_cell_x, int& new_n_cell_y) const
{
  float x01, y01, sigmax1, sigmay1;
  float x02, y02, sigmax2, sigmay2;
  float dx, dy, xmoy, ymoy, sigmax, sigmay;
  float theta0, gam1, gam2, gamma;
  float npart1, npart2, nbpart;
  float Dx, Dy;
  float betax, betay;
  float offset;
  float emittx1, emitty1, emittx2, emitty2;
  float z01, sigmaz1;
  float z02, sigmaz2;
  float sigmaz;

  bff1.beamXyRms(x01, y01, sigmax1, sigmay1);
  bff2.beamXyRms(x02,y02, sigmax2, sigmay2);
  bff1.emittances(emittx1, emitty1);
  bff2.emittances(emittx2, emitty2);
  gam1 = bff1.gamma();
  gam2 = bff2.gamma();

  // emittx = 0.5*(emittx1+emittx2);
  // emitty = 0.5*(emitty1+emitty2);

  npart1  = beam_parameters1_.n_particles();
  npart2  = beam_parameters2_.n_particles();

  xmoy = 0.5*( x01 + x02 );
  ymoy = 0.5*( y01 + y02 );
  dx = 0.5*fabs(x02 - x01);
  dy = 0.5*fabs(y02 - y01);
  //
  sigmax = (sigmax2 > sigmax1) ? sigmax2 : sigmax1; 
  sigmay = (sigmay2 > sigmay1) ? sigmay2 : sigmay1;

  gamma = 0.5*(gam1 + gam2);
  //  std::cout << " gamma = " << gamma << std::endl;
  nbpart = 0.5*(npart1 + npart2);


  // the sigmas are in nanometres
  theta0 = 2.*nbpart*RE/(gamma*(sigmax+sigmay)*1.0e-9);

  bff1.beamZRms(z01, sigmaz1);
  bff2.beamZRms(z02, sigmaz2);

  sigmaz = (sigmaz2 > sigmaz1) ? sigmaz2 : sigmaz1;


  Dy = theta0*sigmaz/sigmay;
  Dx = theta0*sigmaz/sigmax;

  betay = 0.5*(beam_parameters1_.beta_y()+beam_parameters2_.beta_y());
  betax = 0.5*(beam_parameters1_.beta_x()+beam_parameters2_.beta_x());

  test_size_due_to_cdm(gamma, sigmax, sigmay, sigmaz, betax, betay);


  offset = size_due_to_cdm_deflection(theta0, sigmaz, sigmax, betax, Dx,dx);
  //  std::cout << " offset in x " << offset << std::endl;
  size_x = fabs(xmoy) + dx + offset + 6.0*sigmax;
  //  std::cout << " size_x = " << size_x << std::endl;
  offset = size_due_to_cdm_deflection(theta0, sigmaz, sigmay, betay, Dy,dy);
  //  std::cout << " offset in y " << offset << std::endl;
  size_y = fabs(ymoy) + dy + offset + 10.0*sigmay;
  //  std::cout << " size_y = " << size_y << std::endl;
  new_n_cell_x = adjust_nb_cells_from_cut(size_x, sigmax);
  new_n_cell_y = adjust_nb_cells_from_cut(size_y, sigmay);
}
float GUINEA::default_zcut_from_beams(const ABSTRACT_PARTICLE_BEAM& bff1,const ABSTRACT_PARTICLE_BEAM& bff2 ) const
{
  float z01, sigmaz1;
  float z02, sigmaz2;
  float sigmaz;
  bff1.beamZRms(z01, sigmaz1);
  bff2.beamZRms(z02, sigmaz2);

  sigmaz = (sigmaz2 > sigmaz1) ? sigmaz2 : sigmaz1;
  return default_zcut(sigmaz);
}

