/*
Load mpi: module load mpi/openmpi-x86_64

Build using: cmake -Bbuild
Compile using: cmake --build build


This lists a possible set of ways to run BEEP. 
Note, the option "-n 20" can be replace by the number of CPU nodes availabl4e.

Simulation:                 
build/bin/BEEP inputfile="examples/EX_A1.toml" mode="sim" 
	
Multiple simulations:       
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="multisim" nsimulation=100
OPTIONS: nsimulation

Prediction (uses posterior samples from inference to predict the future, with potential model modificiations):                 
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="prediction" 
OPTIONS: prediction_start, prediction_end, modification, nsim_per_sample

MAP inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="map"
OPTIONS: nparticle, ngeneration / cpu_time, nsample_final, posterior_particle

ABC-MBP inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="abcmbp" nparticle=50 ngeneration=5 nrun=4
OPTIONS: nparticle, ngeneration / cutoff_final / cpu_time, GR_max, nrun

ABC-CONT inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="abccont" nparticle=50 ngeneration=5 nrun=4
OPTIONS: nparticle, ngeneration / cutoff_final / cpu_time, GR_max, nrun

PAS inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="pas" nparticle=50 ngeneration=5 nrun=4
OPTIONS: nparticle, ngeneration / invT_final / cpu_time, GR_max, nrun

Simple ABC inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="abc" nsample=100 cutoff_frac=0.1 nrun=4
OPTIONS: nsample / cputime, cutoff / cutoff_frac, nrun

ABC-SMC inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="abcsmc" ngeneration=5 cutoff_frac=0.5 nsample=200 nrun=4
OPTIONS: nsample, ngeneration / cutoff_final / cpu_time, cutoff_frac, nrun

PMCMC inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="pmcmc" nparticle=20 nsample=200
OPTIONS: nparticle, nsample / GR_max / ESSmin, invT, nburnin, nthin, nrun

MCMC-MBP inference:
mpirun -n 4 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="mcmcmbp" invT=303 nsample=200 nrun=4
OPTIONS: nsample / GR_max, invT, nburnin, nthin, nrun

MC3 inference:
mpirun -n 20 build/bin/BEEP inputfile="examples/EX_A1.toml" mode="mc3" nchain=20 invT_final=303 nsample=200 nrun=4
OPTIONS: nchain, nsample / GR_max, invT_start, invT_final, nburnin, nquench, nthin, nrun
*/

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <signal.h>

#include "stdlib.h"
#include "time.h"

#include "inputs.hh"
#include "details.hh"
#include "data.hh"
#include "model.hh"
#include "output.hh"
#include "obsmodel.hh"

#include "utils.hh"
#include "timers.hh"

#include "simulate.hh"
#include "abc.hh"
#include "abcmbp.hh"
#include "abccont.hh"
#include "abcsmc.hh"
#include "mc3.hh"
#include "pas.hh"
#include "pmcmc.hh"
#include "ml.hh"
#include "consts.hh"

#ifdef USE_Data_PIPELINE
#include "pybind11/embed.h"
#include "datapipeline.hh"
#endif

using namespace std;

string gitversion();

unsigned int signum;

void term(int sign)
{
  if(false) cout << sign << " " << signum << "Caught!" << endl;
}

int main(int argc, char** argv)
{
	timersinit();
	timer[TIME_TOTAL].start();
	
	struct sigaction action;
  action.sa_handler = term;
  sigaction(SIGTERM, &action, NULL);	
	
	#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  #endif

	if(false){                                                 	// Outputs the git version
		cout << "BEEP version " << gitversion() << endl << endl; 
	}	

	Inputs inputs(argc,argv);                                   // Loads command line arguments and TOML file into inputs

	Details details(inputs);                                    // Loads up various details of the model
	
	Mpi mpi(details);                                           // Stores mpi information (core and ncore)
	
	bool verbose = (mpi.core == 0);                             // Parameter which ensures that only core 0 outputs results
	
	if(verbose) cout << endl;
	
#ifdef USE_Data_PIPELINE                                      // Sets up data
	pybind11::scoped_interpreter guard{};

	using namespace pybind11::literals;

  pybind11::module::import("logging").attr("basicConfig")("level"_a="DEBUG", "format"_a="%(asctime)s %(filename)s:%(lineno)s %(levelname)s - %(message)s");

 	DataPipeline *dp = new DataPipeline(
		"dpconfig.yaml", "https://github.com/ScottishCovidResponse/BEEP",
		gitversion());

	Data data(inputs,details,mpi,dp);   
#else
	Data data(inputs,details,mpi); 
#endif
	
	Model model(inputs,details,data,mpi);                       // Loads up the model
	
	auto seed = inputs.find_integer("seed",0);                  // Sets up the random seed

	switch(details.siminf){
	case SIMULATE: sran(mpi.core*10000+seed+10); break;
	default: sran(mpi.core*10000+seed+100); break;
	}
	
	ObservationModel obsmodel(details,data,model);              // Creates an observation model

	Output output(details,data,model,inputs,obsmodel,mpi);      // Creates an output class

	if(verbose) cout << endl << "Running...." << endl << endl;

	switch(details.mode){
	case SIM:                                                   // Performs a single simulation from the model 
		{		
			Simulate simu(details,data,model,inputs,output,obsmodel,mpi);
			simu.run();
		}
		break;
	
	case MULTISIM:                                              // Performs multiple simulations from the model
		{
			Simulate simu(details,data,model,inputs,output,obsmodel,mpi);
			simu.multisim();
		}
		break;
			
	case PREDICTION:                                            // Performs prediction from the model using posterior samples
		{
			Simulate simu(details,data,model,inputs,output,obsmodel,mpi);
			simu.model_modification();
		}
		break;
			
	case ABC_SIMPLE:                                            // Performs inference using a simple ABC rejection algorithm
		{	
			ABC abc(details,data,model,inputs,output,obsmodel,mpi);
			abc.run();
		}
		break;
		
	case ABC_SMC:                                               // Performs inference using the ABC-SMC algorithm
		{	
			ABCSMC abcsmc(details,data,model,inputs,output,obsmodel,mpi);
			abcsmc.run();
		}
		break;
		
	case ABC_MBP:                                               // Peforms inference using the ABC-MBP algorithm
		{	
			ABCMBP abcmbp(details,data,model,inputs,output,obsmodel,mpi);
			abcmbp.run();
		}
		break;
		
	case ABC_CONT:                                               // Peforms inference using the ABC-DA algorithm
		{	
			ABCCONT abccont(details,data,model,inputs,output,obsmodel,mpi);
			abccont.run();
		}
		break;
		
	case MC3_INF:                                               // Peforms inference using the MC3 algorithm
		{	
			MC3 mc3(details,data,model,inputs,output,obsmodel,mpi);
			mc3.run();
		}
		break;
		
	case MCMC_MBP:                                              // Peforms inference using the MCMC-MBP algorithm
		{	
			MC3 mc3(details,data,model,inputs,output,obsmodel,mpi);
			mc3.run();
		}
		break;
		
	case PAS_INF:                                              // Peforms inference using the PAS algorithm
		{	
			PAS pas(details,data,model,inputs,output,obsmodel,mpi);
			pas.run();
		}
		break;
		
	case PMCMC_INF:                                             // Peforms inference using the PMCMC algorithm
		{
			PMCMC pmcmc(details,data,model,inputs,output,obsmodel,mpi);
			pmcmc.run();
		}
		break;

	case ML_INF:                                               // Peforms inference using maximum likelihood approach
		{	
			ML ml(details,data,model,inputs,output,obsmodel,mpi);
			ml.run();
		}
		break;
		
	case DATAONLY:
		{
			vector <Particle> particle_store;
			auto param = model.sample_from_prior();   
			State state(details,data,model,obsmodel);
			state.simulate(param);
			particle_store.push_back(state.create_particle(0));
			output.generate_graphs(particle_store,1); 
		}
		break;
		
	default: emsgroot("Mode not recognised"); break;
 	}
	
	timer[TIME_TOTAL].stop();
	
	if(details.siminf == INFERENCE){
		output_timers(details.output_directory+"/Diagnostics/CPU_timings.txt",mpi);
	}
	
	auto time_av = mpi.average(timer[TIME_TOTAL].val);
	auto time_sum = mpi.sum(timer[TIME_TOTAL].val);
	auto alg_time_sum = mpi.sum(timer[TIME_ALG].val);
	auto wait_time_sum = mpi.sum(timer[TIME_WAIT].val);
	
	if(verbose){
		inputs.print_commands_not_used();
		
		if(true){
			cout <<  "Total time: " << prec(double(time_av)/(60.0*CLOCKS_PER_SEC),3) << " minutes." << endl;
		}
		else{
			cout <<  "Total time: " << prec(double(time_sum)/(60.0*CLOCKS_PER_SEC),3) << " minutes." << endl;
			cout <<  "Algorithm time: " << double(alg_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
			cout << "Mpi wait time: "  <<double(wait_time_sum)/(60.0*CLOCKS_PER_SEC) << " minutes." << endl;
		}
	}
	
#ifdef USE_Data_PIPELINE
	delete dp;
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif
}
