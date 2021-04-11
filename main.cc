/*
Load mpi: module load mpi/openmpi-x86_64

Compile using: make

Simulation:  ./beepmbp inputfile="examples/SEIR.toml" mode="sim" 

Inference:    mpirun -n 20 ./beepmbp inputfile="examples/inf.toml" nchain=20

ABC-SMC: ./beepmbp inputfile="examples/inf.toml" 

ABC-MBP: ./beepmbp inputfile="examples/inf.toml" mode="abcmbp" 

		 
// SUCCESSFULLY RUN
  ./beepmbp inputfile="examples/SEIR.toml" mode="sim"  
	mpirun -n 20  ./beepmbp inputfile="examples/SEIR.toml" mode="multisim"  
	mpirun -n 20  ./beepmbp inputfile="examples/SEIR.toml" mode="ppc" nsample=200 ppc_start="2020-4-01"	
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="counter"  	
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="abc" nsample=200 cutoff_frac=0.01
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="abcsmc" ngeneration=5 cutoff_frac=0.5 nsample=200	 
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="pmcmc" nparticle=100 nsample=200
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="abcmbp" nparticle=200 ngeneration=10 
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="pais" nparticle=200 ngeneration=50
	
	./beepmbp inputfile="examples/SEIR.toml" mode="mc3" nchain=80 invT_final=303 nsample=200
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="mc3" nchain=80 invT_final=303 nsample=200
	
	./beepmbp inputfile="examples/SEIR.toml" mode="abcmbp_gr" nparticle=200 ngeneration=50 ngroup=4
	mpirun -n 20 ./beepmbp inputfile="examples/SEIR.toml" mode="abcmbp_gr" nparticle=200  ngeneration=50 ngroup=4
	
	
	mpirun -n 2 gdb --batch --quiet -ex "run" -ex "bt" -ex "quit" --args ./beepmbp inputfile="examples/SEIR.toml" mode="abc" nsample=200 cutoff_frac=0.01
	mpirun -n 2  gdb --batch --quiet -ex "run" -ex "bt" -ex "quit" --args  ./beepmbp inputfile="examples/SEIR.toml" mode="pmcmc"  	
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
#include "areatree.hh"
#include "model.hh"
#include "output.hh"
#include "obsmodel.hh"

#include "utils.hh"
#include "timers.hh"

#include "simulate.hh"
#include "abc.hh"
#include "abcmbp.hh"
#include "abcmbp_gr.hh"
#include "abcsmc.hh"
#include "mc3.hh"
#include "pais.hh"
#include "pmcmc.hh"
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
  if(false) cout << sign << " " << signum << "Caught!\n";
}

int main(int argc, char** argv)
{
	struct sigaction action;
  action.sa_handler = term;
  sigaction(SIGTERM, &action, NULL);	
	
	#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  #endif

	Mpi mpi;                                                    // Stores mpi information (core and ncore)
	
	bool verbose = (mpi.core == 0);                             // Parameter which ensures that only core 0 outputs results
		
	if(false && verbose){                                      	// Outputs the git version
		cout << "BEEPmbp version " << gitversion() << endl << endl; 
	}	

	Inputs inputs(argc,argv);                                   // Loads command line arguments and TOML file into inputs

	Details details(inputs);                                    // Loads up various details of the model

#ifdef USE_Data_PIPELINE                                      // Sets up data
	pybind11::scoped_interpreter guard{};

	using namespace pybind11::literals;

  pybind11::module::import("logging").attr("basicConfig")("level"_a="DEBUG", "format"_a="%(asctime)s %(filename)s:%(lineno)s %(levelname)s - %(message)s");

 	DataPipeline *dp = new DataPipeline(
		"dpconfig.yaml", "https://github.com/ScottishCovidResponse/BEEPmbp",
		gitversion());

	Data data(inputs,details,mpi,dp);   
#else
	Data data(inputs,details,mpi); 
#endif
	
	Model model(inputs,details,data);                           // Loads up the model
	
	AreaTree areatree(data);                                    // Initialises areatree (used later for sampling infections)

	auto seed = inputs.find_integer("seed",0);                  // Sets up the random seed

	switch(details.siminf){
	case SIMULATE: sran(mpi.core*10000+seed+10); break;
	default: sran(mpi.core*10000+seed+100); break;
	}
	
	ObservationModel obsmodel(details,data,model);              // Creates an observation model

	Output output(details,data,model,inputs,obsmodel,mpi);      // Creates an output class

	if(verbose) cout << endl << "Running...." << endl << endl;
	
	timersinit();
	timer[TIME_TOTAL].start();
	
	switch(details.mode){
	case SIM:                                                   // Performs a single simulation from the model 
		{		
			Simulate simu(details,data,model,areatree,inputs,output,obsmodel,mpi);
			simu.run();
		}
		break;
	
	case MULTISIM:                                              // Performs multiple simulations from the model
		{
			Simulate simu(details,data,model,areatree,inputs,output,obsmodel,mpi);
			simu.multirun();
		}
		break;
			
	case COUNTER: case PPC:                                     // Performs counterfactual simulations from the model
		{
			Simulate simu(details,data,model,areatree,inputs,output,obsmodel,mpi);
			simu.counter();
		}
		break;
			
	case ABC_SIMPLE:                                            // Performs a simple ABC algorithm
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
			ABCMBP abcmbp(details,data,model,areatree,inputs,output,obsmodel,mpi);
			abcmbp.run();
		}
		break;
		
	case ABC_MBP_GR:                                            // Peforms inference using the ABC-MBP algorithm with Gelman-Rubin
		{	
			ABCMBP_GR abcmbp_gr(details,data,model,areatree,inputs,output,obsmodel,mpi);
			abcmbp_gr.run();
		}
		break;
		
	case MC3_INF:                                               // Peforms the MC3 algorithm
		{	
			MC3 mc3(details,data,model,areatree,inputs,output,obsmodel,mpi);
			mc3.run();
		}
		break;
		
	case PAIS_INF:                                              // Peforms inference using the PAIS algorithm
		{	
			PAIS pais(details,data,model,areatree,inputs,output,obsmodel,mpi);
			pais.run();
		}
		break;
		
	case PMCMC_INF:                                             // Peforms inference using the PMCMC algorithm
		{
			PMCMC pmcmc(details,data,model,areatree,inputs,output,obsmodel,mpi);
			pmcmc.run();
		}
		break;

	default: emsgroot("Mode not recognised"); break;
 	}
	
	timer[TIME_TOTAL].stop();
	
	output_timers(details.output_directory+"/Diagnostics/CPU_timings.txt",mpi);
	
	auto time_av = mpi.average(timer[TIME_TOTAL].val);
	if(verbose){
		cout <<  "Total time: " << prec(double(time_av)/(60.0*CLOCKS_PER_SEC),3) << " minutes." << endl;
	}
	
#ifdef USE_Data_PIPELINE
	delete dp;
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif
}
