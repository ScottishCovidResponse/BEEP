/*
Load mpi: module load mpi/openmpi-x86_64

Compile using: make

Simulation:  ./beepmbp inputfile="examples/sim.toml"        

Inference:    mpirun -n 20 ./beepmbp inputfile="examples/inf.toml" nchain=20

ABC-SMC: ./beepmbp inputfile="examples/inf.toml" mode="abcsmc"  
*/

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>

#include "stdlib.h"
#include "time.h"

#include "inputs.hh"
#include "details.hh"
#include "data.hh"
#include "poptree.hh"
#include "model.hh"
#include "output.hh"
#include "obsmodel.hh"

#include "utils.hh"
#include "timers.hh"

#include "simulate.hh"
#include "mcmc.hh"
#include "consts.hh"

#include "combinetrace.hh"

#ifdef USE_DATA_PIPELINE
#include "pybind11/embed.h"
#include "datapipeline.hh"
#endif

using namespace std;

string gitversion();

int main(int argc, char** argv)
{
	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  #endif

	Mpi mpi;                                                                 // Stores mpi information (core and ncore)
	
	bool verbose = (mpi.core == 0);                                          // Paramater which ensures that only core 0 outputs results
		
	if(verbose) cout << "BEEPmbp version " << gitversion() << endl << endl;  // Outputs the git version

	Inputs inputs(argc,argv,verbose);                                        // Loads command line arguments and TOML file into inputs
		
	Details details(inputs);                                                 // Loads up various details of the model
	
#ifdef USE_DATA_PIPELINE                                                   // Sets up data
	pybind11::scoped_interpreter guard{};

	using namespace pybind11::literals;

  pybind11::module::import("logging").attr("basicConfig")("level"_a="DEBUG", "format"_a="%(asctime)s %(filename)s:%(lineno)s %(levelname)s - %(message)s");

	// TODO: get the path to the config file from the options
 	DataPipeline *dp = new DataPipeline(
		"dpconfig.yaml", "https://github.com/ScottishCovidResponse/BEEPmbp",
		gitversion());

	DATA data(inputs,details,mpi,dp);   
#else
	DATA data(inputs,details,mpi); 
#endif

	if(details.mode == combinetrace){                                        // If in 'combinetrace' mode then do this and exit
		combine_trace(data,inputs);
		return 0;
	}
	
	MODEL model(inputs,details,data);                                        // Loads up the model
	
	POPTREE poptree(data);                                                   // Initialises poptree

	unsigned int seed = inputs.find_int("seed",0);                           // Sets up the random seed
	sran(mpi.core*10000+seed);
	
	if(verbose){
		data.print_to_terminal();                                              // Summarises the data and model to the terminal
		model.print_to_terminal();
		cout << "Running...." << endl;
	}
	
	Obsmodel obsmodel(details,data,model);                                   // Generates an observation model
	
	Output output(details,data,model,obsmodel);                              // Generates an output class
	
	timersinit();
	timers.timetot = -clock();
	
	switch(details.mode){
	case sim:                                                                // Performs a single simulation from the model 
		{		
			Simulate simu(details,data,model,poptree,mpi,inputs,output,obsmodel);
			simu.run();
		}
		break;
	
	case multisim:                                                           // Performs multiple simulations from the model
		{
			Simulate simu(details,data,model,poptree,mpi,inputs,output,obsmodel);
			simu.multirun();
		}
		break;
			
	case abcsmc:
		{	
			Simulate simu(details,data,model,poptree,mpi,inputs,output,obsmodel);
			simu.abcsmc();
		}
		break;
		
	case inf:                                                                // Performs inference on actual data
		{
			Mcmc mcmc(details,data,model,poptree,mpi,inputs,output,obsmodel);
			mcmc.run();
		}
		break;

	default: emsgroot("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(verbose) cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
	
#ifdef USE_DATA_PIPELINE
	delete dp;
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif
}
