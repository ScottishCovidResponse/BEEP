/*
Load mpi: module load mpi/openmpi-x86_64

Compile using: make

Simulation:  ./beepmbp inputfile="examples/sim.toml"        

Inference:    mpirun -n 20 ./beepmbp inputfile="examples/inf.toml" nchain=20

ABC-SMC: ./beepmbp inputfile="examples/inf.toml" mode="abcsmc"  

ABC-MBP: ./beepmbp inputfile="examples/inf.toml" mode="abcmbp"  
 mpirun -n 2 ./beepmbp inputfile="examples/inf.toml" mode="abcmbp"  
 mpirun -n 20 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="abcmbp"  
 mpirun -n 20 ./beepmbp inputfile="examples/infMSOAnewsim.toml" mode="abcmbp"  
  ./beepmbp inputfile="examples/infMSOAnewsim.toml" mode="sim"  output_directory="Outputnew"

 mpirun -n 2 gdb --batch --quiet -ex "run" -ex "bt" -ex "quit" --args  ./beepmbp inputfile="examples/inf.toml" mode="abcmbp"  

  ./beepmbp inputfile="examples/inftest.toml" mode="sim" output_directory="OutputTest"
	
 mpirun -n 20	 ./beepmbp inputfile="examples/inftest.toml" mode="abcmbp" 
  mpirun -n 20	 ./beepmbp inputfile="examples/inftest.toml" mode="abcsmc"  output_directory="OutputTest4" 
	
	 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="sim"
	 nohup  mpirun -n 10	 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="abcmbp"  output_directory="OutputMSOATest" &
		nohup  mpirun -n 10	 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="abcsmc"  output_directory="OutputMSOATest1" &
		
 nohup  mpirun -n 20	 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="abcmbp"  output_directory="OutputMSOATest" &     beta=2
		
		
		 mpirun -n 20	 ./beepmbp inputfile="examples/infMSOAtest.toml" mode="abcmbp"  
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
#include "areatree.hh"
#include "model.hh"
#include "output.hh"
#include "obsmodel.hh"

#include "utils.hh"
#include "timers.hh"

#include "simulate.hh"
#include "abc.hh"
#include "mcmc.hh"
#include "consts.hh"

#include "combinetrace.hh"

#ifdef USE_Data_PIPELINE
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

	Mpi mpi;                                                         // Stores mpi information (core and ncore)
	
	bool verbose = (mpi.core == 0);                                  // Paramater which ensures that only core 0 outputs results
		
	if(verbose) cout << "BEEPmbp version " << gitversion() << endl << endl;  // Outputs the git version

	Inputs inputs(argc,argv);                                        // Loads command line arguments and TOML file into inputs

	Details details(inputs);                                                 // Loads up various details of the model

#ifdef USE_Data_PIPELINE                                                   // Sets up data
	pybind11::scoped_interpreter guard{};

	using namespace pybind11::literals;

  pybind11::module::import("logging").attr("basicConfig")("level"_a="DEBUG", "format"_a="%(asctime)s %(filename)s:%(lineno)s %(levelname)s - %(message)s");

	// TODO: get the path to the config file from the options
 	DataPipeline *dp = new DataPipeline(
		"dpconfig.yaml", "https://github.com/ScottishCovidResponse/BEEPmbp",
		gitversion());

	Data data(inputs,details,mpi,dp);   
#else
	Data data(inputs,details,mpi); 
#endif

	if(details.mode == COMBINE_TRACE){                                        // If in 'combinetrace' mode then do this and exit
		combine_trace(data,inputs);
		return 0;
	}
	
	Model model(inputs,details,data);                                        // Loads up the model
	
	AreaTree areatree(data);                                                   // Initialises areatree

	unsigned int seed = inputs.find_integer("seed",0);                           // Sets up the random seed
	sran(mpi.core*10000+seed);
	
	if(verbose){
		data.print_to_terminal();                                              // Summarises the data and model to the terminal
		model.print_to_terminal();
		cout << "Running...." << endl;
	}
	
	ObservationModel obsmodel(details,data,model);                           // Generates an observation model
	
	Output output(details,data,model,obsmodel);                              // Generates an output class
	
	timersinit();
	timers.timetot = -clock();
	
	switch(details.mode){
	case SIM:                                                                // Performs a single simulation from the model 
		{		
			Simulate simu(details,data,model,areatree,mpi,inputs,output,obsmodel);
			simu.run();
		}
		break;
	
	case MULTISIM:                                                           // Performs multiple simulations from the model
		{
			Simulate simu(details,data,model,areatree,mpi,inputs,output,obsmodel);
			simu.multirun();
		}
		break;
			
	case ABC_SMC:                                                             // Performs the ABC-SMC algorithm
		{	
			ABC abc(details,data,model,areatree,mpi,inputs,output,obsmodel);
			abc.smc();
		}
		break;
		
	case ABC_MBP:                                                             // Peforms the ABC-MBP algorithm
		{	
			ABC abc(details,data,model,areatree,mpi,inputs,output,obsmodel);
			abc.mbp();
		}
		break;
		
	case MCMCMC:                                                                // Performs inference on actual data
		{
			Mcmc mcmc(details,data,model,areatree,mpi,inputs,output,obsmodel);
			mcmc.run();
		}
		break;

	default: emsgroot("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(verbose) cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
	
#ifdef USE_Data_PIPELINE
	delete dp;
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif
}
