/*
Load mpi: module load mpi/openmpi-x86_64
Compile using: make

Simulate:         
mpirun -n 1 ./run mode=sim model=irish simtype=smallsim area=16384 seed=0 period=16 transdata=I,H,reg,trans.txt housedata=house.txt

MBP Inference:    
mpirun -n 20 ./run mode=mbp model=irish area=16384 nchain=20 nsamp=1000 period=16 transdata=I,H,reg,trans.txt housedata=house.txt

PMCMC Inference:  
mpirun -n 20 ./run mode=pmcmc model=irish area=16384 npart=20 nsamp=1000

The flag -n 1 sets the number of cores (set to 1 for simulation or more for inference)

Here is a description of the various inputs:

mode - Defines how the code operates:
		"sim" generates simulated data.
		"pmcmc" performs inference using particle MCMC.
		"mbp" performs inference using multi-temperature MBP MCMC.

model - This defines the compartmental model being used:
		"original" defines the original SEAIHRD model used.
		"irish" defines a SEIR model for asymtopmatic individuals and SEIRHD for symptomatic.
	
area - Determines	the number of areas into which the houses are divided (should be a power of 4)

npart - The total number of particles used when performing PMCMC (should be a multiple of the number of cores).

nchain - The total number of chains used when performing multi-temperature MBP MCMC (should be a multiple of the number of cores).

nsamp - The number of samples used for inference (note, burnin is assumed to be a quater this value).

seed - Sets the random seed when performing inference (this is set to zero by default)

outputdir - Gives the name of the output directory (optional).

transdata - Transition data. Gives the "from" then "to" compartments, the type ("reg" means regional data and "all" means global) and then the file name. More than one set of transition data can be added for an analysis.

housedata - House data. Gives the file name for a file giving the positions of houses.

period - The period of time over which simulation/inference is performed (e.g. measured in weeks).

simtype - Determines the system on which simulation is performed (simulation mode only):
		"smallsim" a small system with 1024 houses, 10000 individuals and a 2x2 grid of data regions (used for testing)
		"scotsim" a Scotland-like system with 1.5 million houses, 5.5 million individuals and a 4x4 grid of data regions (used for testing)
		"uksim" a UK-like system with 20 million houses, 68 million individuals and a 10x10 grid of data regions (used for testing)
	
 
OUTPUT

simulation


inference
In simulation mode the code generates the files "cases_system.txt" which gives case numbers (where "system" is taken from the options above) and "houses_system.txt" which gives information about the houses. The simulation parameters are currently defined in definemodel() in model.cc.

In inference mode the code reads in "cases_system.txt" and "houses_system.txt" and outputs posterior estimates.
*/

/*
 Commands for running on DiRAC:
 ./submit-csd3 --groups 65536 --seed 1 --samples 10 --nprocs 32 --nnode 1  --walltime 1:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/test
 dos2unix submit-csd3
 squeue -u dc-pool1
 scancel <jobid>
 mybalance
*/

#include <iostream>
#include <sstream>
#include <math.h>

#include "stdlib.h"
#include "time.h"

#include "utils.hh"
#include "timers.hh"
#include "poptree.hh"
#include "model.hh"
#include "data.hh"

#include "simulate.hh"
#include "PMCMC.hh"
//#include "gitversion.hh"
#include "MBP.hh"
#include "consts.hh"

using namespace std;
 
int main(int argc, char** argv)
{
	//cout << "CoronaPMCMC version " << GIT_VERSION << endl;

	int ncore, core;                          // Stores the number of cores and the core of the current process
	int nsamp=-1;                             // The number of samples for inference
	int mode=-1;                              // Sets the mode of operation (sim/PMCMC/MBP)
  int modelsel=-1;                          // Sets the model used
	int npart=-1;                             // The number of particles per core (PMCMC only)
	int nchain=-1;                            // The number of chains per core (MBP only)
	int narea=-1;                             // The number of areas into which houses are seperated
	int period=-1;                            // The period of time for simulation / inference
	int seed=0;                               // Sets the random seed for simulation
	
	int op, j, jst, jmax, flag;
	TRANSDATA transdata;
	string str, command, value;

	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&ncore);
  MPI_Comm_rank(MPI_COMM_WORLD,&core);
  #endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif

	DATA data; data.housefile=""; data.outputdir="Output";                // The default output directory
		
	for(op = 1; op < argc; op++){                                           // Goes the various input options
		str = string(argv[op]);
		j = 0; jmax = str.length(); while(j < jmax && str.substr(j,1) != "=") j++;
		if(j == jmax){
			stringstream ss; ss << "Cannot understand " << str; 
			emsg(ss.str());
		}
		
		command = str.substr(0,j);
		value = str.substr(j+1,jmax-(j+1));
		
		flag = 0;
		
		if(command == "mode"){
			flag = 1;
			if(value == "sim"){ flag = 2; mode = MODE_SIM;}
			if(value == "pmcmc"){ flag = 2; mode = MODE_PMCMC;}
			if(value == "mbp"){ flag = 2; mode = MODE_MBP;}
		}
		
		if(command == "model"){
			flag = 1;
			if(value == "original"){ flag = 2; modelsel = MOD_OLD;}
			if(value == "irish"){ flag = 2; modelsel = MOD_IRISH;}
		}
		
		if(command == "simtype"){
			flag = 1;
			if(value == "smallsim"){ flag = 2; data.simtype = "smallsim";}
			if(value == "scotsim"){ flag = 2; data.simtype = "scotsim";}
			if(value == "uksim"){ flag = 2; data.simtype = "uksim";}
		}
		
		if(command == "area"){
			flag = 2;
			narea = atoi(value.c_str()); 
			if(isnan(narea)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
		}	
		
		if(command == "npart"){
			flag = 2;
			int nparttot = atoi(value.c_str()); 
			if(isnan(nparttot)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
			
			if(nparttot%ncore != 0) emsg("The number of particles must be a multiple of the number of cores");
			npart = nparttot/ncore;
		}	
		
		if(command == "nchain"){
			flag = 2;
			int nchaintot = atoi(value.c_str()); 
			if(isnan(nchaintot)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
			if(nchaintot%ncore != 0) emsg("The number of chains must be a multiple of the number of cores");
			nchain = nchaintot/ncore;
		}	
		
		if(command == "nsamp"){
			flag = 2;
			nsamp = atoi(value.c_str()); 
			if(isnan(nsamp)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
		}	
		
		if(command == "period"){
			flag = 2;
			period = atoi(value.c_str()); 
			if(isnan(period)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
		}	
		
		if(command == "seed"){
			flag = 2;
			seed = atoi(value.c_str()); 
			if(isnan(seed)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
		}	
		
		if(command == "transdata"){
			flag = 2;
			j = 0; jmax = value.length();
			while(j < jmax && value.substr(j,1) != ",") j++;
			if(j == jmax) emsg("Problem with transition data");
			transdata.from = value.substr(0,j);
			j++;
			
			jst = j; while(j < jmax && value.substr(j,1) != ",") j++;
			if(j == jmax) emsg("Problem with transition data");
			transdata.to = value.substr(jst,j-jst);
			j++;
			
			jst = j; while(j < jmax && value.substr(j,1) != ",") j++;
			if(j == jmax) emsg("Problem with transition data");
			transdata.type = value.substr(jst,j-jst);
			if(transdata.type != "reg" && transdata.type != "all") emsg("Transition data type not recognised"); 
			j++;
			
			transdata.file = value.substr(j,jmax-j);
			data.transdata.push_back(transdata);
		}
		
		if(command == "housedata"){
			flag = 2;
			data.housefile = value;
		}
		
		if(command == "outputdir"){
			flag = 2;
			data.outputdir = value;
		}			
		
		if(flag == 0){
			stringstream ss; ss << "Cannot understand the command '" << command << "'";			
			emsg(ss.str());
		}

		if(flag == 1){		
			stringstream ss; ss << "Cannot understand the argument '" << value << "' for '" << command << "'";			
			emsg(ss.str());
		}
	}
	
	if(mode == MODE_SIM && ncore != 1) emsg("Simulation only requires one core");
	
	if(core == 0) cout << "Initialising...." << endl;
 
	if(mode == -1) emsg("The property 'mode' must be set"); 
 
	MPI_Barrier(MPI_COMM_WORLD);

	if(period==-1) emsg("The property 'period' must be set");
	if(data.housefile=="") emsg("The property 'housedata' must be set");
	
	data.readdata(core,mode,period);

	if(narea==-1) emsg("The property 'narea' must be set"); 
	
	POPTREE poptree;
	poptree.init(data,core,narea);	
		
	MODEL model;
	model.definemodel(core,data.period,data.popsize,modelsel);
	
	poptree.setsus(model);
	poptree.setinf(model);
	
	if(core == 0) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();

	MPI_Barrier(MPI_COMM_WORLD);
	switch(mode){
	case MODE_SIM:
		srand(seed); 
		simulatedata(data,model,poptree);
		break;
		
	case MODE_PMCMC:
		if(nsamp == -1) emsg("The number of samples must be set");
		if(npart == -1) emsg("The number of particles must be set");
		PMCMC(data,model,poptree,nsamp,core,ncore,npart);
		break;

	case MODE_MBP: 
		if(nsamp == -1) emsg("The number of samples must be set");
		if(nchain == -1) emsg("The number of chains must be set");
		MBP(data,model,poptree,nsamp,core,ncore,nchain);
		break;

	default: emsg("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(core == 0){
		cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time" << endl;
		if(mode != MODE_MBP) cout << double(timers.timesim)/CLOCKS_PER_SEC << " Simulation time" << endl;
		if(mode == MODE_PMCMC) cout << double(timers.timeboot)/CLOCKS_PER_SEC << " Bootstrap time" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time" << endl;
	}
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}
