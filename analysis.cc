/*
Load mpi: module load mpi/openmpi-x86_64
Compile using: make

INPUTS:

Simulation:         
mpirun -n 1 ./run mode=sim model=irish simtype=smallsim seed=0 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

MBP Inference:    
mpirun -n 1 ./run mode=mbp model=irish simtype=smallsim nchain=1 nsamp=1000 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

PMCMC Inference:  
mpirun -n 20 ./run mode=pmcmc model=irish area=1024 npart=20 nsamp=1000 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt housedata=house.txt

The flag -n 1 sets the number of cores (set to 1 for simulation or more for inference)

Here is a description of the various inputs:

mode - Defines how the code operates:
		"sim" generates simulated data.
		"pmcmc" performs inference using particle MCMC.
		"mbp" performs inference using multi-temperature MBP MCMC.

model - This defines the compartmental model being used:
		"irish" defines a SEAPIRHD model for asymtopmatic / presymptomatic / symptomatic individuals.
	
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
	
 
OUTPUTS:

Simulation - This creates the specified 'transdata' and 'housedata' files along with output directory containing:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives time variation in R0.
3) "parameter.txt", which gives the parameter values used in the simulation.

Inference - The output directory contains postior information (with means and 90% credible intervals) for:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives posterior plotstime variation in R0.
3) "parameter.txt", which gives information about parameters.
4) "trace.txt", which gives trace plots for different models.
5) "traceLi.txt", which gives trace plots for the likelihoods on different chains (MBPs only).
6) "MCMCdiagnostic.txt", which gives diagnostic information on the MCMC algorthm.
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

	unsigned int ncore, core;                 // Stores the number of cores and the core of the current process
	unsigned int nsamp=UNSET;                 // The number of samples for inference
	unsigned int nchain=UNSET;                // The number of chains per core (MBP only)
	unsigned int npart=UNSET;                 // The number of particles per core (PMCMC only)
	
	unsigned int mode=UNSET;                  // Sets the mode of operation (sim/PMCMC/MBP)
  unsigned int modelsel=UNSET;              // Sets the model used
	unsigned int period=UNSET;                // The period of time for simulation / inference
	
	int seed=0;                               // Sets the random seed for simulation
	
	unsigned int j, jst, jmax, flag;
	int op;
	TRANSDATA transdata;
	string str, command, value;

	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore =num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core =num;
  #endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif

	DATA data; 
	data.areadatafile=""; data.democatfile="";
	data.outputdir="Output";                // The default output directory
		
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
			if(value == "irish"){ flag = 2; modelsel = MOD_IRISH;}
		}
		
		if(command == "simtype"){
			flag = 1;
			if(value == "smallsim"){ flag = 2; data.simtype = "smallsim";}
			if(value == "scotsim"){ flag = 2; data.simtype = "scotsim";}
			if(value == "uksim"){ flag = 2; data.simtype = "uksim";}
		}
		
		/*
		if(command == "area"){
			flag = 2;
			narea = atoi(value.c_str()); 
			if(isnan(narea)){
				stringstream ss; ss << "Value '" << value << "' is not a number";
				emsg(ss.str());
			}
		}	
		*/
		
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
			unsigned int nchaintot = atoi(value.c_str()); 
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
		
		/*
		if(command == "housedata"){
			flag = 2;
			data.housefile = value;
		}
		*/
		
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
 
	if(mode == UNSET) emsg("The property 'mode' must be set"); 
 
	MPI_Barrier(MPI_COMM_WORLD);

	if(period == UNSET) emsg("The property 'period' must be set");
	//if(data.housefile=="") emsg("The property 'housedata' must be set");
	
	data.readdata(core,ncore,mode,period);

	//if(narea==-1) emsg("The property 'narea' must be set"); 
	
	POPTREE poptree;
	poptree.init(data,core);	
		
	MODEL model;
	model.definemodel(data,core,data.period,data.popsize,modelsel);
	model.addQ(data);
	
	model.checktransdata(data);
	
	model.setsus(data);
	
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
		if(nsamp == UNSET) emsg("The number of samples must be set");
		if(npart == UNSET) emsg("The number of particles must be set");
		PMCMC(data,model,poptree,nsamp,core,ncore,npart);
		break;

	case MODE_MBP: 
		if(nsamp == UNSET) emsg("The number of samples must be set");
		if(nchain == UNSET) emsg("The number of chains must be set");
		MBP(data,model,poptree,nsamp,core,ncore,nchain);
		break;

	default: emsg("Mode not recognised"); break;
 	}
	
	timers.timetot += clock();
	
	if(core == 0){
		cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time (seconds)" << endl;
		if(mode != MODE_MBP) cout << double(timers.timesim)/CLOCKS_PER_SEC << " Simulation time (seconds)" << endl;
		if(mode == MODE_PMCMC) cout << double(timers.timeboot)/CLOCKS_PER_SEC << " Bootstrap time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timewait)/CLOCKS_PER_SEC << " MBP waiting time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembp)/CLOCKS_PER_SEC << " MBP time (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpinit)/CLOCKS_PER_SEC << " MBP init (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpQmap)/CLOCKS_PER_SEC << " MBP Qmap (seconds)" << endl;
		if(mode == MODE_MBP) cout << double(timers.timembpprop)/CLOCKS_PER_SEC << " MBP prop (seconds)" << endl;
	}
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}
