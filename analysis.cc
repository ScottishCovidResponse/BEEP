
/*
Load mpi: module load mpi/openmpi-x86_64
Compile using:    make

Simulate using:   mpirun -n 1 ./analysis 65536 1

-n 1 sets the number of cores (set to 1 for simulation)
The first number gives the number of areas into which the houses are divided (should be a power of 4)
The second number gives the random seed

Inference using:        mpirun -n 20 ./analysis 65536 500 1

-n 20 gives the number cores
The first number gives the number of areas into which the houses are divided (should be a power of 4)
The second number gives the number of MCMC iterations
The third number gives the number of particles per core
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

#include "stdlib.h"
#include "time.h"

#include "utils.hh"
#include "timers.hh"
#include "poptree.hh"
#include "model.hh"
#include "data.hh"

#include "simulate.hh"
#include "PMCMC.hh"
#include "PMBP.hh"

using namespace std;
 
int main(int argc, char** argv)
{
	int ncore, core, npart;
	POPTREE poptree;
	long nsamp;       // The number of PMCMC samples
	short siminf;     // Set to 1 for simulation and 0 for inference

	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ncore);
  MPI_Comm_rank(MPI_COMM_WORLD, &core);
  #endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif
	
	switch(argc){
	case 3:   // Simulation mode
		siminf = 1;
		poptree.areamax = atoi(argv[1]);   
		npart = 1;
		srand(atoi(argv[2]));
		
		if(ncore != 1){
			cout << "Simulation only requires one core" <<  endl;
			#ifdef USE_MPI
			MPI_Finalize();
			#endif
			return 0;
		}
		break;
		
	case 4:   // Inference mode
		siminf = 0;
		poptree.areamax = atoi(argv[1]);   
		nsamp = atoi(argv[2]); 
		npart = atoi(argv[3]);
		srand(104);
		break;
		
	default:
		emsg("Wrong number of input parameters");
		break;
	}
	
	if(core == 0) cout << "Initialising...." << endl;
 
	MPI_Barrier(MPI_COMM_WORLD);
	DATA data;
	data.readdata(core, siminf);

	poptree.init(data,core);	
		
	MODEL model;
	model.definemodel(core,data.tmax,data.popsize);
	
	poptree.setsus(model);
	poptree.setinf(model);
	
	if(core == 0) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();

	MPI_Barrier(MPI_COMM_WORLD);
			
	if(siminf == 1) simulatedata(data,model,poptree);
	else{
		PMCMC(data,model,poptree,nsamp,core,ncore,npart);   // This is selected to perform particle MCMC inference 
		//PMBP(data,model,poptree,nsamp,core,ncore,npart);  // This is selected to perform PMBP inference
 	}
	
	timers.timetot += clock();
	
	if(core == 0){
		cout << double(timers.timetot)/CLOCKS_PER_SEC << " Total time" << endl;
		cout << double(timers.timesim)/CLOCKS_PER_SEC << " Simulation time" << endl;
		cout << double(timers.timeboot)/CLOCKS_PER_SEC << " Bootstrap time" << endl;
	}
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}
