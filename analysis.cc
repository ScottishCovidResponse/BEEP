
/*
Load mpi: module load mpi/openmpi-x86_64
Compile using:    mpicxx -O3   analysis.cc timers.cc utils.cc model.cc simulate.cc PMCMC.cc PART.cc poptree.cc -o analysis

Simulate using:        mpirun -n 1 ./analysis 65536 0

-n 4 gives the number nodes (i.e. MCMC chains)
The first number gives the number of areas into which the houses are divided (should be a power of 4)
The second number changes the random seed

Inference using:        mpirun -n 1 ./analysis 65536 10000 0

10000 is the number of MCMC iterations
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

using namespace std;
 
int main(int argc, char** argv)
{
	int ncore, core;

	#ifdef USE_MPI
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ncore);
  MPI_Comm_rank(MPI_COMM_WORLD, &core);
  #endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif
	
	POPTREE poptree;
	long nsamp;       // The number of PMCMC samples
	short siminf;     // Set to 1 for simulation and 0 for inference

	switch(argc){
	case 3:   // Simulation mode
		siminf = 1;
		poptree.areamax = atoi(argv[1]);   
		srand(atoi(argv[2]));
		
		if(ncore != 1){
			cout << "Simulation only requires one core\n";
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
		srand(atoi(argv[3]));
		break;
		
	default:
		emsg("Wrong number of input parameters");
		break;
	}
	
	if(core == 0) cout << "Initialising...." << endl;

	DATA data;
	data.readdata(siminf);
	
	poptree.init(data,core);	
	
	MODEL model;

	model.definemodel(core,data.tmax,data.popsize);

	poptree.setsus(model);
	poptree.setinf(model);
	
	if(core == 0) cout << "Running...." << endl;

	timersinit();
	timers.timetot = -clock();

	if(siminf == 1) simulatedata(data,model,poptree);
	else PMCMC(data,model,poptree,nsamp,core,ncore);

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



