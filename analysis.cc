// Compile with g++ analysis.cc -o analysis -O3
// To run in simulation mode: ./analysis 1024 0

// The first number gives the number of areas into which the houses are divided (should be a power of 4)
// The second number changes the random seed

// To run in inference mode: ./analysis 1024 1000 0
// Here the second number gives the number of PMCMC samples

#include <iostream>

#include "stdlib.h"
#include "time.h"

#include "functions.hh"
#include "var.hh"        // Stores all the global variables (Sorry Ian, I know you are shaking your head)
#include "poptree.hh"
#include "model.hh"
#include "simulate.hh"
#include "PMCMC.hh"

using namespace std;

int main(int argc, char** argv)
{
	POPTREE poptree;

	switch(argc){
	case 3:   // Simulation mode
		siminf = 1;
		poptree.areamax = atoi(argv[1]);   
		srand(atoi(argv[2]));
		break;
		
	case 4:   // Inference mode
		siminf = 0;
		poptree.areamax = atoi(argv[1]);   
		nsamp = atoi(argv[2]);               
		srand(atoi(argv[3]));
		break;
		
	default:
		emsg("Wrong number of input parameters\n");
		break;
	}
	
	cout << "Initialising....\n";

	poptree.init();	
	
	MODEL model;

	model.definemodel();

	cout << "Running....\n";

	timetot = -clock();

	if(siminf == 1) simulatedata(model,poptree);
	else PMCMC(model,poptree);

	timetot += clock();
	
	cout << double(timetot)/CLOCKS_PER_SEC << " Total time\n";
	cout << double(timesim)/CLOCKS_PER_SEC << " Simulation time\n";
	cout << double(timeboot)/CLOCKS_PER_SEC << " Bootstrap time\n";
}

