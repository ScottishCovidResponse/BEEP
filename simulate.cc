
// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include "assert.h"
#include "math.h"

#include "utils.hh"
#include "timers.hh"
#include "PART.hh"
#include "data.hh"
#include "output.hh"

using namespace std;

// Generates weekly case data (similar to the actual data we currently have from Covid-19)
void simulatedata(DATA &data, MODEL &model, POPTREE &poptree)
{
	PART *part;
	vector <SAMPLE> opsamp;  
		
	part = new PART(data,model,poptree);
	
	part->partinit(0);
	
	timers.timesim -= clock();
	part->gillespie(0,data.tmax, 1 /* simulating */);
	timers.timesim += clock();
	
	outputsimulateddata(data,model,poptree,part->fev);
	
	opsamp.push_back(outputsamp(1,0,0,data,model,poptree,part->fev));	
	outputresults(data,model,opsamp,1,0);
	outputeventsample(part->fev,data,model,poptree);
}
