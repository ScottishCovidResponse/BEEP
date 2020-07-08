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

/// Generates transition data (e.g. hostipalisation data similar to the actual data we currently have from Covid-19)
void simulatedata(DATA &data, MODEL &model, POPTREE &poptree)
{
	PART *part;
	vector <SAMPLE> opsamp;  
	unsigned int ninf, i;
		
	part = new PART(data,model,poptree);
	
	part->partinit(0);
	
	if(model.setup(model.paramval) == 1) emsg("Simulate: EC9");
		
	timers.timesim -= clock();
	part->gillespie(0,data.period, 1 /* simulating */);
	timers.timesim += clock();
	
	ninf = 0; for(i = 0; i < part->indev.size(); i++) if(part->indev[i].size() > 0) ninf++;
	cout << "# Infected individuals: "<< ninf << endl;
	
	outputsimulateddata(data,model,poptree,part->fev);
	
	opsamp.push_back(outputsamp(1,0,0,data,model,poptree,model.paramval,part->fev));	
	outputresults(data,model,opsamp);
	//outputeventsample(part->fev,data,model,poptree);
}
