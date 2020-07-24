// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include <assert.h>
#include <math.h>

#include "simulate.hh"

#include "utils.hh"
#include "timers.hh"
#include "MBPCHAIN.hh"
#include "data.hh"
#include "output.hh"
#include "obsmodel.hh"

using namespace std;

void proportions(DATA &data, MODEL &model, vector< vector <FEV> > &indev);

/// Simulates data using the MBP algorithm
void simulatedata(DATA &data, MODEL &model, POPTREE &poptree, unsigned int nsamp)
{
	unsigned int s;
	vector <SAMPLE> opsamp; 
    vector <PARAMSAMP> psamp;
    
	model.infmax = large;
		
	switch(data.mode){
	case MODE_SIM:       // Performs a single simulation 
		{
			MBPCHAIN mbpchain(data,model,poptree,1,0);
			proportions(data,model,mbpchain.indevi);
			outputsimulateddata(data,model,poptree,mbpchain.trevi,mbpchain.indevi,data.outputdir);
		}
		break;
		
	case MODE_MULTISIM:  // Performs multiple simulations and plots the distribution of results
		SAMPLE sample;
		PARAMSAMP paramsamp;
		
		for(s = 0; s < nsamp; s++){
			cout << "Simulating sample " << (s+1) << endl;
            MBPCHAIN mbpchain(data,model,poptree,1,0);
			
			sample.meas = getmeas(data,model,poptree,mbpchain.trevi,mbpchain.indevi);
			model.setup(mbpchain.paramval);
			sample.R0 = model.R0calc();
			paramsamp.paramval =  mbpchain.paramval;
			opsamp.push_back(sample);
		}
		outputresults(data,model,psamp,opsamp);
		break;
	}
}

/// Works out the proportion of individuals which visit different compartments
void proportions(DATA &data, MODEL &model, vector< vector <FEV> > &indev)
{
	unsigned int ninf, i, c, e, emax, dp;
	vector <unsigned int> visit;
	vector <unsigned int> demo;
		
	for(c = 0; c < model.comp.size(); c++) visit.push_back(0);
	
	for(dp = 0; dp < data.ndemocatpos; dp++) demo.push_back(0);
	
	ninf = 0;
	for(i = 0; i < indev.size(); i++){
		emax = indev[i].size();
		if(emax > 0){ 
			ninf++;
			demo[data.ind[i].dp]++;
			 
			c = 0;
			visit[c]++;
			for(e = 0; e < emax; e++){
				if(model.trans[indev[i][e].trans].from != model.trans[indev[i][e].trans].to){
					visit[model.trans[indev[i][e].trans].to]++;
				}
			}
		}
	}
	
	cout << endl << "# Infected individuals: "<< ninf << endl << endl;
	
	cout << "Proportion of infected individuals visiting different compartments:" << endl;
	
	for(c = 0; c < model.comp.size(); c++){
		cout << "  " << model.comp[c].name << ": " << double(100*visit[c])/ninf << "%" << endl;
	}
	cout << endl;
}
