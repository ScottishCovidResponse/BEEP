// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include <assert.h>
#include <math.h>

#include "simulate.hh"

#include "utils.hh"
#include "timers.hh"
#include "chain.hh"
#include "data.hh"
#include "output.hh"
#include "obsmodel.hh"

using namespace std;


Simulate::Simulate(DATA &data, MODEL &model, POPTREE &poptree, Mpi &mpi, Inputs &inputs, Mode mode, bool verbose) : data(data), model(model), poptree(poptree), mpi(mpi)
{	
	nsamp = inputs.find("nsamp",verbose,UNSET);                                             // Sets the number of samples for inference
	if(mode == multisim){
		if(nsamp == UNSET) emsgroot("The number of samples must be set");
	}
	
	model.infmax = large;
}

/// Simulates data using the MBP algorithm
void Simulate::run()
{
	Chain chain(data,model,poptree,0);
	proportions(chain.indevi);
	outputsimulateddata(data,model,poptree,chain.trevi,chain.indevi,data.outputdir);
}

void Simulate::multirun()
{		
	unsigned int s;
	SAMPLE sample;
	PARAMSAMP paramsamp;		
	vector <SAMPLE> opsamp; 
  vector <PARAMSAMP> psamp;
	
	for(s = 0; s < nsamp; s++){
		cout << "Simulating sample " << (s+1) << endl;
		Chain chain(data,model,poptree,0);
		
		sample.meas = getmeas(data,model,chain.trevi,chain.indevi);
		model.setup(chain.paramval);
		sample.R0 = model.R0calc();
		paramsamp.paramval = chain.paramval;
		opsamp.push_back(sample);
	}
	outputresults(data,model,psamp,opsamp);
}

/// Works out the proportion of individuals which visit different compartments
void Simulate::proportions(const vector< vector <FEV> > &indev)
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
