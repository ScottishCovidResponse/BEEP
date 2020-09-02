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

/// Initilaises the simulation
Simulate::Simulate(const Details &details, DATA &data, MODEL &model, const POPTREE &poptree, const Mpi &mpi, Inputs &inputs, Output &output, Obsmodel &obsmodel) : details(details), data(data), model(model), poptree(poptree), mpi(mpi), output(output), obsmodel(obsmodel)
{	
	if(details.mode != inf && mpi.ncore != 1) emsgroot("Simulation only requires one core");

	nsamp = inputs.find_int("nsamp",UNSET);                                             // Sets the number of samples for inference
	if(details.mode == multisim){
		if(nsamp == UNSET) emsgroot("The number of samples must be set");
	}
	
	model.infmax = large;
}

/// Performs a simulation
void Simulate::run()
{
	Chain chain(details,data,model,poptree,obsmodel,0);
	proportions(chain.indevp);
	
	output.simulateddata(chain.trevp,chain.indevp,details.outputdir);
	SAMPLE sample;   // TEMP
	PARAMSAMP paramsamp;		
	vector <SAMPLE> opsamp; 
  vector <PARAMSAMP> psamp;
	
	sample.meas = obsmodel.getmeas(chain.trevp,chain.indevp);
	model.setup(chain.paramval);
	sample.R0 = model.R0calc();
	sample.phi = model.phi; 
	paramsamp.paramval = chain.paramval;
	opsamp.push_back(sample);
	psamp.push_back(paramsamp);
	
	output.results(psamp,opsamp);
}

/// Runs multiple simulations
void Simulate::multirun()
{		
	unsigned int s;
	SAMPLE sample;
	PARAMSAMP paramsamp;		
	vector <SAMPLE> opsamp; 
  vector <PARAMSAMP> psamp;
	
	for(s = 0; s < nsamp; s++){
		cout << "Simulating sample " << (s+1) << endl;
		Chain chain(details,data,model,poptree,obsmodel,0);
		
		sample.meas = obsmodel.getmeas(chain.trevp,chain.indevp);
		model.setup(chain.paramval);
		sample.R0 = model.R0calc();
		sample.phi = model.phi; 
		paramsamp.paramval = chain.paramval;
		opsamp.push_back(sample);
		psamp.push_back(paramsamp);
	}
	output.results(psamp,opsamp);
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
