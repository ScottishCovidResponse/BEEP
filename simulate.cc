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

bool Liord (Particle i,Particle j) { return (i.Li > j.Li); }

/// This is an implementation of an ABC-SMC algorithm, which is used to compare against the MBP-MCMC approach 
void Simulate::abcsmc()
{
	Chain chain(details,data,model,poptree,obsmodel,0);
		
	const int N = 100;
	const int Ngather = 500;
	const double fac = 0.5;
	const int G=20;
	vector <Generation> generation; 
	
	double sum; 
	vector <double> sumst;
	vector <double> sigma;
	
	auto nparam = model.param.size();
	
	for(auto g = 0; g < G; g++){
		if(g > 0){
			sumst.resize(N);       // Generate particle sampler
			sum = 0; for(auto i = 0; i < N; i++){ sum += generation[g-1].particle[i].w; sumst[i] = sum;}
			
			sigma.resize(nparam);
			for(auto th = 0u; th < nparam; th++){
				double av = 0, av2 = 0;
				for(auto i = 0u; i < N; i++){
					double val = generation[g-1].particle[i].paramval[th];
					av += val;
					av2 += val*val;
				}
				double var = av2/N - (av/N)*(av/N);
				if(var < tiny) var = 0;
				sigma[th] = sqrt(var);
			}
		}
			
		Generation gen;
		for(auto i = 0u; i < Ngather; i++){ 
			//cout << g << " " << i << " Gathering particles\n";		
			if(g == 0){ // For the first generation 
				chain.sample_from_prior();
			}
			else{
				unsigned int fl;
				do{
					fl = 0;
					double z = ran()*sum; int k = 0; while(k < N && z > sumst[k]) k++;
					if(k == N) emsg("Problem");
			
					PARAMSAMP paramsamp;
					paramsamp.paramval = generation[g-1].particle[k].paramval;
				
					model.setup(paramsamp.paramval);
				
					for(auto th = 0u; th < nparam; th++){
						if(sigma[th] != 0){
							paramsamp.paramval[th] += normal(0,fac*sigma[th]); 
							if(paramsamp.paramval[th] < model.param[th].min || paramsamp.paramval[th] > model.param[th].max) fl = 1;
						}
					}
				
					if(fl == 0) fl = chain.simulate(paramsamp.paramval);
				}while(fl == 1);
			}
			
			Particle part;
			part.paramval = chain.paramval;
			part.Li = obsmodel.Lobs(chain.trevp,chain.indevp);
			part.w = UNSET;
			gen.particle.push_back(part);
		}
		
		sort(gen.particle.begin(),gen.particle.end(),Liord);
		
		cout << "Generation " << g <<  " cuttoff: " << gen.particle[N].Li << endl;
		gen.particle.resize(N);
		
		if(g == 0){
			for(auto i = 0u; i < N; i++) gen.particle[i].w = 1;
		}
		else{
			for(auto i = 0u; i < N; i++){
				double sum = 0;
				for(auto j = 0u; j < N; j++){
					double logsum = 0;
					for(auto th = 0u; th < nparam; th++){
						if(sigma[th] != 0){
							logsum += normalprob(gen.particle[i].paramval[th],generation[g-1].particle[j].paramval[th],sigma[th]*sigma[th]);
						}
					}
					sum += generation[g-1].particle[j].w*exp(logsum);
				}
				gen.particle[i].w = 1.0/sum;
			}
		}
		
		output.abcsmc_output("Posterior_distributions_"+to_string(g)+".txt",gen);
		
		generation.push_back(gen);
	}
}


