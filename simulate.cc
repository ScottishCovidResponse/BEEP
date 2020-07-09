// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include "assert.h"
#include "math.h"

#include "utils.hh"
#include "timers.hh"
#include "MBPCHAIN.hh"
#include "data.hh"
#include "output.hh"

using namespace std;

void proportions(DATA &data, MODEL &model, vector< vector <FEV> > &indev);

/// Simulates data using the MBP algorithm
void simulatedata(DATA &data, MODEL &model, POPTREE &poptree)
{
	
	MBPCHAIN *mbpchain;
	vector <SAMPLE> opsamp; 
		
	model.infmax = large;
		
	mbpchain = new MBPCHAIN(data,model,poptree);
	mbpchain->init(data,model,poptree,1,0);
		
	proportions(data,model,mbpchain->indevi);
	
	outputsimulateddata(data,model,poptree,mbpchain->trevi,mbpchain->indevi);
	
	opsamp.push_back(outputsamp(0,0,0,0,data,model,poptree,mbpchain->paramval,0,mbpchain->trevi,mbpchain->indevi));
	outputresults(data,model,opsamp);
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
	
	cout << "# Infected individuals: "<< ninf << endl << endl;
	
	cout << "Proportion of individuals visiting different compartments:" << endl;
	
	for(c = 0; c < model.comp.size(); c++){
		cout << model.comp[c].name << ": " << double(100*visit[c])/ninf << "%" << endl;
	}
	cout << endl;
}
