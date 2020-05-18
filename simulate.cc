
// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include "assert.h"
#include "math.h"

#include "types.hh"
#include "functions.hh"
#include "var.hh"
#include "PART.hh"

using namespace std;

// Generates weekly case data (similar to the actual data we currently have from Covid-19)
void simulatedata(MODEL &model, POPTREE &poptree)
{
	long week, r;
	vector <long> num;
	
	part[0] = new PART(model,poptree);
  npart = 1;
	
	part[0]->partinit(0);
	
	timesim -= clock();
	assert(siminf == 1);
	part[0]->gillespie(0,tmax, 1 /* simulating */);
	timesim += clock();
		
	num = part[0]->getnumtrans("I","H",0,tmax);
	cout << "\nTotal number of hospitalised cases:\n";
	for(r = 0; r < nregion; r++) cout <<	"Region " <<  r << ": " << num[r] << "\n";
	cout << "\n";
	
	ofstream regplot("Weekly case data.txt");
	//ofstream regplot("Weekly case data.txt");
	for(week = 0; week < tmax/7; week++){
		regplot << week << " ";
		num = part[0]->getnumtrans("I","H",week*7,week*7+7);
		for(r = 0; r < nregion; r++) regplot << num[r] << " "; regplot << "\n";
	}
}
