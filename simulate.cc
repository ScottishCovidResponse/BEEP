
// Implements a modified Gillespie algorithm to simulate from the model

#include <iostream>
#include <fstream>
#include <algorithm>

#include "assert.h"
#include "math.h"

#include "utils.hh"
#include "timers.hh"
#include "PART.hh"

using namespace std;

// Generates weekly case data (similar to the actual data we currently have from Covid-19)
void simulatedata(MODEL &model, POPTREE &poptree)
{
	long week, r;
	vector <long> num;
	PART *part;
	
	part = new PART(model,poptree);
	
	part->partinit(0);
	
	timers.timesim -= clock();
	part->gillespie(0,tmax, 1 /* simulating */);
	timers.timesim += clock();
		
	num = part->getnumtrans("I","H",0,tmax);
	cout << endl << "Total number of hospitalised cases:" << endl;
	for(r = 0; r < nregion; r++) cout <<	"Region " <<  r << ": " << num[r] << endl;
	cout << endl;
	
	ofstream regplot("Weekly case data.txt");
	//ofstream regplot("Weekly case data.txt");
	for(week = 0; week < tmax/7; week++){
		regplot << week << " ";
		num = part->getnumtrans("I","H",week*7,week*7+7);
		for(r = 0; r < nregion; r++) regplot << num[r] << " "; regplot << endl;
	}
}
