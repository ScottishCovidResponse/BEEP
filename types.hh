#pragma once

#include <string>
#include <vector>

using namespace std;

struct PARAM{                              // Store information about a model parameter
 	string name;                             // Its name
 	double val;                              // The simulation value or starting value for inference
	double sim;                              // The simulation value
	double min;                              // The minimum value (assuming a uniform prior) 
	double max;                              // The maximum value (assuming a uniform prior)
	double jump;                             // The size of proposed changes in PMCMC
	long ntr, nac;                           // Store the number of proposals tried and accepted	
};

struct COMP{                               // Stores information about a compartment in the model
	string name;                             // Its name
	double infectivity;                      // How infectious that compartment is
	vector <long> trans;                     // The transitions leaving that compartment
};

struct TRANS{                              // Stores information about a compartmental model transition
	short from;                              // Which compartment the individual is coming from
	short to;                                // Which compartment the individual is going to
	short type;                              // The type of transition (exponential or gamma)
	short param1;                            // First characteristic parameter (e.g. rate)
	short param2;                            // Second characteristic parameter (e.g. standard deviation in the case of gamma) 
};

struct HOUSE {                             // Defines a house
 	double x, y;                             // Position
	vector <long> ind;                       // Individuals which belong to the house
};

struct NODE {                              // Provides information about a node
	vector <long> houseref;                  // References the list of houses within the node
	long parent;                             // The parent node
	vector <long> child;                     // The child nodes
	vector <long> fine;                      // The child nodes on the fine scale
	long popu;                               // The total population in the node
	double x, y;                             // The position of the node (given by the average of all the houses)
	short done;                              // A flag used to determine if this node has been analysed or not
};

struct LEVEL {                             // Stores information about different levels 
	vector <NODE> node;                      // The nodes at a given level
 	vector <long> donelist;                  // List of nodes which have been processed
	vector <double> add;                     // Used when adding up the tree
};

struct IND {                               // Provides information about an individual
	long noderef;                            // The node on the finescale to which the individual belongs
	long houseref;                           // The house to which the individual belongs
	short region;                            // The region to which the individual belongs
};

struct NEV {                               // Information about the immediate next events
  short type; double t;
};

struct FEV {                               // Stores information about a compartmental transition
  long trans;                              // References the transition type
	long ind;                                // The individual on which the transition happens
	double t;                                // The time of the transition
	short done;                              // Set to 1 if that transition is in the past 
};
