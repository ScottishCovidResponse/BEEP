#pragma once

#include <vector>
#include "model.hh"
#include "data.hh"

using namespace std;

struct NODE {                             // Provides information about a node
	vector <int> houseref;                  // References the list of houses within the node
	int parent;                             // The parent node
	vector <int> child;                     // The child nodes
	vector <int> fine;                      // The child nodes on the fine scale
	int population;                         // The total population in the node
	double sussum;                          // The sum of the relative susceptibility of the population in the node
	double x, y;                            // The position of the node (given by the average of all the houses)
	int done;                               // A flag used to determine if this node has been analysed or not
};

struct LEVEL {                            // Stores information about different levels 
	vector <NODE> node;                     // The nodes at a given level
 	vector <int> donelist;                  // List of nodes which have been processed
	vector <double> add;                    // Used when adding up the tree
};

struct IND {                              // Provides information about an individual
	int noderef;                            // The node on the finescale to which the individual beints
	int houseref;                           // The house to which the individual beints
  float sus;                              // The relative susceptibility of an individual
	float inf;                              // The relative infectivity of an individual
	vector <float> X;                       // The design matrix for fixed effects
};

class POPTREE
{
 	public:

	void init(DATA &data, int core, int areama);
	void setsus(MODEL &model);              // Sets the relative susceptibility of individuals
	void setinf(MODEL &model);              // Sets the relative infectivity of individuals

	vector <LEVEL> lev;                     // Stores information about different levels on the tree
	int level;                              // The number of levels of scale in the model
	vector <vector <int> > subpop;          // List of all individuals in node on the fine scale
	int Cfine;                              // Number of nodes on the fine scale
	vector <IND> ind;                       // The individuals in the system
	int areamax;                            // The maximum number of areas 
	int **nMval;                            // These are used to store the matrix M of interactions between individuals
	float ***Mval;
	int ***Mnoderef;
	int **naddnoderef;
	int ***addnoderef;
	
	vector <int> nMfineval;                 // These are used to store the matrix M on fine scale
	vector <vector <float> > Mfineval;
	vector <vector <int> > Mfinenoderef;
	vector< vector <vector <int> > > Mfineadd;
};
