#pragma once

#include <vector>
#include "model.hh"
#include "data.hh"

using namespace std;

struct NODE {                             // Provides information about a node
	vector <int> arearef;                   // References the list of areas within the node
	int parent;                             // The parent node
	vector <int> child;                     // The child nodes
};

struct LEVEL {                            // Stores information about different levels 
	vector <NODE> node;                     // The nodes at a given level
 	vector <int> donelist;                  // List of nodes which have been processed
	vector <double> add;                    // Used when adding up the tree
};

class POPTREE
{
 	public:

	void init(DATA &data, int core);
	
	vector <LEVEL> lev;                     // Stores information about different levels on the tree
	int level;                              // The number of levels of scale in the model
};
