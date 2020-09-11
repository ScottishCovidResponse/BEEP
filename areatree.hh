#ifndef BEEPMBP__AreaTree_HH
#define BEEPMBP__AreaTree_HH

#include <vector>
#include "model.hh"
#include "data.hh"

using namespace std;

struct TreeNode {                             // Provides information about a node
	vector <unsigned int> arearef;          // References the list of areas within the node
	unsigned int parent;                    // The parent node
	vector <unsigned int> child;            // The child nodes
};

struct Level {                            // Stores information about different levels 
	vector <TreeNode> node;                     // The nodes at a given level
 	vector <unsigned int> donelist;         // List of nodes which have been processed
	vector <double> add;                    // Used when adding up the tree
};

class AreaTree
{
 	public:

	AreaTree(Data &data);
	
	vector <Level> lev;                     // Stores information about different levels on the tree
	unsigned int level;                     // The number of levels of scale in the model
};
#endif
