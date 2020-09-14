#ifndef BEEPMBP__AreaTree_HH
#define BEEPMBP__AreaTree_HH

#include <vector>
#include "model.hh"
#include "data.hh"

using namespace std;

struct TreeNode {                         // Provides information about a node
	vector <unsigned int> arearef;          // References the list of areas within the node
	unsigned int parent;                    // The parent node
	vector <unsigned int> child;            // The child nodes
};

struct Level {                            // Stores information about different levels 
	vector <TreeNode> node;                 // The nodes at a given level
};

class AreaTree
{
 	public:
		AreaTree(const Data &data);
		
		vector <Level> lev;                   // Stores information about different levels on the tree
		unsigned int level;                   // The number of levels of scale in the model
		
	private:
		void geo_sort(vector <unsigned int> &vec, Dir dir);
		
	const Data &data; 
};
#endif
