#ifndef BEEPMBP__AreaTree_HH
#define BEEPMBP__AreaTree_HH

using namespace std;

#include "struct.hh"
#include "data.hh"

struct TreeNode {                         // Provides information about a node
	vector <unsigned int> arearef;          // References the list of areas within the node
	unsigned int parent;                    // The parent node
	vector <unsigned int> child;            // The child nodes
	
	vector <Simu_or_mbp> mbp_sim;           // Used in fixedtree to determine which areas are simulated 
};

struct Level {                            // Stores information about different levels 
	vector <TreeNode> node;                 // The nodes at a given level
};

struct AreaSort{                          // Used for time ordering event references	
	unsigned int area;                      // The area    
	double num;                             // The number used to sort (e.g. eastings)
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
