
#pragma once

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


struct HouseRefComparatorX
{
	public:
	HouseRefComparatorX(vector <HOUSE> &house) : house(house)
	{
	}

  bool operator()(const long &a, const long &b) const
	{
		return house[a].x < house[b].x;
	}

private:
	vector <HOUSE> &house;
};

struct HouseRefComparatorY
{
	public:
	HouseRefComparatorY(vector <HOUSE> &house) : house(house)
	{
	}

  bool operator()(const long &a, const long &b) const
	{
		return house[a].y < house[b].y;
	}

private:
	vector <HOUSE> &house;
};




class POPTREE
{
 	public:
	POPTREE() : compX(house), compY(house)
	{
	}

	void init();
	// bool compX(long lhs, long rhs) { return house[lhs].x < house[rhs].x; }
	// bool compY(long lhs, long rhs) { return house[lhs].y < house[rhs].y; }

	vector <HOUSE> house;                      // List of all the houses
	vector <LEVEL> lev;
	short level;                               // The number of levels of scale in the model
	vector <vector <long> > subpop;            // List of all individuals in node on the fine scale
	long Cfine;                                // Number of nodes on the fine scale
	vector <IND> ind;
	long areamax;                              // The maximum number of areas 
	long **nMval;                              // These are used to store the matrix M of interactions between individuals
	float ***Mval;
	long ***Mnoderef;
	long **naddnoderef;
	long ***addnoderef;

	HouseRefComparatorX compX;
	HouseRefComparatorY compY;
};
