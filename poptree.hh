
#pragma once

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
