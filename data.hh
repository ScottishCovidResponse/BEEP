#pragma once

struct TRANSDATA{
	string from;                             // The "from" compartment
 	string to;                               // The "to" compartment 
	string type;                             // The type (either "reg" for regional or "all" for global)
 	string file;                             // The file name of the data 
	vector <vector <int> > num;              // A table giving the number of that transition type
};

struct HOUSE {                             // Defines a house
 	double x, y;                             // Position
	vector <int> ind;                       // Individuals which beint to the house
	int region;                            // The region to which an individual beints
	float density;
};

struct HouseRefComparatorX
{
	public:
	HouseRefComparatorX(vector <HOUSE> &house) : house(house)
	{
	}

  bool operator()(const int &a, const int &b) const
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

  bool operator()(const int &a, const int &b) const
	{
		return house[a].y < house[b].y;
	}

private:
	vector <HOUSE> &house;
};

class DATA
{
	public:
	
	DATA() : compX(house), compY(house)
	{
	}

	int mode;                                // Stores if doing simulation/mbp/pmcmc
	string outputdir;                                // The output directory
	int fediv;                                       // The number of divisions into which the global timeline is divided

	vector <TRANSDATA> transdata;            // Store information about transition data
	string simtype;                          // The system on which simulation is performed
	
	int period;                            // The time over which simulation/inference is performed (e.g. in weeks)
	
	string housefile;                        // The name of the house file
	
  int nregion;                           // Number of data regions
	vector <string> regionname;              // The names of the regions
	int popsize;                             // The total population size 
  int nhouse;                              // The total number of houses
	vector <HOUSE> house;                    // List of all the houses
	
	void sortX(vector <int> &vec);	         // Used for sorting houses by x and y location
	void sortY(vector <int> &vec);	
	HouseRefComparatorX compX;
	HouseRefComparatorY compY;
	void readdata(int core, int mod, int per); 
	
	private:
	void housedensity();
};
