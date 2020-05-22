#pragma once

struct HOUSE {                             // Defines a house
 	double x, y;                             // Position
	vector <long> ind;                       // Individuals which belong to the house
	short region;                            // The region to which an individual belongs
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

class DATA
{
	public:
	
	DATA() : compX(house), compY(house)
	{
	}

	short tmax;                              // The time in days over which inference is performed
	short nweek;                             // The number of weeks
  short nregion;                           // Number of data regions
	vector <string> regionname;              // The names of the regions
	long popsize;                            // The total population size 
  long nhouse;                             // The total number of houses
	vector <HOUSE> house;                    // List of all the houses
	vector <vector <long> > ncase;           // Hospitalised case data
	
	void sortX(vector <long> &vec);	         // Used for sorting houses by x and y location
	void sortY(vector <long> &vec);	
	HouseRefComparatorX compX;
	HouseRefComparatorY compY;
	
	void readdata(short siminf);             // Reads in the data
};
