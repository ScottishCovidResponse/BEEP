#pragma once

//#include "model.hh"

struct TRANSDATA{
	string from;                             // The "from" compartment
 	string to;                               // The "to" compartment 
	string type;                             // The type (either "reg" for regional or "all" for global)
 	string file;                             // The file name of the data 
	vector <vector <unsigned int> > num;     // A table giving the number of that transition type
};

struct DEMOCAT {                           // Stores demographic categories
	string name;                             // The name of the category
	vector <string> value;                   // The postential values it can take
};

struct REGION {                            // Provides information a data region
	string name;														 // The name for the region
	string code;														 // The code for the region
};

struct AREA {                              // Provides information about an area
	string code;                             // The code for the area
	unsigned int region;                     // The data region it belongs to
 	double x, y;                             // The geographic position
	double density;                          // The density of population 
	vector <unsigned int> agepop;            // The populations in different age groups
  vector <unsigned int> pop;               // The population in different demographic categories          
	vector <vector <unsigned int> > ind;     // The individuals in different demographic categories
};

struct IND {                               // Provides information about an individual
	unsigned int area;                       // The area
	unsigned int dp;                         // The demographic category
};

struct AreaRefComparatorX
{
	public:
	AreaRefComparatorX(vector <AREA> &area) : area(area)
	{
	}

  bool operator()(const unsigned int &a, const unsigned int &b) const
	{
		return area[a].x < area[b].x;
	}

private:
	vector <AREA> &area;
};

struct AreaRefComparatorY
{
	public:
	AreaRefComparatorY(vector <AREA> &area) : area(area)
	{
	}

  bool operator()(const unsigned int &a, const unsigned int &b) const
	{
		return area[a].y < area[b].y;
	}

private:
	vector <AREA> &area;
};

class DATA
{
	public:
	
	DATA() : compX(area), compY(area)
	{
	}

	unsigned int mode;                       // Stores if doing simulation/mbp/pmcmc
	string outputdir;                        // The output directory
	unsigned int fediv;                      // The number of divisions into which the global timeline is divided

	vector <TRANSDATA> transdata;            // Store information about transition data
	string simtype;                          // The system on which simulation is performed
	
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
	string democatfile;                      // File giving demographic classifications
	string regiondatafile;                   // File giving information about data regions
	string areadatafile;                     // File giving information about areas
	string Mdatafile;                        // File giving the spatial mixing matrix
	string Ndatafile;                        // File giving the mixing between age classes
	
 	unsigned int popsize;                    // The total population size 
 
	unsigned int nregion;                    // Number of data regions
	vector <REGION> region;                  // The names of the data regions

	unsigned int narea;                      // The number of areas
	vector <AREA> area;                      // List of all the areas
	
	vector <IND> ind;                        // The individuals in the system
		
	unsigned int ndemocat;                   // The number of demographic categories
	vector <DEMOCAT> democat;                // Stores the demographic categories
	
	unsigned int ndemocatpos;                // The number of demographic possibilities
	vector < vector<unsigned int> > democatpos; // Stores all the posible combinations of demographic categories
	
	vector <unsigned int> nM;                // Stores the mixing matrix between areas
	vector< vector <unsigned int> > Mto;
	vector< vector <double> > Mval;

	vector< vector <double> > N;             // The maxtrix giving mixing between age groups 

	unsigned int ntimeperiod;                // The number of different time periods (2: before and after lockdown)
	vector <double> timeperiod;              // The timings of changes to Q;
	
	unsigned int Qnum;                       // The number of Q matrices (for different compartments and time variation)
	vector <string> Qcomp;                   // The compartment when matrix acts
 	vector <unsigned int> Qtimeperiod;       // The time period when matrix acts
	vector <vector <unsigned int> > nQ;      // Stores the mixing matrix between areas and ages at different times
	vector <vector< vector <unsigned int> > > Qto;
	vector <vector <vector< vector <double> > > > Qval;
		
	unsigned int nage;                       // The number of age categories
	unsigned int narage;                     // #area * #age
	unsigned int nardp;                      // #area * #ndemocatpos
	unsigned int ndemocatposperage;          // Demographic states per age group

	void sortX(vector <unsigned int> &vec);	         // Used for sorting houses by x and y location
	void sortY(vector <unsigned int> &vec);	
	AreaRefComparatorX compX;
	AreaRefComparatorY compY;
	
	void readdata(unsigned int core, unsigned int ncore, unsigned int mod, unsigned int per); 
	
	private:
	string strip(string line);
	void copydata(unsigned int core);
};
