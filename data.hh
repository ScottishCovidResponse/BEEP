#pragma once

//#include "model.hh"

struct TABLE {                             // Loads a table
	unsigned int ncol;                       // The number of columns
	unsigned int nrow;                       // The number if rows
	vector <string> heading;                 // The headings for the columns
	vector <vector <string> > ele;        	 // The elements of the table
};

struct TIMEP { 														 // Stores a time period
	string name;														 // The name of the time period
	double tend;														 // The end time
};

struct QTENSOR {                           // Stores information about a Q tensor
	string comp;                             // The compartment on which the tensor acts
	unsigned int timep;                      // The time period over which the tensor acts
	string file;  													 // The name of the file
	vector <vector <unsigned int> > to;      // Stores the mixing matrix between areas and ages at different times
	vector <vector< vector <double> > > val;	
};

struct TRANSDATA{
	string from;                             // The "from" compartment
 	string to;                               // The "to" compartment 
	string type;                             // The type (either "reg" for regional or "all" for global)
 	string file;                             // The file name of the data 
	unsigned int start;                      // The start time for the data
	unsigned int units;                      // The units used (e.g. 1=days, 7=weeks) 
	vector <vector <unsigned int> > num;     // A table giving the number of that transition type
	unsigned int rows;                       // The number of rows of data
};

struct DEMOCAT {                           // Stores demographic categories
	string name;                             // The name of the category
	vector <string> value;                   // The postential values it can take
	vector <string> param;                   // The parameters used for the susceptibility
	vector <unsigned int> col;               // The columns in the table
};

struct COVAR {                             // Stores the  covariate for the area
	string name;                             // The name of the covariate (i.e. the column in the area data file)
	string param;                            // The parameters used
	string func;                             // The functional transformation
	unsigned int col;                        // The column in the table
};

struct REGION {                            // Provides information a data region
	string name;														 // The name for the region
	string code;														 // The code for the region
};

struct AREA {                              // Provides information about an area
	string code;                             // The code for the area
	unsigned int region;                     // The data region it belongs to
 	double x, y;                             // The geographic position
	vector <unsigned int> agepop;            // The populations in different age groups
  vector <unsigned int> pop;               // The population in different demographic categories          
	vector <double> covar;                   // The covariates for that area

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
	unsigned int fediv;                      // # Divisions into which the global timeline is divided for events
	unsigned int fepertime;                  // # fediv per nsettime
	unsigned int settpertime;                // # nsettime per unit time

	unsigned int nsettime;                   // # Divisions into which the global timeline is divided for update of Q
	vector <double> settime;                 // The timings at which beta changes
	
	vector <TRANSDATA> transdata;            // Store information about transition data
	
	unsigned int period;                     // The time over which simulation/inference is performed (e.g. in weeks)
	
	string datadir; 												 // The data directory
	string regiondatafile;                   // File giving information about data regions
	string areadatafile;                     // File giving information about areas
	//string Mdatafile;                        // File giving the spatial mixing matrix
	//string Ndatafile;                        // File giving the mixing between age classes
	
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
	
	unsigned int ncovar;                     // The number of covariates for area  
	vector <COVAR> covar;                    // Covariates for area
	
	vector <TIMEP> timeperiod;               // The timings of changes to Q;
	
	vector <QTENSOR> Q;                      // Stores the list of Q tensors
	
	unsigned int nage;                       // The number of age categories
	unsigned int narage;                     // #area * #age
	unsigned int nardp;                      // #area * #ndemocatpos
	unsigned int ndemocatposperage;          // Demographic states per age group
	unsigned int nsettardp;                  // #sett * #area * #ndemocatpos

	void sortX(vector <unsigned int> &vec);	         // Used for sorting houses by x and y location
	void sortY(vector <unsigned int> &vec);	
	AreaRefComparatorX compX;
	AreaRefComparatorY compY;
	
	void readdata(unsigned int core, unsigned int ncore, unsigned int mod, unsigned int per); 
	void adddemocat(string name, vector <string> &st, vector <string> &params);
	void addcovar(string name, string param, string func);
	void addtimep(string name, double tend);
	void addQtensor(string timep, string comp, string file);
	
	private:
	string strip(string line);
	void copydata(unsigned int core);
	TABLE loadtable(string file);
	unsigned int findcol(TABLE &tab, string name);
	void normaliseQ(unsigned int q);
};
