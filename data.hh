#ifndef BEEPMBP__DATA_HH
#define BEEPMBP__DATA_HH

#include <string>

using namespace std;

#include "consts.hh"
#include "details.hh"
#include "utils.hh"

struct QTENSOR {                           // Stores information about a Q tensor
	string comp;                             // The compartment on which the tensor acts
	unsigned int timep;                      // The time period over which the tensor acts
	string name;     			 									 // The name of the file
	unsigned int Qtenref;                    // References the actual tensor information in genQ
};

struct SPARSETENSOR{                       // Stores the Q tensor in a sparse way
	string name;                             // The reference name
	vector <unsigned short> ntof;                 	  // Stores the mixing matrix between areas and ages at different times
	unsigned short **tof;
	float ***valf;
};

struct GENQ{
	string Nall;                             // The age matrix of all interations
	string Nhome;                            // The age matrix of home interations
	string Nother;                           // The age matrix of other interations
	string Nschool;                          // The age matrix of school interations
	string Nwork;                            // The age matrix of work interations
	string M;                                // The geographic mixing matrix
	string localhome;                        // The Q matrix for someone at home
	string flowall;                          // The Q matrix for general daily life

	vector <SPARSETENSOR> Qten;              // Stores the actual tensors
};

struct TABLE {                             // Loads a table
	string file; 														 // The file from which the tables was loaded
	unsigned int ncol;                       // The number of columns
	unsigned int nrow;                       // The number if rows
	vector <string> heading;                 // The headings for the columns
	vector <vector <string> > ele;        	 // The elements of the table
};

struct TIMEP { 														 // Stores a time period
	string name;														 // The name of the time period
	double tend;														 // The end time
};

struct TRANSDATA{                          // Stores data about transitions
	string fromstr;                          // The "from" compartment
 	string tostr;                            // The "to" compartment 
	unsigned int trans;                      // The transition number
	string type;                             // The type (either "reg" for regional or "all" for global)
 	string file;                             // The file name of the data 
	unsigned int start;                      // The start time for the data
	unsigned int units;                      // The units used (e.g. 1=days, 7=weeks) 
	vector <vector <unsigned int> > num;     // A table giving the number of that transition type
	unsigned int rows;                       // The number of rows of data
};

struct POPDATA{                            // Stores population data
	string compstr;                          // The compartment string
 	unsigned int comp;                       // The comparmtent number
	string type;                             // The type (either "reg" for regional or "all" for global)
 	string file;                             // The file name of the data 
	unsigned int start;                      // The start time for the data
	unsigned int units;                      // The units used (e.g. 1=days, 7=weeks) 
	vector <vector <unsigned int> > num;     // A table giving population numbers 
	unsigned int rows;                       // The number of rows of data
};

struct MARGDATA{                           // Stores data about marginal distributions
	string fromstr;                          // The "from" compartment
 	string tostr;                            // The "to" compartment 
	unsigned int trans;                      // The transition number
	string type;                             // The type (either "reg" for regional or "all" for global)
	unsigned int democat;                    // The demographic category
 	string file;                             // The file name of the data 
	vector < vector <double> > percent;      // A table giving the percentage of transition in different demographic values / regions
};

struct MEAS{                               // Stores values for all the measurements made on 
	vector < vector <vector <unsigned int> > > transnum;
	vector < vector <vector <unsigned int> > > popnum;
	vector < vector <vector <unsigned int> > > margnum;	
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
	unsigned int dp;                         // The demographic category possibility
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

class DataPipeline;

class MODEL;
class Inputs;

class DATA
{
public:
	
	DATA(const Inputs &inputs, const Details &details, const Mpi &mpi, DataPipeline *dp=0);

	DataPipeline *datapipeline;              // DataPipeline object
	
	string datadir; 												 // The data directory

	unsigned int threshold;                  // The limit under which numbers cannot be specified exactly 
	double thres_h;                          // The height of the threshold observation model

	unsigned int ndemocat;                   // The number of demographic categories
	vector <DEMOCAT> democat;                // Stores the demographic categories
	
	unsigned int ncovar;                     // The number of covariates for area  
	vector <COVAR> covar;                    // Covariates for area

	unsigned int nage;                       // The number of age categories
	unsigned int narage;                     // #area * #age
	unsigned int nardp;                      // #area * #ndemocatpos
	unsigned int ndemocatposperage;          // Demographic states per age group
	unsigned int nsettardp;                  // #sett * #area * #ndemocatpos

	vector <TIMEP> timeperiod;               // The timings of changes to Q;

	GENQ genQ; 															 // Stores information about generating the Q matrix
	vector <QTENSOR> Q;                      // Stores the list of Q tensors
	
	unsigned int nregion;                    // Number of data regions
	vector <REGION> region;                  // The names of the data regions

	unsigned int narea;                      // The number of areas
	vector <AREA> area;                      // List of all the areas
	
	vector <IND> ind;                        // The individuals in the system
		
	unsigned int popsize;                    // The total population size 
 
	unsigned int ndemocatpos;                // The number of demographic possibilities
	vector < vector<unsigned int> > democatpos; // Stores all the posible combinations of demographic categories

	vector <TRANSDATA> transdata;            // Store information about transition data
	
	vector <POPDATA> popdata;                // Store information about population data
	
	vector <MARGDATA> margdata;              // Store information about marginalised distribution data
	
	vector <double> agedist; 								 // Gives the overall age distribution
	
	void sortX(vector <unsigned int> &vec);	 // Used for sorting houses by x and y location
	void sortY(vector <unsigned int> &vec);	
	AreaRefComparatorX compX;
	AreaRefComparatorY compY;
	
	void print_to_terminal() const;
	
private:
	void calc_democatpos();
	void read_data_files(const Inputs &inputs, const Mpi &mpi);
	void load_region_file(const Inputs &inputs);
	
	string strip(string line) const;
	void copydata(unsigned int core);
	TABLE loadtable(string file, string dir="") const;
	TABLE loadtablefromdatapipeline(string file) const;
	TABLE loadtablefromfile(string file, string dir) const;

	void table_createcol(string head,vector <unsigned int> cols, TABLE &tab) const;
	void table_selectdates(unsigned int t, unsigned int units, TABLE &tab, string type) const;
	unsigned int findcol(const TABLE &tab, string name) const;
	unsigned int getint(const string& st, const string& file) const;
	
	void plotrawdata();   // These are temporary function used for analysis (they will be removed later)
	void generatedeathdata();
	void convertOAtoM();
	void convertRegion_M();
	
	const Details &details;
};

#endif
