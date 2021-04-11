#ifndef BEEPMBP__Data_HH
#define BEEPMBP__Data_HH

#include <string>

using namespace std;

#include "struct.hh"
#include "details.hh"
#include "inputs.hh"
#include "mpi.hh"

class DataPipeline;

class Data
{
	public:
		Data(const Inputs &inputs, const Details &details, Mpi &mpi, DataPipeline *dp=0);

		DataPipeline *datapipeline;              // DataPipeline object
		
		string data_directory; 						    	 // The data directory

		unsigned int threshold;                  // The limit under which numbers cannot be specified exactly 
	
		unsigned int ndemocat;                   // The number of demographic categories
		vector <DemographicCategory> democat;    // Stores the demographic categories
		
		vector <unsigned int> age_from, age_to;  // Time ranges for different age classifications
		
		unsigned int ncovar;                     // The number of covariates for area  
		vector <Covariate> covar;                // Coariate information

		vector <GeographicMap> region_effect;    // If a regional effect is added this stores the map
		ParamSpec sigma;                         // This gives the parameter specification for the standard deviation sigma
		
		unsigned int nage;                       // The number of age categories
		unsigned int narage;                     // #area * #age
		unsigned int nardp;                      // #area * #ndemocatpos
		unsigned int ndemocatposperage;          // Demographic states per age group
		unsigned int nsettardp;                  // #sett * #area * #ndemocatpos
	
		GenerateQ genQ; 								   			 // Stores information about age / geographical mixing
	
		unsigned int narea;                      // The number of areas
		vector <Area> area;                      // List of all the areas
		
		vector <Individual> ind;                 // The individuals in the system
			
		unsigned int popsize;                    // The total population size 
	 
		unsigned int ndemocatpos;                // The number of demographic possibilities
		vector < vector<unsigned int> > democatpos; // Stores all the posible combinations of demographic categories

		vector <DataTable> datatable;            // Stores information about data tables
		
		unsigned int nobs;
		vector <Observation> obs;                // A list of all the individual observations
		
		vector <Graph> graph;                    // A list of all the graphs which need to be plotted
	 	
		vector <double> agedist; 				  	  	 // Gives the overall age distribution
		
		vector < vector <double> > democat_dist; // Gives distributions in different demographic categories
		
		vector <double> democatpos_dist; 				 // Gives the overall demographic possibility distribution
		
		vector <CounterFact> counterfact;        // Used to implement counterfactuals
		
		string print() const;
	
	private:
		void remove_empty_rows(Table& tab) const;
		void calc_democatpos();
		void read_data_files(const Inputs &inputs, Mpi &mpi);
		void copy_data(unsigned int core);
		Table load_table(string file, string dir="", bool heading=true) const;
		Table load_table_from_datapipeline(string file) const;
		Table load_table_from_file(string file, string dir, bool heading, char sep) const;
		void read_areas(Table &tab, string file);
		void load_counterfactual(const Inputs &inputs, const Table &tabarea);
		
		void filter_areas(Table &tab);
		void filter_table(const string st, Table &tab) const;
		void load_datatable(const Table &tabarea);
		void check_datatable();
		void load_timeseries_datatable(const Table &tabarea, unsigned int i, bool sim);
		void load_marginal_datatable(const Table &tabarea, unsigned int i, bool sim);
		vector <DataFilter> get_datafilter(const Table &tabarea, string geo, string democats);
		void load_region_effect(const Inputs &inputs, const Table& tab);
		void table_create_column(string head, const vector <unsigned int> &cols, Table &tab) const;
		void table_select_dates(unsigned int t, unsigned int units, Table &tab, string type, unsigned int shift) const;
		unsigned int find_column(const Table &tab, string name) const;
		unsigned int get_integer(const string& st, const string& file) const;
		void table_add_age(string name, unsigned int ti, unsigned int tf, Table &tab);
		void generate_matrices();
		void plotmat(const Matrix& mat, const string& title);
		void geo_normalise(SparseMatrix &mat);
		void agematrix_normalise(Matrix &mat);
		Matrix matfromtable(const Table& tab, unsigned int N);
		SparseMatrix loadsparse(const string& file, unsigned int N, GenerateQ &genQ);
		vector <GeographicMap> create_geomap(const Table &tab, string geo) const;
		vector <unsigned int> create_dp_sel(string dp_str) const;
		void set_datatable_weights();
		void coarse_grain(const Table &tab, string coarse);
		vector <unsigned int> convert_areas(string type, vector <unsigned int> &vec, const vector <GeographicMap> &map) const;
		bool vector_contains(const vector <unsigned int> &vec, unsigned int num) const;
		bool vector_contains(const vector <string> &vec, string num) const;
		void vector_remove(vector <unsigned int> &vec, unsigned int num) const;
		void print_obs() const;
		
		/* Used for raw data analysis */
		void raw();
		void IFR();
		void case_age_distribution();
		void death_age_distribution();
		string remove_comma(string st);
		void calculate_IR();
		void fractional_change(string file);
		void fractional_change2(string file);
		void plotrawdata(); 
		void cases_age();
		void generatehospdata();		
		void generatedeathdata_scotland();
		void generatedeathdata_england_wales();
		void convert_Mdata();
		void generate_admission_age();
		void split_deaths_age_data();
		void deaths_age();
		void deaths_hospital_england();
	
		const Details &details;
};

#endif
