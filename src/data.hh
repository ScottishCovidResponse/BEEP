#ifndef BEEPMBP__Data_HH
#define BEEPMBP__Data_HH

#include <string>
#include <vector>
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

#include "struct.hh"
#include "details.hh"
#include "inputs.hh"
#include "mpi.hh"

class DataPipeline;

class Data
{
	public:
		Data(Inputs &inputs, const Details &details, Mpi &mpi, DataPipeline *dp=0);

		DataPipeline *datapipeline;              // DataPipeline object
		
		string data_directory; 						    	 // The data directory

		string areas_file;                       // The files which gives information about areas
		
		string threshold_str;                    // The string used for threshold
	
		string nodata_str;                       // The string used to represent no data
	
		string init_pop;                         // File containing initial population information (optional) 
	
		unsigned int ndemocat;                   // The number of demographic categories
		vector <DemographicCategory> democat;    // Stores the demographic categories
		
		unsigned int nstrain;                    // The number of disease strains
		vector <Strain> strain;                  // Information about the strains
		
		LevelEffect level_effect;                // Stores a level effect
		
		AreaEffect area_effect;                  // Stores the area effect
		
		unsigned int ncovar;                     // The number of covariates
		vector <Covariate> covar;                // Covariate information

		unsigned int nage;                       // The number of age categories
		unsigned int narage;                     // #area * #age
		unsigned int ndemocatpos_per_age;        // Demographic states per age group
		unsigned int ndemocatpos_per_strain;     // Demographic states per strain
	
		GenerateQ genQ; 								   			 // Stores information about age / geographical mixing
	
		unsigned int narea;                      // The number of areas
		vector <Area> area;                      // List of all the areas
			
		unsigned int popsize;                    // The total population size 
	 
		unsigned int ndemocatpos;                // The number of demographic possibilities
		vector < vector<unsigned int> > democatpos;// Stores all the posible combinations of demographic categories
		vector <string> democatpos_name;         // Stores as a name  

		vector <DataTable> datatable;            // Stores information about data tables
		
		vector <DemocatChange> democat_change;   // Stores information about demographic changes
		
		unsigned int nobs;                       // The number of observations
		vector <Observation> obs;                // A list of all the individual observations
		
		vector <Graph> graph;                    // A list of all the graphs which need to be plotted
	 	
		vector <double> agedist; 				  	  	 // Gives the overall age distribution
		
		vector < vector <double> > democat_dist; // Gives distributions in different demographic categories
		
		vector <double> democatpos_dist; 				 // Gives the overall demographic possibility distribution
		
		vector <Modification> modification;      // Used to implement modification to the model
		
		void chop_dir(string &file, const string dir) const;
		Table get_table(const string file, const string dir) const;
		vector <string> get_table_column(const string col_name, const Table &tab) const;
		vector <double> get_table_column(const unsigned int col, const Table &tab) const;
		vector <string> get_table_column_str(const unsigned int col, const Table &tab) const;
		void generate_file_from_table(const string file_data, const string file, const vector <string> &cols) const;
		string get_array_JSON(const string file, const string dir) const;
		string get_table_cols_JSON(const string file, const string dir, const vector <unsigned int> cols) const;
		string get_table_JSON(const string file, const string dir) const;
		Matrix get_matrix(const string file, const string dir) const;
		void make_table_with_time(const string file_in, const string file_out, const string col) const;
		void load_boundaries(const string file, vector < vector < vector <Coord> > > &bound) const;
		void create_boundaries(string x, string y, vector < vector < vector <Coord> > > &bound) const;
		void rescale_boundary(vector < vector < vector <Coord> > > &bound) const;
		void make_circle_boundary(const string xcol, vector < vector < vector <Coord> > > &bound) const;
		string observation_description(const DataType type, const string obs, const unsigned int timestep = 1) const;
		string print() const;
	
	private:
		void calc_democatpos();
		void read_data_files(Inputs &inputs, Mpi &mpi);
		void copy_data(unsigned int core);
		Table load_table(const string file, const string dir="", const bool heading=true, const bool supop=false) const;
		Table load_table_from_datapipeline(const string file) const;
		Table load_table_from_file(const string file, const string dir, const bool heading, const bool supop, const char sep) const;
		void read_covars();
		void read_level_effect();
		void check_or_create_column(Table &tab, string head, unsigned int d) const;
		void read_initial_population(Table &tab, Inputs &inputs);
		void read_init_pop_file(const vector <string> comps);
		void read_initial_population_areas(const unsigned int co_sus, Table &tab);
		void load_modification(Inputs &inputs, const Table &tabarea);
		
		void filter_areas(Table &tab);
		void filter_table(const string st, Table &tab) const;
		vector < vector <double> > get_demo_frac(const unsigned int d, const Table &tab, const vector <double> &pop) const;
		void load_democat_change(const Table &tabarea);
		void load_datatable(const Table &tabarea);
		void normal_approximation();
		void check_datatable();
		void load_timeseries_datatable(const Table &tabarea, const unsigned int i, const bool sim);
		void load_marginal_datatable(const Table &tabarea, const unsigned int i, const bool sim);
		vector <DataFilter> get_datafilter(const Table &tabarea, const string geo_dep, const string democats_dep, const string geo_filt, const string democats_filt, const DataType type, const string observation, const string datafile) const;
		void table_create_column(const string head, const vector <unsigned int> &cols, Table &tab) const;
		void table_select_dates(int &start, int &end, unsigned int &timestep, Table &tab, const DataType type, const unsigned int shift, vector <int> &times) const;
		unsigned int find_column(const Table &tab, const string name) const;
		unsigned int find_column_noerror(const Table &tab, const string name) const;
		void generate_matrices();
		void plotmat(const Matrix &mat, const string &title);
		void geo_normalise(SparseMatrix &mat);
		void agematrix_normalise(Matrix &mat);
		Matrix age_mixing_matrix(const Table &tab) const;
		SparseMatrix load_geo_mixing_matrix(const Table &tab) const;
		SparseMatrix load_geo_mixing_matrix_sparse() const;
		vector <GeographicMap> create_geomap(const Table &tab, const string geo) const;
		vector <unsigned int> create_dp_sel(const string dp_str, const string em) const;
		vector <unsigned int> create_area_sel(const Table &tabarea, const string str, const string em) const;
		void set_datatable_weights();
		bool vector_contains(const vector <unsigned int> &vec, const unsigned int num) const;
		bool vector_contains(const vector <string> &vec, const string num) const;
		void vector_remove(vector <unsigned int> &vec, const unsigned int num) const;
		vector <TreeNode> area_split(SparseMatrix M) const;
		void split_in_two(const vector <unsigned int> &arearef, vector <unsigned int> &ch1,  vector <unsigned int> &ch2, const vector <vector <double> > &G) const;
		double get_split_fit(const vector <unsigned int> &area_list, const vector <bool> &gr, const vector <vector <double> > &G) const;
		double get_split_fit_dif(const vector <unsigned int> &area_list, vector <bool> &gr, const vector <vector <double> > &G, unsigned int isel) const;
		void print_obs() const;
		
		/* Used in data_boundary */
		FileType filetype(const string file) const;
		void load_KML(const string file, vector < vector < vector <Coord> > > &bound) const;
		void placemark_KML(XMLNode* root, vector < vector < vector <Coord> > > &bound) const;
		unsigned int get_area_name_KML(XMLNode* child) const;
		void get_polygon_KML(XMLNode* child, const unsigned int c, vector < vector < vector <Coord> > > &bound) const;
		void get_coordinate_KML(XMLNode* child, const unsigned int c, vector < vector < vector <Coord> > > &bound) const;
		void load_geojson(const string file, vector < vector < vector <Coord> > > &bound) const;	
			
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
		void generate_initpop();
		void generate_tvcovar();
		void generate_level_effect();
		void reducesize_geojson(const string file);
		void generate_age_mixing_perturb();
		void generate_geo_mixing_matrix();
	
		const Details &details;
		Inputs &inputs;
};

#endif
