#ifndef BEEP__OUTPUT_HH
#define BEEP__OUTPUT_HH

#include <fstream>

using namespace std;

#include "utils.hh"

struct OutputLine{                                            // Specifies a line within an output plot
	string name;                                                // The name of the line 
	string file;                                                // The file providing the raw data
	unsigned int xcol, ycol;                                    // Which columns in the file provide x and y information
	unsigned int EB;                                            // The column which gives the error bar information (optional)
	LineType style;                                             // The style of the line
};


struct OutputPlot{                                            // Specifies a plot to appear in the final pdf
	OutputPlot(OutputPlotType type_, string title_, string fulldesc_, string tab_, string tab2_, string tab3_, string tab4_, string xaxis_, string yaxis_, double min_, double max_);

	OutputPlotType type;                                        // The type of the output
	string fulldesc;                                            // A full description of the output
	string tab, tab2, tab3, tab4;                               // Used for the menu
	string title;                                               // The title of the plot
	string xaxis;                                               // The name of the x axis
	string yaxis;                                               // The name of the y axis
	double min, max;                                            // The minimum and maximum values in the plot
	double xmin, xmax;                                          // The minimum and maximum x values in the plot (optional)
	vector <OutputLine> line;                                   // The lines within the plot_distribution
	
	vector <string> label;                                      // Stores label (in the case of marginals)
	bool legend_labelson;                                       // Determines if labels are put into description
	
	string spline_param_JSON;                                   // Describes parameters on splines
	
	string source;                                              // Stores information about the source files.
	
	void addline(const string name, const string file, const unsigned int xcol, const unsigned int ycol, const LineType style, const unsigned int EB = UNSET);
};

struct Statistics{                                            // Stores statistical information
	double mean;                                                // The mean
	double CImin, CImax;                                        // The minimum and maximum of the 90% credible interval
	double sd;                                                  // The standard deviation
};

struct Distribution{                                          // Stores a probability distribution
	bool variation;                                             // Determines if there is any variation in the quantity
	vector <string> value;	                                    // Stores the x values for the distribution
	vector <string> prob;                                       // Stores the probability values for the distributions
};

struct WeightedPoint{                                         // A weighted point
	double val;                                                 // The value of the point
	double w;                                                   // The weight attached to the point
};

class Output
{
	public:
		Output(const Details &details, const Data &data, const Model &model, Inputs &inputs, const ObservationModel &obsmodel, Mpi &mpi);
		
		void generate_graphs(vector <Particle> &particle_store, const double invT) const;
		void final_model_evidence(const vector <double> &ME_list, const double invT_final, const double cutoff_final) const;
		void trace_plot_inititialise(const string name, ofstream &trace) const;
		void trace_plot(const unsigned int samp, const double Li, const vector <double> &paramval, ofstream &trace) const;
		void simulated_data(const vector <double>& obs_value, const string dir) const;
		void set_generation_time(Generation &gen) const; 
		void generation_results(const vector <Generation> &generation) const;
		vector <unsigned int> get_effective_sample_size(const vector <ParamSample> &psamp) const;
		vector <double> get_correlation(const vector <ParamSample> &before, const vector <ParamSample> &after) const;
		vector <double>	get_Gelman_Rubin_statistic(const vector <ParamSample> &psamp, const vector <double> &w, const unsigned int nrun) const;
		vector <double>	get_Gelman_Rubin_statistic(const vector <ParamSample> &psamp) const;
		vector <double> get_Gelman_Rubin_statistic(const vector < vector < vector <double> > > &param_GR) const;
		Statistics get_statistic(const vector <double> &vec) const;
		Statistics get_statistic_75_percent(const vector <double> &vec) const;
		void print_percentage(const double s, const unsigned int nsamp, unsigned int &percentage) const;
		void ensure_directory(const string &path) const; // temporarily made public
		
	private:
		void EF_datatable_plot(const string file, const vector <Generation> &generation) const;
		void generate_graphs(vector <ParamSample> &psamp, const vector <Sample> &opsamp, const double invT) const;
		void generation_plot(const string file, const vector <Generation> &generation) const;
		void EF_dist(const vector <ParamSample> &psamp) const;
		void set_graph_source(vector <OutputPlot> &op) const;
		void generate_pdf(const string file, const string desc) const;
		void generate_visualisation(const vector <OutputPlot> &op, const string grfile) const;
		void generate_pdf_description(vector <OutputPlot> &op, const string grfile) const;
		void spline_plots(const vector <Sample> &opsamp, vector <OutputPlot> &opplot) const;
		void datatable_maps(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void graph_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op, const double invT) const;
		void get_line_colours(vector <LineColour> line_colour, vector <LineType> &lt, vector <LineType> &lt2) const;
		void posterior_parameter_estimates(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void susceptibility_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		vector <DemoDep> get_demodep(const string dep) const;
		void area_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void level_effect_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void mean_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void branch_prob_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void age_mixing_perturb_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void derived_parameter_distributions(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void posterior_parameter_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void spatial_R_map(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void age_mixing_matrix(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void covar_data(vector <OutputPlot> &op) const;
		void level_data(vector <OutputPlot> &op) const;
		void spatial_mixing_map(vector <OutputPlot> &op) const;
		void compartmental_model(vector <OutputPlot> &op) const;
		void add_generation_plots(vector <OutputPlot> &op) const;
		LineType get_linestyle(const unsigned int i) const;
		void add_trace_plots(vector <OutputPlot> &op) const;
		void add_democat_change(vector <OutputPlot> &op) const;
		void set_labelon(vector <OutputPlot> &op) const;
		void add_pred_timeplot(const string name, unsigned int time, vector <TimePlot> &pred_timeplot) const;
		vector <TimePlot> get_pred_timeplot() const;
		void generate_gnuplot_file(const vector <OutputPlot> &op, const string grfile) const;
		string load_boundaries() const;
	
		Statistics get_statistic_with_weight(vector <WeightedPoint> vec) const;	
		Distribution get_distribution(const vector <double> &vec, const double priormin, const double priormax) const;
		void ensure_directories();
		void print_model_specification() const;
		string rus(string st) const;
		string reference_label(const unsigned int i) const;
		string label(string lab) const;

		void readme() const;
		
		bool plot_param_values;                                 // Determines if simulated parameter values are put on plots
		
		OutputProp prop;                                        // Defines properties of the output
		unsigned int nrun;                                      // Stores the number of runs (used for trace plots)
		
		StateUncertainty stateuncer;                            // Determines if state uncertaity uses EB / curve
		
		string post_dir;                                        // The output directory depending on mode
		
		string boundaries;                                      // Loads up boundary information
		
		AreaPlot area_plot;                                     // Stores information about plotting areas
		
		const Inputs &inputs;
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
		Mpi &mpi;
};
#endif
