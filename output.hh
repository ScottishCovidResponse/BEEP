#ifndef BEEPMBP__OUTPUT_HH
#define BEEPMBP__OUTPUT_HH

#include <fstream>

using namespace std;

#include "struct.hh"

struct OutputLine{                                            // Specifies a line within an output plot
	string name;                                                // The name of the line 
	string file;                                                // The file providing the raw data
	unsigned int xcol, ycol;                                    // Which columns in the file provide x and y information
	LineType style;                                             // The style of the line
};


struct OutputPlot{                                            // Specifies a plot to appear in the final pdf
	OutputPlot(OutputPlotType type_, string title_, string xaxis_, string yaxis_, double min_, double max_);

	OutputPlotType type;                                        // The type of the output
	string title;                                               // The title of the plot
	string xaxis;                                               // The name of the x axis
	string yaxis;                                               // The name of the y axis
	double min, max;                                            // The minimum and maximum values in the plot
	vector <OutputLine> line;                                   // The lines within the plot_distribution
	
	void addline(string name, string file, unsigned int xcol, unsigned int ycol, LineType style);
};

struct Statistics{                                            // Stores statistical information
	string mean;                                                // The mean
	string CImin, CImax;                                        // The minimum and maximum of the 90% credible interval
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
		Output(const Details &details, const Data &data, const Model &model, const Inputs &inputs, const ObservationModel &obsmodel, Mpi &mpi);
		
		void generate_graphs(vector <Particle> &particle_store) const;
		void trace_plot_inititialise(string name, ofstream &trace) const;
		void trace_plot(unsigned int samp, double Li, const vector <double> &paramval, ofstream &trace) const;
		void simulated_data(const vector <double>& obs_value, string dir) const;
		void generation_results(const vector <Generation> &generation) const;
		Statistics get_statistic(const vector <double> &vec) const;
		void print_percentage(const double s, const unsigned int nsamp, unsigned int &percentage) const;
		
	private:
		void EF_datatable_plot(string file, const vector <Generation> &generation) const;
		void generate_graphs(vector <ParamSample> &psamp, const vector <Sample> &opsamp) const;
		void generation_plot(string file, const vector <Generation> &generation) const;
		void generate_pdf() const;
		void generate_pdf_description(const vector <OutputPlot> &op) const;
		void spline_plots(const vector <Sample> &opsamp, vector <OutputPlot> &opplot) const;
		void graph_plots(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void susceptibility_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void branch_prob_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void derived_parameter_distributions(const vector <Sample> &opsamp, vector <OutputPlot> &op) const;
		void posterior_parameter_distributions(const vector <ParamSample> &psamp, vector <OutputPlot> &op) const;
		void add_generation_plots(vector <OutputPlot> &op) const;
		void generate_gnuplot_file(const vector <OutputPlot> &op) const;
		
		void posterior_parameter_estimates(vector <ParamSample> &psamp) const;
		void spatial_R_map(const vector <Sample> &opsamp) const;
	
		void dirichlet_normalisation(vector <ParamSample> &psamp) const;
		Statistics get_statistic_with_weight(vector <WeightedPoint> vec) const;	
		Distribution get_distribution(const vector <double> &vec, const double priormin, const double priormax) const;
		void ensure_directories() const;
		void ensure_directory(const string &path) const;
		void print_model_specification() const;
		string rus(string st) const;

		void readme() const;
		
		OutputProp prop;                                        // Defines properties of the output
		
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
		Mpi &mpi;
};
#endif
