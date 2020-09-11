#ifndef BEEPMBP__OUTPUT_HH
#define BEEPMBP__OUTPUT_HH

#include <fstream>

using namespace std;

#include "areatree.hh"
#include "simulate.hh"
#include "model.hh"

struct Sample{                                        // Stores information about a sample from the posterior
	Measurements meas;                                          // Stores measurements corresponding to the data file
	vector <double> R0;	                                // Time variation in R0
	vector <double> phi;	                              // Time variation in phi
};

struct ParamSample{                                     // Stores information about a sample from the posterior
	vector <double> paramval;                           // A parameter sample
};

struct Statistics{                                          // Stores statistical information
	string mean;                                        // The mean
	string CImin, CImax;                                // The minimum and maximum of the 90% credible interval
	string ESS;                                         // The estimated sample size
};

struct Distribution{                                          // Stores a probability distribution
	vector <string> value;
	vector <string> prob;
};

struct WeightedPoint{                                            // A weighted point
	double val;
	double w;
};

class ObservationModel;
struct Generation;

class Output
{
	public:
		Output(const Details &details, const Data &data, const Model &model, const ObservationModel &obsmodel);
		
		void print_populations(unsigned int sett, const vector <int> &N) const;
		void trace_plot_inititialise();
		void trace_plot(unsigned int samp, double Li, double Pri, unsigned int ninf, const vector <double> &paramval);
		void L_trace_plot_inititialise(unsigned int nchain_total);
		void L_trace_plot(unsigned int samp, const vector <double> &Litot);
		void final_results(const vector <ParamSample> &psamp, const vector <Sample> &opsamp) const;
		void save_event_sample(const vector < vector <Event> > &fev) const;
		void simulated_data(const vector < vector <EventRef> > &transev, const vector < vector <Event> > &indev, string dir) const;
		void combined_trace_file(const vector <string> &paramname, const vector < vector < vector <double> > > &vals, 
											 string file, string distfile, unsigned int burnin) const;
		void plot_distribution(string file, const Generation &gen) const;
		void generation_plot(string file, const vector <Generation> generation) const;
		double model_evidence_plot(string file, const vector <Generation> &generation) const;

	private:
		Statistics get_statastic(const vector <double> &vec) const;
		Statistics get_statistic_with_weight(vector <WeightedPoint> vec) const;	
		Distribution get_distribution(const vector <double> &vec) const;
		void posterior_plot(const vector <Sample> &opsamp, unsigned int d, unsigned int r, unsigned int type) const;
		void ensure_directory(const string &path) const;

		ofstream trace, traceLi;
		
		const Details &details;
		const Data &data;
		const Model &model;
		const ObservationModel &obsmodel;
};
#endif
