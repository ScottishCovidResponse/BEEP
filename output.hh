#ifndef BEEPMBP__OUTPUT_HH
#define BEEPMBP__OUTPUT_HH

#include <fstream>

using namespace std;

#include "poptree.hh"
#include "simulate.hh"
#include "model.hh"

struct SAMPLE{                                        // Stores information about a sample from the posterior
	MEAS meas;                                          // Stores measurements corresponding to the data file
	vector <double> R0;	                                // Time variation in R0
	vector <double> phi;	                              // Time variation in phi
};

struct PARAMSAMP{                                     // Stores information about a sample from the posterior
	vector <double> paramval;                           // A parameter sample
};

struct STAT{                                          // Stores statistical information
	string mean;                                        // The mean
	string CImin, CImax;                                // The minimum and maximum of the 90% credible interval
	string ESS;                                         // The estimated sample size
};

struct DIST{                                          // Stores a probability distribution
	vector <string> value;
	vector <string> prob;
};

struct PW{                                            // A weighted point
	double val;
	double w;
};

class Obsmodel;
struct Generation;

class Output
{
public:
	Output(const Details &details, const DATA &data, MODEL &model, Obsmodel &obsmodel);
	
	void trace_plot_init();
	void trace_plot(unsigned int samp, double Li, double Pri, unsigned int ninf, const vector <double> &paramval);
	void Li_trace_plot_init(unsigned int nchaintot);
	void Li_trace_plot(unsigned int samp, const vector <double> &Litot);
	void plot(string file, const vector < vector <FEV> > &xi, double tmin, double period) const;

	void results(const vector <PARAMSAMP> &psamp, const vector <SAMPLE> &opsamp) const;
	void eventsample(const vector < vector <FEV> > &fev) const;
	void simulateddata(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, string dir) const;
	void combinedtrace(const vector <string> &paramname, const vector < vector < vector <double> > > &vals, 
	                   string file, string distfile, unsigned int burnin) const;
	void plot_distribution(string file, const Generation &gen) const;
	void generation_plot(string file, const vector <Generation> generation) const;
	double model_evidence_plot(string file, const vector <Generation> &generation) const;

private:
	STAT getstat(const vector <double> &vec) const;
	STAT getstat_with_w(vector <PW> vec) const;	
	DIST getdist(const vector <double> &vec) const;
	void posterior_plot(const vector <SAMPLE> &opsamp, unsigned int d, unsigned int r, unsigned int type) const;
	void ensuredirectory(const string &path);

	ofstream trace, traceLi;
	
	const Details &details;
	const DATA &data;
	MODEL &model;
	Obsmodel &obsmodel;
};
#endif
