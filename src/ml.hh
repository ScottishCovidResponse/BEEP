#ifndef BEEP__ML_HH
#define BEEP__ML_HH

#include "struct.hh"
#include "state.hh"

class ML
{
public:	
	ML(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi);
	void run();

private:
	vector <ObsSlice> obs_slice;             // The observation at different time points

	Accuracy ac;                             // The accuracy of the caluculation (DOUBLE or FLOAT)
	
	double invT;                             // The inverse temperature for the power posterior
	
	unsigned int G;                          // The total number of generations 
	
	double cpu_time;                         // Maximum cpu for execution

	vector <Generation> generation;          // Stores information about each generation (cmaes)

	unsigned int nvar;                       // The number of parameters which are not fixed
	
	vector <unsigned int> var;               // The free parameters which can be changed
	
	unsigned int param_per_core;             // The number of parameters caluculated per core
	
	MLAlg algorithm;                         // The type of algorithm used
	
	unsigned int P;                          // The nunber of particles used to obtain the posterior samples
	
	unsigned int nsample_final;              // The numner of final posterior samples
	
	unsigned int npart;                      // The number of particles used for CMAES
	unsigned int npart_per_core;             // The number of particles per core 
			
	State state;                             // Stores the state of the system
	
	void find_parameters();
	void gradient_decent();
	void cmaes();
	void generate_samples(vector <ParamSample> &ps_per_core, const vector <double> &mean, const vector < vector <double> > &C, const double sigma, Generation &gen);
	vector < vector <double> > sample_param(const vector <double> &mean, const vector < vector <double> > &C, const double sigma, const unsigned int num) const;
	bool terminate_generation(unsigned int g, const vector <double> &EFbest_store) const;
	PointInfo calculate_point_info(vector <double> param);
	void calculate_mimimum();
	vector < vector <double> > calculate_covariance_martrix(const vector <double> &mean);
	double scale_covariance_martrix(const vector <double> &mean, const vector < vector <double> > &C);
	vector < vector <double> > calculate_hessian(const vector <double> &mean, MatrixType mattype);
	void calculate_heatmap();
	void BFGS(vector <double> &param);
	void find_jump_distance(double &eta, const vector <double> &param, vector <double> &param_new, const double L, double &L_new, const vector <double> &grad, const vector <double> &p);
	void BHHH(vector <double> &param);
	void newton_raphson(vector <double> &param);
	void model_evidence(const vector <double> &mean, const vector < vector <double> > &C, const double sigma);
	
	const Details &details;
	const Data &data;
	const Model &model;
	const Output &output;
	const ObservationModel &obsmodel;	
	Mpi &mpi;
};

#endif
