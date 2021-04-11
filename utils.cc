// Stores functions for sampling from different distributions and deals with error messages

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <stdexcept>
#include "stdlib.h"
#include "math.h"
#include <sys/stat.h>
#include <cstring>

using namespace std;

#include "utils.hh"
#include "consts.hh"

static std::mt19937 mt(0);

default_random_engine generator;


/// Sets the seed for the random number generator
void sran(int seed)
{
	generator.seed(seed);
#ifdef OLD_RAND
	srand(seed);
#else
	mt = std::mt19937(seed);
#endif
}


/// Draws a random number between 0 and 1
double ran()
{
#ifdef OLD_RAND
	if(RAND_MAX == 32767) {
		double v = (double(rand())*32766.0+rand())/(32767.0*RAND_MAX);
		if(v == 0 || v == 1) {
			return 0.1;
		}
		else {
			return v;
		}
	}
	else {
		return double(0.999999999*rand())/RAND_MAX;
	}
#else
	std::uniform_real_distribution<> dist(0.0000000001,0.999999999);
	return dist(mt);
#endif
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(float mu, float sd)
{
	normal_distribution<double> distribution(mu,sd);
	return distribution(generator);
}


/// Samples from the exponential distribution using a rate
double exp_sample(double rate)
{
	if(rate <= 0){ 
		if(rate < -TINY){ cout << rate << "rate\n"; emsgEC("Utils",8);}
		return LARGE;
	}
	exponential_distribution<double> distribution(rate);
	return distribution(generator);
}


/// Samples from the exponential distribution using a mean time
double exp_sample_time(double time)
{
	return -log(ran())*time;
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gamma_sample(double a, double b)
{
	gamma_distribution<double> distribution(a,1.0/b);
	return distribution(generator);
}


/// The log of the probability from the normal distribution
double normal_probability(double x, double mean, double var)
{
  if(var <= 0) emsgEC("Utils",2);
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}


/// The log of the probability from the gamma distribution
double gamma_probability(double x, double a, double b)
{
	// Scale parameter is 1/b in several sources
  if(x <= 0 || a <= 0 || b <= 0) emsgEC("Utils",3);
  return (a-1)*log(x) - b*x + a*log(b) - lgamma(a);
}


/// The log of the probability from the negative binomial distribution
double negative_binomial_probability(unsigned int val, double m, double shape)
{
	return lgamma(shape+val) - lgamma(val+1)- lgamma(shape) + shape*log(shape/(shape+m)) + val*log(m/(shape+m));
}


/// The log of the lognormal probability distribution
double lognormal_probability(double x, double mean, double var)
{
	if(x <= 0) emsgEC("Utils",4);
	if(var <= 0) emsgEC("Utils",5); 
	auto val = log(x);
	return -0.5*log(2*M_PI*var) - (mean-val)*(mean-val)/(2*var) - val;
}


///Generates a sample from the log-lornormal probability distribution
double lognormal_sample(double mean, double sd)
{
	return exp(normal_sample(mean,sd));
}


/// Generates a sample from the Poisson distribution
int poisson_sample(double lam)
{
	if(lam > LARGE) emsgEC("Utils",10);
  poisson_distribution<int> distribution(lam);
	return distribution(generator);
}


/// The log probability of the poisson distribution
double poisson_probability(int i, double lam)
{
	if(lam < 0) emsg("Poisson rate negative");
  if(lam == 0){
    if(i == 0) return 0;
    else return -LARGE;
  }
  else{
    if(i < 0) return -LARGE;
    else{
			return i*log(lam)-lam - lgamma(i+1);
		}
  }
}


/// A sample from the binomial distribution
unsigned int binomial_sample(double p, unsigned int n)
{
	if(n == 0) return 0;
		
	binomial_distribution<int> distribution(n,p);
	return distribution(generator);
}


/// Checks the binomial sampler
void binomial_check()
{
	const auto loopmax = 100000u;
	
	auto n = 100u;
	auto p = 0.0;
	for(p = 0; p < 1; p += 0.1){
		auto av = 0.0, av2 = 0.0;
		for(auto loop = 0u; loop < loopmax; loop++){
			auto X = binomial_sample(p,n);
			av += X; av2 += X*X;
		} 
	
		cout << "p: " << p << " " << av/loopmax << " " << p*n << " "
    		<< av2/loopmax - (av/loopmax)*(av/loopmax) << " " << n*p*(1-p) <<  "\n";
	}
	
	p = 0.05;
	for(n = 0; n < 1000; n++){
		auto av = 0.0, av2 = 0.0;
		for(auto loop = 0u; loop < loopmax; loop++){
			auto X = binomial_sample(p,n);
			av += X; av2 += X*X;
		} 
	
		cout << "n: " << n << " " << av/loopmax << " " << p*n << " "
    		<< av2/loopmax - (av/loopmax)*(av/loopmax) << " " << n*p*(1-p) <<  "\n";
	}
}


/// Split up a string at a specified delimiter
vector<string> split(const string& s, char delimiter)                                                               
{                              
  std::vector<std::string> splits;                       
  std::string split;                                      
  std::istringstream ss(s);                               
  while (std::getline(ss, split, delimiter)) splits.push_back(split);                                                  
  return splits;                                           
}


/// Removes '\r' and quotations from a string
string strip(string line) 
{
	unsigned int len;
	
	len = line.length();
	if(len > 0 && line.substr(len-1,1) == "\r") line = line.substr(0,len-1);
	len = line.length();
	if(len > 0 && line.substr(0,1) == "\"") line = line.substr(1,len-1);
	len = line.length();
	if(len > 0 && line.substr(len-1,1) == "\"") line = line.substr(0,len-1);
	
	return line;
}	

/// @cond EMSG

// Alter error behaviour to throwing exception, which can be caught in unit
// tests
bool emsg_throws = false;

void emsgfile(string st)
{
	ofstream ems("error.txt");
	ems << st << " msg\n";
}

/// Displays an error message
void emsg(const string& msg)
{
	emsgfile(msg);
	
	if(emsg_throws)
		throw(std::runtime_error(msg));
	cout << msg << endl;
	exit (EXIT_FAILURE);
}


/// Displays an internal error message
void emsgEC(const string& section, unsigned int ec)
{
	std::ostringstream oss;
	oss << "Unfortunately BEEPmbp has generated an internal error. We are very sorry about this!" << endl;
	oss << "The error occurred in '" << section << "' with code '" << ec << "'";
	emsg(oss.str());
}


/// Displays an error message on the root core
void emsgroot(const string& msg)
{
	int core;
	MPI_Comm_rank(MPI_COMM_WORLD,&core);
	if(core == 0) emsg(msg);
	MPI_Barrier(MPI_COMM_WORLD);
	exit (EXIT_FAILURE);
}
/// @endcond


string filebasename(const string &path)
{
	auto filecomps = split(path,'/');
	return filecomps.at(filecomps.size()-1);
}


bool stringhasending (std::string const &fullString, std::string const &ending)
{
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}


/// Gets a positive integer from a string
unsigned int get_integer(
	const std::string& st,
	unsigned int threshold
	)
{
	char* endptr;
	
	long longi = strtol(&st[0],&endptr,10);
	ptrdiff_t j = endptr-&st[0];
	unsigned int i;
	
	if (st.length()-j == 0 && longi >= 0
			&& longi <= std::numeric_limits<unsigned int>::max()) {
		i = static_cast<unsigned int>(longi);
	} else if(st == "NA" || st == "missing") {
		i = UNKNOWN;
	} else if (st == "*") {
		i = THRESH;
		if (threshold == UNSET) {
			throw(std::runtime_error(
							"since '*' is used "
							"there must be a threshold set with the 'threshold' "
							"command in the input TOML file."));
		}
	}	else {
		throw (std::runtime_error("the quantity '"+st+"' is not a number"));
	}
	return i;
}


/// Inverts a matrix
vector <vector <double> > invert_matrix(const vector <vector <double> > &mat)   
{
	unsigned int nvar = mat.size();
	vector <vector <double> > inv_M;
	
	double A2[nvar][nvar];

	inv_M.resize(nvar);
  for(auto i = 0u; i < nvar; i++){
		inv_M[i].resize(nvar);
    for(auto j = 0u; j < nvar; j++){
      A2[i][j] = mat[i][j];
      if(i == j) inv_M[i][j] = 1; else inv_M[i][j] = 0;
    }
  }

  for(auto ii = 0u; ii < nvar; ii++){
    double r = A2[ii][ii];
    for(auto i = 0u; i < nvar; i++){
      A2[ii][i] /= r; inv_M[ii][i] /= r; 
    }

    for(auto jj = ii+1; jj < nvar; jj++){
      double r = A2[jj][ii];
			if(r != 0){
				for(auto i = 0u; i < nvar; i++){ 
					A2[jj][i] -= r*A2[ii][i];
					inv_M[jj][i] -= r*inv_M[ii][i];
				}
			}
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
    for(int jj = ii-1; jj >= 0; jj--){
      double r = A2[jj][ii];
			if(r != 0){
				for(auto i = 0u; i < nvar; i++){ 
					A2[jj][i] -= r*A2[ii][i];
					inv_M[jj][i] -= r*inv_M[ii][i];
				}
			}
    }
  }

	if(false){ // checks inverse
		for(auto j = 0u; j < nvar; j++){
			for(auto i = 0u; i < nvar; i++){
				double sum = 0; for(auto ii = 0u; ii < nvar; ii++) sum += mat[j][ii]*inv_M[ii][i];
				
				if(i != j){ if(sum < -TINY || sum > TINY) emsg("prob");}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsg("prob");}		
			}
		}
	}
	
	return inv_M;
}

vector <double> initial_guess;

/// Determines the largest eigenvalue and vector from a matrix
double largest_eigenvalue(const vector < vector <double> > &M, vector <double> &eigenvector)
{
	auto N = M.size();
	
	if(false){
		cout << "mat\n";
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << M[i][j] << " ";
			cout << " M\n";
		}
	}
	
	if(initial_guess.size() == 0){
		for(auto i = 0u; i < N; i++) initial_guess.push_back(ran());
	}
	
	vector <double> vec = initial_guess;
	
	auto ev = 0.0;
	auto loop = 0u;
	do{
		vector <double> vec2(N);
		for(auto i = 0u; i < N; i++){
			auto sum = 0.0; for(auto ii = 0u; ii < N; ii++) sum += M[i][ii]*vec[ii];
			vec2[i] = sum;
		}
		
		auto sum = 0.0; for(auto i = 0u; i < N; i++) sum += vec2[i];
		for(auto i = 0u; i < N; i++) vec[i] = vec2[i]/sum;
			
		auto ev_new = sum;
		
		sum = 0.0; for(auto i = 0u; i < N; i++) sum += vec[i];
		if(ev_new-ev > -TINY && ev_new-ev < TINY) break;
		ev = ev_new;
		loop++;
		if(loop > 10000) emsg("Convergence problem"); 
	}while(1 == 1);

	eigenvector = vec;

	initial_guess = vec;

	return ev;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <unsigned int> vec, const unsigned int val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <string> vec, const string val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Outputs a number to a given precision
string prec(double num, unsigned int pre)
{
	stringstream ss; ss << std::fixed; ss.precision(pre); 
	ss << num;
	return ss.str();
}


/// Represents a number as a percentage
string per(string per)
{
	return prec(100*atof(per.c_str()),0)+"%";
}

/// Represents a number as a percentage
string per(double per)
{
	return prec(100*per,0)+"%";
}
