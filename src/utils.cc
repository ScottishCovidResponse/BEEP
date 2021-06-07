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
#include <signal.h>

using namespace std;

const bool debug = false;

#include "utils.hh"
#include "consts.hh"

static std::mt19937 mt(0);

default_random_engine generator;


/// Sets the seed for the random number generator
void sran(const int seed)
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
double normal_sample(const double mu, const double sd)
{
	normal_distribution<double> distribution(mu,sd);
	return distribution(generator);
}


/// Samples from the exponential distribution using a rate
double exp_sample(const double rate)
{
	if(rate <= 0){ 
		if(rate < -TINY) emsgEC("Utils",1);
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
double gamma_sample(const double a, const double b)
{
	gamma_distribution<double> distribution(a,1.0/b);
	return distribution(generator);
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double var)
{
  if(var <= 0) emsgEC("Utils",2);
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double a, const double b)
{
	// Scale parameter is 1/b in several sources
  if(x <= 0 || a <= 0 || b <= 0) emsgEC("Utils",3);
  return (a-1)*log(x) - b*x + a*log(b) - lgamma(a);
}


/// The log of the probability from the negative binomial distribution
double negative_binomial_probability(const unsigned int val, const double m, const double shape)
{
	return lgamma(shape+val) - lgamma(val+1)- lgamma(shape) + shape*log(shape/(shape+m)) + val*log(m/(shape+m));
}


/// The log of the lognormal probability distribution
double lognormal_probability(const double x, const double mean, const double var)
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
int poisson_sample(const double lam)
{
	if(lam > LARGE) emsgEC("Utils",6);
  poisson_distribution<int> distribution(lam);
	return distribution(generator);
}


/// The log probability of the poisson distribution
double poisson_probability(const int i, const double lam)
{
	if(lam < 0) emsgEC("Utils",7);
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
unsigned int binomial_sample(const double p, const unsigned int n)
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
    		<< av2/loopmax - (av/loopmax)*(av/loopmax) << " " << n*p*(1-p) << endl;
	}
	
	p = 0.05;
	for(n = 0; n < 1000; n++){
		auto av = 0.0, av2 = 0.0;
		for(auto loop = 0u; loop < loopmax; loop++){
			auto X = binomial_sample(p,n);
			av += X; av2 += X*X;
		} 
	
		cout << "n: " << n << " " << av/loopmax << " " << p*n << " "
    		<< av2/loopmax - (av/loopmax)*(av/loopmax) << " " << n*p*(1-p) << endl;
	}
}


/// Split up a string at a specified delimiter
vector<string> split(const string &s, char delimiter)                                                               
{                              
  std::vector<std::string> splits;                       
  std::string split;                                      
  std::istringstream ss(s);                               
  while (std::getline(ss, split, delimiter)) splits.push_back(split);   

	for(auto &spl : splits) strip(spl);
	
  return splits;                                           
}


/// Coverts a string to lower case
string toLower(string st)
{	
	transform(st.begin(), st.end(), st.begin(), ::tolower);
	return st;
}

/// Removes space, '\r' and quotations from a string
void strip(string &line) 
{
	while(line.length() > 0 && (line.substr(0,1) == "\r" || line.substr(0,1) == "\"" || line.substr(0,1) == " ")){
		line = line.substr(1,line.length()-1);
	}
	
	while(line.length() > 0 &&  (line.substr(line.length()-1,1) == "\r" || line.substr(line.length()-1,1) == "\"" 
	                           || line.substr(line.length()-1,1) == " ")){
		line = line.substr(0,line.length()-1);
	}
}	

/// @cond EMSG

// Alter error behaviour to throwing exception, which can be caught in unit
// tests
bool emsg_throws = false;

/// Displays an error message
void emsg(const string& msg)
{
	if(emsg_throws) throw(std::runtime_error(msg));
	cout << msg;
	if(msg.length() > 0 && msg.substr(msg.length()-1,1) != ".") cout << ".";
	cout << endl;
	
	if(debug == true) raise(SIGABRT);
	exit (EXIT_FAILURE);
}


/// Displays an internal error message
void emsgEC(const string &section, unsigned int ec)
{
	std::ostringstream oss;
	oss << "Unfortunately BEEPmbp has generated an internal error. We are very sorry about this!" << endl;
	oss << "The error occurred in '" << section << "' with code '" << ec << "'";
	emsg(oss.str());
}


/// Displays an error message on the root core
void emsgroot(const string &msg)
{
	int core;
	MPI_Comm_rank(MPI_COMM_WORLD,&core);
	if(core == 0) emsg(msg);
	MPI_Barrier(MPI_COMM_WORLD);
	if(debug == true) raise(SIGABRT);
	exit (EXIT_FAILURE);
}
/// @endcond

/// Displays a warning message on the root core
void warning(const string &msg)
{
	int core;
	MPI_Comm_rank(MPI_COMM_WORLD,&core);
	if(core == 0) cout << "WARNING: " << msg << endl;
}

string filebasename(const string &path)
{
	auto filecomps = split(path,'/');
	return filecomps.at(filecomps.size()-1);
}


bool stringhasending(std::string const &fullString, std::string const &ending)
{
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
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
				
				if(i != j){ if(sum < -TINY || sum > TINY) emsgEC("Utils",8);}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsgEC("Utils",9);}		
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
		cout << "mat" << endl;
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++) cout << M[i][j] << " ";
			cout << " M" << endl;
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
		if(loop > 10000) emsg("Eigen-vector convergence problem"); 
	}while(1 == 1);

	eigenvector = vec;

	initial_guess = vec;

	return ev;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <unsigned int> &vec, const unsigned int val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <string> &vec, const string &val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds if a given character exists in a string
unsigned int find_char(const string st, const string char_in)
{
	for(auto i = 0u; i < st.length(); i++){
		auto j = 0u; while(j < char_in.length() && char_in.substr(j,1) != st.substr(i,1)) j++;
		if(j < char_in.length()) return i;
	}
	return UNSET;
}


/// Finds the maximum value in a vector
double vec_max(const vector <double> &vec)
{
	double num = -LARGE;
	for(auto val : vec){ if(val != UNSET && val > num) num = val;}
	if(num == -LARGE) return UNSET;
	else return num;
}	


/// Outputs a number to a given precision
string prec(const double num, const unsigned int pre)
{
	stringstream ss; ss << std::fixed; ss.precision(pre); 
	ss << num;
	return ss.str();
}


/// Represents a number as a percentage
string per(const string per)
{
	return prec(100*atof(per.c_str()),0)+"%";
}


/// Represents a number as a percentage
string per(const double per)
{
	return prec(100*per,0)+"%";
}


/// Checks a string has a certain set of characters
bool allow_string(const string st, const string ok_char)
{
	for(auto i = 0u; i < st.length(); i++){
		auto j = 0u; while(j == ok_char.length() && st.substr(i,1) != ok_char.substr(j,1)) j++;
		if(j == ok_char.length()) return false;
	}
	return true;
}


/// Return a double from a string and allows for the possibility of unset and to be set
double get_double_with_unset(string st, const string em)
{
	strip(st);
	if(st == "") return UNSET;
	else return get_double(st,em);
}


/// Return a double from a string and allows for the possibility of unset and to be set
double get_double_with_tobeset(string st, const string em)
{
	strip(st);
	if(st == "*") return TOBESET;
	else return get_double(st,em);
}


/// Return a double from a string and ensures it is positive
double get_double_positive(const string st, const string em)
{
	auto val = get_double(st,em);
	if(val < 0) emsgroot(em+" the expression '"+st+"' is not positive.");
	return val;
}


/// Return a double from a string
double get_double(string st, string em)
{
	strip(st);
	if(allow_string(st,"-0123456789.") == false) emsgroot(em+" the expression '"+st+"' is not a number.");
	else{
		char* endptr;
		double val = strtod(&st[0],&endptr);
		ptrdiff_t j = endptr-&st[0];
		if(j != (int)st.length()) emsgroot(em+" the expression '"+st+"' is not a number.");
		return val;
	}
}


/// Return an unsigned int from a string
unsigned int get_int(string st, const string em)
{
	strip(st);
	
	if(allow_string(st,"-0123456789") == false) emsgroot(em+" the expression '"+st+"' is not an integer.");

	char* endptr;
	int val = strtol(&st[0],&endptr,10);
	ptrdiff_t j = endptr-&st[0];
	if(j != (int)st.length()) emsgroot(em+" the expression '"+st+"' is not an integer.");
	if(val < 0) emsgroot(em+" the expression '"+st+"' cannot be negative.");
	
	return val;
}


/// Return an unsigned int from a stringand returns UNSET if not possible
unsigned int get_int_error(string st)
{
	strip(st);
	
	if(allow_string(st,"-0123456789") == false) return UNSET;
	int val = atoi(st.c_str());
	if(val == 0 && st != "0") return UNSET;
	return val;
}


/// Gets a data value
double get_data(string val, const string em, const string &threshold_str, const string &nodata_str)
{
	strip(val);
	if(val == threshold_str) return THRESH;
	else{
		if(val == nodata_str) return UNKNOWN;
		else return get_double(val,em);
	}
}	


/// Replace one string with another
string replace(string st, const string s1, const string s2)
{
	auto i = 0u;
	while(i + s1.length() < st.length()){
		if(st.substr(i,s1.length()) == s1){
			st = st.substr(0,i)+s2+st.substr(i+s1.length());
			i += s2.length();
		}
		else i++;
	}
	
	return st;
}


/// Return a parameter fixed to one
ParamSpec ps_one()
{
	ParamSpec ps; ps.name = "one"; ps.value = "1"; ps.prior = "Fixed(1)"; ps.smooth = UNSET; ps.factor = 1;
	return ps;
}
	

/// Return a parameter fixed to one
ParamSpec ps_zero()
{
	ParamSpec ps; ps.name = "zero"; ps.value = "0"; ps.prior = "Fixed(0)"; ps.smooth = UNSET; ps.factor = 1;
	return ps;
}

