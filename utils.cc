// Stores functions for sampling from different distributions and deals with error messages

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <stdexcept>
#include "stdlib.h"
#include "math.h"
#include <sys/stat.h>
#include <cstring>

#include "utils.hh"
#include "consts.hh"

using namespace std;

static std::mt19937 mt(0);

void sran(int seed)
{
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
	//std::uniform_real_distribution<> dist(0.,1.);
	std::uniform_real_distribution<> dist(0.0000000001,0.999999999);
	return dist(mt);
#endif
}

/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(float mu, float sd)
{
	return mu + sd*sqrt(-2*log(ran()))*cos(2*M_PI*ran());
}

/// Samples from the exponential distribution using a rate
double exp_sample(double rate)
{
	if(rate < 0) emsgEC("Utils",8);
	if(rate == 0) return LARGE;
	return - log(ran())/rate;
}

/// Samples from the exponential distribution using a mean time
double exp_sample_time(double time)
{
	return - log(ran())*time;
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gamma_sample(double a, double b)
{
  if(a < 0 || b < 0) emsgEC("Utils",1);

  if(a < 1){
    double u = ran();
    return gamma_sample(1.0+a,b)*pow(u,1.0/a);
  }
  else{
    double d = a - 1.0/3.0;
    double c = (1.0/3.0)/sqrt(d);
		double v, x;
    while(true){
      do{
        x = sqrt(-2*log(ran()))*cos(2*M_PI*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      double u = ran();

      if(u < 1 - 0.0331*x*x*x*x) break;

      if(log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
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

/// The log of the lognormal probability distribution
double lognormal_probability(double x, double mean, double var)
{
	if(x <= 0) emsgEC("Utils",4);
	if(var <= 0) emsgEC("Utils",5); 
	auto val = log(x);
	return -0.5*log(2*M_PI*var) - (mean-val)*(mean-val)/(2*var) - val;
}

double lognormal_sample(double mean, double sd)
{
	return exp(normal_sample(mean,sd));
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

/// @cond EMSG

// Alter error behaviour to throwing exception, which can be caught in unit
// tests
bool emsg_throws = false;

/// Displays an error message
void emsg(const string& msg)
{
	if(emsg_throws)
		throw(std::runtime_error(msg));
	cout << msg << endl;
	MPI_Barrier(MPI_COMM_WORLD);
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
	} else if(st == "NA") {
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
