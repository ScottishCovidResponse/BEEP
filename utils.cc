// Stores functions for sampling from different distributions and deals with error messages

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "stdlib.h"
#include "math.h"

#include "utils.hh"
#include "consts.hh"

using namespace std;

/// Draws a random number between 0 and 1
double ran()
{
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
}

/// Draws a normally distributed number with mean mu and standard deviation sd
double normal(float mu, float sd)
{
	return mu + sd*sqrt(-2*log(ran()))*cos(2*M_PI*ran());
}

/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gammasamp(double a, double b)
{
  if(a < 0 || b < 0) emsg("Part: EC2");

  if(a < 1){
    double u = ran();
    return gammasamp(1.0 + a, b) * pow (u, 1.0 / a);
  }
  else{
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while(1 == 1){
      do{
        x = sqrt(-2*log(ran()))*cos(2*M_PI*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      u = ran();

      if (u < 1 - 0.0331*x*x*x*x) break;

      if (log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
}


/// Splits up a string
vector<string> split(const string& s, char delimiter)                                                               
{                                 
   std::vector<std::string> splits;                       
   std::string split;                                      
   std::istringstream ss(s);                               
   while (std::getline(ss, split, delimiter)) splits.push_back(split);                                                  
   return splits;                                           
}


/// Displays an error message
void emsg(string msg)
{
	cout << msg << endl;
	MPI_Finalize();
	exit (EXIT_FAILURE);
}

/// Displays an error message on the root core
void emsgroot(string msg)
{
	int core;
	MPI_Comm_rank(MPI_COMM_WORLD,&core);
	if(core == 0) cout << msg << endl;
	MPI_Finalize();
	exit (EXIT_FAILURE);
}
