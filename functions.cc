
#include <iostream>
#include <string>
#include <vector>
#include "stdlib.h"
#include "math.h"

#include "functions.hh"
#include "consts.hh"

using namespace std;

/// Draws a random number between 0 and 1
double ran(){
	if(RAND_MAX == 32767){
		double v = (double(rand())*32766.0+rand())/(32767.0*RAND_MAX); if(v == 0 || v == 1) return 0.1;
		else return v;
	}
	else return double(0.999999999*rand())/RAND_MAX;
}

/// Displays any error messages
void emsg(string msg){ cout << msg << endl; exit (EXIT_FAILURE);}
