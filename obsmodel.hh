#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>


using namespace std;

#include "data.hh"
#include "poptree.hh"
#include "PART.hh"

vector <int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, int ti, int tf);
double Lobstot(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, double invT);
double Lobs(DATA &data, MODEL &model, POPTREE &poptree, int w, vector < vector <FEV> > &fev, double invT);

	