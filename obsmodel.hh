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

vector <unsigned int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v);

double Lobs(DATA &data, MODEL &model, POPTREE &poptree,  vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev);

MEAS getmeas(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev);