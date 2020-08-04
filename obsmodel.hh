#ifndef BEEPMBP__OBSMODEL_HH
#define BEEPMBP__OBSMODEL_HH

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

vector <unsigned int> getnumtrans(DATA &data, const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v);

double Lobs(DATA &data, MODEL &model, const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev);

MEAS getmeas(DATA &data, MODEL &model, const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev);
#endif
