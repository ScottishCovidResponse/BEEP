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

vector <unsigned int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, unsigned int ti, unsigned int tf);
double Lobstot(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, double invT);
double Lobs(DATA &data, MODEL &model, POPTREE &poptree, unsigned int w, vector < vector <FEV> > &fev, double invT);

double Lobs_mbp(DATA &data, MODEL &model, POPTREE &poptree,  vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev);
vector <unsigned int> getnumtrans_mbp(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string from, string to, unsigned int ti, unsigned int tf);
