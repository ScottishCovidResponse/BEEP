#pragma once

using namespace std;

#include "PART.hh"

void packinit();
int packsize();
double * packbuffer();

void pack(int num);
void pack(double num);
void pack(vector <int> &vec);
void pack(vector <double> &vec);
void pack(vector< vector <int> > &vec);
void pack(vector< vector <double> > &vec);
void pack(vector <string> &vec);
void pack(vector< vector <FEV> > &vec, int fedivmin, int fedivmax);
void pack(vector <HOUSE> &vec, int min, int max);

void unpack(int &num);
void unpack(double &num);
void unpack(vector <int> &vec);
void unpack(vector <double> &vec);
void unpack(vector< vector <int> > &vec);
void unpack(vector< vector <double> > &vec);
void unpack(vector <string> &vec);
void unpack(vector< vector <FEV> > &vec, int fedivmin, int fedivmax);
void unpack(vector <HOUSE> &vec, int min, int max);
	
	