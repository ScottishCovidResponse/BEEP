#pragma once

using namespace std;

#include "model.hh"

void packinit();
int packsize();
double * packbuffer();

void pack(unsigned int num);
void pack(double num);
void pack(string &vec);
void pack(vector <unsigned int> &vec);
void pack(vector <int> &vec);
void pack(vector <double> &vec);
void pack(vector< vector <unsigned int> > &vec);
void pack(vector< vector <double> > &vec);
void pack(vector< vector< vector <double> > > &vec);
void pack(vector <string> &vec);
void pack(vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax);
void pack(vector <AREA> &vec);
void pack(vector <REGION> &vec);
void pack(vector <DEMOCAT> &vec);
void pack(vector <vector <EVREF> > &vec);

void unpack(unsigned int &num);
void unpack(double &num);
void unpack(string &vec);
void unpack(vector <unsigned int> &vec);
void unpack(vector <int> &vec);
void unpack(vector <double> &vec);
void unpack(vector< vector <unsigned int> > &vec);
void unpack(vector< vector <double> > &vec);
void unpack(vector< vector< vector <double> > > &vec);
void unpack(vector <string> &vec);
void unpack(vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax);
void unpack(vector <AREA> &vec);
void unpack(vector <REGION> &vec);
void unpack(vector <DEMOCAT> &vec);
void unpack(vector <vector <EVREF> > &vec);

	
	