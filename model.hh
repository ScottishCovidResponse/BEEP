#pragma once

class MODEL
{
public:
	void definemodel();
	void betaspline();

	double settime[nsettime];
	double beta[nsettime];
	short nspline;                             // The spline points which are parameters in the model
	vector <double> splinet;
	vector <PARAM> param;
	vector <TRANS> trans;
	vector <COMP> comp;	

private:
	void addcomp(string name, double infectivity);
	void addparam(string name, double val, double min, double max);
	void addtrans(string from, string to, short type, string param1, string param2);
};
