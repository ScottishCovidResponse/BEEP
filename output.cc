// Outputs various graphs and statistics

#include <iostream>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

#include "utils.hh"
#include "model.hh"
#include "PART.hh"
#include "output.hh"

struct STAT{
	double mean;
	double CImin, CImax;
	double ESS;
};

void ensuredirectory(const string &path);
STAT getstat(vector <double> &vec);

ofstream trace;

/// Initialises trace plot for parameters
void outputinit(MODEL &model, long nparttot)
{
	long r, p;
	
	ensuredirectory("Output");
	
	stringstream ss; ss << "Output/trace_" << nparttot << ".txt";

	trace.open(ss.str().c_str());		
	trace << "state";
	for(p = 0; p < model.param.size(); p++) trace << "\t" << model.param[p].name; 
	trace << "\tLi"; 
	trace << "\tinvT"; 
	trace << endl;
}

/// Outputs trace plot form parameters and store state data for plotting later
SAMPLE outputsamp(double invT, long samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev)
{
	SAMPLE sa;
	long p, np, r, week, st;
	vector <long> num;
	double rAI, rAR, rIH, rIR, rID, tinfav;
	
	np = model.param.size();
	
	trace << samp; 
	for(p = 0; p < np; p++) trace << "\t" << model.param[p].val; 
	trace << "\t" << Li; 
	trace << "\t" << invT; 
	trace << endl;
	
	sa.paramval.resize(np);
	for(p = 0; p < np; p++) sa.paramval[p] = model.param[p].val;
	
	sa.ncase.resize(data.nregion); for(r = 0; r < data.nregion; r++) sa.ncase[r].resize(data.nweek);
	
	for(week = 0; week < data.nweek; week++){
		num = getnumtrans(data,model,poptree,fev,"I","H",week*timestep,(week+1)*timestep);
		for(r = 0; r < data.nregion; r++) sa.ncase[r][week] = num[r];
	}

	switch(modelsel){
	case MOD_IRISH:
		double afrac, aI, tIaR, tIH;
		afrac = model.param[model.afracparam].val;
		aI = model.param[model.aIparam].val;
		tIaR = 1.0/model.getrate("Ia","R"); tIH = 1.0/model.getrate("I","H"); 
		tinfav = afrac*aI*tIaR + (1-afrac)*tIH;
		break;
		
	case MOD_OLD:
		rAI = model.getrate("A","I"); rAR = model.getrate("A","R");
		rIH = model.getrate("I","H");	rIR = model.getrate("I","R");	rID = model.getrate("I","D");
		tinfav = 0.2/(rAI + rAR) + 1*(rAI/(rAI+rAR))*(1.0/(rIH + rIR + rID));
		break;
	}
	
	sa.R0.resize(nsettime);
	for(st = 0; st < nsettime; st++) sa.R0[st] = model.beta[st]*tinfav;
	
	return sa;
}

/// Returns the number of transitions for individuals going from compartment "from" to compartment "to" 
/// in different regions over the time range ti - tf
vector <long> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, short ti, short tf)
{
	long d, k, r, tra;
	FEV fe;
	vector <long> num;
	
	tra = 0; 
	while(tra < model.trans.size() && 
	     !(model.comp[model.trans[tra].from].name == from && model.comp[model.trans[tra].to].name == to)) tra++;
	if(tra == model.trans.size()) emsg("Finescale: Cannot find transition");
	
	for(r = 0; r < data.nregion; r++) num.push_back(0);

	if(fev.size() < fediv) emsg("fevsize");
	for(d = long(fediv*double(ti)/data.tmax); d <= long(fediv*double(tf)/data.tmax); d++){
		if(d < fediv){
			for(k = 0; k < fev[d].size(); k++){
				fe = fev[d][k];
				if(fe.t > tf) break;
				if(fe.t > ti && fe.trans == tra) num[data.house[poptree.ind[fe.ind].houseref].region]++;
			}
		}
	}
	
	return num;
}

/// Generates posterior plots for cases in different regions, variation in R0 over time, parameter statistices and MCMC diagnoistics 
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp, short siminf, long nparttot)
{      
	long p, r, s, st, nopsamp, week;
	vector <double> vec;
	STAT stat;
	
	nopsamp = opsamp.size();
	
	for(r = 0; r < data.nregion; r++){
		stringstream ss; ss << "Output/" << data.regionname[r] << "_data_" <<  nparttot << ".txt";
		ofstream dataout(ss.str().c_str());
		for(week = 0; week < data.nweek; week++){
			vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].ncase[r][week]);
			stat = getstat(vec);
			
			dataout << week+0.5 << " ";
			if(siminf == 0) dataout << data.ncase[r][week]; else dataout << stat.mean;
			dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
		}
	}
	
	stringstream sst; sst << "Output/R0_" << nparttot << ".txt";
	ofstream R0out(sst.str().c_str());
	for(st = 0; st < nsettime; st++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		stat = getstat(vec);
		
		R0out << (st+0.5)*data.tmax/nsettime << " " 
		      << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}
	
	stringstream ss; ss << "Output/params_" << nparttot << ".txt";
	ofstream paramout(ss.str().c_str());
	for(p = 0; p < model.param.size(); p++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].paramval[p]);
		stat = getstat(vec);
			
		paramout << model.param[p].name  <<" " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << endl; 
	}
	
	// This gives the acceptance rates for different MCMC proposals on different parameters
	cout << "MCMC diagnostics:" << endl;
	cout << "Base acceptance rate " << double(model.nac)/model.ntr << endl;
	for(p = 0; p < model.param.size(); p++){
		cout << model.param[p].name << ": ";
		if(model.param[p].ntr == 0) cout << "Fixed" << endl;
		else cout << "Acceptance rate " << double(model.param[p].nac)/model.param[p].ntr << endl;
	}
}
	
/// Outputs a given siumlation event sample
void outputeventsample(vector < vector <FEV> > &fev, DATA &data, MODEL &model, POPTREE &poptree)
{
	long d, i, j, e, h, nind;
	vector< vector <FEV> > indev;
	TRANS tr;
	
	nind = poptree.ind.size();
	indev.resize(nind);
	for(d = 0; d < fev.size(); d++){
		for(j = 0; j < fev[d].size(); j++) indev[fev[d][j].ind].push_back(fev[d][j]);
	}
	
	stringstream sst; sst << "events_" << type << ".txt";
	ofstream evsamp(sst.str().c_str());
	
	for(i = 0; i < nind; i++){
		if(indev[i].size() > 0){
			h = poptree.ind[i].houseref;
			evsamp << i << "\t" << data.house[h].x << "\t" << data.house[h].y << "\t" << indev[i].size() << "\t";
			for(e = 0; e < indev[i].size(); e++){
				tr = model.trans[indev[i][e].trans];
				if(e == 0) evsamp << model.comp[tr.from].name << "\t";
				evsamp << indev[i][e].t << "\t" << model.comp[tr.to].name << "\t";
			}
			evsamp << endl;
		}
	}
}

/// Generates case data based on a simulation
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev)
{
	long week, r, h, tot;
	vector <long> num;
	
	num = getnumtrans(data,model,poptree,fev,"I","H",0,data.tmax);
	tot = 0;
	cout << endl << "Total number of hospitalised cases:" << endl;
	for(r = 0; r < data.nregion; r++){
		cout <<	"Region " <<  r << ": " << num[r] << endl;
		tot += num[r];
	}
	cout << "Total: " << tot << endl; 
	cout << endl;
	
	stringstream sst; sst << "cases_" << type << ".txt";
	ofstream regplot(sst.str().c_str());
	
	regplot << "Week"; 
	for(r = 0; r < data.nregion; r++){ regplot << "\t" << data.regionname[r];} regplot << endl;
	for(week = 0; week < data.nweek; week++){
		regplot << week;
		num = getnumtrans(data,model,poptree,fev,"I","H",week*timestep,(week+1)*timestep);
		for(r = 0; r < data.nregion; r++){ regplot <<  "\t" << num[r];} regplot << endl;
	}
	
	if(type.substr(0,4) != "real"){
		stringstream ssw; ssw << "houses_" << type << ".txt";
		ofstream houseout(ssw.str().c_str());
		houseout << data.popsize << " " << data.nhouse << " "  << data.nregion << " " <<  endl;
		for(r = 0; r < data.nregion; r++) houseout << data.regionname[r] << endl;
		for(h = 0; h < data.nhouse; h++){
			houseout << data.house[h].x << " " << data.house[h].y << " " << data.house[h].region << endl;
		}
	}
}

/// Calculates diagnostic statistics
STAT getstat(vector <double> &vec)                           
{
	long n, i, d;
	double sum, sum2, sd, a, cor, f;
	STAT stat;
	
	n = vec.size();
	
	sum = 0; sum2 = 0; for(i = 0; i < n; i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
	sum /= n; sum2 /= n;
	stat.mean = sum; 
	
	sort(vec.begin(),vec.end());
			
	i = long((n-1)*0.05); f = (n-1)*0.05 - i;
	stat.CImin = vec[i]*(1-f) + vec[i+1]*f;
			
	i = long((n-1)*0.95); f = (n-1)*0.95 - i;
	stat.CImax = vec[i]*(1-f) + vec[i+1]*f;

	sd = sqrt(sum2 - sum*sum);
	for(i = 0; i < n; i++) vec[i] = (vec[i]-sum)/sd;
		
	sum = 1;
	for(d = 0; d < n/2; d++){             // calculates the effective sample size
		a = 0; for(i = 0; i < n-d; i++) a += vec[i]*vec[i+d]; 
		cor = a/(n-d);
		if(cor < 0) break;
		sum += 0.5*cor;			
	}
	stat.ESS = n/sum;
	
	return stat;
}

// Create a directory if it doesn't already exist
void ensuredirectory(const string &path) {
	struct stat st = {0};
	if (stat(path.c_str(), &st) == -1)
	{
		// Directory not found
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (ret == -1)
			emsg("Error creating directory "+path);
	}
}
