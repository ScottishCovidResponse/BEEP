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
#include "obsmodel.hh"
#include "consts.hh"

struct STAT{                                           // Stores statistical information
	double mean;                                         // The mean
	double CImin, CImax;                                 // The minimum and maximum of the 90% credible interval
	double ESS;                                          // The estimated sample size
};

static void ensuredirectory(const string &path);
STAT getstat(vector <double> &vec);

ofstream trace, traceLi;

/// Initialises trace plot for parameters
void outputinit(DATA &data, MODEL &model)
{
	unsigned int p;
	
	ensuredirectory(data.outputdir);
	stringstream ss; ss << data.outputdir << "/trace.txt";

	trace.open(ss.str().c_str());		
	trace << "state";
	for(p = 0; p < model.param.size(); p++) trace << "\t" << model.param[p].name; 
	trace << "\tLi"; 
	trace << "\tinvT"; 
	trace << endl;
}

/// Initialises trace plot for likelihoods on difference chains (MBP only)
void outputLiinit(DATA &data, unsigned int nchaintot)
{
	unsigned int p;
	
	ensuredirectory(data.outputdir);
	stringstream ss; ss << data.outputdir << "/traceLi.txt";

	traceLi.open(ss.str().c_str());		
	traceLi << "state";
	for(p = 0; p < nchaintot; p++) traceLi << "\tchain" <<  p; 
	traceLi << endl;
}

/// Outputs trace plot for likelihoods on difference chains (MBP only)
void outputLi(unsigned int samp, unsigned int nparttot, double *Litot)
{
	unsigned int p;
	
	traceLi << samp;
	for(p = 0; p < nparttot; p++) traceLi << "\t" << Litot[p]; 
	traceLi << endl;
}

/// Outputs trace plot for parameters and store state data for plotting later
SAMPLE outputsamp(double invT, unsigned int samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector <double> &paramval, vector < vector <FEV> > &fev)
{
	SAMPLE sa;
	unsigned int p, np, r, t, st, sum, td;
	vector <unsigned int> num;
	double tinfav=UNSET, probA, Ainf, Pinf, tA, tP, tI;
	
	np = paramval.size();
	
	trace << samp; 
	for(p = 0; p < np; p++) trace << "\t" << paramval[p]; 
	trace << "\t" << Li; 
	trace << "\t" << invT; 
	trace << endl;
	
	sa.paramval = paramval;
	
	sa.transnum.resize(data.transdata.size());
	for(td = 0; td < data.transdata.size(); td++){
		if(data.transdata[td].type == "reg"){
			sa.transnum[td].resize(data.nregion); for(r = 0; r < data.nregion; r++) sa.transnum[td][r].resize(data.period);
			
			for(t = 0; t < data.period; t++){
				num = getnumtrans(data,model,poptree,fev,data.transdata[td].from,data.transdata[td].to,t,t+1);
				for(r = 0; r < data.nregion; r++) sa.transnum[td][r][t] = num[r];
			}
		}
		
		if(data.transdata[td].type == "all"){
			sa.transnum[td].resize(1); for(r = 0; r < 1; r++) sa.transnum[td][r].resize(data.period);
			
			for(t = 0; t < data.period; t++){
				num = getnumtrans(data,model,poptree,fev,data.transdata[td].from,data.transdata[td].to,t,t+1);
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
				sa.transnum[td][0][t] = sum;
			}
		}
	}

	model.paramval = paramval; model.betaspline(data);

	switch(model.modelsel){
	case MOD_IRISH:
		probA = model.getparam("probA"); Ainf = model.getparam("Ainf"); Pinf = model.getparam("Pinf");
		tA = 1.0/model.getparam("rA"); tP = 1.0/model.getparam("rP"); tI = 1.0/model.getparam("rI");

		tinfav = probA*Ainf*tA + (1-probA)*(Pinf*tP+tI); 
		break;
		
	default: emsg("Output: EC10"); break;
	}

	sa.R0.resize(data.nsettime);
	for(st = 0; st < data.nsettime; st++) sa.R0[st] = model.beta[st]*tinfav;
	
	return sa;
}

/// Outputs trace plot for parameters and store state data for plotting later
SAMPLE outputsamp_mbp(double invT, unsigned int samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector <double> &paramval, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev)
{
	SAMPLE sa;
	unsigned int p, np, r, t, st, sum, td;
	vector <unsigned int> num;
	double tinfav=UNSET, probA, Ainf, Pinf, tA, tP, tI;
	
	np = paramval.size();
	
	trace << samp; 
	for(p = 0; p < np; p++) trace << "\t" << paramval[p]; 
	trace << "\t" << Li; 
	trace << "\t" << invT; 
	trace << endl;
	
	sa.paramval = paramval;
	
	sa.transnum.resize(data.transdata.size());
	for(td = 0; td < data.transdata.size(); td++){
		if(data.transdata[td].type == "reg"){
			sa.transnum[td].resize(data.nregion); for(r = 0; r < data.nregion; r++) sa.transnum[td][r].resize(data.period);
			
			for(t = 0; t < data.period; t++){
				num = getnumtrans_mbp(data,model,poptree,trev,indev,data.transdata[td].from,data.transdata[td].to,t,t+1);
				for(r = 0; r < data.nregion; r++) sa.transnum[td][r][t] = num[r];
			}
		}
		
		if(data.transdata[td].type == "all"){
			sa.transnum[td].resize(1); for(r = 0; r < 1; r++) sa.transnum[td][r].resize(data.period);
			
			for(t = 0; t < data.period; t++){
				num = getnumtrans_mbp(data,model,poptree,trev,indev,data.transdata[td].from,data.transdata[td].to,t,t+1);
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
				sa.transnum[td][0][t] = sum;
			}
		}
	}

	model.paramval = paramval; model.betaspline(data);

	switch(model.modelsel){
	case MOD_IRISH:
		probA = model.getparam("probA"); Ainf = model.getparam("Ainf"); Pinf = model.getparam("Pinf");
		tA = 1.0/model.getparam("rA"); tP = 1.0/model.getparam("rP"); tI = 1.0/model.getparam("rI");

		tinfav = probA*Ainf*tA + (1-probA)*(Pinf*tP+tI); 
		break;
		
	default: emsg("Output: EC10"); break;
	}

	sa.R0.resize(data.nsettime);
	for(st = 0; st < data.nsettime; st++) sa.R0[st] = model.beta[st]*tinfav;
	
	return sa;
}

/// Generates posterior plots for transitions, variation in R0 over time, parameter statistics and MCMC diagnostics 
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp)
{      
	unsigned int p, r, s, st, nopsamp, t, td, j, jmax;
	vector <double> vec;
	STAT stat;
	string name;
	
	ensuredirectory(data.outputdir);
		
	nopsamp = opsamp.size();
	
	cout << endl;
	if(data.mode == MODE_SIM) cout << "Outputs:" << endl;
	else cout << "Posterior Outputs:" << endl;
	
	for(td = 0; td < data.transdata.size(); td++){
		name = data.transdata[td].file;
	
		j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
		name = name.substr(0,j);
	
		if(data.transdata[td].type == "reg"){
			for(r = 0; r < data.nregion; r++){
				stringstream ss; ss << data.outputdir << "/" << name << "_" << data.region[r].name << ".txt";
				ofstream dataout(ss.str().c_str());
				
				cout << "'" << ss.str() << "' gives numbers of " << data.transdata[td].from << "→" << data.transdata[td].to << " transitions for region '" << data.region[r].name << "'." << endl;
		
				for(t = 0; t < data.period; t++){
					vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].transnum[td][r][t]);
					stat = getstat(vec);
					
					dataout << t+0.5 << " ";
					if(data.mode != MODE_SIM) dataout << data.transdata[td].num[r][t]; else dataout << stat.mean;
					dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
				}
			}
		}
		
		if(data.transdata[td].type == "all"){
			stringstream ss; ss << data.outputdir << "/" << name << ".txt";
			ofstream dataout(ss.str().c_str());
			
			cout << "'" << ss.str() << "' gives numbers of " << data.transdata[td].from << "→" << data.transdata[td].to << " transitions." << endl;
	
			for(t = 0; t < data.period; t++){
				vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].transnum[td][0][t]);
				stat = getstat(vec);
				
				dataout << t+0.5 << " ";
				if(data.mode != MODE_SIM) dataout << data.transdata[td].num[0][t]; else dataout << stat.mean;
				dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
			}
		}
	}
	
	stringstream sst; sst << data.outputdir << "/R0" << ".txt";
	ofstream R0out(sst.str().c_str());
	
	cout << "'" << sst.str() << "' gives the time variation in R0." << endl;
	
	for(st = 0; st < data.nsettime; st++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		stat = getstat(vec);
		
		R0out << (st+0.5)*data.period/data.nsettime << " " 
		      << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}
	
	stringstream ss; ss << data.outputdir << "/parameters" << ".txt";
	ofstream paramout(ss.str().c_str());
	
	cout << "'" << ss.str() << "' gives the model parameters." << endl;
	
	for(p = 0; p < model.param.size(); p++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].paramval[p]);
		stat = getstat(vec);
			
		paramout << model.param[p].name  <<" " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << endl; 
	}
	
	if(data.mode != MODE_SIM){
		cout << "'" << data.outputdir << "/trace.txt' gives trace plots for model parameters." << endl;
	}
	
	if(data.mode == MODE_MBP){
		cout << "'" << data.outputdir << "/traceLi.txt' gives trace plots for the observation likelihoods on different chains." << endl;
	}

	// This gives the acceptance rates for different MCMC proposals on different parameters
	
	if(data.mode != MODE_SIM){
		stringstream ss; ss << data.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
	
		cout << "'" << ss.str() << "' gives MCMC diagnostics." << endl;
	
		diag << "MCMC diagnostics:" << endl;

		if(data.mode == MODE_PMCMC) diag << "Base acceptance rate " << double(model.nac)/model.ntr << endl;
		
		for(p = 0; p < model.param.size(); p++){
			diag << model.param[p].name << ": ";
			if(model.param[p].ntr == 0) diag << "Fixed" << endl;
			else diag << "Acceptance rate " << double(model.param[p].nac)/model.param[p].ntr << endl;
		}
	}
	cout << endl;
}
	
/// Outputs an event sample fev
void outputeventsample(vector < vector <FEV> > &fev, DATA &data, MODEL &model, POPTREE &poptree)
{
	unsigned int d, j, nind;
	vector< vector <FEV> > indev;
	TRANS tr;
	
	nind = data.ind.size();
	indev.resize(nind);
	for(d = 0; d < fev.size(); d++){
		for(j = 0; j < fev[d].size(); j++) indev[fev[d][j].ind].push_back(fev[d][j]);
	}
	
	ensuredirectory(data.outputdir);
	stringstream sst; sst << data.outputdir << "/events.txt";
	ofstream evsamp(sst.str().c_str());
	
	/*
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
	*/
}

/// Outputs a population plot for event sequence xi
void outputplot(string file, DATA &data, MODEL &model,  vector < vector <FEV> > &xi, double tmin, double period)
{
	unsigned int c, tra, td, tdf;
	double t;
	vector <int> N;
	TRANS tr;
	
	N.resize(model.comp.size()); for(c = 0; c < model.comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	td = 0; tdf = 0; while(td < data.fediv && xi[td].size()==0) td++;
	
	ofstream plot(file.c_str());
	for(t = tmin; t < period; t += (period-tmin)/100){
		while(td < data.fediv && xi[td][tdf].t < t){
			tra = xi[td][tdf].trans;
			tr = model.trans[tra];
			N[tr.from]--; N[tr.to]++;
			
			tdf++;
			if(tdf == xi[td].size()){
				td++; tdf = 0; 
				while(td < data.fediv && xi[td].size() == 0) td++;
			}
		}
		
		plot << t << " ";
		for(c = 0; c < model.comp.size(); c++) plot << N[c] << " ";
		plot << "\n";
	}
}

/// Generates case data based on a simulation
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev)
{
	unsigned int t, r, tot, td, sum;
	vector <unsigned int> num;
	
	for(td = 0; td < data.transdata.size(); td++){
		num = getnumtrans(data,model,poptree,fev,data.transdata[td].from,data.transdata[td].to,0,data.period);
		
		tot = 0;
		cout << endl << "The following number of " << data.transdata[td].from << "→" << data.transdata[td].to << " transitions were observed:" << endl;
		for(r = 0; r < data.nregion; r++){
			cout <<	data.region[r].name << ": " << num[r] << endl;
			tot += num[r];
		}
		cout << "Total: " << tot << endl; 
		cout << endl;
	}
	
	cout << "Simulated Data:" << endl;
	for(td = 0; td < data.transdata.size(); td++){
		ofstream transout(data.transdata[td].file.c_str());
		
		cout << "'" << data.transdata[td].file.c_str() << "' gives the observed weekly number of " << data.transdata[td].from << "→" << data.transdata[td].to << " transitions";
		if(data.transdata[td].type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		transout << "Week"; 
		for(r = 0; r < data.nregion; r++){ transout << "\t" << data.region[r].name;} transout << endl;
		for(t = 0; t < data.period; t++){
			transout << t;
			num = getnumtrans(data,model,poptree,fev,data.transdata[td].from,data.transdata[td].to,t,t+1);
			if(data.transdata[td].type == "reg"){				
				for(r = 0; r < data.nregion; r++){ transout <<  "\t" << num[r];} transout << endl;
			}
			
			if(data.transdata[td].type == "all"){
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r]; 
				transout <<  "\t" << sum << endl;
			}
		}
	}
	
	cout << endl;
}

/// Calculates diagnostic statistics
STAT getstat(vector <double> &vec)                           
{
	unsigned int n, i, d;
	double sum, sum2, sd, a, cor, f;
	STAT stat;
	
	n = vec.size();
	
	sum = 0; sum2 = 0; for(i = 0; i < n; i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
	sum /= n; sum2 /= n;
	stat.mean = sum; 
	
	sort(vec.begin(),vec.end());
			
	i = (unsigned int)((n-1)*0.05); f = (n-1)*0.05 - i;
	stat.CImin = vec[i]*(1-f) + vec[i+1]*f;
			
	i = (unsigned int)((n-1)*0.95); f = (n-1)*0.95 - i;
	stat.CImax = vec[i]*(1-f) + vec[i+1]*f;

	sd = sqrt(sum2 - sum*sum);
	if(sd == 0) stat.ESS = 0;
	else{
		for(i = 0; i < n; i++) vec[i] = (vec[i]-sum)/sd;
			
		sum = 1;
		for(d = 0; d < n/2; d++){             // calculates the effective sample size
			a = 0; for(i = 0; i < n-d; i++) a += vec[i]*vec[i+d]; 
			cor = a/(n-d);
			if(cor < 0) break;
			sum += 0.5*cor;			
		}
		stat.ESS = n/sum;
	}
		
	return stat;
}

/// Create a directory if it doesn't already exist
static void ensuredirectory(const string &path) 
{
	struct stat st = {0};
	if (stat(path.c_str(), &st) == -1)
	{
		// Directory not found
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (ret == -1)
			emsg("Error creating directory "+path);
	}
}
