// Outputs various graphs and statistics

#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>

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

void outputinit(MODEL &model)
{
	long r, p;
	
	trace.open("trace.txt");		
	trace << "state";
	for(p = 0; p < model.param.size(); p++) trace << "\t" << model.param[p].name; 
	trace << "\tLi"; 
	trace << endl;
	
	ensuredirectory("Output");
}

SAMPLE outputsamp(long samp, double Li, DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev)
{
	SAMPLE sa;
	long p, np, r, week;
	vector <long> num;
	
	np = model.param.size();
	trace << samp; 
	for(p = 0; p < np; p++) trace << "\t" << model.param[p].val; 
	trace << "\t" << Li; 
	trace << endl;
	
	sa.paramval.resize(np);
	for(p = 0; p < np; p++) sa.paramval[p] = model.param[p].val;
	
	sa.ncase.resize(data.nregion); for(r = 0; r < data.nregion; r++) sa.ncase[r].resize(data.nweek);
	
	for(week = 0; week < data.nweek; week++){
		num = getnumtrans(data,model,poptree,fev,"I","H",week*7,week*7+7);
		for(r = 0; r < data.nregion; r++) sa.ncase[r][week] = num[r];
	}
	
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
		
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp)
{      
	long p, r, s, nopsamp, week;
	vector <double> vec;
	STAT stat;
	
	nopsamp = opsamp.size();
	
	for(r = 0; r < data.nregion; r++){
		stringstream ss; ss << "Output/" << data.regionname[r] << "_data.txt";
		ofstream dataout(ss.str().c_str());
		for(week = 0; week < data.nweek; week++){
			vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].ncase[r][week]);
			stat = getstat(vec);
			
			dataout << week << " " <<data.ncase[r][week] << " " 
			        << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << "\n"; 
		}
	}
	
	stringstream ss; ss << "Output/params.txt";
	ofstream paramout(ss.str().c_str());
	for(p = 0; p < model.param.size(); p++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].paramval[p]);
		stat = getstat(vec);
			
		paramout << model.param[p].name  <<" " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << "\n"; 
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
	
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev)
{
	long week, r, h;
	vector <long> num;
	
	num = getnumtrans(data,model,poptree,fev,"I","H",0,data.tmax);
	cout << endl << "Total number of hospitalised cases:" << endl;
	for(r = 0; r < data.nregion; r++) cout <<	"Region " <<  r << ": " << num[r] << endl;
	cout << endl;
	
	stringstream sst; sst << "cases_" << type << ".txt";
	ofstream regplot(sst.str().c_str());
	
	regplot << "Week"; for(r = 0; r < data.nregion; r++){ regplot << "\t" << data.regionname[r];} regplot << endl;
	for(week = 0; week < data.nweek; week++){
		regplot << week;
		num = getnumtrans(data,model,poptree,fev,"I","H",week*7,week*7+7);
		for(r = 0; r < data.nregion; r++){ regplot <<  "\t" << num[r];} regplot << endl;
	}
	
	stringstream ssw; ssw << "houses_" << type << ".txt";
	ofstream houseout(ssw.str().c_str());
	houseout << data.popsize << " " << data.nhouse << "\n";
	for(h = 0; h < data.nhouse; h++){
		houseout << data.house[h].x << " " << data.house[h].y << " " << data.house[h].region << "\n";
	}
}

STAT getstat(vector <double> &vec)                                          // Calculates diagnostic statistics
{
	long n, i, d;
	double sum, a, cor, f;
	STAT stat;
	
	n = vec.size();
	
	sum = 0; for(i = 0; i < n; i++) sum += vec[i];
	stat.mean = sum/n;
	
	sum = 1;
	for(d = 0; d < n/2; d++){             // calculates the effective sample size
		a = 0; for(i = 0; i < n-d; i++) a += vec[i]*vec[i+d]; 
		cor = a/(n-d);
		if(cor < 0) break;
		sum += 0.5*cor;			
	}
	stat.ESS = n/sum;
	
	sort(vec.begin(),vec.end());
			
	i = long((n-1)*0.05); f = (n-1)*0.05 - i;
	stat.CImin = vec[i]*(1-f) + vec[i+1]*f;
			
	i = long((n-1)*0.95); f = (n-1)*0.95 - i;
	stat.CImax = vec[i]*(1-f) + vec[i+1]*f;
			
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
