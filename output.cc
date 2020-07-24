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
#include "output.hh"
#include "obsmodel.hh"
#include "consts.hh"
#include "data.hh"

struct STAT{                                           // Stores statistical information
	string mean;                                         // The mean
	string CImin, CImax;                                 // The minimum and maximum of the 90% credible interval
	string ESS;                                          // The estimated sample size
};

STAT getstat(vector <double> &vec);

ofstream trace, traceLi;

/// Initialises trace plot for parameters
void outputinit(DATA &data, MODEL &model)
{
	unsigned int p, pc;
	
	ensuredirectory(data.outputdir);
	stringstream ss; ss << data.outputdir << "/trace.txt";

	trace.open(ss.str().c_str());		
	trace << "state";
	for(p = 0; p < model.param.size(); p++) trace << "\t" << model.param[p].name; 
	for(pc = 0; pc < model.priorcomps.size(); pc++) trace << "\tProb in " << model.comp[model.priorcomps[pc].comp].name;
	trace << "\tLi"; 
	trace << "\tPri"; 
	trace << "\tinvT"; 
	trace << "\tninf";
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
void outputLi(unsigned int samp, unsigned int nchaintot, double *Litot)
{
	unsigned int p;
	
	traceLi << samp;
	for(p = 0; p < nchaintot; p++) traceLi << "\t" << Litot[p]; 
	traceLi << endl;
}

/// Outputs trace plot for parameters and store state data for plotting later

SAMPLE outputsamp(double invT, unsigned int samp, double Li, double Pri, DATA &data, MODEL &model, POPTREE &poptree, vector <double> &paramval, unsigned int ninf, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev)
{
	SAMPLE sa;
	unsigned int p, np, pc, c, a;
	double pr;
	vector <unsigned int> num, numtot;
	
	np = paramval.size();

	trace << samp; 
	for(p = 0; p < np; p++) trace << "\t" << paramval[p]; 
	
	if(model.priorcomps.size() > 0){
		model.calcprobin();
		for(pc = 0; pc < model.priorcomps.size(); pc++){
			c = model.priorcomps[pc].comp;
			pr = 0; for(a = 0; a < data.nage; a++) pr += data.agedist[a]*model.comp[c].probin[a];
			trace << "\t" << pr;
		}
	}

	trace << "\t" << Li; 
	trace << "\t" << Pri; 
	trace << "\t" << invT; 
	trace << "\t" << ninf; 
	trace << endl;
	
	sa.paramval = paramval;
	
	sa.meas = getmeas(data,model,poptree,trev,indev);

	sa.R0 = model.R0calc();
	
	return sa;
}


/// Outputs a posterior graph
void outputplot(DATA &data, vector <SAMPLE> &opsamp, unsigned int d, unsigned int r, unsigned int type)
{
	unsigned int j, jmax, row, s, t=0, nopsamp, opsampmin, rr, rrmin, rrmax, nrow=0, dc=0, jj;
	double sum, valsum;
	string name, file, filefull;
	vector <double> vec;
	STAT stat;
	
	nopsamp = opsamp.size();
	opsampmin = nopsamp/4;
	
	switch(type){
	case POP_DATA: name = data.popdata[d].file; break;
	case TRANS_DATA: name = data.transdata[d].file; break;
	case MARG_DATA: name = data.margdata[d].file; break;
	default: emsg("Output: EC10"); break;
	}
	
	j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
	name = name.substr(0,j);
	
	switch(r){
	case UNSET:	file = "Posterior_"+name+".txt"; break;
	case ADD:	file = "Posterior_"+name+"_sum.txt"; break;
	default: file = "Posterior_"+name+"_"+data.region[r].code+".txt"; break;
	}
	
	filefull = data.outputdir+"/"+file;
	ofstream dataout(filefull.c_str());
	if(!dataout) emsg("Cannot output the file '"+filefull+"'");
						
	switch(type){
	case POP_DATA: 
		cout << "'" << file << "' gives the population in '" << data.popdata[d].compstr << "'";
		dataout << "# The population in '" << data.popdata[d].compstr << "'"; 
		nrow = data.popdata[d].rows;
		break;
		
	case TRANS_DATA:
		cout << "'" << file << "' gives numbers of " << data.transdata[d].fromstr << "→" << data.transdata[d].tostr << " transitions";
		dataout << "# Population in " << data.transdata[d].fromstr << "→" << data.transdata[d].tostr << " transitions";
		nrow = data.transdata[d].rows;
		break;
	
	case MARG_DATA:
		dc = data.margdata[d].democat;
		nrow = data.democat[dc].value.size();
		cout << "'" << file << "' gives " << data.margdata[d].fromstr << "→" << data.margdata[d].tostr << " transitions stratified by '" << data.democat[dc].name << "'";
		dataout << "# The " << data.margdata[d].fromstr << "→" << data.margdata[d].tostr << " transitions stratified by '" << data.democat[dc].name << "'";
		break;
	}	
	
	if(r != UNSET && r != ADD){
		cout << " for region '" << data.region[r].name << "'." << endl;
		dataout << " for region '" << data.region[r].name << "'." << endl;
	}
	else{
		cout << "." << endl;
		dataout << "." << endl;
	}
		
	switch(r){
	case UNSET: rrmin = 0; rrmax = 1; break;
	case ADD: rrmin = 0; rrmax = data.nregion; if(type == MARG_DATA) emsg("Output: EC64"); break;
	default: rrmin = r; rrmax = r+1; break;
	}

	switch(type){
	case TRANS_DATA: case POP_DATA: dataout << "# Time from start, " << data.tformat; break;
	case MARG_DATA: dataout << "category"; break;
	}	
	dataout << ", data, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;
	for(row = 0; row < nrow; row++){
		vec.clear(); 
		for(s = opsampmin; s < nopsamp; s++){
			valsum = 0;
			for(rr = rrmin; rr < rrmax; rr++){
				switch(type){
				case POP_DATA: valsum += opsamp[s].meas.popnum[d][row][rr]; break;
				case TRANS_DATA: valsum += opsamp[s].meas.transnum[d][row][rr]; break;
				case MARG_DATA:
					sum = 0; for(jj = 0; jj < data.democat[dc].value.size(); jj++) sum += opsamp[s].meas.margnum[d][jj][rr];
					valsum += 100.0*opsamp[s].meas.margnum[d][row][rr]/sum;
					break;
				}	
			}
			vec.push_back(valsum);
		}
		
		stat = getstat(vec);
		switch(type){
		case POP_DATA: t = data.popdata[d].start+row*data.popdata[d].units; break;
		case TRANS_DATA: t = data.transdata[d].start+row*data.transdata[d].units; break;
		}
		if(type == MARG_DATA) dataout << data.democat[dc].value[row] << " ";
		else dataout << t << " " << data.getdate(t) << " ";
		
		if(data.mode != MODE_SIM){
			valsum = 0;
			for(rr = rrmin; rr < rrmax; rr++){
				if(type == MARG_DATA) valsum += data.margdata[d].percent[rr][row]; 
				else{
					switch(type){
					case POP_DATA: jj = data.popdata[d].num[rr][row]; break;
					case TRANS_DATA: jj = data.transdata[d].num[rr][row]; break;
					default: jj = 0; break;
					}
					switch(jj){
					case UNKNOWN: break;
					case THRESH: valsum += data.threshold; break;
					default: valsum += jj; break;
					}
				}
			}
			dataout << valsum;
		}
		else dataout << stat.mean;
		dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}
}

/// Generates posterior plots for transitions, variation in R0 over time, parameter statistics and MCMC diagnostics 
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp)
{      
	unsigned int p, r, s, st, nopsamp, opsampmin, d;
	string file, filefull;
	vector <double> vec;
	STAT stat;
	string name;

	ensuredirectory(data.outputdir);
		
	nopsamp = opsamp.size();
	opsampmin = nopsamp/4;
	
	cout << endl;
	if(data.mode == MODE_SIM) cout << "Outputs in directory '" << data.outputdir << "':" << endl;
	else cout << "Posterior outputs in directory '" << data.outputdir << "':" << endl;
	
	for(d = 0; d < data.transdata.size(); d++){
		if(data.transdata[d].type == "reg"){
			for(r = 0; r < data.nregion; r++) outputplot(data,opsamp,d,r,TRANS_DATA);
			outputplot(data,opsamp,d,ADD,TRANS_DATA);
		}
		if(data.transdata[d].type == "all") outputplot(data,opsamp,d,UNSET,TRANS_DATA);
	}
	
	for(d = 0; d < data.popdata.size(); d++){
		if(data.popdata[d].type == "reg"){
			for(r = 0; r < data.nregion; r++) outputplot(data,opsamp,d,r,POP_DATA);
			outputplot(data,opsamp,d,ADD,POP_DATA);
		}
		if(data.popdata[d].type == "all") outputplot(data,opsamp,d,UNSET,POP_DATA);
	}
	
	for(d = 0; d < data.margdata.size(); d++){
		if(data.margdata[d].type == "reg"){ for(r = 0; r < data.nregion; r++) outputplot(data,opsamp,d,r,MARG_DATA);}
		if(data.margdata[d].type == "all") outputplot(data,opsamp,d,UNSET,MARG_DATA);
	}
		
	file = "Posterior_R0.txt";
	filefull = data.outputdir+"/"+file;
	ofstream R0out(filefull.c_str());
	if(!R0out) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time variation in R0." << endl;
	R0out << "# Gives the time variation in R0." << endl;	
	R0out << "# Time from start, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	for(st = 0; st < data.nsettime; st++){
		vec.clear(); for(s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		stat = getstat(vec);
		
		R0out << (st+0.5)*data.period/data.nsettime << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}
	
	file = "Posterior_parameters.txt";
	filefull = data.outputdir+"/"+file;
	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the model parameters." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200." << endl;
	paramout << endl;
	
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	for(p = 0; p < model.param.size()-1; p++){
		vec.clear(); for(s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].paramval[p]);
		stat = getstat(vec);
			
		paramout << model.param[p].name  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << endl; 
	}
	
	if(data.mode != MODE_SIM){
		cout << "'" << data.outputdir << "/trace.txt' gives trace plots for model parameters." << endl;
	}
	
	if(data.mode == MODE_INF){
		cout << "'" << data.outputdir << "/traceLi.txt' gives trace plots for the observation likelihoods on different chains." << endl;
	}
}

void outputcombinedtrace(vector <string> &paramname, vector < vector < vector <double> > > &vals, string file)
{
	unsigned int p, inp, s, nopsamp, opsampmin, N, M, nmax, nmin;
	double W, B, valav, varr, muav;
	string GR;
	vector <double> vec, mu, vari;
	STAT stat;
		
	ofstream paramout(file.c_str());
	if(!paramout) emsg("Cannot output the file '"+file+"'");
	
	cout << "'" << file << "' gives the model parameters combining the input trace files." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200 and PSRF should be between 0.9 and 1.1." << endl;
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, effective sample size (ESS), potential scale reduction factor (PSRF) otherwise known as the Gelman-Rubin convergence diagnostic." << endl;
	paramout << endl;

	M = vals.size();
	nmax = large;
	for(inp = 0; inp < M; inp++) if(vals[inp][0].size() < nmax) nmax = vals[inp][0].size();
	nmin = nmax/4;
	N = nmax-nmin; 
		
	mu.resize(M); vari.resize(M);
	for(p = 0; p < paramname.size(); p++){
		GR = "---";
		if(M > 1){
			muav = 0;
			for(inp = 0; inp < M; inp++){
				valav = 0; for(s = nmin; s < nmax; s++) valav += vals[inp][p][s]/N;
				varr = 0; for(s = nmin; s < nmax; s++) varr += (vals[inp][p][s]-valav)*(vals[inp][p][s]-valav)/(N-1);
				mu[inp] = valav;
				vari[inp] = varr;
				muav += mu[inp]/M;
			}
					
			W = 0; for(inp = 0; inp < M; inp++) W += vari[inp]/M;
			B = 0; for(inp = 0; inp < M; inp++) B += (mu[inp]-muav)*(mu[inp]-muav)*N/(M-1);
			if(W > tiny) GR = to_string(sqrt(((1-1.0/N)*W+B/N)/W));
		}
	
		vec.clear(); 
		for(inp = 0; inp < M; inp++){
			nopsamp = vals[inp][p].size();
			opsampmin = nopsamp/4;
			for(s = opsampmin; s < nopsamp; s++) vec.push_back(vals[inp][p][s]);
		}
		stat = getstat(vec);
			
		paramout << paramname[p]  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << " " << GR << endl; 
	}
}
	
/// Outputs an event sample fev
void outputeventsample(vector < vector <FEV> > &fev, DATA &data, MODEL & /*model*/, POPTREE & /*poptree*/)
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
	if(!evsamp) emsg("Cannot output the file '"+sst.str()+"'");
	
	/*
	for(i = 0; i < nind; i++){
		if(indev[i].size() > 0){
			h = poptree.ind[i].houseref;
			evsamp << i << "\t" << data.house[h].x << "\t" << data.house[h].y << "\t" << indev[i].size() << "\t";
			for(e = 0; e < indev[i].size(); e++){
				tr = model.trans[indev[i][e].trans];
				if(e == 0) evsamp << model.comp[tr.from].name << "\t";
				evsamp << indev[i][e].t << "\t" << model.comp[tr.tostr].name << "\t";
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
	if(!plot) emsg("Cannot output the file '"+file+"'");
	
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
		plot << endl;
	}
}

/// Generates case data based on a simulation using the MBP algorithm
void outputsimulateddata(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string dir)
{
	unsigned int row, r, td, pd, md, d, j, jj;
	string file, filefull;
	double sum;
	MEAS meas;
	
	ensuredirectory(dir);
		
	meas = getmeas(data,model,poptree,trev,indev);
	
	cout << "Simulated data in directory '" << dir <<"':" << endl;
	for(td = 0; td < data.transdata.size(); td++){
		file = data.transdata[td].file;
		filefull =  dir+"/"+file;
		ofstream transout(filefull);
		if(!transout) emsg("Cannot output the file '"+filefull+"'");
		
		cout << "  '" << file << "' gives the observed weekly number of " << data.transdata[td].fromstr << "→" << data.transdata[td].tostr << " transitions";
		if(data.transdata[td].type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		switch(data.tform){
		case TFORM_NUM: transout << "time"; break;
		case TFORM_YMD: transout << "date"; break;
		}

		if(data.transdata[td].type == "reg"){ for(r = 0; r < data.nregion; r++){ transout << "\t" << data.region[r].code;}}
		else transout << "\t" << "all";
		transout << endl;
		
		for(row = 0; row < data.transdata[td].rows; row++){
			transout << data.getdate(data.transdata[td].start + row*data.transdata[td].units);
			for(r = 0; r < meas.transnum[td][row].size(); r++){ transout <<  "\t" << meas.transnum[td][row][r];} transout << endl;
		}
	}
	
	for(pd = 0; pd < data.popdata.size(); pd++){
		file = data.popdata[pd].file;
		filefull = dir+"/"+file;
		ofstream popout(filefull);
		if(!popout) emsg("Cannot output the file '"+filefull+"'");
		
		cout << "  '" << file << "' gives the numbers in population '" << data.popdata[pd].compstr << "'";
		if(data.popdata[pd].type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		switch(data.tform){
		case TFORM_NUM: popout << "time"; break;
		case TFORM_YMD: popout << "date"; break;
		}
		
		if(data.popdata[pd].type == "reg"){ for(r = 0; r < data.nregion; r++){ popout << "\t" << data.region[r].code;}}
		else popout << "\t" << "all";
		popout << endl;
		
		for(row = 0; row < data.popdata[pd].rows; row++){
			popout << data.getdate(data.popdata[pd].start + row*data.popdata[pd].units);
			for(r = 0; r <  meas.popnum[pd][row].size(); r++){ popout <<  "\t" << meas.popnum[pd][row][r];} popout << endl;
		}
	}
	
	for(md = 0; md < data.margdata.size(); md++){
		d = data.margdata[md].democat;
			
		file = data.margdata[md].file;
		filefull = dir+"/"+file;
		ofstream margout(filefull);
		if(!margout) emsg("Cannot output the file '"+filefull+"'");

		cout << "  '" << file << "' gives the '" << data.democat[d].name << "' stratified number of " << data.margdata[md].fromstr << "→" << data.margdata[md].tostr << " transitions";
		if(data.margdata[md].type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		margout << "Category"; 
		if(data.margdata[md].type == "reg"){ for(r = 0; r < data.nregion; r++){ margout << "\t" << data.region[r].code;}}
		else margout << "\t" << "all";
		margout << endl;
	
		for(j = 0; j < data.democat[d].value.size(); j++){
			margout << data.democat[d].value[j];
			for(r = 0; r < meas.margnum[md][j].size(); r++){ 
				sum = 0; for(jj = 0; jj < data.democat[d].value.size(); jj++) sum += meas.margnum[md][jj][r];
				margout << "\t" << (100.0*meas.margnum[md][j][r])/sum;
			} 
			margout << endl;
		}
	}
	cout << endl;
}

/// Calculates diagnostic statistics
STAT getstat(vector <double> &vec)                           
{
	unsigned int n, i, d;
	double sum, sum2, var, sd, a, cor, f;
	STAT stat;
	vector <double> vec2;
	
	n = vec.size();
	if(n == 0){
		stat.mean = "---"; stat.CImin = "---"; stat.CImax = "---"; stat.ESS = "---"; 
	}
	else{
		sum = 0; sum2 = 0; for(i = 0; i < n; i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
		sum /= n; sum2 /= n;
		stat.mean = to_string(sum); 
		
		vec2 = vec;
		sort(vec2.begin(),vec2.end());
	
		if(n >= 2){
			i = (unsigned int)((n-1)*0.05); f = (n-1)*0.05 - i;
			stat.CImin = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
				
			i = (unsigned int)((n-1)*0.95); f = (n-1)*0.95 - i;
			stat.CImax = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
		}
		else{
			stat.CImin = to_string(vec2[0]);
			stat.CImax = to_string(vec2[0]);
		}

		var = sum2 - sum*sum;
		if(var <= tiny || n <= 2)  stat.ESS = "---";
		else{	
			sd = sqrt(var);
			for(i = 0; i < n; i++) vec[i] = (vec[i]-sum)/sd;
				
			sum = 1;
			for(d = 1; d < n/2; d++){             // calculates the effective sample size
				a = 0; for(i = 0; i < n-d; i++) a += vec[i]*vec[i+d]; 
				cor = a/(n-d);
				if(cor < 0) break;
				sum += 2*cor;			
			}
			stat.ESS = to_string(int(n/sum));
		}
	}
	
	return stat;
}
