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

struct STAT{                                           // Stores statistical information
	double mean;                                         // The mean
	double CImin, CImax;                                 // The minimum and maximum of the 90% credible interval
	double ESS;                                          // The estimated sample size
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

/// Generates posterior plots for transitions, variation in R0 over time, parameter statistics and MCMC diagnostics 
void outputresults(DATA &data, MODEL &model, vector <SAMPLE> &opsamp)
{      
	unsigned int p, r, s, st, nopsamp, row, td, pd, md, j, jj, jmax, d;
	double sum;
	string file, filefull;
	vector <double> vec;
	STAT stat;
	string name;

	ensuredirectory(data.outputdir);
		
	nopsamp = opsamp.size();
	
	cout << endl;
	if(data.mode == MODE_SIM) cout << "Outputs in directory '" << data.outputdir << "':" << endl;
	else cout << "Posterior outputs in directory '" << data.outputdir << "':" << endl;
	
	for(td = 0; td < data.transdata.size(); td++){
		name = data.transdata[td].file;
	
		j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
		name = name.substr(0,j);
	
		if(data.transdata[td].type == "reg"){
			for(r = 0; r < data.nregion; r++){
				file = "Posterior_"+name+"_"+data.region[r].code+".txt";
				filefull = data.outputdir+"/"+file;
				ofstream dataout(filefull.c_str());
				if(!dataout) emsg("Cannot output the file '"+filefull+"'");
						
				cout << "'" << file << "' gives numbers of " << data.transdata[td].fromstr << "→" << data.transdata[td].tostr << " transitions for region '" << data.region[r].name << "'." << endl;
		
				for(row = 0; row < data.transdata[td].rows; row++){
					vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].meas.transnum[td][row][r]);
					stat = getstat(vec);
					
					dataout << data.transdata[td].start + (row+0.5)*data.transdata[td].units << " ";
					if(data.mode != MODE_SIM) dataout << data.transdata[td].num[r][row]; else dataout << stat.mean;
					dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
				}
			}
		}
		
		if(data.transdata[td].type == "all"){
			file = "Posterior_"+name+".txt";
			filefull = data.outputdir+"/"+file;
			ofstream dataout(filefull.c_str());
			if(!dataout) emsg("Cannot output the file '"+filefull+"'");
			
			cout << "'" << file << "' gives numbers of " << data.transdata[td].fromstr << "→" << data.transdata[td].tostr << " transitions." << endl;
	
			for(row = 0; row < data.transdata[td].rows; row++){
				vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].meas.transnum[td][row][0]);
				stat = getstat(vec);
				
				dataout << data.transdata[td].start + (row+0.5)*data.transdata[td].units << " ";
				if(data.mode != MODE_SIM) dataout << data.transdata[td].num[0][row]; else dataout << stat.mean;
				dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
			}
		}
	}
	
	for(pd = 0; pd < data.popdata.size(); pd++){
		name = data.popdata[pd].file;
	
		j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
		name = name.substr(0,j);
	
		if(data.popdata[pd].type == "reg"){
			for(r = 0; r < data.nregion; r++){
				file = "Posterior_"+name+"_"+data.region[r].code+".txt";
				filefull = data.outputdir+"/"+file;
				ofstream dataout(filefull.c_str());
				if(!dataout) emsg("Cannot output the file '"+filefull+"'");
						
				cout << "'" << file << "' gives the population in '" << data.popdata[pd].compstr << "' for region '" << data.region[r].name << "'." << endl;
		
				for(row = 0; row < data.popdata[pd].rows; row++){
					vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].meas.popnum[pd][row][r]);
					stat = getstat(vec);
					
					dataout << data.popdata[pd].start << " ";
					if(data.mode != MODE_SIM) dataout << data.popdata[pd].num[r][row]; else dataout << stat.mean;
					dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
				}
			}
		}
		
		if(data.popdata[pd].type == "all"){
			file = "Posterior_"+name+".txt";
			filefull = data.outputdir+"/"+file;
			ofstream dataout(filefull.c_str());
			if(!dataout) emsg("Cannot output the file '"+filefull+"'");
			
			cout << "'" << file << "' gives the population in '" << data.popdata[pd].compstr << "'." << endl;
	
			for(row = 0; row < data.popdata[pd].rows; row++){
				vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].meas.popnum[pd][row][0]);
				stat = getstat(vec);
				
				dataout << data.popdata[pd].start << " ";
				if(data.mode != MODE_SIM) dataout << data.popdata[pd].num[0][row]; else dataout << stat.mean;
				dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
			}
		}
	}
	
	for(md = 0; md < data.margdata.size(); md++){
		d = data.margdata[md].democat;
			
		name = data.margdata[md].file;

		j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
		name = name.substr(0,j);
	
		if(data.margdata[md].type == "reg"){
			for(r = 0; r < data.nregion; r++){
				file = "Posterior_"+name+"_"+data.region[r].code+".txt";
				filefull = data.outputdir+"/"+file;
				ofstream dataout(filefull.c_str());
				if(!dataout) emsg("Cannot output the file '"+filefull+"'");
						
				cout << "'" << file << "' gives " << data.margdata[md].fromstr << "→" << data.margdata[md].tostr << " transitions stratified by '" << data.democat[d].name << "' for region '" << data.region[r].name << "'." << endl;
		
				for(j = 0; j < data.democat[d].value.size(); j++){
					vec.clear();
					for(s = 0; s < nopsamp; s++){
						sum = 0; for(jj = 0; jj < data.democat[d].value.size(); jj++) sum += opsamp[s].meas.margnum[md][jj][r];
						vec.push_back(opsamp[s].meas.margnum[md][j][r]/sum);
					}
					stat = getstat(vec);
					
					dataout << data.democat[d].value[j] << " ";
					if(data.mode != MODE_SIM) dataout << data.margdata[md].percent[r][j]; else dataout << stat.mean;
					dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
				}
			}
		}
		
		if(data.margdata[md].type == "all"){
			file = "Posterior_"+name+".txt";
			filefull = data.outputdir+"/"+file;
			ofstream dataout(filefull.c_str());
			if(!dataout) emsg("Cannot output the file '"+filefull+"'");
			
			cout << "'" << file << "' gives the percentages in the " << data.margdata[md].fromstr << "→" << data.margdata[md].tostr << " transition stratified by '" << data.democat[d].name << "'." << endl;
		
			for(j = 0; j < data.democat[d].value.size(); j++){
				vec.clear(); 
				for(s = 0; s < nopsamp; s++){
					sum = 0; for(jj = 0; jj < data.democat[d].value.size(); jj++) sum += opsamp[s].meas.margnum[md][jj][0];
					vec.push_back(opsamp[s].meas.margnum[md][j][0]/sum);
				}
				stat = getstat(vec);
				
				dataout << data.democat[d].value[j] << " ";
				if(data.mode != MODE_SIM) dataout << data.margdata[md].percent[0][j]; else dataout << stat.mean;
				dataout << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
			}
		}
	}
	
	file = "Posterior_R0.txt";
	filefull = data.outputdir+"/"+file;
	ofstream R0out(filefull.c_str());
	if(!R0out) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time variation in R0." << endl;
	
	for(st = 0; st < data.nsettime; st++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		stat = getstat(vec);
		
		R0out << (st+0.5)*data.period/data.nsettime << " " 
		      << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}
	
	file = "Posterior_parameters.txt";
	filefull = data.outputdir+"/"+file;
	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the model parameters." << endl;
	
	for(p = 0; p < model.param.size(); p++){
		vec.clear(); for(s = 0; s < nopsamp; s++) vec.push_back(opsamp[s].paramval[p]);
		stat = getstat(vec);
			
		paramout << model.param[p].name  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << endl; 
	}
	
	if(data.mode != MODE_SIM){
		cout << "'" << data.outputdir << "/trace.txt' gives trace plots for model parameters." << endl;
	}
	
	if(data.mode == MODE_INF){
		cout << "'" << data.outputdir << "/traceLi.txt' gives trace plots for the observation likelihoods on different chains." << endl;
	}

	// This gives the acceptance rates for different MCMC proposals on different parameters
	
	if(data.mode != MODE_SIM){
		stringstream ss; ss << data.outputdir << "/MCMCdiagnostic.txt";
		ofstream diag(ss.str().c_str()); 
		if(!diag) emsg("Cannot output the file '"+ss.str()+"'");
	
		cout << "'" << ss.str() << "' gives MCMC diagnostics." << endl;
	
		diag << "MCMC diagnostics:" << endl;

		for(p = 0; p < model.param.size(); p++){
			diag << model.param[p].name << ": ";
			if(model.param[p].ntr == 0) diag << "Fixed" << endl;
			else diag << "Acceptance rate " << double(model.param[p].nac)/model.param[p].ntr << endl;
		}
	}
	cout << endl;
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
		
		transout << "Week"; 
		if(data.transdata[td].type == "reg"){ for(r = 0; r < data.nregion; r++){ transout << "\t" << data.region[r].code;}}
		else transout << "\t" << "all";
		transout << endl;
		
		for(row = 0; row < data.transdata[td].rows; row++){
			transout << row;
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
		
		popout << "Week"; 
		if(data.popdata[pd].type == "reg"){ for(r = 0; r < data.nregion; r++){ popout << "\t" << data.region[r].code;}}
		else popout << "\t" << "all";
		popout << endl;
		
		for(row = 0; row < data.popdata[pd].rows; row++){
			popout << row;
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
	double sum, sum2, sd, a, cor, f;
	STAT stat;
	
	n = vec.size();
	
	sum = 0; sum2 = 0; for(i = 0; i < n; i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
	sum /= n; sum2 /= n;
	stat.mean = sum; 
	
	sort(vec.begin(),vec.end());
	
	if (n >= 2)
	{
		i = (unsigned int)((n-1)*0.05); f = (n-1)*0.05 - i;
		stat.CImin = vec[i]*(1-f) + vec[i+1]*f;
			
		i = (unsigned int)((n-1)*0.95); f = (n-1)*0.95 - i;
		stat.CImax = vec[i]*(1-f) + vec[i+1]*f;
	}
	else
	{
		stat.CImin = vec[0];
		stat.CImax = vec[0];
	}

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
