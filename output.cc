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

Output::Output(const Details &details, const DATA &data, MODEL &model, Obsmodel &obsmodel) :  details(details), data(data), model(model), obsmodel(obsmodel)
{
}
	
/// Initialises trace plot for parameters
void Output::trace_plot_init()
{
	unsigned int p, pc;
	
	ensuredirectory(details.outputdir);
	stringstream ss; ss << details.outputdir << "/trace.txt";

	trace.open(ss.str().c_str());		
	trace << "state";
	for(p = 0; p < model.param.size(); p++) trace << "\t" << model.param[p].name; 
	for(pc = 0; pc < model.priorcomps.size(); pc++) trace << "\tProb in " << model.comp[model.priorcomps[pc].comp].name;
	trace << "\tLi"; 
	trace << "\tPri"; 
	trace << "\tninf";
	trace << endl;
}

/// Output trace plot
void  Output::trace_plot(unsigned int samp, double Li, double Pri, unsigned int ninf, const vector <double> &paramval)
{
	unsigned int p, c, pc, a;
	double pr;
	
	trace << samp; 
	for(p = 0; p < paramval.size(); p++) trace << "\t" << paramval[p]; 
	
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
	trace << "\t" << ninf; 
	trace << endl;
}

/// Initialises trace plot for likelihoods on difference chains (MBP only)
void Output::Li_trace_plot_init(unsigned int nchaintot)
{
	unsigned int p;
	
	ensuredirectory(details.outputdir);
	stringstream ss; ss << details.outputdir << "/traceLi.txt";

	traceLi.open(ss.str().c_str());		
	traceLi << "state";
	for(p = 0; p < nchaintot; p++) traceLi << "\tchain" <<  p; 
	traceLi << endl;
}

/// Outputs trace plot for likelihoods on difference chains (MBP only)
void Output::Li_trace_plot(unsigned int samp, unsigned int nchaintot, const vector <double> &Litot)
{
	unsigned int p;
	
	traceLi << samp;
	for(p = 0; p < nchaintot; p++) traceLi << "\t" << Litot[p]; 
	traceLi << endl;
}

/// Outputs a posterior graph
void Output::posterior_plot(const vector <SAMPLE> &opsamp, unsigned int d, unsigned int r, unsigned int type) const
{
	unsigned int j, jmax, row, s, t=0, nopsamp, opsampmin, rr, rrmin, rrmax, nrow=0, dc=0, jj;
	double sum, valsum;
	string name, file, filefull;
	vector <double> vec;
	STAT stat;
	
	nopsamp = opsamp.size();
	opsampmin = nopsamp/4;
	
	switch(type){
	case pop_data: name = data.popdata[d].file; break;
	case trans_data: name = data.transdata[d].file; break;
	case marg_data: name = data.margdata[d].file; break;
	default: emsgEC("Output",1); break;
	}
	
	name = filebasename(name); // TODO: think about this

	j = 0; jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
	name = name.substr(0,j);
	
	switch(r){
	case UNSET:	file = "Posterior_"+name+".txt"; break;
	case ADD:	file = "Posterior_"+name+"_sum.txt"; break;
	default: file = "Posterior_"+name+"_"+data.region[r].code+".txt"; break;
	}
	
	filefull = details.outputdir+"/"+file;
	ofstream dataout(filefull.c_str());
	if(!dataout) emsg("Cannot output the file '"+filefull+"'");
						
	switch(type){
	case pop_data: 
		cout << "'" << file << "' gives the population in '" << data.popdata[d].compstr << "'";
		dataout << "# The population in '" << data.popdata[d].compstr << "'"; 
		nrow = data.popdata[d].rows;
		break;
		
	case trans_data:
		cout << "'" << file << "' gives numbers of " << data.transdata[d].fromstr << "→" << data.transdata[d].tostr << " transitions";
		dataout << "# Population in " << data.transdata[d].fromstr << "→" << data.transdata[d].tostr << " transitions";
		nrow = data.transdata[d].rows;
		break;
	
	case marg_data:
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
	case ADD: rrmin = 0; rrmax = data.nregion; if(type == marg_data) emsgEC("Output",2); break;
	default: rrmin = r; rrmax = r+1; break;
	}

	switch(type){
	case trans_data: case pop_data: dataout << "# Time from start, " << details.tformat; break;
	case marg_data: dataout << "category"; break;
	}	
	dataout << ", data, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;
	for(row = 0; row < nrow; row++){
		vec.clear(); 
		for(s = opsampmin; s < nopsamp; s++){
			valsum = 0;
			for(rr = rrmin; rr < rrmax; rr++){
				switch(type){
				case pop_data: valsum += opsamp[s].meas.popnum[d][row][rr]; break;
				case trans_data: valsum += opsamp[s].meas.transnum[d][row][rr]; break;
				case marg_data:
					sum = 0; for(jj = 0; jj < data.democat[dc].value.size(); jj++) sum += opsamp[s].meas.margnum[d][jj][rr];
					valsum += 100.0*opsamp[s].meas.margnum[d][row][rr]/sum;
					break;
				}	
			}
			vec.push_back(valsum);
		}
		
		stat = getstat(vec);
		switch(type){
		case pop_data: t = data.popdata[d].start+row*data.popdata[d].units; break;
		case trans_data: t = data.transdata[d].start+row*data.transdata[d].units; break;
		}
		if(type == marg_data) dataout << data.democat[dc].value[row] << " ";
		else dataout << t << " " << details.getdate(t) << " ";
		
		if(details.mode != sim){
			valsum = 0;
			for(rr = rrmin; rr < rrmax; rr++){
				if(type == marg_data) valsum += data.margdata[d].percent[rr][row]; 
				else{
					switch(type){
					case pop_data: jj = data.popdata[d].num[rr][row]; break;
					case trans_data: jj = data.transdata[d].num[rr][row]; break;
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
		dataout << " " << stat.mean << " " << stat.CImin << " " << stat.CImax << " " << stat.ESS << endl; 
	}
}

/// Generates posterior plots for transitions, variation in R0 over time, parameter statistics and MCMC diagnostics 
void Output::results(const vector <PARAMSAMP> &psamp, const vector <SAMPLE> &opsamp) const
{      
	unsigned int p, r, s, st, nopsamp, opsampmin, npsamp, psampmin, d, b, c;
	double areaav;
	string file, filefull;
	vector <double> vec;
	STAT stat;
	string name;
	vector <DIST> paramdist;
	vector <double> paramav, Rav;

	ensuredirectory(details.outputdir);
		
	nopsamp = opsamp.size();
	opsampmin = nopsamp/4;
	
	cout << endl;
	if(details.mode == sim || details.mode == multisim) cout << "Outputs in directory '" << details.outputdir << "':" << endl;
	else cout << "Posterior outputs in directory '" << details.outputdir << "':" << endl;
	
	for(d = 0; d < data.transdata.size(); d++){
		if(data.transdata[d].type == "reg"){
			for(r = 0; r < data.nregion; r++) posterior_plot(opsamp,d,r,trans_data);
			posterior_plot(opsamp,d,ADD,trans_data);
		}
		if(data.transdata[d].type == "all") posterior_plot(opsamp,d,UNSET,trans_data);
	}
	
	for(d = 0; d < data.popdata.size(); d++){
		if(data.popdata[d].type == "reg"){
			for(r = 0; r < data.nregion; r++) posterior_plot(opsamp,d,r,pop_data);
			posterior_plot(opsamp,d,ADD,pop_data);
		}
		if(data.popdata[d].type == "all") posterior_plot(opsamp,d,UNSET,pop_data);
	}
	
	for(d = 0; d < data.margdata.size(); d++){
		if(data.margdata[d].type == "reg"){ for(r = 0; r < data.nregion; r++) posterior_plot(opsamp,d,r,marg_data);}
		if(data.margdata[d].type == "all") posterior_plot(opsamp,d,UNSET,marg_data);
	}

	file = "Posterior_R0.txt";
	filefull = details.outputdir+"/"+file;
	ofstream R0out(filefull.c_str());
	if(!R0out) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time variation in R0." << endl;
	R0out << "# Gives the time variation in R0." << endl;	
	R0out << "# Time from start, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	Rav.resize(details.nsettime);
	for(st = 0; st < details.nsettime; st++){
		vec.clear(); for(s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		stat = getstat(vec);
		Rav[st] = atof(stat.mean.c_str());
		
		R0out << (st+0.5)*details.period/details.nsettime << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}

	file = "Posterior_phi.txt";
	filefull = details.outputdir+"/"+file;
	ofstream phiout(filefull.c_str());
	if(!phiout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time variation in phi." << endl;
	phiout << "# Gives the time variation in phi (expressed as the number of infected per 1000000 individuals)." << endl;	
	phiout << "# Time from start, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	for(st = 0; st < details.nsettime; st++){
		vec.clear(); for(s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].phi[st]*1000000);
		stat = getstat(vec);
	
		phiout << (st+0.5)*details.period/details.nsettime << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}

	npsamp = psamp.size();
	psampmin = npsamp/4;
	
	file = "Posterior_parameters.txt";
	filefull = details.outputdir+"/"+file;
	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the model parameters." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200." << endl;
	paramout << endl;
	
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	paramav.resize(model.param.size());
	for(p = 0; p < model.param.size(); p++){
		vec.clear(); for(s = psampmin; s < npsamp; s++) vec.push_back(psamp[s].paramval[p]);
		stat = getstat(vec);
		paramav[p] = atof(stat.mean.c_str());
		if(p < model.param.size()-1){
			paramout << model.param[p].name  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << endl; 
		}
	}
	
	file = "Posterior_distributions.txt";
	filefull = details.outputdir+"/"+file;
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the probability distributions for parameters." << endl;
	
  distout << "# Posterior probability distributions for model parameters." << endl;
	distout << endl;

	paramdist.resize(model.param.size()-1);
	for(p = 0; p < model.param.size()-1; p++){
		vec.clear(); for(s = psampmin; s < npsamp; s++) vec.push_back(psamp[s].paramval[p]);
		paramdist[p] = getdist(vec);
	}
	
	for(p = 0; p < model.param.size()-1; p++){
		distout << model.param[p].name << "\t" << "Probability\t";
	}	
	distout << endl;
	
	for(b = 0; b < BIN; b++){
		for(p = 0; p < model.param.size()-1; p++){
			distout << paramdist[p].value[b] << "\t" << paramdist[p].prob[b] << "\t";
		}
		distout << endl;
	}
	
	file = "Posterior_Rmap.txt";
	filefull = details.outputdir+"/"+file;
	ofstream Rmapout(filefull.c_str());
	if(!Rmapout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time and spatial variation in R0." << endl;
	Rmapout << "# Gives the time and spatial variation in R." << endl;	
	Rmapout << "# Area, Day: ";
	for(st = 0; st < details.nsettime; st++){ Rmapout << (st+1); if(st < details.nsettime-1) Rmapout << ", ";}
	Rmapout << endl;
	
	model.setup(paramav);
	areaav = 0; for(c = 0; c < data.narea; c++)	areaav +=  model.areafac[c]/ data.narea;

	for(c = 0; c < data.narea; c++){
		Rmapout << data.area[c].code;
		for(st = 0; st < details.nsettime; st++) Rmapout << "\t" << (model.areafac[c]/areaav)*Rav[st];
		Rmapout << endl;
	}
	
	if(details.mode != sim){
		cout << "'trace.txt' gives trace plots for model parameters." << endl;
	}
	
	if(details.mode == inf){
		cout << "'traceLi.txt' gives trace plots for the observation likelihoods on different chains." << endl;
	}
}

void Output::combinedtrace(const vector <string> &paramname, const vector < vector < vector <double> > > &vals, 
                           string file, string distfile, unsigned int burnin) const
{
	unsigned int p, inp, s, npsamp, psampmin, N, M, nmax, nmin, b;
	double W, B, valav, varr, muav;
	string GR;
	vector <double> vec, mu, vari;
	STAT stat;
	vector <DIST> paramdist;
			
	ofstream paramout(file.c_str());
	if(!paramout) emsg("Cannot output the file '"+file+"'");
	
	cout << endl << "'" << file << "' gives the model parameters combining the input trace files." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200 and PSRF should be between 0.9 and 1.1." << endl;
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, effective sample size (ESS), potential scale reduction factor (PSRF) otherwise known as the Gelman-Rubin convergence diagnostic." << endl;
	paramout << endl;

	M = vals.size();
	nmax = large;
	for(inp = 0; inp < M; inp++) if(vals[inp][0].size() < nmax) nmax = vals[inp][0].size();

	if(burnin != UNSET){
		if(burnin >= nmax) emsg("The 'burnin' period must be greater than the number of samples.");
		nmin = burnin;
	}
	else nmin = nmax/4;
	
	N = nmax-nmin; 
		
	paramdist.resize(paramname.size());
			
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
			npsamp = vals[inp][p].size();
			if(burnin != UNSET){
				if(burnin >= npsamp) emsg("The 'burnin' period must be greater than the number of samples.");
				psampmin = burnin;
			}
			else psampmin = npsamp/4;
			for(s = psampmin; s < npsamp; s++) vec.push_back(vals[inp][p][s]);
		}
		paramdist[p] = getdist(vec);
		stat = getstat(vec);
				
		paramout << paramname[p]  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << " " << GR << endl; 
	}
	
	if(distfile != ""){
		ofstream distout(distfile.c_str());
		if(!distout) emsg("Cannot output the file '"+distfile+"'");
		
		cout << "'" << distfile << "' gives the probability distributions for parameters." << endl;
		
		distout << "# Posterior probability distributions for model parameters." << endl;
		distout << endl;

		for(p = 0; p < paramname.size(); p++){
			distout << paramname[p] << "\t" << "Probability\t";
		}	
		distout << endl;
		
		for(b = 0; b < BIN; b++){
			for(p = 0; p < paramname.size(); p++){
			distout << paramdist[p].value[b] << "\t" << paramdist[p].prob[b] << "\t";
			}
			distout << endl;
		}				
	}
}
	
/// Outputs an event sample fev
void Output::eventsample(const vector < vector <FEV> > &fev) const
{
	unsigned int d, j, nind;
	vector< vector <FEV> > indev;
	TRANS tr;
	
	nind = data.ind.size();
	indev.resize(nind);
	for(d = 0; d < fev.size(); d++){
		for(j = 0; j < fev[d].size(); j++) indev[fev[d][j].ind].push_back(fev[d][j]);
	}
	
	ensuredirectory(details.outputdir);
	stringstream sst; sst << details.outputdir << "/events.txt";
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
void Output::plot(string file, const vector < vector <FEV> > &xi, double tmin, double period) const
{
	unsigned int c, tra, td, tdf;
	double t;
	vector <int> N;
	TRANS tr;
	
	N.resize(model.comp.size()); for(c = 0; c < model.comp.size(); c++) N[c] = 0;
	N[0] = data.popsize;
		
	td = 0; tdf = 0; while(td < details.fediv && xi[td].size()==0) td++;
	
	ofstream plot(file.c_str());
	if(!plot) emsg("Cannot output the file '"+file+"'");
	
	for(t = tmin; t < period; t += (period-tmin)/100){
		while(td < details.fediv && xi[td][tdf].t < t){
			tra = xi[td][tdf].trans;
			tr = model.trans[tra];
			N[tr.from]--; N[tr.to]++;
			
			tdf++;
			if(tdf == xi[td].size()){
				td++; tdf = 0; 
				while(td < details.fediv && xi[td].size() == 0) td++;
			}
		}
		
		plot << t << " ";
		for(c = 0; c < model.comp.size(); c++) plot << N[c] << " ";
		plot << endl;
	}
}

/// Generates case data based on a simulation using the MBP algorithm
void Output::simulateddata(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, string dir) const
{
	unsigned int row, r, td, pd, md, d, j, jj;
	string file, filefull;
	double sum;
	MEAS meas;
	
	ensuredirectory(dir);
		
	meas = obsmodel.getmeas(trev,indev);
	
	cout << "Simulated data in directory '" << dir <<"':" << endl;
	for(td = 0; td < data.transdata.size(); td++){
		file = data.transdata[td].file;
		filefull =  dir+"/"+file;
		ofstream transout(filefull);
		if(!transout) emsg("Cannot output the file '"+filefull+"'");
		
		cout << "  '" << file << "' gives the observed ";
		switch(data.transdata[td].units){
		case 1: cout << "daily"; break;
		case 7: cout << "weekly"; break;
		default: emsg("Problem with units"); break;
		}
		
		cout << " number of " << data.transdata[td].fromstr << "→" << data.transdata[td].tostr << " transitions";
		if(data.transdata[td].type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		switch(details.tform){
		case tform_num: transout << "time"; break;
		case tform_ymd: transout << "date"; break;
		}

		if(data.transdata[td].type == "reg"){ for(r = 0; r < data.nregion; r++){ transout << "\t" << data.region[r].code;}}
		else transout << "\t" << "all";
		transout << endl;
		
		for(row = 0; row < data.transdata[td].rows; row++){
			transout << details.getdate(data.transdata[td].start + row*data.transdata[td].units);
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
		
		switch(details.tform){
		case tform_num: popout << "time"; break;
		case tform_ymd: popout << "date"; break;
		}
		
		if(data.popdata[pd].type == "reg"){ for(r = 0; r < data.nregion; r++){ popout << "\t" << data.region[r].code;}}
		else popout << "\t" << "all";
		popout << endl;
		
		for(row = 0; row < data.popdata[pd].rows; row++){
			popout << details.getdate(data.popdata[pd].start + row*data.popdata[pd].units);
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

bool PW_ord (PW p1,PW p2) { return (p1.val < p2.val); }

/// Gets mean and credible interval for a series of weighted samples
STAT Output::getstat_with_w(vector <PW> vec) const  
{
	STAT stat;
	
	auto n = vec.size();
	double sum = 0, sumw = 0; 
	for(auto i = 0u; i < n; i++){ sum += vec[i].val*vec[i].w; sumw += vec[i].w;}
	
	stat.mean = to_string(sum/sumw); 
	
	sort(vec.begin(),vec.end(),PW_ord);

	if(n >= 2){
		double sumtail = 0;
		auto i = 0u; while(i < n && sumtail < 0.05*sumw){ sumtail += vec[i].w; i++;}
		i--;
		double valw = vec[i].w;
		sumtail -= valw;
		double f = (0.05*sumw-sumtail)/valw;

		stat.CImin = to_string(vec[i].val*(1-f) + vec[i+1].val*f);
		
		sumtail = 0;
		i = 0; while(i < n && sumtail < 0.95*sumw){ sumtail += vec[i].w; i++;}
		i--;
		valw = vec[i].w;
		sumtail -= valw;
		f = (0.95*sumw-sumtail)/valw;

		stat.CImax = to_string(vec[i].val*(1-f) + vec[i+1].val*f);
	}
	else{
		stat.CImin = to_string(vec[0].val);
		stat.CImax = to_string(vec[0].val);
	}
	
	return stat;
}
		
/// Calculates diagnostic statistics
STAT Output::getstat(const vector <double> &vec) const                       
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

		vec2 = vec;
		var = sum2 - sum*sum;
		if(var <= tiny || n <= 2)  stat.ESS = "---";
		else{	
			sd = sqrt(var);
			for(i = 0; i < n; i++) vec2[i] = (vec2[i]-sum)/sd;
				
			sum = 1;
			for(d = 1; d < n/2; d++){             // calculates the effective sample size
				a = 0; for(i = 0; i < n-d; i++) a += vec2[i]*vec2[i+d]; 
				cor = a/(n-d);
				if(cor < 0) break;
				sum += 2*cor;			
			}
			stat.ESS = to_string(int(n/sum));
		}
	}
	
	return stat;
}


/// Gets the probability distributions for a given set of samples
DIST Output::getdist(const vector <double> &vec) const
{
	unsigned int i, b;
	double min, max, d;
	DIST dist;
	vector <unsigned int> bin;
	
	min = large; max = -large;	
	for(i = 0; i < vec.size(); i++){
		if(vec[i] < min) min = vec[i];
		if(vec[i] > max) max = vec[i];
	}
	
	dist.value.resize(BIN);
	dist.prob.resize(BIN);
	
	if(min == max){
		for(b = 0; b < BIN; b++){ dist.value[b] = "---"; dist.prob[b] = "---";} 
	}
	else{
		d = max-min;
		min -= 0.2*d; max += 0.2*d; 
	
		bin.resize(BIN);
		for(b = 0; b < BIN; b++) bin[b] = 0;
		
		for(i = 0; i < vec.size(); i++){
			b = (unsigned int)(BIN*(vec[i]-min)/(max-min+tiny)); if(b >= BIN) emsgEC("Output",3);
			bin[b]++;
		}

		for(b = 0; b < BIN; b++){ 
			dist.value[b] = to_string(min+(b+0.5)*(max-min)/BIN); 
			dist.prob[b] = to_string(double(bin[b])/vec.size());
		} 
	}
	
	return dist;
}

/// Outputs the probability distributions generated by 
void Output::plot_distribution(string file, const Generation &gen) const 
{
	const int smax = 100000;
	
	ensuredirectory(details.outputdir);
		
	if(details.mode != abcsmc) emsg("The mode should be ANC-SMC");
	
	string filefull = details.outputdir+"/"+file;
	ofstream distout(filefull.c_str());
	if(!distout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the probability distributions for parameters." << endl;
	
  distout << "# Posterior probability distributions for model parameters." << endl;
	distout << endl;

	vector <DIST> paramdist;
	paramdist.resize(model.param.size()-1);
	
	auto N = gen.param_samp.size();
	for(auto th = 0u; th < model.param.size()-1; th++){	
		vector <double> vec;
		
		double sumst[N];    // Generate particle sampler
		double sum = 0; for(auto i = 0u; i < N; i++){ sum += gen.w[i]; sumst[i] = sum;}
	
		for(auto s = 0u; s < smax; s++){
			double z = ran()*sum; auto k = 0u; while(k < N && z > sumst[k]) k++;
			if(k == N) emsg("Problem");
		
			vec.push_back(gen.param_samp[k][th]);
		}
	
		paramdist[th] = getdist(vec);
	}
	
	for(auto th = 0u;th < model.param.size()-1; th++){
		distout << model.param[th].name << "\t" << "Probability\t";
	}	
	distout << endl;
	
	for(auto b = 0u; b < BIN; b++){
		for(auto th = 0u; th < model.param.size()-1; th++){
			distout << paramdist[th].value[b] << "\t" << paramdist[th].prob[b] << "\t";
		}
		distout << endl;
	}
}

/// Outputs how statistics improve as a function of generation
void Output::generation_plot(string file, const vector <Generation> generation) const
{
	auto nparam = model.param.size(); 
	
	ensuredirectory(details.outputdir);
		
	string filefull = details.outputdir+"/"+file;
	ofstream genout(filefull.c_str());
	if(!genout) emsg("Cannot output the file '"+filefull+"'");

	genout << "# Genetarion"; 
	for(auto th = 0u; th < nparam; th++) genout << model.param[th].name << " (mean,CImin,CImax) ";
	genout << endl;
	
	long timestart = generation[0].time;
	for(auto g = 0u; g < generation.size(); g++){
		const Generation &gen=generation[g];
		genout << g << " " << (gen.time - timestart)/(60.0*CLOCKS_PER_SEC) << " " << gen.EFcut;
		
		for(auto th = 0u; th < nparam; th++){
			vector <PW> vec;
			for(auto i = 0u; i < gen.param_samp.size(); i++){
				PW pw;
				pw.val = gen.param_samp[i][th];
				switch(details.mode){
					case abcsmc: pw.w = gen.w[i]; break;
					case abcmbp: pw.w = 1; break;
					default: emsgEC("Output",11); break;
				}
				vec.push_back(pw);
			}

			STAT stat = getstat_with_w(vec);
			
			genout << " " << stat.mean << " " << stat.CImin << " " << stat.CImax;
		}
		genout << endl;
	}	
}

/// Plots the log of the model evidence based as a function of likelihood
double Output::model_evidence_plot(string file, const vector<Generation> &generation) const 
{
	cout << "'" << file << "' gives the log model evidence as a function of the cutoff in the log of the observation probability." << endl;

	ensuredirectory(details.outputdir);
		
	string filefull = details.outputdir+"/"+file;
	ofstream MEout(filefull.c_str());
	if(!MEout) emsg("Cannot output the file '"+filefull+"'");
	
	double ME = 0;
	for(auto g = 0u; g < generation.size()-1; g++){
		auto num = 0u;
		auto num_cut = 0u;
		auto EFcut = generation[g+1].EFcut;
		auto EF_upper_limit = generation[g].EFcut;  
		
		for(auto gg = 0u; gg <= g; gg++){
			for(auto EF : generation[gg].EF_samp){
				if(gg == g && EF > EF_upper_limit){ cout << g << " " << EF << " " << EF_upper_limit << " y\n";  emsgEC("Output",10);}
				if(EF < EF_upper_limit){
					num++;
					if(EF < EFcut) num_cut++;
				}
			}
		}
		ME += log(double(num_cut)/num);

		MEout << EFcut << " " << ME << endl;
	}
	
	return ME;
}
