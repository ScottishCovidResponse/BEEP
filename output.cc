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

Output::Output(const Details &details, const DATA &data, const MODEL &model, const Obsmodel &obsmodel) :  details(details), data(data), model(model), obsmodel(obsmodel)
{
	ensuredirectory(details.outputdir);
}
	
/// Create a directory if it doesn't already exist
void Output::ensuredirectory(const string &path) const
{
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	// Directory not found
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(ret == -1) emsg("Error creating directory '"+path+"'");
	}
}

/// Prints the populations to the terminal	
void Output::populations(unsigned int sett, const vector <int> &N) const
{
	cout  << "  Time: " << details.settime[sett];
	for(auto c = 0u; c < model.comp.size(); c++) cout << "  " << model.comp[c].name << ":"	<< N[c];
	cout << endl;	
}
		
/// Initialises trace plot for parameters
void Output::trace_plot_init()
{
	stringstream ss; ss << details.outputdir << "/trace.txt";

	trace.open(ss.str().c_str());		
	trace << "state";
	for(const auto& par : model.param) trace << "\t" << par.name; 
	for(const auto& pri : model.priorcomps) trace << "\tProb in " << model.comp[pri.comp].name;
	trace << "\tLi"; 
	trace << "\tPri"; 
	trace << "\tninf";
	trace << endl;
}

/// Output trace plot
void  Output::trace_plot(unsigned int samp, double Li, double Pri, unsigned int ninf, const vector <double> &paramval)
{
	trace << samp; 
	for(const auto& pval : paramval) trace << "\t" << pval; 
	
	if(model.priorcomps.size() > 0){
		vector <CompTrans> comptrans;
		model.create_comptrans(comptrans,paramval);
		vector <CompProb> compprob = model.create_compprob(comptrans);
	
		for(const auto& pri : model.priorcomps){
			auto c = pri.comp;
			auto pr = 0.0; for(auto a = 0u; a < data.nage; a++) pr += data.agedist[a]*compprob[c].value[a];
			trace << "\t" << pr;
		}
	}

	trace << "\t" << Li; 
	trace << "\t" << Pri; 
	trace << "\t" << ninf; 
	trace << endl;
}

/// Initialises trace plot for likelihoods on difference chains (MBP only)
void Output::L_trace_plot_init(unsigned int nchaintot)
{
	stringstream ss; ss << details.outputdir << "/traceLi.txt";

	traceLi.open(ss.str().c_str());		
	traceLi << "state";
	for(auto p = 0u; p < nchaintot; p++) traceLi << "\tchain" <<  p; 
	traceLi << endl;
}

/// Outputs trace plot for likelihoods on difference chains (MBP only)
void Output::L_trace_plot(unsigned int samp, const vector <double> &Litot)
{
	traceLi << samp;
	for(const auto& Lito : Litot) traceLi << "\t" << Lito; 
	traceLi << endl;
}

/// Outputs a posterior graph
void Output::posterior_plot(const vector <SAMPLE> &opsamp, unsigned int d, unsigned int r, unsigned int type) const
{
	auto nopsamp = opsamp.size();
	auto opsampmin = nopsamp/4;
	
	string name;
	switch(type){
	case pop_data: name = data.popdata[d].file; break;
	case trans_data: name = data.transdata[d].file; break;
	case marg_data: name = data.margdata[d].file; break;
	default: emsgEC("Output",1); break;
	}
	
	name = filebasename(name); // TODO: think about this

	auto j = 0u; auto jmax = name.length(); while(j < jmax && name.substr(j,1) != ".") j++;
	name = name.substr(0,j);
	
	string file;
	switch(r){
	case UNSET:	file = "Posterior_"+name+".txt"; break;
	case ADD:	file = "Posterior_"+name+"_sum.txt"; break;
	default: file = "Posterior_"+name+"_"+data.region[r].code+".txt"; break;
	}
	
	auto filefull = details.outputdir+"/"+file;
	ofstream dataout(filefull.c_str());
	if(!dataout) emsg("Cannot output the file '"+filefull+"'");
				
	unsigned int nrow = UNSET;
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
		auto dc = data.margdata[d].democat;
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
	
	unsigned int rrmin, rrmax;
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
	for(auto row = 0u; row < nrow; row++){
		vector <double> vec;
		for(auto s = opsampmin; s < nopsamp; s++){
			auto valsum = 0.0;
			for(auto rr = rrmin; rr < rrmax; rr++){
				switch(type){
				case pop_data: valsum += opsamp[s].meas.popnum[d][row][rr]; break;
				case trans_data: valsum += opsamp[s].meas.transnum[d][row][rr]; break;
				case marg_data:
					auto dc = data.margdata[d].democat;
					auto sum = 0.0; for(auto jj = 0u; jj < data.democat[dc].value.size(); jj++) sum += opsamp[s].meas.margnum[d][jj][rr];
					valsum += 100.0*opsamp[s].meas.margnum[d][row][rr]/sum;
					break;
				}	
			}
			vec.push_back(valsum);
		}
		
		auto stat = getstat(vec);
		
		unsigned int t=0;
		switch(type){
		case pop_data: t = data.popdata[d].start+row*data.popdata[d].units; break;
		case trans_data: t = data.transdata[d].start+row*data.transdata[d].units; break;
		}
		if(type == marg_data){
			auto dc = data.margdata[d].democat;
			dataout << data.democat[dc].value[row] << " ";
		}
		else dataout << t << " " << details.getdate(t) << " ";
		
		if(details.mode != sim){
			auto valsum = 0.0;
			for(auto rr = rrmin; rr < rrmax; rr++){
				if(type == marg_data) valsum += data.margdata[d].percent[rr][row]; 
				else{
					unsigned int jj = 0;
					switch(type){
					case pop_data: jj = data.popdata[d].num[rr][row]; break;
					case trans_data: jj = data.transdata[d].num[rr][row]; break;
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
	auto nopsamp = opsamp.size();
	auto opsampmin = nopsamp/4;
	
	cout << endl;
	if(details.mode == sim || details.mode == multisim) cout << "Outputs in directory '" << details.outputdir << "':" << endl;
	else cout << "Posterior outputs in directory '" << details.outputdir << "':" << endl;
	
	for(auto d = 0u; d < data.transdata.size(); d++){
		if(data.transdata[d].type == "reg"){
			for(auto r = 0u; r < data.nregion; r++) posterior_plot(opsamp,d,r,trans_data);
			posterior_plot(opsamp,d,ADD,trans_data);
		}
		if(data.transdata[d].type == "all") posterior_plot(opsamp,d,UNSET,trans_data);
	}
	
	for(auto d = 0u; d < data.popdata.size(); d++){
		if(data.popdata[d].type == "reg"){
			for(auto r = 0u; r < data.nregion; r++) posterior_plot(opsamp,d,r,pop_data);
			posterior_plot(opsamp,d,ADD,pop_data);
		}
		if(data.popdata[d].type == "all") posterior_plot(opsamp,d,UNSET,pop_data);
	}
	
	for(auto d = 0u; d < data.margdata.size(); d++){
		if(data.margdata[d].type == "reg"){ for(auto r = 0u; r < data.nregion; r++) posterior_plot(opsamp,d,r,marg_data);}
		if(data.margdata[d].type == "all") posterior_plot(opsamp,d,UNSET,marg_data);
	}

	auto file = "Posterior_R0.txt";
	auto filefull = details.outputdir+"/"+file;
	ofstream R0out(filefull.c_str());
	if(!R0out) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the time variation in R0." << endl;
	R0out << "# Gives the time variation in R0." << endl;	
	R0out << "# Time from start, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	vector <double> Rav(details.nsettime);
	for(auto st = 0u; st < details.nsettime; st++){
		vector <double> vec; for(auto s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].R0[st]);
		auto stat = getstat(vec);
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

	for(auto st = 0u; st < details.nsettime; st++){
		vector <double> vec; for(auto s = opsampmin; s < nopsamp; s++) vec.push_back(opsamp[s].phi[st]*1000000);
		auto stat = getstat(vec);
	
		phiout << (st+0.5)*details.period/details.nsettime << " " << stat.mean << " " << stat.CImin << " "<< stat.CImax << " " << stat.ESS << endl; 
	}

	auto npsamp = psamp.size();
	auto psampmin = npsamp/4;
	
	file = "Posterior_parameters.txt";
	filefull = details.outputdir+"/"+file;
	ofstream paramout(filefull.c_str());
	if(!paramout) emsg("Cannot output the file '"+filefull+"'");
	
	cout << "'" << file << "' gives the model parameters." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200." << endl;
	paramout << endl;
	
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, estimated sample size" << endl;

	vector <double> paramav(model.param.size());
	for(auto p = 0u; p < model.param.size(); p++){
		vector <double> vec; for(auto s = psampmin; s < npsamp; s++) vec.push_back(psamp[s].paramval[p]);
		auto stat = getstat(vec);
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

	vector <DIST> paramdist(model.param.size()-1);
	for(auto p = 0u; p < model.param.size()-1; p++){
		vector <double> vec; for(auto s = psampmin; s < npsamp; s++) vec.push_back(psamp[s].paramval[p]);
		paramdist[p] = getdist(vec);
	}
	
	for(auto p = 0u; p < model.param.size()-1; p++){
		distout << model.param[p].name << "\t" << "Probability\t";
	}	
	distout << endl;
	
	for(auto b = 0u; b < BIN; b++){
		for(auto p = 0u; p < model.param.size()-1; p++){
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
	for(auto st = 0u; st < details.nsettime; st++){ Rmapout << (st+1); if(st < details.nsettime-1) Rmapout << ", ";}
	Rmapout << endl;
	
	vector <double> areafac = model.create_areafac(paramav);
	auto areaav = 0.0; for(const auto& areafa : areafac) areaav += areafa/data.narea;

	for(auto c = 0u; c < data.narea; c++){
		Rmapout << data.area[c].code;
		for(auto st = 0u; st < details.nsettime; st++) Rmapout << "\t" << (areafac[c]/areaav)*Rav[st];
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
	ofstream paramout(file.c_str());
	if(!paramout) emsg("Cannot output the file '"+file+"'");
	
	cout << endl << "'" << file << "' gives the model parameters combining the input trace files." << endl;
	paramout << "# Posterior distributions for model parameters." << endl;
	paramout << "# For convergence ESS should be greater than 200 and PSRF should be between 0.9 and 1.1." << endl;
	paramout << "# Name, mean, minimum of 95% credible interval, maximum of 95% credible interval, effective sample size (ESS), potential scale reduction factor (PSRF) otherwise known as the Gelman-Rubin convergence diagnostic." << endl;
	paramout << endl;

	auto M = vals.size();
	auto nmax = large;
	for(const auto& vs : vals){ if(vs[0].size() < nmax) nmax = vs[0].size();}

	auto nmin = nmax/4;
	if(burnin != UNSET){
		if(burnin >= nmax) emsg("The 'burnin' period must be greater than the number of samples.");
		nmin = burnin;
	}
	
	auto N = nmax-nmin; 
		
	vector <DIST> paramdist(paramname.size());
		
	vector <double> mu(M), vari(M);
	for(auto p = 0u; p < paramname.size(); p++){
		string GR = "---";
		if(M > 1){
			auto muav = 0.0;
			for(auto inp = 0u; inp < M; inp++){
				auto valav = 0.0; for(auto s = nmin; s < nmax; s++) valav += vals[inp][p][s]/N;
				auto varr = 0.0; for(auto s = nmin; s < nmax; s++) varr += (vals[inp][p][s]-valav)*(vals[inp][p][s]-valav)/(N-1);
				mu[inp] = valav;
				vari[inp] = varr;
				muav += mu[inp]/M;
			}
					
			auto W = 0.0; for(auto var : vari) W += var/M;
			auto B = 0.0; for(auto muu : mu) B += (muu-muav)*(muu-muav)*N/(M-1);
			if(W > tiny) GR = to_string(sqrt(((1-1.0/N)*W+B/N)/W));
		}
	
		vector <double> vec;
		for(const auto& vs : vals){
			auto npsamp = vs[p].size();
			auto psampmin = npsamp/4;
			if(burnin != UNSET){
				if(burnin >= npsamp) emsg("The 'burnin' period must be greater than the number of samples.");
				psampmin = burnin;
			}
			for(auto s = psampmin; s < npsamp; s++) vec.push_back(vs[p][s]);
		}
		paramdist[p] = getdist(vec);
		STAT stat = getstat(vec);
				
		paramout << paramname[p]  << " " <<  stat.mean << " (" << stat.CImin << " - "<< stat.CImax << ") " << stat.ESS << " " << GR << endl; 
	}
	
	if(distfile != ""){
		ofstream distout(distfile.c_str());
		if(!distout) emsg("Cannot output the file '"+distfile+"'");
		
		cout << "'" << distfile << "' gives the probability distributions for parameters." << endl;
		
		distout << "# Posterior probability distributions for model parameters." << endl;
		distout << endl;

		for(auto parname : paramname) distout << parname << "\t" << "Probability\t";
		distout << endl;
		
		for(auto b = 0u; b < BIN; b++){
			for(const auto& pardist : paramdist) distout << pardist.value[b] << "\t" << pardist.prob[b] << "\t";
			distout << endl;
		}				
	}
}
	
/// Outputs an event sample fev
void Output::eventsample(const vector < vector <FEV> > &fev) const
{
	auto nind = data.ind.size();
	vector< vector <FEV> > indev;
	indev.resize(nind);
	for(auto& fe: fev){
		for(auto& ev : fe) indev[ev.ind].push_back(ev);
	}
	
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
	vector <int> N(model.comp.size());
	for(auto& NN : N) NN = 0;
	N[0] = data.popsize;
		
	auto td = 0u, tdf = 0u; while(td < details.fediv && xi[td].size()==0) td++;
	
	ofstream plot(file.c_str());
	if(!plot) emsg("Cannot output the file '"+file+"'");
	
	for(auto t = tmin; t < period; t += (period-tmin)/100){
		while(td < details.fediv && xi[td][tdf].t < t){
			auto tra = xi[td][tdf].trans;
			TRANS tr = model.trans[tra];
			N[tr.from]--; N[tr.to]++;
			
			tdf++;
			if(tdf == xi[td].size()){
				td++; tdf = 0; 
				while(td < details.fediv && xi[td].size() == 0) td++;
			}
		}
		
		plot << t << " ";
		for(auto NN : N) plot << NN << " ";
		plot << endl;
	}
}

/// Generates case data based on a simulation using the MBP algorithm
void Output::simulateddata(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, string dir) const
{	
	MEAS meas = obsmodel.getmeas(trev,indev);
	
	cout << "Simulated data in directory '" << dir <<"':" << endl;
	for(auto td = 0u; td < data.transdata.size(); td++){
		auto tdata = data.transdata[td];
		
		auto file = tdata.file;
		auto filefull =  dir+"/"+file;
		ofstream transout(filefull);
		if(!transout) emsg("Cannot output the file '"+filefull+"'");
		
		cout << "  '" << file << "' gives the observed ";
		switch(tdata.units){
		case 1: cout << "daily"; break;
		case 7: cout << "weekly"; break;
		default: emsg("Problem with units"); break;
		}
		
		cout << " number of " << tdata.fromstr << "→" << tdata.tostr << " transitions";
		if(tdata.type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		switch(details.tform){
		case tform_num: transout << "time"; break;
		case tform_ymd: transout << "date"; break;
		}

		if(tdata.type == "reg"){ for(const auto& reg : data.region){ transout << "\t" << reg.code;}}
		else transout << "\t" << "all";
		transout << endl;
		
		for(auto row = 0u; row < tdata.rows; row++){
			transout << details.getdate(tdata.start + row*tdata.units);
			for(const auto& trnum : meas.transnum[td][row]){ transout <<  "\t" << trnum;} transout << endl;
		}
	}
	
	for(auto pd = 0u; pd < data.popdata.size(); pd++){
		auto pdata = data.popdata[pd];
		
		auto file = pdata.file;
		auto filefull = dir+"/"+file;
		ofstream popout(filefull);
		if(!popout) emsg("Cannot output the file '"+filefull+"'");
		
		cout << "  '" << file << "' gives the numbers in population '" << pdata.compstr << "'";
		if(pdata.type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		switch(details.tform){
		case tform_num: popout << "time"; break;
		case tform_ymd: popout << "date"; break;
		}
		
		if(pdata.type == "reg"){ for(const auto& reg : data.region){ popout << "\t" << reg.code;}}
		else popout << "\t" << "all";
		popout << endl;
		
		for(auto row = 0u; row < pdata.rows; row++){
			popout << details.getdate(pdata.start + row*pdata.units);
			for(auto popnum : meas.popnum[pd][row]){ popout <<  "\t" << popnum;} popout << endl;
		}
	}
	
	for(auto md = 0u; md < data.margdata.size(); md++){
		auto mdata = data.margdata[md];
		
		auto d = mdata.democat;
			
		auto file = mdata.file;
		auto filefull = dir+"/"+file;
		ofstream margout(filefull);
		if(!margout) emsg("Cannot output the file '"+filefull+"'");

		cout << "  '" << file << "' gives the '" << data.democat[d].name << "' stratified number of " << mdata.fromstr << "→" << mdata.tostr << " transitions";
		if(mdata.type == "reg") cout << " for different regions." << endl;
		else cout << "." << endl;
		
		margout << "Category"; 
		if(mdata.type == "reg"){ for(const auto& reg : data.region){ margout << "\t" << reg.code;}}
		else margout << "\t" << "all";
		margout << endl;
	
		for(auto j = 0u; j < data.democat[d].value.size(); j++){
			margout << data.democat[d].value[j];
			for(auto r = 0u; r < meas.margnum[md][j].size(); r++){ 
				auto sum = 0.0; for(auto jj = 0u; jj < data.democat[d].value.size(); jj++) sum += meas.margnum[md][jj][r];
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
	
	double sum = 0, sumw = 0; 
	for(auto v : vec){ sum += v.val*v.w; sumw += v.w;}
	
	stat.mean = to_string(sum/sumw); 
	
	sort(vec.begin(),vec.end(),PW_ord);

	auto n = vec.size();
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
	STAT stat;
	
	auto n = vec.size();
	if(n == 0){
		stat.mean = "---"; stat.CImin = "---"; stat.CImax = "---"; stat.ESS = "---"; 
	}
	else{
		auto sum = 0.0, sum2 = 0.0; for(auto v : vec){ sum += v; sum2 += v*v;}
		sum /= n; sum2 /= n;
		stat.mean = to_string(sum); 
		
		auto vec2 = vec;
		sort(vec2.begin(),vec2.end());
	
		if(n >= 2){
			auto i = (unsigned int)((n-1)*0.05); auto f = (n-1)*0.05 - i;
			stat.CImin = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
				
			i = (unsigned int)((n-1)*0.95); f = (n-1)*0.95 - i;
			stat.CImax = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
		}
		else{
			stat.CImin = to_string(vec2[0]);
			stat.CImax = to_string(vec2[0]);
		}

		vec2 = vec;
		auto var = sum2 - sum*sum;
		if(var <= tiny || n <= 2)  stat.ESS = "---";
		else{	
			auto sd = sqrt(var);
			for(auto& v : vec2) v = (v-sum)/sd;
				
			auto sum = 1.0;
			for(auto d = 1u; d < n/2; d++){             // calculates the effective sample size
				auto a = 0.0; for(auto i = 0u; i < n-d; i++) a += vec2[i]*vec2[i+d]; 
				auto cor = a/(n-d);
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
	DIST dist;
	
	auto min = large, max = -large;	
	for(auto v : vec){
		if(v < min) min = v;
		if(v > max) max = v;
	}
	
	dist.value.resize(BIN);
	dist.prob.resize(BIN);
	
	if(min == max){
		for(auto b = 0u; b < BIN; b++){ dist.value[b] = "---"; dist.prob[b] = "---";} 
	}
	else{
		auto d = max-min;
		min -= 0.2*d; max += 0.2*d; 
	
		vector <unsigned int> bin(BIN);
		for(auto b = 0u; b < BIN; b++) bin[b] = 0;
		
		for(auto v : vec){
			auto b = (unsigned int)(BIN*(v-min)/(max-min+tiny)); if(b >= BIN) emsgEC("Output",3);
			bin[b]++;
		}

		for(auto b = 0u; b < BIN; b++){ 
			dist.value[b] = to_string(min+(b+0.5)*(max-min)/BIN); 
			dist.prob[b] = to_string(double(bin[b])/vec.size());
		} 
	}
	
	return dist;
}

/// Outputs the probability distributions generated by abc
void Output::plot_distribution(string file, const Generation &gen) const 
{
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
	
		const int smax = 100000;
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

	string filefull = details.outputdir+"/"+file;
	ofstream genout(filefull.c_str());
	if(!genout) emsg("Cannot output the file '"+filefull+"'");

	genout << "# Genetarion"; 
	for(const auto& param : model.param) genout << param.name << " (mean,CImin,CImax) ";
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
				if(gg == g && EF > EF_upper_limit) emsgEC("Output",10);
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
