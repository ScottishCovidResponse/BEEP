// This describes the observation model. This prescribes how likely the data is given a true event state.

using namespace std;

#include "obsmodel.hh"
#include "utils.hh"

vector <unsigned int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, unsigned int ti, unsigned int tf);

/// The total observation probability for a complete set of events fev
double Lobstot(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, double invT)
{
	unsigned int t;
	double Ltot;
	
	Ltot = 0; for(t = 0; t < data.period; t++) Ltot += Lobs(data,model,poptree,t,fev,invT);
	return Ltot;
}

/// Measures how well the particle agrees with the observations for a given time range t to t+1
/// (e.g. weekly hospitalised case data)
double Lobs(DATA &data, MODEL &model, POPTREE &poptree, unsigned int t, vector < vector <FEV> > &fev, double invT)
{
	unsigned int r, td, sum;
	double mean, var, L;
	vector <unsigned int> num;
	
	L = 0;	
	for(td = 0; td < data.transdata.size(); td++){
		num = getnumtrans(data,model,poptree,fev,data.transdata[td].from,data.transdata[td].to,t,t+1);
		
		if(data.transdata[td].type == "reg"){
			for(r = 0; r < data.nregion; r++){
				mean = data.transdata[td].num[r][t];
				var = mean; if(var < 5) var = 5;
				var *= varfac;
				L += -0.5*log(2*M_PI*var) - (mean-num[r])*(mean-num[r])/(2*var);
			}
		}
		
		if(data.transdata[td].type == "all"){
			sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
	
			mean = data.transdata[td].num[0][t];
			var = mean; if(var < 5) var = 5;
			var *= varfac;
			L += -0.5*log(2*M_PI*var) - (mean-sum)*(mean-sum)/(2*var);
		}
	}
		
	return invT*L;
}

/// Returns the number of transitions for individuals going from compartment "from" to compartment "to" 
/// in different regions over the time range ti - tf
vector <unsigned int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, unsigned int ti, unsigned int tf)
{
	unsigned int d, dmax, k, r, tra;
	FEV fe;
	vector <unsigned int> num;
	
	tra = 0; 
	while(tra < model.trans.size() && 
	     !(model.comp[model.trans[tra].from].name == from && model.comp[model.trans[tra].to].name == to)) tra++;
	if(tra == model.trans.size()) emsg("Obsmodel: Cannot find transition");
	
	for(r = 0; r < data.nregion; r++) num.push_back(0);

	dmax = fev.size();
	if(dmax != data.fediv) emsg("fevsize");
	for(d = (unsigned int)(dmax*double(ti)/data.period); d <= (unsigned int)(dmax*double(tf)/data.period); d++){
		if(d < dmax){
			for(k = 0; k < fev[d].size(); k++){
				fe = fev[d][k];
				if(fe.t > tf) break;
				if(fe.t > ti && fe.trans == tra) num[data.area[data.ind[fe.ind].area].region]++;
			}
		}
	}
	
	return num;
}

/// Measures how well the particle agrees with the observations for a given time range t to t+1
/// (e.g. weekly hospitalised case data)
double Lobs_mbp(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev)
{
	unsigned int r, td, sum, t;
	double mean, var, L;
	vector <unsigned int> num;
	
	L = 0;	
	for(t = 0; t < data.period; t++){
		for(td = 0; td < data.transdata.size(); td++){
			num = getnumtrans_mbp(data,model,poptree,trev,indev,data.transdata[td].from,data.transdata[td].to,t,t+1);
			
			if(data.transdata[td].type == "reg"){
				for(r = 0; r < data.nregion; r++){
					mean = data.transdata[td].num[r][t];
					var = mean; if(var < 5) var = 5;
					var *= varfac;
					L += -0.5*log(2*M_PI*var) - (mean-num[r])*(mean-num[r])/(2*var);
				}
			}
			
			if(data.transdata[td].type == "all"){
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
		
				mean = data.transdata[td].num[0][t];
				var = mean; if(var < 5) var = 5;
				var *= varfac;
				L += -0.5*log(2*M_PI*var) - (mean-sum)*(mean-sum)/(2*var);
			}
		}
	}
	
	return L;
}

/// Returns the number of transitions for individuals going from compartment "from" to compartment "to" 
/// in different regions over the time range ti - tf
vector <unsigned int> getnumtrans_mbp(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, string from, string to, unsigned int ti, unsigned int tf)
{
	unsigned int d, dmax, k, r, tra, sett, i, j;
	FEV fe;
	vector <unsigned int> num;
	
	tra = 0; 
	while(tra < model.trans.size() && 
	     !(model.comp[model.trans[tra].from].name == from && model.comp[model.trans[tra].to].name == to)) tra++;
	if(tra == model.trans.size()) emsg("Obsmodel: Cannot find transition");
	
	for(r = 0; r < data.nregion; r++) num.push_back(0);

	for(sett = data.settpertime*ti; sett < data.settpertime*tf; sett++){
		for(j = 0; j < trev[sett].size(); j++){
			i = trev[sett][j].ind;
			if(tra == indev[i][trev[sett][j].e].trans) num[data.area[data.ind[i].area].region]++;
		}
	}
	
	return num;
}
