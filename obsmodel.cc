
using namespace std;

#include "obsmodel.hh"
#include "utils.hh"

vector <int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, int ti, int tf);

double Lobstot(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, double invT)
{
	int t;
	double Ltot;
	
	Ltot = 0; for(t = 0; t < data.period; t++) Ltot += Lobs(data,model,poptree,t,fev,invT);
	return Ltot;
}

/// Measures how well the particle agrees with the observations for a given time range t to t+1
/// (e.g. weekly hospitalised case data)
double Lobs(DATA &data, MODEL &model, POPTREE &poptree, int t, vector < vector <FEV> > &fev, double invT)
{
	int r, td, sum;
	double mean, var, L;
	vector <int> num;
	
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
vector <int> getnumtrans(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <FEV> > &fev, string from, string to, int ti, int tf)
{
	int d, dmax, k, r, tra;
	FEV fe;
	vector <int> num;
	
	tra = 0; 
	while(tra < model.trans.size() && 
	     !(model.comp[model.trans[tra].from].name == from && model.comp[model.trans[tra].to].name == to)) tra++;
	if(tra == model.trans.size()) emsg("Finescale: Cannot find transition");
	
	for(r = 0; r < data.nregion; r++) num.push_back(0);

	dmax = fev.size();
	if(dmax != data.fediv) emsg("fevsize");
	for(d = int(dmax*double(ti)/data.period); d <= int(dmax*double(tf)/data.period); d++){
		if(d < dmax){
			for(k = 0; k < fev[d].size(); k++){
				fe = fev[d][k];
				if(fe.t > tf) break;
				if(fe.t > ti && fe.trans == tra) num[data.house[poptree.ind[fe.ind].houseref].region]++;
			}
		}
	}
	
	return num;
}
