// This describes the observation model. This prescribes how likely the data is given a true event state.

#include "utils.hh"
#include "model.hh"
#include "poptree.hh"
#include "obsmodel.hh"

using namespace std;

/// Gets all measured quantities
MEAS getmeas(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev)
{
	unsigned int td, pd, md, ti, tf, r, row, sett, c, j, i, tra, d;
	double sum;
	vector <unsigned int> num;
	MEAS meas;
	
	meas.transnum.resize(data.transdata.size());
	for(td = 0; td < data.transdata.size(); td++){                                   // Incorporates transition observations
		meas.transnum[td].resize(data.transdata[td].rows);
		for(row = 0; row < data.transdata[td].rows; row++){
			ti = data.transdata[td].start + row*data.transdata[td].units;
			tf = ti + data.transdata[td].units;

			num = getnumtrans(data,model,poptree,trev,indev,data.transdata[td].trans,ti,tf,UNSET,UNSET);
			
			if(data.transdata[td].type == "reg"){
				meas.transnum[td][row].resize(data.nregion);
				for(r = 0; r < data.nregion; r++) meas.transnum[td][row][r] = num[r];
			}
			
			if(data.transdata[td].type == "all"){
				meas.transnum[td][row].resize(1);
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
				meas.transnum[td][row][0] = sum;
			}
		}
	}
	
	meas.popnum.resize(data.popdata.size());
	for(pd = 0; pd < data.popdata.size(); pd++){                                   // Incorporates population observations
		meas.popnum[pd].resize(data.popdata[pd].rows);
		
		num.resize(data.nregion); for(r = 0; r < data.nregion; r++) num[r] = 0;
		c = data.popdata[pd].comp;
		
		sett = 0;
		for(row = 0; row < data.popdata[pd].rows; row++){
			ti = data.popdata[pd].start + row*data.popdata[pd].units;
				
			while(data.settime[sett] < ti){				
				for(j = 0; j < trev[sett].size(); j++){
					i = trev[sett][j].ind;
					tra = indev[i][trev[sett][j].e].trans;
					if(model.trans[tra].from == c) num[data.area[data.ind[i].area].region]--;
					if(model.trans[tra].to == c) num[data.area[data.ind[i].area].region]++;
				}
				sett++;
			}
	
			if(data.popdata[pd].type == "reg"){
				meas.popnum[pd][row].resize(data.nregion);
				for(r = 0; r < data.nregion; r++) meas.popnum[pd][row][r] = num[r];
			}
			
			if(data.popdata[pd].type == "all"){
				meas.popnum[pd][row].resize(1);
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
				meas.popnum[pd][row][0] = sum;
			}
		}
	}

	meas.margnum.resize(data.margdata.size());
	for(md = 0; md < data.margdata.size(); md++){                                     // Incorporates marginal distributions
		d = data.margdata[md].democat;
		
		meas.margnum[md].resize(data.democat[d].value.size());
	
		if(data.margdata[md].type == "reg"){
			for(j = 0; j < data.democat[d].value.size(); j++){
				num = getnumtrans(data,model,poptree,trev,indev,data.margdata[md].trans,0,data.period,d,j);
				
				meas.margnum[md][j].resize(data.nregion);
				for(r = 0; r < data.nregion; r++) meas.margnum[md][j][r] = num[r];
			}
		}
		
		if(data.margdata[md].type == "all"){
			for(j = 0; j < data.democat[d].value.size(); j++){
				num = getnumtrans(data,model,poptree,trev,indev,data.margdata[md].trans,0,data.period,d,j);
				sum = 0; for(r = 0; r < data.nregion; r++) sum += num[r];
				
				meas.margnum[md][j].resize(1);
				meas.margnum[md][j][0] = sum;
			}
		}
	}
	
	return meas;
}

/// Measures how well the particle agrees with the observations for a given time range t to t+1
/// (e.g. weekly hospitalised case data)
double Lobs(DATA &data, MODEL &model, POPTREE &poptree, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev)
{
	unsigned int td, pd, md, row, r, j;
	double mean, val, var, L, sum;
	MEAS meas;
	
	meas = getmeas(data,model,poptree,trev,indev);
	
	L = 0;	
	for(td = 0; td < meas.transnum.size(); td++){                                   // Incorporates transition observations
		for(row = 0; row < meas.transnum[td].size(); row++){
			for(r = 0; r < meas.transnum[td][row].size(); r++){
				val = meas.transnum[td][row][r];
				mean = data.transdata[td].num[r][row];
				var = mean; if(var < 5) var = 5;
				var *= varfac;
				L += normalprob(val,mean,var);
			}
		}
	}
	
	for(pd = 0; pd < meas.popnum.size(); pd++){                                   // Incorporates population observations
		for(row = 0; row < meas.popnum[pd].size(); row++){
			for(r = 0; r < meas.popnum[pd][row].size(); r++){
				val = meas.popnum[pd][row][r];
				mean = data.popdata[pd].num[r][row];
				var = mean; if(var < 5) var = 5;
				var *= varfac;
				L += normalprob(val,mean,var);
			}
		}
	}
	
	for(md = 0; md < meas.margnum.size(); md++){                                   // Incorporates population observations
		for(row = 0; row < meas.margnum[md].size(); row++){
			for(r = 0; r < meas.margnum[md][row].size(); r++){
				sum = 0; for(j = 0; j < meas.margnum[md].size(); j++) sum += meas.margnum[md][j][r];
				
				val = meas.margnum[md][row][r];
				mean = data.margdata[md].percent[r][row]*sum/100.0;
				var = mean; if(var < 5) var = 5;
				var *= varfac;
				L += normalprob(val,mean,var);
			}
		}
	}

	return L;
}

/// Returns the number of transitions for individuals going down a transition
/// in different regions over the time range ti - tf
/// If the demographic catergoty d is set then it must have the value v
vector <unsigned int> getnumtrans(DATA &data, MODEL & /*model*/, POPTREE &/*poptree*/, vector < vector <EVREF> > &trev, vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v)
{
	unsigned int r, sett, i, j;
	vector <unsigned int> num;

	for(r = 0; r < data.nregion; r++) num.push_back(0);

	for(sett = data.settpertime*ti; sett < data.settpertime*tf; sett++){
		for(j = 0; j < trev[sett].size(); j++){
			i = trev[sett][j].ind;
			if(d == UNSET || data.democatpos[data.ind[i].dp][d] == v){		
				if(tra == indev[i][trev[sett][j].e].trans) num[data.area[data.ind[i].area].region]++;
			}
		}
	}
	
	return num;
}
