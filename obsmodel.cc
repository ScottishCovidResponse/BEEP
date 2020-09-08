// This describes the observation model. This prescribes how likely the data is given a true event state.

#include "utils.hh"
#include "model.hh"
#include "obsmodel.hh"

using namespace std;

Obsmodel::Obsmodel(const Details &details, const DATA &data, const MODEL &model) : details(details), data(data), model(model)
{
}

/// Gets all measured quantities
MEAS Obsmodel::getmeas(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev) const
{
	MEAS meas;
	meas.transnum.resize(data.transdata.size());
	for(auto td = 0u; td < data.transdata.size(); td++){                                   // Incorporates transition observations
		meas.transnum[td].resize(data.transdata[td].rows);
		for(auto row = 0u; row < data.transdata[td].rows; row++){
			auto ti = data.transdata[td].start + row*data.transdata[td].units;
			auto tf = ti + data.transdata[td].units;

			auto num = getnumtrans(trev,indev,data.transdata[td].trans,ti,tf,UNSET,UNSET);
			
			if(data.transdata[td].type == "reg"){
				meas.transnum[td][row].resize(data.nregion);
				for(auto r = 0u; r < data.nregion; r++) meas.transnum[td][row][r] = num[r];
			}
			
			if(data.transdata[td].type == "all"){
				meas.transnum[td][row].resize(1);
				auto sum = 0.0; for(auto r = 0u; r < data.nregion; r++) sum += num[r];
				meas.transnum[td][row][0] = sum;
			}
		}
	}
	
	meas.popnum.resize(data.popdata.size());
	for(auto pd = 0u; pd < data.popdata.size(); pd++){                                   // Incorporates population observations
		meas.popnum[pd].resize(data.popdata[pd].rows);
		
		vector <unsigned int> num(data.nregion);
		for(auto r = 0u; r < data.nregion; r++) num[r] = 0;
		
		auto c = data.popdata[pd].comp;
		
		auto sett = 0u;
		for(auto row = 0u; row < data.popdata[pd].rows; row++){
			auto ti = data.popdata[pd].start + row*data.popdata[pd].units;
				
			while(details.settime[sett] < ti){				
				for(auto j = 0u; j < trev[sett].size(); j++){
					auto i = trev[sett][j].ind;
					auto tra = indev[i][trev[sett][j].e].trans;
					if(model.trans[tra].from == c) num[data.area[data.ind[i].area].region]--;
					if(model.trans[tra].to == c) num[data.area[data.ind[i].area].region]++;
				}
				sett++;
			}
	
			if(data.popdata[pd].type == "reg"){
				meas.popnum[pd][row].resize(data.nregion);
				for(auto r = 0u; r < data.nregion; r++) meas.popnum[pd][row][r] = num[r];
			}
			
			if(data.popdata[pd].type == "all"){
				meas.popnum[pd][row].resize(1);
				auto sum = 0.0; for(auto r = 0u; r < data.nregion; r++) sum += num[r];
				meas.popnum[pd][row][0] = sum;
			}
		}
	}

	meas.margnum.resize(data.margdata.size());
	for(auto md = 0u; md < data.margdata.size(); md++){                                     // Incorporates marginal distributions
		auto d = data.margdata[md].democat;
		
		meas.margnum[md].resize(data.democat[d].value.size());
	
		if(data.margdata[md].type == "reg"){
			for(auto j = 0u; j < data.democat[d].value.size(); j++){
				auto num = getnumtrans(trev,indev,data.margdata[md].trans,0,details.period,d,j);
				
				meas.margnum[md][j].resize(data.nregion);
				for(auto r = 0u; r < data.nregion; r++) meas.margnum[md][j][r] = num[r];
			}
		}
		
		if(data.margdata[md].type == "all"){
			for(auto j = 0u; j < data.democat[d].value.size(); j++){
				auto num = getnumtrans(trev,indev,data.margdata[md].trans,0,details.period,d,j);
				auto sum = 0.0; for(auto r = 0u; r < data.nregion; r++) sum += num[r];
				
				meas.margnum[md][j].resize(1);
				meas.margnum[md][j][0] = sum;
			}
		}
	}
	
	return meas;
}

/// The contribution from a single measurement 
double Obsmodel::singobs(unsigned int mean, unsigned int val) const
{
	switch(mean){
	case UNKNOWN: return 0;     // The data value is unknown
	case THRESH:                // The case in which there is a threshold applied to the observation
		if(val <= data.threshold){
			if(details.mode == abcmbp || details.mode == abcsmc) return 0;
			else return data.thres_h;
		}
		else{
			auto var = minvar;
			if(details.mode == abcmbp || details.mode == abcsmc) return (val-data.threshold)*(val-data.threshold)/var;
			else return data.thres_h - (val-data.threshold)*(val-data.threshold)/(2*var);
		}
	default:                    // A measurement is made
		auto var = mean; if(var < minvar) var = minvar;
		if(details.mode == abcmbp || details.mode == abcsmc) return (val-mean)*(val-mean)/var;
		return normalprob(val,mean,var);
	}
}

/// Measures how well the particle agrees with the observations for a given time range t to t+1
/// (e.g. weekly hospitalised case data)
double Obsmodel::Lobs(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev) const 
{
	MEAS meas = getmeas(trev,indev);
	
	auto L = 0.0;	
	for(auto td = 0u; td < meas.transnum.size(); td++){                                   // Incorporates transition observations
		for(auto row = 0u; row < meas.transnum[td].size(); row++){
			for(auto r = 0u; r < meas.transnum[td][row].size(); r++){
				L += singobs(data.transdata[td].num[r][row],meas.transnum[td][row][r]);
			}
		}
	}

	for(auto pd = 0u; pd < meas.popnum.size(); pd++){                                   // Incorporates population observations
		for(auto row = 0u; row < meas.popnum[pd].size(); row++){
			for(auto r = 0u; r < meas.popnum[pd][row].size(); r++){
				L += singobs(data.popdata[pd].num[r][row],meas.popnum[pd][row][r]);
			}
		}
	}
	
	for(auto md = 0u; md < meas.margnum.size(); md++){                                   // Incorporates marginal observations
		for(auto row = 0u; row < meas.margnum[md].size(); row++){
			for(auto r = 0u; r < meas.margnum[md][row].size(); r++){
				auto sum = 0.0; for(auto j = 0u; j < meas.margnum[md].size(); j++) sum += meas.margnum[md][j][r];
				
				auto val = meas.margnum[md][row][r];
				auto mean = data.margdata[md].percent[r][row]*sum/100.0;
				auto var = mean; if(var < minvar) var = minvar;
				if(details.mode == abcmbp || details.mode == abcsmc) L += (val-mean)*(val-mean)/var;
				else L += normalprob(val,mean,var);
			}
		}
	}

	return L;
}

/// Returns the number of transitions for individuals going down a transition
/// in different regions over the time range ti - tf
/// If the demographic catergoty d is set then it must have the value v
vector <unsigned int> Obsmodel::getnumtrans(const vector < vector <EVREF> > &trev, const vector < vector <FEV> > &indev, unsigned int tra, unsigned int ti, unsigned int tf, unsigned int d, unsigned int v) const
{
	vector <unsigned int> num(data.nregion);
	for(auto r = 0u; r < data.nregion; r++) num[r] = 0;

	for(auto sett = details.settpertime*ti; sett < details.settpertime*tf; sett++){
		for(auto j = 0u; j < trev[sett].size(); j++){
			auto i = trev[sett][j].ind;
			if(d == UNSET || data.democatpos[data.ind[i].dp][d] == v){		
				if(tra == indev[i][trev[sett][j].e].trans) num[data.area[data.ind[i].area].region]++;
			}
		}
	}
	
	return num;
}
