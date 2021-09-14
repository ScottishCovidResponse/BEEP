/// This provides all the functions relating to generating MCMC proposals in parameter space

#include <assert.h>
#include <math.h>
#include <sstream>
#include <algorithm>   

using namespace std;

#include "param_prop.hh"
#include "mvn.hh"
#include "model.hh"
#include "output.hh"

ParamProp::ParamProp(const Details &details, const Data &data, const Model &model, const Output &output, Mpi &mpi) : details(details), data(data), model(model), output(output), mpi(mpi)
{
	if(sim_only == true){
		add_single();
	}
	else{
		if(details.mcmc_update.full_mvn == true){
			MVN mv("Full MVN",model.param_not_fixed,0.3,ALL_PARAM,MULTIPLE);
			mvn.push_back(mv);
		}
		
		if(details.mcmc_update.mvn == true){
			add_mvn("Distribution value",DISTVAL_PARAM,0.5);
			add_mvn("Branching probability",BRANCHPROB_PARAM,0.5);
			add_mvn("Infectivity",INF_PARAM,0.5);
			add_mvn("Regional effect",RE_PARAM,0.5);
			add_mvn("Sigma",SIGMA_PARAM,0.3);
			add_mvn("Observation parameters",OBS_PARAM,0.5);
			add_mvn("Covariates",COVAR_PARAM,0.5);
			add_mvn("R / efoi",R_EFOI_PARAM,0.3);
			add_mvn("Geomix",GEOMIX_PARAM,0.5);
			add_mvn("Modify param",MODIFY_PARAM,0.5);
			add_mvn("Susceptibility",SUSCEPTIBILITY_PARAM,0.5);
			add_mvn("Demo Specific",DEMO_SPECIFIC_PARAM,0.5);
		}
		
		if(details.mcmc_update.single == true) add_single();

		//if(details.mcmc_update.demo_spec == true) add_demographic_specific();
		
		if(details.mcmc_update.mean_time == true) mean_time_init();
		
		if(details.mcmc_update.neighbour == true) neighbour_init();
		
		if(details.mcmc_update.joint == true) joint_init();
		
		covar_area_init();
		
		if(details.mode == ABC_MBP || details.mode == MC3_INF || details.mode == MCMC_MBP || details.mode == PAIS_INF){
			FixedTree ft; ft.n = 0; ft.sim_frac = 1;
			fixedtree.push_back(ft);
		
			SliceTime st; st.sett_i = 0; st.sett_f = details.ndivision; st.sim_frac = 1;
			slicetime.push_back(st);
		}
	}
	
	zero_ntr_nac();
}


/// Adds MVN proposals which just affect a single age group
void ParamProp::add_demographic_specific()
{
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		vector <unsigned int> var;
	
		for(auto c = 0u; c < data.ndemocat; c++){
			if(data.democat[c].sus_vari == true){
				auto th = data.democat[c].sus_param[data.democatpos[dp][c]];
				if(model.param[th].priortype != FIXED_PRIOR) var.push_back(th); 
			}				
		}
		
		for(auto c = 0u; c < model.comp.size(); c++){
			auto kmax = model.comp[c].trans.size();
			if(kmax > 1){
				for(auto k = 0u; k < kmax; k++){
					auto th = model.trans[model.comp[c].trans[k]].param_prob[dp];
					if(model.param[th].priortype != FIXED_PRIOR) var.push_back(th);
				}
			}
		}
		
		if(var.size() > 1){
			stringstream ss;
			ss << "Demogaphic: ";
			for(auto c = 0u; c < data.ndemocat; c++){
				if(c > 0) ss << ",";
				ss << data.democat[c].value[data.democatpos[dp][c]];
			}		
			MVN mv(ss.str(),var,1.0,DEMO_SPECIFIC_PARAM,MULTIPLE);
			mvn.push_back(mv);
		}
	}
}


/// Adds a multivariate normal proposal distribution
void ParamProp::add_mvn(const string name, const ParamType type, const double size)
{
	auto var = model.parameter_type[type];
	if(var.size() > 0){
		MVN mv(name,var,size,type,MULTIPLE);
		mvn.push_back(mv);
	}
}


/// Adds a multivariate normal proposal distribution
void ParamProp::add_single()
{
	for(auto th = 0u; th < model.param.size(); th++){
		auto type = model.param[th].type;
		if(model.param[th].priortype != FIXED_PRIOR){
			vector <unsigned int> vec; vec.push_back(th);
			MVN mv("Univariate: "+model.param[th].name,vec,1,type,SINGLE);
			mvn.push_back(mv);
		}
	}
}


/// Joint proposals on covariate and area effects 
void ParamProp::covar_area_init()
{
	if(data.area_effect.on == false) return;
	
	for(auto loop = 0u; loop < 3; loop++){
		for(auto i = 0u; i < model.covariate_param.size(); i++){
			CovarArea ca;
			ca.covar_ref = i;
				
			ca.size = 0.05;
			covar_area.push_back(ca);
		}
	}
}


/// Initialises mean time proposals
void ParamProp::mean_time_init()
{
	for(auto &co : model.comp){
		if(co.trans.size() > 0){
			
			bool flag = true;
			
			vector <unsigned int> list;
				for(auto th : co.param_mean){
				if(th == UNSET || model.param[th].priortype == FIXED_PRIOR) flag = false;
				else{
					if(find(list.begin(),list.end(),th) == list.end()) list.push_back(th);
				}
			}
		
			if(list.size() > 0 && flag == true){
				vector <unsigned int> list_rev;
				for(auto trto : co.trans){
					for(auto th : model.comp[model.trans[trto].to].param_mean){
						if(th == UNSET || model.param[th].priortype == FIXED_PRIOR) flag = false;
						else{
							if(find(list_rev.begin(),list_rev.end(),th) == list_rev.end()) list_rev.push_back(th);
						}
					}
				}
				
				if(list_rev.size() > 0 && flag == true){
					for(auto li : list){  // Makes sure parameters not on both
						for(auto lirev : list_rev){
							if(li == lirev) flag = false;
						}
					}
					
					if(flag == true){
						MeanTime mt;
						mt.param_mean = list;
						mt.param_mean_rev = list_rev;
			
						mt.size = 1;
						mean_time.push_back(mt);
					}
				}
			}
		}
	}	
	
	if(false &&  mean_time.size() > 0){
		for(auto &mt : mean_time){
			for(auto th : mt.param_mean) cout << model.param[th].name << ", ";
			cout << ": ";
			for(auto th : mt.param_mean_rev) cout << model.param[th].name << ", ";
			cout << "Mean Time proposal" << endl;
		}
		emsg("Done");
	}
}


/// Initialises spline neighbour proposals
void ParamProp::neighbour_init()
{
	for(const auto &spl : model.spline){
		auto np = spl.p.size();
	
		for(auto i = 0u; i < np-1; i++){	
			if(spl.p[i].t != spl.p[i+1].t){
				auto param1 = spl.p[i].param;
				auto param2 = spl.p[i+1].param;
				if(param1 != param2){
					if(model.param[param1].priortype != FIXED_PRIOR && model.param[param2].priortype != FIXED_PRIOR){
						Neighbour rn;
						rn.param1 = param1;
						rn.param2 = param2;
						rn.size = 0.1;
						neighbour.push_back(rn);
					}
				}
			}
		}	
	}
	
	if(false && neighbour.size() > 0){
		for(auto &rn : neighbour){
			cout << "Param: " << model.param[rn.param1].name << ", " << model.param[rn.param2].name << endl;
			cout << "Neighbour proposal" << endl;
		}
		emsg("Done");
	}
}


/// Makes proposals jointly on many variables on splines (e.g. sinusoidal variation)
void ParamProp::joint_init()
{
	for(const auto &spl : model.spline){
		auto np = spl.p.size();
	
		vector <unsigned int> var_list;
		for(auto i = 0u; i < np; i++){	
			auto th = spl.p[i].param;
			if(model.param[th].priortype != FIXED_PRIOR){
				auto j = 0u; while(j < var_list.size() && th != var_list[j]) j++;
				if(j == var_list.size()) var_list.push_back(th);
			}
		}
		
		if(var_list.size() > 1){
			Joint rn;
			rn.var_list = var_list;
			rn.size = 0.1;
			rn.type = UP_DOWN;
			joint.push_back(rn);
			
			unsigned int max = var_list.size()/3;
			for(auto j = 0u; j < max; j++){
				rn.type = SINE;
				rn.sinenum = j;
				joint.push_back(rn);
			}
		}
	}	
	
	if(false && joint.size() > 0){
		for(auto &rn : joint){
			cout << "Param: ";
			for(auto th : rn.var_list) cout << model.param[th].name << ", ";
			cout << endl;
			cout << "Joint proposal" << endl;
		}
		emsg("Done");
	}
}


/// Zeros quanties relating to acceptance probability
void ParamProp::zero_ntr_nac()
{
	self.ntr = 0; self.nac = 0; self.nbo = 0;
	
	for(auto &mv : mvn){ mv.ntr = 0; mv.nac = 0; mv.nbo = 0;}
	
	for(auto &mt : mean_time){ mt.ntr = 0; mt.nac = 0; mt.nbo = 0;}
	
	for(auto &rn : neighbour){ rn.ntr = 0; rn.nac = 0; rn.nbo = 0;}
	
	for(auto &jo : joint){ jo.ntr = 0; jo.nac = 0; jo.nbo = 0;}
	
	for(auto &ca : covar_area){ ca.ntr = 0; ca.nac = 0;}
	
	for(auto &ft : fixedtree){ ft.ntr = 0; ft.nac = 0;}
	
	for(auto &st : slicetime){ st.ntr = 0; st.nac = 0;}
}


/// Updates proposal sizes based on acceptance probability
void ParamProp::update_proposals()
{	
	for(auto &mv : mvn){
		auto a_rate = mpi.get_acrate(mv.nbo,mv.ntr); 
		mv.bo_rate = a_rate;
		
		a_rate = mpi.get_acrate(mv.nac,mv.ntr); 
		mv.ac_rate = a_rate;
		
		update(mv.size,a_rate);
	}
	
	for(auto &mt : mean_time){
		auto a_rate = mpi.get_acrate(mt.nbo,mt.ntr); 
		mt.bo_rate = a_rate;
		
		a_rate = mpi.get_acrate(mt.nac,mt.ntr); 
		mt.ac_rate = a_rate;
		
		update(mt.size,a_rate);
	}
	
	for(auto &rn : neighbour){
		auto a_rate = mpi.get_acrate(rn.nbo,rn.ntr); 
		rn.bo_rate = a_rate;
		
		a_rate = mpi.get_acrate(rn.nac,rn.ntr); 
		rn.ac_rate = a_rate;
		
		update(rn.size,a_rate);
	}
	
	for(auto &rn : joint){
		auto a_rate = mpi.get_acrate(rn.nbo,rn.ntr); 
		rn.bo_rate = a_rate;
		
		a_rate = mpi.get_acrate(rn.nac,rn.ntr); 
		rn.ac_rate = a_rate;
		
		update(rn.size,a_rate);
	}
	
	for(auto &ca : covar_area){
		auto a_rate = mpi.get_acrate(ca.nac,ca.ntr); 
		ca.ac_rate = a_rate;
		
		update(ca.size,a_rate);
	}
	
	for(auto &ft : fixedtree){
		auto a_rate = mpi.get_acrate(ft.nac,ft.ntr); 
		ft.ac_rate = a_rate;
		
		update_high(ft.sim_frac,a_rate);
		if(ft.sim_frac > 1) ft.sim_frac = 1;
	}
	
	for(auto &st : slicetime){
		auto a_rate = mpi.get_acrate(st.nac,st.ntr); 
		st.ac_rate = a_rate;
		
		update_high(st.sim_frac,a_rate);
		if(st.sim_frac > 1) st.sim_frac = 1;
	}
	
	if(mpi.core == 0 && diagnotic_output == true) cout << print_proposal_information(true);
	
	update_fixedtree();
	
	update_splicetime();
}


/// Sets the acceptance rate (PMCMC)
void ParamProp::set_ac_rate()
{	
	self.bo_rate = double(self.nbo)/self.ntr;
	self.ac_rate = double(self.nac)/self.ntr;
		
	for(auto &mv : mvn){
		mv.bo_rate = double(mv.nbo)/mv.ntr;
		mv.ac_rate = double(mv.nac)/mv.ntr;
	}
	
	for(auto &mt : mean_time){
		mt.bo_rate = double(mt.nbo)/mt.ntr;
		mt.ac_rate = double(mt.nac)/mt.ntr;
	}
	
	for(auto &rn : neighbour){
		rn.bo_rate = double(rn.nbo)/rn.ntr;
		rn.ac_rate = double(rn.nac)/rn.ntr;
	}
	
	for(auto &rn : joint){
		rn.bo_rate = double(rn.nbo)/rn.ntr;
		rn.ac_rate = double(rn.nac)/rn.ntr;
	}
	
	for(auto &ca : covar_area){
		ca.ac_rate = double(ca.nac)/ca.ntr;
	}
}


/// Calculates the number of times a proposal needs to be made for simulation-only proposals 
void ParamProp::update_sim_proposals()
{
	for(auto &mv : mvn){
		auto a_rate = mpi.get_acrate(mv.nac,mv.ntr); 
		mv.ac_rate = a_rate;
		
		mv.number = (unsigned int)(1.0/a_rate + 0.5); if(mv.number < 1) mv.number = 1;
		if(mpi.core == 0) cout << mv.name << " " <<  mv.number << "nu" << endl;
	}
}	

		
/// Updates fixed tree proposals
void ParamProp::update_fixedtree()
{
	auto fti = 0u;
	auto ftimax = fixedtree.size();
	while(fti < ftimax){
		auto &ft = fixedtree[fti];
		
		auto flag = 0u;
		if(ft.sim_frac < 0.5){
			auto n = ft.n;
			auto &no = data.genQ.treenode[n];
			if(no.child.size() > 0){
				for(auto nn : no.child){
					FixedTree ft; ft.n = nn; ft.sim_frac = 0.8; ft.ac_rate = 0;
					fixedtree.push_back(ft);
				}
				flag = 1;
			}
		}
		
		if(flag == 1){ fixedtree.erase(fixedtree.begin()+fti); ftimax--;}
		else fti++;
	}
}


/// Updates slice time proposals
void ParamProp::update_splicetime()
{
	auto sti = 0u;
	auto stimax = slicetime.size();
	while(sti < stimax){
		auto &st = slicetime[sti];
		
		auto flag = 0u;
		if(st.sim_frac < 0.5){
			auto ti = st.sett_i, tf = st.sett_f;
				
			if(tf-ti > details.ndivision*0.05){
				flag = 1;
				unsigned int tmid = (ti+tf)/2;
				unsigned int tleft = ti*0.75 + 0.25*tf;
				unsigned int tright = ti*0.25 + 0.75*tf;
			
				SliceTime st; st.sim_frac = 0.6; st.ac_rate = 0;
				
				if(tmid > ti){
					st.sett_i = ti; st.sett_f = tmid; 
					slicetime.push_back(st);
				}
				
				if(tf > tmid){
					st.sett_i = tmid; st.sett_f = tf; 
					slicetime.push_back(st);
				}
				
				if(tright > tleft){
					st.sett_i = tleft; st.sett_f = tright; 
					slicetime.push_back(st);
				}
			}
		}
		
		if(flag == 1){ slicetime.erase(slicetime.begin()+sti); stimax--;}
		else sti++;
	}
}

	
/// Updates the size of parameter proposals
void ParamProp::update(double &val, double acrate)
{	
	if(acrate > 0.8) val *= 2;
	else{
		if(acrate > 0.6) val *= 1.5;
		else{
			if(acrate > 0.5) val *= 1.2;
			else{
				if(acrate < 0.4) val *= 0.8;
				else{
					if(acrate < 0.3) val *= 0.7;
					else{
						if(acrate < 0.2) val *= 0.5;
					}
				}
			}
		}
	}
}


/// Updates the size of parameter proposals (a larger change)
void ParamProp::update_high(double &val, double acrate)
{	
	if(acrate > 0.8) val *= 1.5;
	else{
		if(acrate > 0.6) val *= 1.2;
		else{
			if(acrate > 0.5) val *= 0.8;
			else{
				if(acrate < 0.4) val *= 0.7;
				else{
					if(acrate < 0.3) val *= 0.5;
					else{
						if(acrate < 0.2) val *= 0.3;
					}
				}
			}
		}
	}
}


/// Returns a list of all the proposals to be performed 
vector <Proposal> ParamProp::get_proposal_list(const vector <ParamSample> &param_samp)
{				
	for(auto &mv : mvn) mv.setup(param_samp);
	zero_ntr_nac();

	vector <Proposal> prop_list;
	
	if(details.mode == PMCMC_INF){
		Proposal prop; prop.type = SELF_PROP; prop.num = UNSET;
		prop_list.push_back(prop);
	}
	
	for(auto i = 0u; i < mvn.size(); i++){
		auto num = 1u;
		if(details.mcmc_update.mvn_multiple == true){
			auto numf = (details.mcmc_update.multiple_factor/(mvn[i].size*mvn[i].size));
			num = (unsigned int)(numf+0.5);
			
			switch(mvn[i].mvntype){
				case MULTIPLE: if(num > 50) num = 50; break;
				case SINGLE: if(num > 10) num = 10; break;
			}
		}
		if(num < 3) num = 3; 
		//if(num < 1) num = 1; 
		//if(num < 6) num = 6; 
		
		Proposal prop; prop.type = MVN_PROP; prop.num = i;
		for(auto j = 0u; j < num; j++) prop_list.push_back(prop);
	}
		
	for(auto i = 0u; i < mean_time.size(); i++){
		Proposal prop; prop.type = MEAN_TIME_PROP; prop.num = i;
		prop_list.push_back(prop);
	}
	
	for(auto i = 0u; i < neighbour.size(); i++){
		Proposal prop; prop.type = NEIGHBOUR_PROP; prop.num = i;
		prop_list.push_back(prop);
	}	
	
	for(auto i = 0u; i < joint.size(); i++){
		Proposal prop; prop.type = JOINT_PROP; prop.num = i;
		prop_list.push_back(prop);
	}	
	
	for(auto i = 0u; i < covar_area.size(); i++){
		Proposal prop; prop.type = COVAR_AREA_PROP; prop.num = i;
		prop_list.push_back(prop);
	}	
	
	if(details.mode == ABC_MBP || details.mode == MC3_INF || details.mode == MCMC_MBP || details.mode == PAIS_INF){
		for(auto fti = 0u; fti < fixedtree.size(); fti++){		
			Proposal prop; prop.type = FIXEDTREE_PROP; prop.num = fti;
			prop_list.push_back(prop);
		}
	
		for(auto sti = 0u; sti < slicetime.size(); sti++){		
			Proposal prop; prop.type = SLICETIME_PROP; prop.num = sti;
			prop_list.push_back(prop);
		}
	}
	
	randomise(prop_list);
		
	return prop_list;
}
	

/// Prints a list of proposals
void ParamProp::print_prop_list(const vector <Proposal> &prop_list) const
{
	for(const auto &prop : prop_list){
		auto num = prop.num;
		switch(prop.type){
			case MVN_PROP: cout << "MVN: " << mvn[num].name << " " << mvn[num].size << endl; break;
			case MEAN_TIME_PROP: cout << "Mean Time: " << num << " " << mean_time[num].size << endl; break;
			case NEIGHBOUR_PROP: cout << "Neighbour: " << num << " " << neighbour[num].size << endl; break;
			case JOINT_PROP: cout << "Joint: " << num << " " << joint[num].size << endl; break;
			case COVAR_AREA_PROP: cout << "Covar Area: " << num << " " << covar_area[num].size << endl; break;
			case SELF_PROP: cout << "Self: " << endl; break;
			case FIXEDTREE_PROP: cout << "Fixed tree: " << num << endl; break;
			case SLICETIME_PROP: cout << "Slice time: " << num << endl; break;
		}
	}
}
	
	
/// Randomises a vector of proposals
void ParamProp::randomise(vector <Proposal> &vlist)
{
	vector <Proposal> vlist_random; 
	auto n = vlist.size();
	for(auto v = 0u; v < n; v++){
		auto i = (unsigned int)(ran()*vlist.size());
		if(vlist.size() == 0) emsgEC("ParamProp",1);
		if(i >= vlist.size()) emsgEC("ParamProp",2);
		
		vlist_random.push_back(vlist[i]);
		
		if(vlist.size() > 1) vlist[i] = vlist[vlist.size()-1];
		vlist.pop_back();
	}
	if(vlist.size() != 0) emsgEC("ParamProp",3);
	
	vlist = vlist_random;
}


/// Generates a covariance matrix from a sets of parameter samples
vector <double> ParamProp::variance_vector(const vector <ParamSample> &param_samp) const 
{
	vector <double> vec;
	auto N = param_samp.size();                             // Generates the covariance matrix

	auto nvar = model.param_not_fixed.size();
	vec.resize(nvar);
	for(auto i = 0u; i < nvar; i++){
		auto th = model.param_not_fixed[i]; 
	
		double av = 0;
		for(auto k = 0u; k < N; k++) av += param_samp[k].paramval[th];
		av /= N;
		
		double av2 = 0;
		for(auto k = 0u; k < N; k++){
			double val = param_samp[k].paramval[th] - av;
			av2 += val*val;
		}		
		vec[i] = av2/N;
	}
	
	return vec;
}


/// Prints current proposal infomation 
string ParamProp::print_proposal_information(const bool brief) const
{
	stringstream ss;
	ss.precision(2);
	
	ss << "Diagnostics for MCMC proposals:" << endl;
	if(details.mode == PMCMC_INF){
		ss << "Self  -   Acceptance: "  << per(self.ac_rate) << endl;
	}
	
	for(auto &mv : mvn){
		ss << mv.name << " -   Acceptance: "  << per(mv.ac_rate) <<  "    Size: " << mv.size << endl;
	}
	
	for(auto &mt : mean_time){
		ss << "Mean Time: ";
		for(auto th : mt.param_mean) ss << model.param[th].name << " "; 
		ss << " <> ";
		for(auto th : mt.param_mean_rev) ss << model.param[th].name << " ";
		ss << " -   Acceptance: " <<  per(mt.ac_rate) << "    Size: " << mt.size << endl;
	}
	
	for(auto &rn : neighbour){
		ss << "Spline Neighbour: ";
		ss << model.param[rn.param1].name << " " <<  model.param[rn.param2].name; 
		ss << " -   Acceptance: " <<  per(rn.ac_rate) << "    Size: " << rn.size << endl;
	}
	
	for(auto &rn : joint){
		ss << "Spline Joint: ";
		for(auto i = 0u; i < rn.var_list.size(); i++){
			ss << model.param[rn.var_list[i]].name << " ";
			if(i > 3){ ss << "..."; break;}
		}
		ss << " ";
		switch(rn.type){
			case UP_DOWN: ss << "UpDown "; break;
			case SINE: ss << "Sine " << rn.sinenum; break;
		}
		ss << " -   Acceptance: " <<  per(rn.ac_rate) << "    Size: " << rn.size << endl;
	}
	
	for(auto &ca : covar_area){
		ss << "Covar Area: ";
		ss << model.param[model.covariate_param[ca.covar_ref]].name; 
		ss << " -   Acceptance: " <<  per(ca.ac_rate) << "    Size: " << ca.size << endl;
	}

	if(brief == true){
		cout << "Fixed Tree: #" << fixedtree.size();
		
		vector <double> vec; for(auto &ft : fixedtree) vec.push_back(ft.ac_rate);
		auto stat = output.get_statistic(vec);
		cout << " -   Acceptance: " << per(stat.mean) << " (" << per(stat.CImin) << " - " << per(stat.CImax) << ") "; 
		
		vector <double> vecb; for(auto &ft : fixedtree) vecb.push_back(ft.sim_frac);
		stat = output.get_statistic(vecb);
		cout << "   Sim Fraction: " <<  per(stat.mean) << " (" << per(stat.CImin) << " - " << per(stat.CImax) << ")" << endl;
		
		cout << "Slice Time: #" << slicetime.size();
		
		vector <double> vecc; for(auto &st : slicetime) vecc.push_back(st.ac_rate);
		stat = output.get_statistic(vecc);
		cout << " -   Acceptance: " << per(stat.mean) << " (" << per(stat.CImin) << " - " << per(stat.CImax) << ") "; 
		
		vector <double> vecd; for(auto &st : slicetime) vecd.push_back(st.sim_frac);
		stat = output.get_statistic(vecd);
		cout << "   Sim Fraction: " <<  per(stat.mean) << " (" << per(stat.CImin) << " - " << per(stat.CImax) << ")" << endl;
	}
	else{
		for(auto &ft : fixedtree){
			ss << "Fixed Tree: n=" << ft.n;
			ss << " -   Acceptance: " << per(ft.ac_rate) << "    Sim Fraction: " << per(ft.sim_frac) << endl;
		}
	
		for(auto &st : slicetime){
			ss << "Slice Time: ti=" << st.sett_i << " tf=" << st.sett_f;
			ss << " -   Acceptance " << per(st.ac_rate) << "    Sim Fraction: " << per(st.sim_frac) << endl;
		}
	}
	
	ss << endl;
	
	return ss.str();
}


/// The Metropolis-Hastings proposal
Status Self::MH(double al, double &invT, ParamUpdate pup)
{
	if(al == -1) nbo++;
	ntr++;
	if(ran() < al){
		if(pup == FAST_UPDATE) invT *= fac_up_invT_fast;
		if(pup == SLOW_UPDATE) invT *= fac_up_invT;
		//if(invT > 1) invT = 1;
		nac++; 
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) invT *= fac_down_invT_fast;
	if(pup == SLOW_UPDATE) invT *= fac_down_invT;
	return FAIL;
}


/// The Metropolis-Hastings proposal
Status MeanTime::MH(double al, ParamUpdate pup)
{
	if(al == -1) nbo++; 
	ntr++;
	if(ran() < al){
		if(pup == FAST_UPDATE) size *= fac_up_fast;
		if(pup == SLOW_UPDATE) size *= fac_up;
		nac++;
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) size *= fac_down_fast;
	if(pup == SLOW_UPDATE) size *= fac_down;
	return FAIL;
}


/// The mean time proposal
Status  MeanTime::propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	param_prop = paramval;
	
	auto dt = normal_sample(0,size);

	for(auto th : param_mean) param_prop[th] += dt;
	for(auto th : param_mean_rev) param_prop[th] -= dt;
	
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// The Metropolis-Hastings proposal
Status Neighbour::MH(double al, ParamUpdate pup)
{
	if(al == -1) nbo++;
	ntr++;
	if(ran() < al){
		if(pup == FAST_UPDATE) size *= fac_up_fast;
		if(pup == SLOW_UPDATE) size *= fac_up;
		nac++; 
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) size *= fac_down_fast;
	if(pup == SLOW_UPDATE) size *= fac_down;

	return FAIL;
}


/// The neighbouring points on a spline proposal
Status Neighbour::propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	param_prop = paramval;
	
	auto fac = exp(normal_sample(0,size));

	param_prop[param1] *= fac;
	
	param_prop[param2] /= fac;

	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}



/// The Metropolis-Hastings proposal
Status Joint::MH(double al, ParamUpdate pup)
{
	if(al == -1) nbo++;
	ntr++;
	if(ran() < al){ 
		if(pup == FAST_UPDATE) size *= fac_up_fast;
		if(pup == SLOW_UPDATE) size *= fac_up;
		nac++; 
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) size *= fac_down_fast;
	if(pup == SLOW_UPDATE) size *= fac_down;

	return FAIL;
}


/// The joint proposals on a spline
Status Joint::propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model)
{
	param_prop = paramval;
	
	auto d = normal_sample(0,size);
	
	switch(type){
	case UP_DOWN:
		for(auto th : var_list) param_prop[th] += d;
		break;
		
	case SINE:
		double num = sinenum;
		auto jmax = var_list.size();
		double phi = 2*M_PI*ran();
		
		if(num == 0) num = 0.5;
		auto wavelenth = (jmax-1)/num;
		
		for(auto j = 0u; j < jmax; j++){
			param_prop[var_list[j]] += d*sin(phi + 2*M_PI*j/wavelenth);
		}
		break;
	}

	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// The Metropolis-Hastings proposal
Status CovarArea::MH(double al, ParamUpdate pup)
{
	ntr++;
	if(ran() < al){ 
		if(pup == FAST_UPDATE) size *= fac_up_fast;
		if(pup == SLOW_UPDATE) size *= fac_up;
		nac++; 
		return SUCCESS;
	}
	if(pup == FAST_UPDATE) size *= fac_down_fast;
	if(pup == SLOW_UPDATE) size *= fac_down;

	return FAIL;
}


/// The joint proposals on a spline
Status CovarArea::propose(vector <double> &param_prop, const vector <double> &paramval, const Model &model, const Data &data)
{
	param_prop = paramval;
	
	auto th = model.covariate_param[covar_ref];
	
	auto param_st = param_prop[th];
	param_prop[th] += normal_sample(0,size);
	
	auto narea = data.narea;
	
	auto sum = 0.0;
	for(auto c = 0u; c < narea; c++){
		auto fac = exp((param_st-param_prop[th])*data.covar[covar_ref].value[c][0]);
		sum += log(fac);
	}

	for(auto c = 0u; c < narea; c++){
		auto fac = exp((param_st-param_prop[th])*data.covar[covar_ref].value[c][0] - sum/narea);
		param_prop[model.area_effect_param_list[c]] *= fac;
	}
	
	if(model.inbounds(param_prop) == false) return FAIL;
	else return SUCCESS;
}


/// Outputs diagnostic information
void ParamProp::diagnostics() const
{		
	if(mpi.core == 0){
		string filefull = details.output_directory+"/Diagnostics/MCMC_proposals.txt";
		ofstream dia(filefull);
		if(!dia) emsg("Cannot open the file '"+filefull+"'");
		
		dia << print_proposal_information(false);
	}	
}
