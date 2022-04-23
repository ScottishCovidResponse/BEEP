/// These are algorithms to check state is working correctly

#include <math.h>  

using namespace std;

#include "state.hh"
#include "matrix.hh"

/// Checks quantities in the state are correct
void State::check(const unsigned int checknum)
{
	if(model.inbounds(paramval) == false) emsg("Parameters not in bounds");

	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto c = 0u; c < data.narea; c++){
			for(auto tr = 0u; tr < model.trans.size(); tr++){
				auto from = model.trans[tr].from;
			
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){	 
					auto num = transnum[sett][c][tr][dp];
					if(num < 0) emsgEC("State_check",15);
					
					auto p = pop[sett][c][from][dp];
					if(p <= 0){ if(num != 0){ cout << p<< " " << checknum << "checknum" << endl; emsgEC("State_check",74);}}
					else{	
						if(model.trans[tr].inf == TRANS_NOTINFECTION){
							if(transrate[tr][dp] == 0){
								if(num != 0 && model.trans[tr].inf == TRANS_NOTINFECTION) emsgEC("State_check",75);
							}
						}
					}
				}
			}
		}
	}		
	
	for(auto sett = 0u; sett < details.ndivision-1; sett++){       // Checks pop is correct
		auto popst = pop[sett+1];
		update_pop(sett);
		democat_change_pop_adjust(sett+1);
			
		for(auto c = 0u; c < data.narea; c++){
			for(auto co = 0u; co < model.comp.size(); co++){
				for(auto dp = 0u; dp < data.ndemocatpos; dp++){
					if(popst[c][co][dp] != pop[sett+1][c][co][dp]){
						emsgEC("State_check",1);
					}
				}
			}
		}
	}	
	
	set_Imap(1);  // Checks to make sure Imap is OK
	
	for(auto sett = 0u; sett < details.ndivision; sett++){   // Check I map consistency with pop
		set_I_from_pop(sett,true);
	}
	
	if(details.mode == ABC_SMC || details.mode == ABC_MBP){ 
		auto EF_temp = -2*obsmodel.calculate(this);
		auto dd = EF - EF_temp;
		if(dd*dd > TINY) emsgEC("State_check",2);
	}
	
	auto dd = Pr - model.prior(paramval); if(sqrt(dd*dd) > TINY) emsgEC("State_check",59);
	
	auto paramv_dir = model.dirichlet_correct(paramval);
	
	auto af = model.create_areafactor(paramv_dir);                   // Checks areafactor
	for(auto t = 0u; t < details.period; t++){
		for(auto c = 0u; c < data.narea; c++){
			if(af[t][c] != areafactor[t][c]){
				emsgEC("State_check",88);
			}
		}
	}
	
	auto sus = model.create_susceptibility(paramv_dir);              // Checks susceptibility
	for(auto dp = 0u; dp < data.ndemocatpos; dp++){
		if(sus[dp] != susceptibility[dp]) emsgEC("State_check",89);
	}

	for(auto sp = 0u; sp < model.spline.size(); sp++){               // Checks disc_spline
		auto spl = model.create_disc_spline(sp,paramval);  
		for(auto se = 0u; se < details.ndivision; se++){
			if(spl[se] != disc_spline[sp][se]) emsgEC("State_check",90);
		}
	}
}

/* The following functions refer to checking approximate methods */


/// Perform the integatated matrix producted (measures how close the distributions are
double State::integrated_matrix_prod(const vector <double> &mu1, const vector < vector <double> > &Minv1, const vector <double> &mu2, const vector < vector <double> > &Minv2) const
{
	auto Minv = matrix_add(Minv1,Minv2);
	auto M = invert_matrix(Minv);
	
	auto Minv_mu1 = matrix_mult(Minv1,mu1);
	auto Minv_mu2 = matrix_mult(Minv2,mu2);
	
	auto Minv_mu = vec_add(Minv_mu1,Minv_mu2);
	auto mu = matrix_mult(M,Minv_mu);
	
	auto mu_Minv_mu1 = vec_mult(mu1,Minv_mu1);
	auto mu_Minv_mu2 = vec_mult(mu2,Minv_mu2);
		
	auto mu_Minv_mu = vec_mult(mu,Minv_mu);

	return 0.5*(determinant_fast(M) + determinant_fast(Minv1) + determinant_fast(Minv2) - (mu_Minv_mu1 + mu_Minv_mu2 - mu_Minv_mu));
}


/// Checks that 
void State::integrated_matrix_prod_check() const
{
	vector <unsigned int> vec; vec.push_back(0); vec.push_back(1); 
	MVN mvn1("MVN",vec,1,ALL_PARAM,MULTIPLE), mvn2("MVN",vec,1,ALL_PARAM,MULTIPLE);
	
	const auto n = vec.size();
	vector <double> mu1(n), mu2(n);
	vector < vector <double> > M1, M2;
	
	M1.resize(n); M2.resize(n);
	for(auto i = 0u; i < n; i++){ M1[i].resize(n); M2[i].resize(n);}
		
	auto file = details.output_directory+"/scan.txt";
	ofstream op(file);
	
	auto sh1 = 0.0, sh2 = 0.0;
	
	for(auto scan = -0.5; scan < 0.5; scan += 0.01){	
		M1[0][0] = 1; M1[1][1] = 2; M1[0][1] = scan; M1[1][0] = M1[0][1]; 
		M2[0][0] = 1; M2[1][1] = 0.7; M2[0][1] = 0.03; M2[1][0] = M2[0][1]; 
	
		mu1[0] = 1; mu1[1] = 0;
		mu2[0] = 1.0; mu2[1] = 0.5;
			
		auto Minv1 = invert_matrix(M1);
		auto Minv2 = invert_matrix(M2);
		
		mvn1.set_covariance_matrix(M1);
		mvn2.set_covariance_matrix(M2);
		
		auto loopmax = 1000000u;
		auto prob = 0.0;
		for(auto loop = 0u; loop < loopmax; loop++){
			auto x = mvn1.propose(mu1);
			prob += exp(mvn2.probability(x,mu2))/loopmax;
		}
		if(sh1 == 0){ sh1 = log(prob); sh2 = integrated_matrix_prod(mu1,Minv1,mu2,Minv2);}
 		op << scan << " " << log(prob)-sh1 << " " << integrated_matrix_prod(mu1,Minv1,mu2,Minv2)-sh2 << endl; 
	}
	
	/*
	const n = 2;
	vector <double> mu1(n), mu2(n);
	vector < vector <double> > Minv1, Minv2;
	
	Minv1.resize(n); Minv1.resize(n);
	*/
	
}

/// Tests the gradients
void State::test_transmean_gradients(const unsigned int sett)
{
	cout << "Testing transmean_gradients" << endl;
	
	set_Imap_sett(sett);
	set_transmean(sett,0);
	set_transmean_gradients(sett);
	
	cout << "start" << endl;
	
	vector < vector <double> > dm_dp;
	dm_dp.resize(ntrans_approx);
	for(auto tr = 0u; tr < ntrans_approx; tr++) dm_dp[tr].resize(ncomp_approx);
	
	for(auto c = 0u; c < data.narea; c++) set_transmean(sett,c);
	auto tmean_store = transmean[sett];
	
	auto dif = 0.1;
	for(auto ta = 0u; ta < ntrans_approx; ta++){
		const auto &tap = trans_approx[ta];
		for(auto ca = 0u; ca < ncomp_approx; ca++){
			const auto &cap = comp_approx[ca];
			
			pop[sett][cap.c][cap.co][cap.dp] += dif;
							
			set_I_from_pop(sett,false);
			set_transmean(sett,tap.c);
							
			dm_dp[ta][ca] = (transmean[sett][tap.c][tap.tr][tap.dp] - tmean_store[tap.c][tap.tr][tap.dp])/dif;
						
			pop[sett][cap.c][cap.co][cap.dp] -= dif;
		}
	}
	
	for(auto ta = 0u; ta < ntrans_approx; ta++){
		cout << model.trans[trans_approx[ta].tr].name << " transition" << endl;
		for(auto ca = 0u; ca < ncomp_approx; ca++){
			if(dm_dp[ta][ca] != 0 || dtransmean_dp[ta][ca] != 0){
				cout << dm_dp[ta][ca] << " " << dtransmean_dp[ta][ca] << " compare" << endl;
			}
		}
		cout << endl;
	}
}


