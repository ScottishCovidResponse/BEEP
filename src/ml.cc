/// Implements a maximum likelihood approach
//REMINDER:m likelihood gives different answers
//convergence

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "ml.hh"
#include "mpi.hh"
#include "matrix.hh"
#include "output.hh"

using namespace std;

/// Initilaises the ML class
ML::ML(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), data(data), model(model), output(output), obsmodel(obsmodel), mpi(mpi)
{
	find_parameters();
	
	inputs.find_algorithm(algorithm,npart,G,cpu_time,P,nsample_final,mpi.core,mpi.ncore,nvar);
	
	invT = inputs.find_double("invT",1);  
	
	if(algorithm == ML_CMAES){
		if(npart%mpi.ncore != 0) emsgroot("'nparticle' must be a multiple of the number of cores.");
		npart_per_core = (unsigned int)(npart/mpi.ncore);
	}
	
	param_per_core = (unsigned int)((nvar+mpi.ncore)/mpi.ncore);
	
	if(mpi.core == 0 && false){
		auto per_error = obsmodel.mean_percentage_error(invT);
		cout << "For invT='" << invT << "', the mean percentage error in the observations is " <<  prec(per_error,3) << "%." << endl;
	}
	
	obs_slice = obsmodel.generate_obs_slice();
	state.initialise_approx();
	
	ac = DOUBLE; 
	//ac = FLOAT;
}


/// Finds all the free parameters (removes all fixed and some non-independent Dirichlet)
void ML::find_parameters()
{
	for(auto th : model.param_not_fixed){
		bool flag = false;
		for(const auto &dir : model.dirichlet){
			if(dir.param[0].th == th) flag = true; 
		}
	
		if(flag == false) var.push_back(th);
	}
	nvar = var.size();
}


/// Runs the inference algorithm
void ML::run()
{	
	timer[TIME_ALG].start();
 
	//calculate_mimimum();  emsg("Done");
	//calculate_heatmap(); emsg("Done");

	switch(algorithm){
		case ML_GD: gradient_decent(); break;
		case ML_CMAES: cmaes(); break;
	}
	timer[TIME_ALG].stop();
}


/// Performs the gradient descent algorithm
void ML::gradient_decent()
{
	auto param = model.sample_from_prior();  
	
	auto pi =	calculate_point_info(param);

	auto eta = LARGE;                                              // Sets the initial value for the step size
	if(mpi.core == 0){
		for(auto i = 0u; i < nvar; i++){
			auto na = model.prior_normal_approx(var[i]);
			
			if(eta*pi.grad[i] > na.mean) eta = na.mean/pi.grad[i];
		}
	}
	
	auto loop = 0u;
	do{
		if(mpi.core == 0) cout << "Iteration: " << loop << "  Likelihood: " << pi.value << "  Jump size: " << eta << endl;
			
		auto pi_store = pi;
		auto param_store = param;
		
		auto fail = false;
		if(mpi.core == 0){
			for(auto i = 0u; i < nvar; i++){
				auto th = var[i];
			
				param[th] += eta*pi.grad[i];
			}
			
			if(model.inbounds(param) != true) fail = true;
		}
		
		mpi.bcast(fail);
		
		if(fail == false){
			pi = calculate_point_info(param);
			if(pi.value < pi_store.value) fail = true;
			mpi.bcast(fail);
		}
	
		if(mpi.core == 0){
			if(fail == true){
				eta *= 0.5; 
				pi = pi_store; 
				param = param_store;
			}
			else{
				eta *= 1.5;
			}
		}
		loop++;
	}while(loop < 100);
	
	vector <Particle> particle_store; 
	
	if(mpi.core == 0){
		state.EF = -0.5*state.likelihood_approx(param,obs_slice,invT,ac);		
		particle_store.push_back(state.create_particle(UNSET));
	}
	
	output.generate_graphs(particle_store,invT);		           // Outputs the results
}


// Used to order particles by EF
bool ParamSample_ord(ParamSample p1,ParamSample p2)                      
{ return (p1.EF < p2.EF); };  


/// Performs the Covariance matrix adaptation evolution strategy algorithm
void ML::cmaes()
{
	const auto c_c_limit = 0.9;              // These values are used to limit the system in the case of few parameters
	const auto c_sigma_limit = 0.9;
	const auto c1_limit = 0.3;
	const auto c_mu_limit = 0.3;
	const auto c_s_limit = 0.3;

	auto sigma = 1.0;
	vector <double> p_sigma(nvar), p_c(nvar);
	for(auto i = 0u; i < nvar; i++){
		p_sigma[i] = 0.0; p_c[i] = 0.0;
	}
	
	auto mean = model.sample_from_prior();
	
	vector < vector <double> > C;
	C.resize(nvar); 
	for(auto j = 0u; j < nvar; j++){
		C[j].resize(nvar);
		for(auto i = 0u; i < nvar; i++) C[j][i] = 0;
	}

	for(auto i = 0u; i < nvar; i++){
		auto th = var[i];
		auto na = model.prior_normal_approx(th);
		mean[th] = na.mean;
		C[i][i] = na.var;
	}

	/* Calculates quantities used later */
	auto mu = (unsigned int)(npart/2);
	vector <double> w(mu);
	//for(auto i = 0u; i < mu; i++) w[i] = mu-i;
	for(auto i = 0u; i < mu; i++) w[i] = log(mu+0.5) - log(i+1);

	auto sum = 0.0; for(auto i = 0u; i < mu; i++) sum += w[i];
	for(auto i = 0u; i < mu; i++) w[i] /= sum;
	
	auto su = 0.0; for(auto i = 0u; i < mu; i++) su += w[i]*w[i];
	auto mu_w = 1.0/su;

	vector <ParamSample> ps_per_core(npart_per_core);
		
	if(mpi.core == 0) cout << endl << "Start Inference..." << endl;
	
	vector <double> EFbest_store;
	
	auto g = 0u;
	do{
		mpi.bcast(mean);
		mpi.bcast(C);
		mpi.bcast(sigma);

		Generation gen;
	
		/* Samples parameters */
		timer[TIME_GENERATE_SAMPLES].start();
		generate_samples(ps_per_core,mean,C,sigma,gen);
		timer[TIME_GENERATE_SAMPLES].stop();
		
		timer[TIME_CMAES].start();
		auto ps = mpi.gather_psamp(ps_per_core);

		if(mpi.core == 0){
			sort(ps.begin(),ps.end(),ParamSample_ord);   
				
			auto EFav = 0.0; for(auto i = 0u; i < mu; i++) EFav += w[i]*ps[i].EF;
	
			//auto Prav = 0.0; for(auto i = 0u; i < mu; i++) Prav += w[i]*model.prior(ps[i].paramval);
	
			cout << "Generation " << g << " - Best log(Post. prob.): " << -0.5*ps[0].EF <<  "   Sigma: " << sigma << endl;
			//cout << "Generation " << g << " - EF: " << ps[0].EF << "    EF average: " << EFav << "   Sigma: " << sigma << endl;
			//cout << "Generation " << g << " - EF: " << ps[0].EF << endl;
			
			EFbest_store.push_back(ps[0].EF);
			
			gen.EFcut = EFav;
		
			auto mean_store = mean;

			/* Updates the mean */
			for(auto i = 0u; i < nvar; i++){
				auto th = var[i];
				auto sum = 0.0; for(auto k = 0u; k < mu; k++) sum += ps[k].paramval[th]*w[k];
				mean[th] = sum; 
			}		
			
			/* Update p_sigma */
			vector <double> dif(nvar);
			for(auto i = 0u; i < nvar; i++){
				auto th = var[i];
				dif[i] = mean[th] - mean_store[th]; 
			}
			
			auto Cinvsqrt = invert_matrix_square_root(C);
			
			auto vec = matrix_mult(Cinvsqrt,dif);
			
			auto c_sigma = 3.0/nvar; if(c_sigma > c_sigma_limit) c_sigma = c_sigma_limit;
			
			vector <double> p_sigma_new(nvar);
			for(auto i = 0u; i < nvar; i++){
				p_sigma_new[i] = (1-c_sigma)*p_sigma[i] + sqrt(1-(1-c_sigma)*(1-c_sigma))*sqrt(mu_w)*vec[i]/sigma;
			}
			p_sigma = p_sigma_new;
			
			/* Update p_c */
			auto alpha = 1.5;
			
			auto p_sigma_mag = 0.0; for(auto i = 0u; i < nvar; i++) p_sigma_mag += p_sigma[i]*p_sigma[i];
			p_sigma_mag = sqrt(p_sigma_mag);
			
			auto ind = 1u; if(p_sigma_mag > alpha*sqrt(nvar)) ind = 0;
				
			auto c_c = 4.0/nvar; if(c_c > c_c_limit) c_c = c_c_limit;
				
			vector <double> p_c_new(nvar);
			for(auto i = 0u; i < nvar; i++){
				p_c_new[i] = (1-c_c)*p_c[i] + ind*sqrt(1-(1-c_c)*(1-c_c))*sqrt(mu_w)*dif[i]/sigma;
			}
			p_c = p_c_new;
			
			/* Update C */
			auto c1 = 2.0/(nvar*nvar); if(c1 > c1_limit) c1 = c1_limit;
			auto c_mu = mu_w/(nvar*nvar); if(c_mu > c_mu_limit) c_mu = c_mu_limit; if(c_mu > 1-c1) c_mu = 1-c1;
			auto c_s = (1-ind)*c1*c_c*(2-c_c); if(c_s > c_s_limit) c_s = c_s_limit;
			
			vector < vector <double> > C_new;
			C_new.resize(nvar);
			for(auto j = 0u; j < nvar; j++){
				auto thj = var[j];
				C_new[j].resize(nvar);
				for(auto i = 0u; i < nvar; i++){
					auto thi = var[i];
					
					auto sum = 0.0;
					for(auto k = 0u; k < mu; k++){
						sum += w[k]*(ps[k].paramval[thi]-mean_store[thi])*(ps[k].paramval[thj]-mean_store[thj])/(sigma*sigma);
					}
					
					C_new[j][i] = (1-c1-c_mu+c_s)*C[j][i] + c1*p_c[i]*p_c[j] + c_mu*sum;
				}
			}
			C = C_new;
			
			/* Update sigma */
			//auto d_sigma = 1;
			auto d_sigma = 1 + 2*max(0,sqrt((mu_w-1.0)/(nvar+1.0))-1) + c_sigma;
		
			auto EN = sqrt(nvar)*(1-(1.0/(4*nvar)) + 1.0/(21*nvar*nvar));
			
			//EN *= 2;
			auto fac = (c_sigma/d_sigma)*((p_sigma_mag/EN) - 1);
			if(fac > 0.2) fac = 0.2;
			sigma *= exp(fac);
		}
		
		mpi.exchange_samples(gen);	                            // Exchanges parameter samples across MPI cores
	
		output.set_generation_time(gen);                        // Sets the CPU time for the generation
		
		generation.push_back(gen);
		timer[TIME_CMAES].stop();
		
		g++;
	}while(terminate_generation(g,EFbest_store) == false);

	if(mpi.core == 0) cout << "Generating posterior samples..." << endl;

	timer[TIME_POSTERIOR_SAMPLE].start();
	vector <Particle> particle_store;
	auto num_per_core = nsample_final/mpi.ncore;
	
	mpi.bcast(mean);
	mpi.bcast(C);

	timer[TIME_SCALE_COVARIANCE].start();
	//auto covar = calculate_covariance_martrix(mean);
	if(mpi.core == 0) cout << "Scale covariance..." << endl;
	sigma = scale_covariance_martrix(mean,C);
	timer[TIME_SCALE_COVARIANCE].stop();
			
	if(mpi.core == 0) cout << "SMC sampler..." << endl;	
	auto psamp = sample_param(mean,C,sigma,num_per_core);
	for(auto i = 0u; i < num_per_core; i++){               // Generates the final posterior samples
		if(details.stochastic == true){
			particle_store.push_back(state.posterior_particle_sample(psamp[i],obs_slice,P,invT));
		}
		else{
			state.simulate(psamp[i]);
			particle_store.push_back(state.create_particle(UNSET));
		}
	}
	timer[TIME_POSTERIOR_SAMPLE].stop();

	if(mpi.core == 0) cout << "Output results..." << endl;

	output.generation_results(generation);                    // Generates pdf of graphs
	
	output.generate_graphs(particle_store,invT);		  
}


/// Generates samples based on a mean and covariance matrix
void ML::generate_samples(vector <ParamSample> &ps_per_core, const vector <double> &mean, const vector < vector <double> > &C, const double sigma, Generation &gen)
{
	auto psamp = sample_param(mean,C,sigma,ps_per_core.size());

	timer[TIME_EF_CALCULATE].start();
	for(auto i = 0u; i < ps_per_core.size(); i++){
		const auto &param = psamp[i];
		ps_per_core[i].paramval = psamp[i];
		
		if(details.stochastic == true){
			auto L = state.likelihood_approx(param,obs_slice,invT,ac);

			auto Pr = model.prior(param);
			ps_per_core[i].EF = -(L + Pr);
 			if(ps_per_core[i].EF < 0) emsgEC("ML",44);
		}
		else{
			state.simulate(param);
			ps_per_core[i].EF = invT*state.EF - 0.5*state.Pr;
		}
		
		gen.param_samp.push_back(state.create_param_sample(0));
		gen.EF_datatable.push_back(obsmodel.get_EF_datatable(&state));
	}
	timer[TIME_EF_CALCULATE].stop();
}


/// Samples a set of parameters from a given set of parameters
vector < vector <double> > ML::sample_param(const vector <double> &mean, const vector < vector <double> > &C, const double sigma, const unsigned int num) const
{
	MVN mv("MVN",var,sigma,ALL_PARAM,MULTIPLE);
	mv.set_covariance_matrix(C);
	
	auto nparam = model.param.size();
	
	vector < vector <double> > psamp;
	
	for(auto i = 0u; i < num; i++){
		vector <double> param;
		do{
			param = mv.propose(mean);
				
			for(auto th = 0u; th < nparam; th++){   // This reflects values to ensure inside prior
				const auto &pa = model.param[th];
				if(pa.priortype == UNIFORM_PRIOR){
					auto val = param[th];
					auto val1 = pa.val1, val2 = pa.val2;
			
					bool reflect;
					do{
						reflect = false;
						if(val > val2){
							val = val2 - (val-val2);
							reflect = true;
						}
						else{
							if(val < val1){
								val = val1 + (val1-val);
								reflect = true;
							}
						}
					}while(reflect == true);
					param[th] = val;
				}
				
				switch(pa.priortype){   // Unsures value is positive
					case DIRICHLET_PRIOR: case DIRICHLET_ALPHA_PRIOR: case MDIRICHLET_PRIOR: case EXP_PRIOR: case DIRICHLET_FLAT_PRIOR:
						{
							auto val = param[th];
							if(val < 0) val = -val;
							param[th] = val;
						}
						break;
						
					default: break;
				}
			}
			
			if(model.inbounds(param) == false) emsgEC("ML",10);
		}while(model.inbounds(param) == false);


		for(const auto &dir : model.dirichlet){                       // Fills in dependent diriclet values
			auto sum = 0.0;				
			for(auto j = 1u; j < dir.param.size(); j++){
				sum += param[dir.param[j].th];
			}
		
			if(sum > 1){
				auto summax = 1 - ran()*0.5/dir.param.size();
				auto fac = summax/sum;
				for(auto j = 1u; j < dir.param.size(); j++){
					param[dir.param[j].th] *= fac;
				}
				sum = summax;						
			}
			
			if(sum > 1) emsgEC("ML",13);
			param[dir.param[0].th] = 1-sum;
		}

		psamp.push_back(param);
	}
	
	return psamp;
}


/// Determines when generations are stopped
bool ML::terminate_generation(unsigned int g, const vector <double> &EFbest_store) const 
{
	if(g == G) return true;  // Iterate until specifies generation
	
	if(cpu_time != UNSET){ 	 // Iterate until time limit
		auto time_av = mpi.average((clock() - details.time_start)/(60.0*CLOCKS_PER_SEC));
 
		if(time_av > cpu_time){
			if(mpi.core == 0)	cout << "Maximum execution time reached." << endl << flush;
			return true;
		}
	}
	
	if(G == ITERATE_GENERATION){          // Iterate until termination conidition 
		auto term = false;
		if(mpi.core == 0){
			auto ng = EFbest_store.size();
			if(ng >= ML_GENERATION_TERM_COND){
				auto av = 0.0;
				for(auto i = ng-ML_GENERATION_TERM_COND; i < ng; i++) av += EFbest_store.size();
				av /= ML_GENERATION_TERM_COND;
				
				auto tol = av/1000;
				auto value = EFbest_store[ng-1];
				auto i = ng-ML_GENERATION_TERM_COND; 
				while(i < ng && EFbest_store[i] > value-tol && EFbest_store[i] < value+tol) i++;
				if(i == ng) term = true;
			}
		}
		mpi.bcast(term);
		
		return term;
	}
	
	return false;
}


/// Calculates the gradient in the likelihood for the non-fixed parameters
PointInfo ML::calculate_point_info(vector <double> param) 
{
	const auto delta = 0.1;
	
	mpi.bcast(param);
	
	vector <double> val(param_per_core), dif(param_per_core);
	for(auto k = 0u; k < param_per_core; k++){
		auto p = mpi.core*param_per_core + k;
	
		if(p < nvar){
			auto th = var[p];
			auto param_st = param[th];
			param[th] += delta; 
			if(model.inbounds(param) == true){
				val[k] = state.likelihood_approx(param,obs_slice,invT,ac) + model.prior(param);
				dif[k] = delta;
			}
			else{
				param[th] -= 2*delta;
				if(model.inbounds(param) == false) emsgEC("ML",1);
				val[k] = state.likelihood_approx(param,obs_slice,invT,ac) + model.prior(param);
				dif[k] = delta;
			}
			param[th] = param_st;
		}
		else{
			if(p == nvar){
				val[k] = state.likelihood_approx(param,obs_slice,invT,ac) + model.prior(param);
				dif[k] = 0;
			}
			else{
				val[k] = UNSET;
				dif[k] = UNSET;
			}
		}
	}
	
	auto val_tot = mpi.gather(val);
	auto dif_tot = mpi.gather(dif);
	
	PointInfo pi;
	if(mpi.core == 0){
		pi.value = val_tot[nvar];
		pi.grad.resize(nvar);
		for(auto i = 0u; i < nvar; i++) pi.grad[i] = (val_tot[i] - val_tot[nvar])/dif_tot[i];
	}
		
	return pi;
}
	
	
/// Calculates the minimum value for the error function
void ML::calculate_mimimum() 
{
	auto param = model.simulated_values();

	double EF;
	if(details.stochastic == true){
		auto time = clock();
		EF = -0.5*state.likelihood_approx(param,obs_slice,invT,ac);
		cout.precision(10);		
		cout << EF << " " << double(clock()-time)/CLOCKS_PER_SEC << "time" << endl;
		cout << double(timer[TIME_DETERMINANT].val)/CLOCKS_PER_SEC << " j" << endl;
		cout << double(timer[TIME_TEMP1].val)/CLOCKS_PER_SEC << " j" << endl;
		cout << double(timer[TIME_TEMP2].val)/CLOCKS_PER_SEC << " j" << endl;
		cout << double(timer[TIME_TEMP3].val)/CLOCKS_PER_SEC << " j" << endl;
		emsg("d");
	}
	else{
		state.simulate(param);
		EF = invT*state.EF;
	}
		
	cout << EF << " Minimum" << endl;

	vector <Particle> particle_store; 
	
	if(mpi.core == 0){
		particle_store.push_back(state.create_particle(UNSET));
	}
	
	output.generate_graphs(particle_store,invT);		           // Outputs the results
}
	

/// Creates a heatmap for two variables (used the heatmap.html too to view)
void ML::calculate_heatmap()
{	
	auto th1 = 2, th2 = 3;

	auto param = model.simulated_values();

	if(mpi.core == 0){
		for(auto th = 0u; th < model.param.size(); th++) cout << th << " " <<  model.param[th].name << " param " << endl;
	}
	
	auto min1 = 1.0, max1 = 10.0;
	auto min2 = 0.5, max2 = 4.0;
	auto MX = 100u, MY = 100u;
	
	auto M = MX*MY;
	auto M_per_core = (unsigned int)((M+mpi.ncore-1)/mpi.ncore);
	
	vector <double> answer(M_per_core);
	for(auto i = 0u; i < M_per_core; i++){
		auto k = mpi.core*M_per_core+i;
		if(k < M){
			param[th1] = min1+(k%MX + 0.5)*(max1-min1)/MX;
			param[th2] = (k/MX + 0.5)*(max2-min2)/MY;
			answer[i] = -(state.likelihood_approx(param,obs_slice,invT,ac) + 0*model.prior(param));
		}
	}
	
	auto answer_tot = mpi.gather(answer);
	
	auto imin = 0u, jmin = 0u;
	auto valmin = 1000000.0;
	if(mpi.core == 0){
		ofstream heatmap(details.output_directory+"/heatmap.txt");
		heatmap << model.param[th1].name << "," << min1 << "," << max1 << "," << model.param[th2].name << "," << min2 << "," << max2 << endl;
		
		for(auto j = 0u; j < MY; j++){
			for(auto i = 0u; i < MX; i++){
				auto val = answer_tot[j*MX+i];
				if(val < valmin){ valmin = val; imin  = i; jmin = j;}
				heatmap << val; if(i < MX-1) heatmap << ",";
			}
			heatmap << endl;
		}
	}
	
	
	vector <Particle> particle_store; 
	
	if(mpi.core == 0){
		
		param[th1] = min1+(imin + 0.5)*(max1-min1)/MX;
		auto GT = param[th1]*0.5; // Generation time
		param[th2] = exp(((jmin + 0.5)*(max2-min2)/MY)*GT);
		
		
		auto Li = state.likelihood_approx(param,obs_slice,invT,ac);
		
		cout << model.param[th1].name << ": " << param[th1] << "  "
					<< model.param[th2].name << ": " << param[th2] << " minimum" << endl;;
	
	
		param[th1] = 6; param[th2] = 2;
		
		auto Li2 = state.likelihood_approx(param,obs_slice,invT,ac);
		
		cout << Li << " " << Li2 << " Li" << endl;
		
		particle_store.push_back(state.create_particle(UNSET));
	}
	
	output.generate_graphs(particle_store,invT);		           // Outputs the results

	mpi.barrier();
}


/// Calculates the Hessian matrix and uses this to estimate the covariance matrix
vector < vector <double> > ML::calculate_covariance_martrix(const vector <double> &mean)
{
	auto H = calculate_hessian(mean,FULL);

	return invert_matrix(H);
}


/// Works out the scaling factor for the covaiance matrix
double ML::scale_covariance_martrix(const vector <double> &mean, const vector < vector <double> > &covar)
{
	auto H_diag = calculate_hessian(mean,DIAG);

	double ratio;
	if(mpi.core == 0){
		auto H = invert_matrix(covar);
	
		if(false){
			print_matrix("H",H);
			print_matrix("H_diag",H_diag);
		
			for(auto i = 0u; i < nvar; i++){
				cout << H[i][i] << " " << H_diag[i][i] << " " << H_diag[i][i]/H[i][i] << " compare" << endl;
			}
		}
		
		auto sum_log_ratio = 0.0;
		for(auto i = 0u; i < nvar; i++){
			auto ra = H[i][i]/H_diag[i][i];
			if(ra < 0.0001) ra = 0.0001;
			if(ra > 10000) ra = 10000;
			sum_log_ratio += log(ra);
		}
		sum_log_ratio /= nvar;
		
		ratio = exp(0.5*sum_log_ratio);
	}
	
	mpi.bcast(ratio);

	return ratio;
}


/// Calculates the Hessian matrix
vector < vector <double> > ML::calculate_hessian(const vector <double> &mean, MatrixType mattype)
{
	vector <double> d(nvar);                             // Determines grid size to use
	
	for(auto i = 0u; i < nvar; i++){  
		auto approx = model.prior_normal_approx(var[i]);
		d[i] = 0.01*sqrt(approx.var);
		if(d[i] == 0) emsgEC("ML",43);
	}
	
	// Should check not overlapping boundaries...
	
	vector <MatrixElement> matrix_element;
	for(auto j = 0u; j < nvar; j++){
		auto imax = nvar;
		if(mattype == DIAG) imax = j+1; 
		for(auto i = j; i < imax; i++){
			MatrixElement el; el.j = j; el.i = i;
			matrix_element.push_back(el);
		}
	}
		
	auto element_per_core = (unsigned int)((matrix_element.size()+mpi.ncore-1)/mpi.ncore);
	
	vector <double> ur, dr, ul, dl;
	
	for(auto k = 0u; k < element_per_core; k++){
		auto l = mpi.core*element_per_core + k;
		if(l < matrix_element.size()){	
			auto &el = matrix_element[l];
	
			auto vari = el.i, varj = el.i;
			
			auto param = mean; param[var[vari]] += d[vari]; param[var[varj]] += d[varj];
			model.calculate_dirichlet_missing(param);
			ur.push_back(state.likelihood_approx(param,obs_slice,invT,DOUBLE) + model.prior(param));
			
			param = mean; param[var[vari]] += d[vari]; param[var[varj]] -= d[varj];
			model.calculate_dirichlet_missing(param);
			dr.push_back(state.likelihood_approx(param,obs_slice,invT,DOUBLE) + model.prior(param));
			
			param = mean; param[var[vari]] -= d[vari]; param[var[varj]] += d[varj];
			model.calculate_dirichlet_missing(param);
			ul.push_back(state.likelihood_approx(param,obs_slice,invT,DOUBLE) + model.prior(param));
			
			param = mean; param[var[vari]] -= d[vari]; param[var[varj]] -= d[varj];
			model.calculate_dirichlet_missing(param);
			dl.push_back(state.likelihood_approx(param,obs_slice,invT,DOUBLE) + model.prior(param));
		}
		else{
			ur.push_back(UNSET); dr.push_back(UNSET); ul.push_back(UNSET); dl.push_back(UNSET);
		}
	}

	auto ur_tot = mpi.gather(ur);
	auto dr_tot = mpi.gather(dr);
	auto ul_tot = mpi.gather(ul);
	auto dl_tot = mpi.gather(dl);
	
	vector < vector <double> > H;
	H.resize(nvar);
	for(auto j = 0u; j < nvar; j++){
		H[j].resize(nvar); for(auto i = 0u; i < nvar; i++) H[j][i] = UNSET;
	}
		
	if(mpi.core == 0){	
		for(auto k = 0u; k < matrix_element.size(); k++){
			const auto &el = matrix_element[k];
			H[el.j][el.i] = -(ur_tot[k] + dl_tot[k] - ul_tot[k] - dr_tot[k])/(4*d[el.i]*d[el.j]);
			H[el.i][el.j] = H[el.j][el.i];
		}
		
		for(auto j = 0u; j < nvar; j++){
			switch(mattype){
				case FULL:
					for(auto i = 0u; i < nvar; i++){
						if(H[j][i] == UNSET || std::isnan(H[j][i])) emsgEC("ML",32);
					}
					break;
				
				case DIAG:
					if(H[j][j] == UNSET || std::isnan(H[j][j])) emsgEC("ML",33);
					break;
			}
		}		
	}		

	mpi.bcast(H);

	return H;
}
