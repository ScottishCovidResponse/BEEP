/// This routines check that Mbp is working correctly

#include <math.h>   

using namespace std;

#include "mbp.hh"

/// Checks that MBPs are correctly working by setting observations to simulated data
void Mbp::check_MBP()
{
	auto nparam = model.param.size();
	
	const unsigned int nsamp = 10;
	const unsigned int burnin = 100;
	vector <double> jump(nparam);
	for(auto th = 0u; th < nparam; th++) jump[th] = initial->paramval[th]/10;
	vector <double> ntr(nparam), nac(nparam);	
	for(auto th = 0u; th < nparam; th++){ ntr[th] = 0; nac[th] = 0;}
	
	EFcut = 0.08;
	initial->set_EF();
	initial->set_Pr();
	
	string file = details.output_directory+"/trace.txt";
	ofstream trace(file.c_str());		
	
	trace << "state";
	for(const auto& par : model.param) trace << "\t" << par.name; 
	trace << "\tEF"; 
	trace << "\tPri"; 
	trace << endl;
	
	vector <Sample> opsamp; 
	vector <ParamSample> psamp;
		
	cout << obsmodel.calculate(initial) << " initial error\n";
		
	for(auto samp = 0u; samp < nsamp; samp++){
		cout << samp << "samp\n";
		
		trace << samp; 
		for(const auto& pval : initial->paramval) trace << "\t" << pval;
		trace << "\t" << initial->EF; 
		trace << "\t" << initial->Pr; 
		trace << endl;
	
		ParamSample paramsamp;		
		paramsamp.paramval = initial->paramval;
		opsamp.push_back(initial->create_sample());
		psamp.push_back(paramsamp);
								
		for(auto th = 0u; th < model.param.size(); th++){
			//if(model.param[th].type ==GEOMIX_PARAM&& model.param[th].min != model.param[th].max){
			if(model.param[th].priortype != FIXED_PRIOR){
				
				auto param_propose = initial->paramval;
				param_propose[th] += normal_sample(0,jump[th]);
				auto probif = 0;
			
				auto al = 0.0;
			
				if(mbp(param_propose,INF_UPDATE) == SUCCESS){
					propose->set_EF();
					propose->set_Pr();
				
					auto probfi = 0;
					if(propose->EF < EFcut) al = exp(propose->Pr - initial->Pr + probfi - probif); 
				}
				
				ntr[th]++;		
				if(ran() < al){
					nac[th]++;
					swap_initial_propose_state();
					if(samp < burnin) jump[th] *= 1.1;
				}
				else{
					if(samp < burnin) jump[th] *= 0.95;
				}
			}
		}	
	}
	
	for(auto th = 0u; th < nparam; th++){
		cout << mpi.core << model.param[th].name << " " << nac[th]/ntr[th] << " " << jump[th] << "\n";
	}
	
	//output.generate_graphs(psamp,opsamp);
}
