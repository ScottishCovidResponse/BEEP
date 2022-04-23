/// Implements an importance sampler 

#include <iostream>
#include <fstream>
#include <algorithm> 
#include <sstream> 

using namespace std;

#include "importance.hh"

/// Initilaises the PAIS class
IMPORTANCE::IMPORTANCE(const Details &details, const Data &data, const Model &model, Inputs &inputs, const Output &output, const ObservationModel &obsmodel, Mpi &mpi) : state(details,data,model,obsmodel), details(details), data(data), model(model), output(output), obsmodel(obsmodel),mpi(mpi)
{	
	inputs.find_nsample(Ntot,1);
	if(Ntot != UNSET && Ntot%mpi.ncore != 0)  emsgroot("'nsample' must be a multiple of the number of cores");
	N = Ntot/mpi.ncore;
	
	inputs.find_invT(invT);
	
	obs_slice = obsmodel.generate_obs_slice();
}

void IMPORTANCE::run()
{
	//state.integrated_matrix_prod_check(); return; 

	auto param = model.sample_from_prior();  
	
	vector < vector <double> > av;
	vector < vector < vector <double> > > M;
	
	auto N = 3u;
	av.resize(details.ndivision); M.resize(details.ndivision);
	for(auto sett = 0u; sett < details.ndivision; sett++){
		av[sett].resize(N); M[sett].resize(N);
		for(auto i = 0u; i < N; i++){
			av[sett][i] = 0;
			M[sett][i].resize(N);
			for(auto ii = 0u; ii < N; ii++){
				M[sett][i][ii] = 0;
			}
		}
	}
	
	auto kmax = 1000u;   
	for(auto k = 0u; k < kmax; k++){
		cout << k <<  "l" << "k\n";
		state.simulate(param);
		
		for(auto sett = 0u; sett < details.ndivision; sett++){
			for(auto i = 0u; i < N; i++){
				av[sett][i] += state.pop[sett][0][i][0]; 
				for(auto ii = 0u; ii < N; ii++){
					M[sett][i][ii] += state.pop[sett][0][i][0]*state.pop[sett][0][ii][0];
				}
			}
		}
	}
	
	for(auto sett = 0u; sett < details.ndivision; sett++){
		for(auto i = 0u; i < N; i++){
			av[sett][i] /= kmax;
			for(auto ii = 0u; ii < N; ii++){
				M[sett][i][ii] /= kmax;
			}
		}
	}
	
	
	auto filesim = details.output_directory+"/sim.txt";
	ofstream opsim(filesim);

	for(auto sett = 0u; sett < details.ndivision; sett++){
		opsim << sett << " " << av[sett][0] << " " << av[sett][1] << " " << av[sett][2] << " "
			<< M[sett][0][0]-av[sett][0]*av[sett][0] << " " 
			<< M[sett][1][1]-av[sett][1]*av[sett][1] << " " 
			<< M[sett][0][1]-av[sett][0]*av[sett][1] << " hello\n";
	}
	
	state.initialise_approx();

//param[3] = 2.4;  

	auto L = state.likelihood_approx(param,obs_slice,invT,DOUBLE);  
		cout << L << "HEREDONE\n";
	//	return;
		
	auto path = details.output_directory+"/Approx";
	output.ensure_directory(path);
	
	for(auto &gr : data.graph){
		const auto &dt = data.datatable[gr.datatable];
		if(dt.optype == DATA){		 
			auto file = path+"/"+gr.file.substr(0,gr.file.length()-4)+".txt";
			ofstream op(file);
			for(auto sett = 0u; sett < details.ndivision; sett++){
				auto var = 0.0;
				auto sum = 0.0;
				for(auto co : dt.complist){
					for(auto c : gr.area){
						for(auto dp : gr.dp_sel){
							sum += state.pop[sett][c][co][dp];
							auto cc = state.comp_approx_ref[c][co][dp];
							//cout << cc << " " << state.pop_covar.size() << " " << state.pop_covar[sett].size() << " cc\n";
							
							var = state.pop_covar[sett][cc][cc];
						}
					}
				}
				op << sett/details.division_per_time << " " << sum << " " << sum - sqrt(var) << " " << sum + sqrt(var) <<  endl;
			}
			
			auto file_data = path+"/"+gr.file.substr(0,gr.file.length()-4)+"_data.txt";
			ofstream opdata(file_data);	
			for(const auto &gp : gr.point){
				opdata << 0.5*(gp.xi+gp.xf) << " " << data.obs[gp.obs].value << endl;
			}
			
//				op << sett << " " << state.pop[sett][0][0][0] << " " << state.pop[sett][0][1][0] << " " << state.pop[sett][0][2][0] << " "
	//				 << state.pop_covar[sett][0][0] << " "	<< state.pop_covar[sett][1][1] << " " << state.pop_covar[sett][0][1] << " \n";
	//		}
	
			cout << gr.name << " " <<  gr.datatable <<  "name\n";	
		}
	}
	
	auto file = details.output_directory+"/compare.txt";
	ofstream op(file);

	for(auto sett = 0u; sett < details.ndivision; sett++){
		op << sett << " " << state.pop[sett][0][0][0] << " " << state.pop[sett][0][1][0] << " " << state.pop[sett][0][2][0] << " "
		<< state.pop_covar[sett][0][0] << " "	<< state.pop_covar[sett][1][1] << " " << state.pop_covar[sett][0][1] << " \n";
	}
	
	state.future_obs_approx(obs_slice,invT,DOUBLE);
	
	for(auto &gr : data.graph){
		const auto &dt = data.datatable[gr.datatable];
		if(dt.optype == DATA){		 
			auto file = path+"/"+gr.file.substr(0,gr.file.length()-4)+"_future.txt";
			ofstream op(file);
			for(auto sett = 0u; sett < details.ndivision; sett++){
				auto var = 0.0;
				auto sum = 0.0;
				for(auto co : dt.complist){
					for(auto c : gr.area){
						for(auto dp : gr.dp_sel){
							sum += state.pop[sett][c][co][dp];
							auto cc = state.comp_approx_ref[c][co][dp];
							
							var = state.obs_covar[sett][cc][cc];
						}
					}
				}
				auto min = sum - sqrt(var); if(min < 0) min = 0;
				auto max =  sum + sqrt(var); if(max > 2*sum) max = 2*sum;
				op << sett/details.division_per_time << " " << sum << " " << min << " " << max <<  endl;
			}
		}
	}
	
	//auto L = state.likelihood_approx(param,obs_slice);  
		
	//state.generate_steer(param);  
	//state.simulate(param);  
}
