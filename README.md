
# CoronaPMCMC

C. M. Pooley† [1] and Glenn Marion [1]

[1] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK 

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

CoronaPMCMC is a code for analysing coronavirus regional level data. This divides the area under study (e.g. Scotland or the UK) into a large number of geographical groupings of houses (e.g. ~65k groups). The model captures short range and long range disease transmission (as well as household effects and the inclusion of schools / workplaces / care homes, which  will be added soon). The data to be analysed with be weekly case numbers at a Health board level along with national mortality data.

Parameter inference is performed using particle MCMC (PMCMC). This relies on being able to generate fast simulations from the model (of the order of seconds). For this purpose we have developed a fast Gillespie algorithm which relies on a multiscale approach. A parallelised version will be created soon to increase speed further.

## Performing analysis

To compile the code use:
make


Simulation:  
       
mpirun -n 1 ./run mode=sim model=irish simtype=smallsim seed=0 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

MBP Inference (expected to be best under most circumstances):   
 
mpirun -n 1 ./run mode=mbp model=irish simtype=smallsim nchain=1 nsamp=100 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt

PMCMC Inference (an alternative):  

(TODO: update this)

mpirun -n 20 ./run mode=pmcmc model=irish area=1024 npart=20 nsamp=1000 period=16 transdata=I,H,reg,cases.txt transdata=H,D,all,deaths.txt housedata=house.txt outputdir=Output

The flag -n 1 sets the number of cores (set to 1 for simulation or more for inference)

INPUTS:

(TODO: update this)

Here is a description of the various inputs:

**mode** - Defines how the code operates:
		"sim" generates simulated data.
		"pmcmc" performs inference using particle MCMC.
		"mbp" performs inference using multi-temperature MBP MCMC.

**model** - This defines the compartmental model being used:
		"irish" defines a SEAPIRHD model for asymtopmatic / presymptomatic / symptomatic individuals.
	
**area** - Determines	the number of areas into which the houses are divided (should be a power of 4)

**npart** - The total number of particles used when performing PMCMC (should be a multiple of the number of cores).

**nchain** - The total number of chains used when performing multi-temperature MBP MCMC (should be a multiple of the number of cores).

**nsamp** - The number of samples used for inference (note, burnin is assumed to be a quater this value).

**seed** - Sets the random seed when performing inference (this is set to zero by default)

**outputdir** - Gives the name of the output directory (optional).

**transdata** - Transition data. Gives the "from" then "to" compartments, the type ("reg" means regional data and "all" means global) and then the file name. More than one set of transition data can be added for an analysis.

**housedata** - House data. Gives the file name for a file giving the positions of houses.

**period** - The period of time over which simulation/inference is performed (e.g. measured in weeks).

**simtype** - Determines the system on which simulation is performed (simulation mode only):
		"smallsim" a small system with 1024 houses, 10000 individuals and a 2x2 grid of data regions (used for testing)
		"scotsim" a Scotland-like system with 1.5 million houses, 5.5 million individuals and a 4x4 grid of data regions (used for testing)
		"uksim" a UK-like system with 20 million houses, 68 million individuals and a 10x10 grid of data regions (used for testing)
	
	
OUTPUTS:

(TODO: update this)

Simulation - This creates the specified 'transdata' and 'housedata' files along with output directory containing:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives time variation in R0.
3) "parameter.txt", which gives the parameter values used in the simulation.

Inference - The output directory contains postior information (with means and 90% credible intervals) for:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives posterior plotstime variation in R0.
3) "parameter.txt", which gives information about parameters.
4) "trace.txt", which gives trace plots for different models.
5) "traceLi.txt", which gives trace plots for the likelihoods on different chains (MBPs only).
6) "MCMCdiagnostic.txt", which gives diagnostic information on the MCMC algorthm.
