
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

Simulate using:   mpirun -n 1 ./analysis 65536 1

-n 1 sets the number of cores (set to 1 for simulation)

The first number gives the number of areas into which the houses are divided (should be a power of 4)

The second number gives the random seed

Simulation generates the files

- houses_small.txt: Randomly generated data about houses; eventually this will come from real data
- cases_small.txt: Simulated case data
- events_small.txt: Simulated events

Inference can then be performed on this data by using:

mpirun -n 20 ./analysis 65536 500 1

-n 20 gives the number cores

The first number gives the number of areas into which the houses are divided (should be a power of 4)

The second number gives the number of MCMC iterations

The third number gives the number of particles per core

Output appears in the Output directory:

- R0_<N>.txt: R0
- params_<N>.txt: Estimated parameters with credible intervals etc
- region_<A>_<B>_data_<N>.txt: Information about what's happening in each region
- trace_<N>.txt: Model parameter samples
