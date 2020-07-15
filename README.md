
# BEEPmbp

| Branch        | Test status   |
| ------------- | ------------- |
| master        | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=master)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |
| dev           | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=dev)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |
| chrispooley   | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=chrispooley)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |

C. M. Pooley† [1] and Glenn Marion [1]

[1] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK 

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEPmbp (Bayesian Estimation of Epidemic Parameters using Model Based Proposals) is a code for analysing coronavirus using regional level data. This analysis is performed by dividing the area under study (e.g. Scotland or the UK) into small geographical groupings, e.g. at medium super output area (MSOA) level or output area (OA) level, and modelling the spread of disease. The model captures short range and long range disease transmission by making use of census flow data and previously published age mixing matrices. The data to be analysed is weekly case numbers at a healthboard level along with national mortality data. The time-varying disease transmission rate and infection rate from abroad are estimated, along with the effects of covariates (e.g. age, sex, and population density) on disease susceptibility. 

Parameter inference is performed using a multi-temperature model-based proposal MCMC (MBP-MCMC). This runs MCMC chains at different "temperatures" spanning from the posterior to the prior. This enables the model evidence to be estimated allowing for reliable comparison between different potential models. 

## Performing analysis

Compilation: make

Simulation:  ./beepmbp inputfile="examples/sim.toml" 

Inference:   mpirun -n 2 ./beepmbp inputfile="examples/inf.toml" nchain=2
(nchain and -n must be the same; they are the number of chains run, and hence the number of processes used)

The input TOML file provides details of simulation or inference and contains all the information BEEPmbp needs to define the compartmental model and provide the filenames for the data. Examples of these files can be found in the "examples" directory, along with an simple test dataset.
 
# INPUTS:

Here is a description of the various commands used in the TOML files:

DETAILS

**mode** - Defines how the code operates:
		"sim" generates simulated data.
		"inf" performs inference using multi-temperature MBP-MCMC.

**period** - The time period of simulation / inference (in weeks).

**seed** - Sets the random seed when performing inference (this is set to zero by default)

**nchain** - The total number of chains used when performing MCMC (should be a multiple of the number of cores).

**nsamp** - The number of samples used for inference (note, burnin is assumed to be a quarter this value).

**outputdir** - Gives the name of the output directory (optional).

THE MODEL

Note, for examples of how these commands are used, see 'inf.toml'.

**comps** - Defines the compartments in the model.

**trans** - Defines transitions between compartments.

**params** - Defines parameters in the model (simulation only).

**priors** - Defines priors in the model (inference only).

**indmax** - The maximum number of infected individuals (placed as a prior).

**betaspline** - Defines a linear spline used to capture time variation in transmission rate beta.

**phispline** - Defines the linear spline used to represent external force of infection phi.

**ages** - The age groups used in the analysis.

**democats** - Used to define other demographic categories.

THE DATA 

**datadir** - The data directory.

**regions** - Filename for a table giving data regions.

**areas** - Filename for a table giving information about areas (e.g. MSOAs or OAs).

**geocont** - Filename for a table informing the matrix of contacts between different areas.
 
**agecont** - Filename for a matrix giving the contact rates between different age groups.

**transdata** - Transition data. Gives the observed numbers of transitions between different compartments. More than one set of transition data can be used in an analysis.

Note, all the single variable quantities in the TOML file can be overridden using equivalent comand line definitions.

For example: mpirun -n 1 ./beepmbp inputfile="inf.toml" nsamp=10000 

generate 10000 samples (irrespective of the definition given in "inf.toml").
	
# OUTPUTS:

Simulation - This creates the specified 'transdata' files along with an output directory containing:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives time variation in R0.
3) "parameter.txt", which gives the parameter values used in the simulation.

Inference - The output directory contains posterior information (with means and 90% credible intervals) for:
1) Plots for the transitions corresponding to the 'transdata' files.
2) "R0.txt", which gives posterior plots time variation in R0.
3) "parameter.txt", which gives information about parameters.
4) "trace.txt", which gives trace plots for different models.
5) "traceLi.txt", which gives trace plots for the likelihoods on different chains.
6) "MCMCdiagnostic.txt", which gives diagnostic information on the MCMC algorithm.

# Development

- [Code documentation](https://projectdata.scrc.uk/coronapmcmc/branches/master/doxygen/html/index.html) -- not yet automatically updated
  - [Call tree](https://projectdata.scrc.uk/coronapmcmc/branches/master/doxygen/html/analysis_8cc.html#a3c04138a5bfe5d72780bb7e82a18e627)
- Currently work on feature branches and then merge into master
- Continuous integration is implemented using Github Actions,
  controlled by a [workflow](.github/workflows/ci.yml) file. Whenever
  a branch is pushed to Github, the workflow is run on that
  branch. The workflow compiles and runs the code in both simulation
  and inference mode. The results of the CI run should be emailed to
  the committer.  They are also visible on the
  [Actions](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions)
  page.
