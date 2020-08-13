# BEEPmbp

| Branch        | Test status   | Codacy grade |
| ------------- | ------------- | ------------ |
| master        | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=master)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=master)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=master) |
| dev           | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=dev)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=dev)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=dev) |
| chrispooley   | [![](https://github.com/ScottishCovidResponse/CoronaPMCMC/workflows/CI/badge.svg?branch=chrispooley)](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions?query=workflow%3ACI) |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=chrispooley)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=chrispooley) |

C. M. Pooley† \[1\], I. Hinder \[2\], R. Bailey \[3\], R. Williams\[4\], S. Catterall \[1\],  A. Doeschl-Wilson \[3\] and Glenn Marion \[1\]

\[1\] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK.

\[2\] The University of Manchester, Oxford Rd, Manchester, M13 9PL, UK.

\[3\] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

\[4\] University of Bristol, Queen's Building, University Walk, Clifton BS8 1TR, UK.

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEPmbp (Bayesian Estimation of Epidemic Parameters using Model Based Proposals) is a code for analysing coronavirus using regional level data. This analysis is performed by dividing the area under study (e.g. Scotland or the UK) into small geographical groupings, e.g. at medium super output area (MSOA) level, and modelling the spread of disease. The model captures short range and long range disease transmission by making use of census flow data and previously published age mixing matrices. The data to be analysed is daily hospitalisations and weekly deaths for Covid-19 patients at a health board level along with national demographic data. The time-varying disease transmission rate and infection rate from abroad are estimated, along with the effects of covariates (e.g. age, sex, and population density) on disease progression. 

Parameter inference is performed using a multi-temperature model-based proposal MCMC (MBP-MCMC) approach. This runs MCMC chains at different "temperatures" spanning from the posterior to the prior. This enables the model evidence to be estimated allowing for reliable comparison between different potential models. 

## Downloading

The code should be downloaded from git using the "--recursive" flag to
ensure that submodules are included. For example:
```sh
git clone --recursive https://github.com/ScottishCovidResponse/BEEPmbp.git
```

If you have downloaded without the --recursive flag, you can update the submodules with
```sh
git submodule update --init
```

## Requirements

You need to have MPI installed to compile and run this code.

## Performing analysis

To compile the code, from the repository directory:
```sh
make
```
The code uses mpicxx to compile, so this must be available on your PATH. You can use multiple make processes, e.g. `make -j 4` to speed up compilation.

To run a simulation using demographic data from examples/Data_example
and a fixed set of epidemic parameters:

```sh
./beepmbp inputfile="examples/sim.toml" outputdir="OutputSim"
```

The output appears in OutputSim. If you have Python 3, MatPlotLib and Pandas installed, you can visualise the number of hospitalisations in region 0 with
```sh
python3 -c 'import pandas as pd; import matplotlib.pyplot as plt; plt.plot(pd.read_csv("OutputSim/H.txt",sep="\t",index_col="time")["r0"]); plt.show()'
```

To run inference on simulated data provided in examples/Data_example:

```sh
mpirun -n 2 ./beepmbp inputfile="examples/inf.toml" nchain=2 outputdir="OutputInf"
```

(Note: if you use the same outputdir for inference as simulation,
inference will be performed on the data produced by the simulation,
not the example data.)

You can plot the histogram of samples for f0, which should be
roughly peaked near to the simulated value (0.3). Increase the number
of samples (e.g. nsamp=1000) on the command line to improve this.

```sh
python3 -c 'import pandas as pd; import matplotlib.pyplot as plt; plt.hist(pd.read_csv("OutputInf/trace.txt",sep="\t")["beta0"][10:],30); plt.show()'
```

The input TOML file provides details of simulation or inference and contains all the information BEEPmbp needs to define the compartmental model and provide the filenames for the data. Examples of these files can be found in the "examples" directory, along with a simple test dataset.

## Inputs

Here is a description of the various parameters used in the input TOML files:

### General

**mode** - Defines how the code operates:
		"sim" generates simulated data.
		"inf" performs inference using multi-temperature MBP-MCMC.

**timeformat** - The time format used (either dates or numerical times).

**start** and **end** - The time period for simulation / inference.

**seed** - Sets the random seed when performing inference (this is set to zero by default).

**nchain** - The total number of chains used when performing MCMC (should be a multiple of the number of cores).

**nsamp** - The number of samples used for inference (note, burnin is assumed to be a quarter this value).

**invTmax** - Sets the inverse temperature of the posterior chain (optional, set to 0.25 by default).  

**invTmin** - Sets the inverse temperature of the prior chain (optional, set to 0.0 by default but can be set to invTmax to speed up inference if model evidence not required). 

**outputdir** - Gives the name of the output directory (optional, set to "Output" by default).

### The model

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

**timep** - Used to define discrete time periods (e.g. before and after lockdown).

### The data

**datadir** - The data input directory.

**regions** - Filename for a table giving data regions.

**areas** - Filename for a table giving information about areas (e.g. MSOAs or OAs).

**genQ** - Sets whether spatial and age mixing matrices are used to derive the Q tensors for the model.

**geomix** - Filename for a table informing the matrix of contacts between different areas.

**agemix** - Filename for a matrix giving the contact rates between different age groups.

**genQoutput** - Filenames for the derived Q tensors.

**Q** - Assigns which Q tensor is applied to which compartment and time period (as defined in 'timep').

**threshold** - Used to define the threshold when data is truncated below this value (represented by '*' in the data files).

**transdata** - Transition data. Gives the observed numbers of transitions between different compartments. 

**popdata** - Population data. Gives the number of individuals in specified compartments at specified time points.

**margdata** - Allows for marginalised distributions to be added (e.g. the percentage of overall cases in different age categories).

Note: all the single variable quantities in the TOML file can be overridden using equivalent command line definitions.

For example:
```sh
mpirun -n 1 ./beepmbp inputfile="inf.toml" nsamp=10000 
```

will generate 10000 samples (irrespective of the definition given in "inf.toml").
	
## Outputs

### Simulation

This creates the specified 'transdata', 'popdata' and/or 'margdata' files in the output directory.

### Inference

The output directory contains posterior information (with means and 90% credible intervals) for:

1.  Plots for the transitions corresponding to the 'transdata', 'popdata' and/or 'margdata' files.
2.  "Posterior_R0.txt" gives posterior plots for the time variation in the basic reproduction number R0.
3.  "Posterior_Rmap.txt'" gives the basic reproduction number R0 as a function of time for different areas.
4.  "Posterior_phi.txt" gives posterior plots for time variation in the external force of infection (in units of infections per 100000 individuals).
5.  "Posterior_parameter.txt" gives information about parameters.
6.  "Posterior_distributions.txt" gives posterior distributions for parameters (using a binning procedure).
7.  "trace.txt" gives trace plots for different models.
8.  "traceLi.txt" gives trace plots for the likelihoods on different chains.
9.  "MCMCdiagnostic.txt" gives diagnostic information on the MCMC algorithm.
10. "MCMCdiagnostic_timings.txt", provides information about CPU times for different parts of the algorithm.

Diagnostic checks: two types of checks can be performed to ensure that the results obtained are reliable:

1.  Estimates for the effective sample size in "Posterior_parameter.txt". These should exceed 200 for all parameters if the number of samples is sufficiently large. If this is not the case it indicates that MCMC should be run with more samples (see the 'nsamp' option in the input TOML file).

2.  Results from different runs can be combined to ensure that they all converge on the same posterior distribution (if the likelihood exhibits significant multimodality then under some circumstances different runs can converge on different solutions rendering the results questionable). This is achieved by running BEEPmbp in 'combinetrace' mode. For example, if two sets of inference results using different seeds have been placed into directories 'OutputA' and 'OutputB', the following command:
    ```sh
    ./beepmbp mode="combinetrace" dirs="OutputA,OutputB" output="parameter_combined.txt" dist="distribution_combined.txt"
    ```
    generates a file combining the two sets of samples along with Gelman–Rubin convergence diagnostic results that test for convergence across runs. Optionally parameter distributions can also be generated by setting the 'dist' property. 

## Development

-   [Code documentation](https://projectdata.scrc.uk/beepmbp/branches/dev/doxygen/index.html) -- not yet automatically updated

    -   [Call tree](https://projectdata.scrc.uk/beepmbp/branches/dev/doxygen/main_8cc.html#a3c04138a5bfe5d72780bb7e82a18e627)

-   Work on feature branches and then create a PR for merge into dev branch

-   Continuous integration is implemented using Github Actions,
  controlled by a [workflow](.github/workflows/ci.yml) file. Whenever
  a branch is pushed to Github, the workflow is run on that
  branch. The workflow compiles and runs the code in both simulation
  and inference mode. The results of the CI run should be emailed to
  the committer.  They are also visible on the
  [Actions](https://github.com/ScottishCovidResponse/CoronaPMCMC/actions)
  page.

-   Regression tests can be run with

```sh
make test
```

This will report whether the code gives the same results as when the reference data in `tests/*/refdata` was
  committed. The output of the tests will be stored in `regression_test_results`.  If the tests fail because you
  have made a change which *should* change the results, run

```sh
make test-update
```
	
    which will store the new results from `regression_test_results` into `tests/*/refdata'.  If you ran the tests
  with an uncommitted version of the code, this will fail; you need to commit the code changes and rerun "make test"
  before storing the results. This ensures that the reference data corresponds to a committed version of the code.
  You can then commit the new results with
  
```sh
git add tests/*/refdata
git commit -m "Regenerate test reference data"
```
	
    The tests pass on every Linux platform we have run on, but runs on macOS give sufficiently different results that they fail. See
  [#628](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/628).

-   Code tests (unit tests, and other tests which might not be considered unit tests) can be run with

```sh
make codetest
```

    There are currently only example tests implemented (1% coverage as of 30-Jul-2020). Coverage is reported in the build logs accessible from the [GitHub Actions Page (dev branch)](https://github.com/ScottishCovidResponse/BEEPmbp/actions?query=branch%3Adev).

-   Running "make" stores intermediate build objects in the "build"
  directory, but the executable (for historical and convenience
  reasons) is written to the current directory.

-   Clean the build with `make clean`
