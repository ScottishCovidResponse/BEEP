# Outputs

This directory contains the outputs from the BEEPmbp analysis.
We give a brief description below of the various files and folders contained within it:

## Parameter_estimates.csv

If inference is performed, this file gives the posterior means and 95% credible intervals for each of the model parameters.

## Graphs.pdf

This visualises outputs from BEEPmbp.

## Graphs_description.txt

Provides information about the source data for the graphs.

## Model_specification.txt

This gives a summary of the model used to perform the analysis.

## Simulation

This folder contains files giving posterior probability distributions for the model parameters.

These are arranged in a number of different ways:

**parameter** - Contains files for the model parameters.

**state** - Contains files which compare the system state with that observed in the actual data files.

**spline** - Contains files giving time variation in splines used within the model.

**susceptibility** - Contains files giving the variation in susceptibility for different demographic classes within the model.

**sample** - Contains files giving raw posterior samples.

**Rmap.csv** - For spatial models this gives the variation in R across different regions.


## Simulated_data

This gives simulated data files corresponding to the specifications provided in the input TOML file.

## Diagnostics

This provides diagnostic information to inform how well the algorithm is performing:

**Generation.csv** - Shows how model parameters move from the prior distribution to an approximation of the posterior distribution as a function of the generation number ('generations' are used within the ABC-MBP and ABC-SMC procedures to filter particles such that they provide a better and better approximation to the posterior).

**MCMC_proposals.txt** - This shows the performance of any MCMC proposals used.
