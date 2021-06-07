# BEEPmbp

C. M. Pooley† \[1\], I. Hinder \[2\], R. Bailey \[3\], R. Williams\[4\], S. Catterall \[1\],  A. Doeschl-Wilson \[3\] and Glenn Marion \[1\]

\[1\] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK.

\[2\] The University of Manchester, Oxford Rd, Manchester, M13 9PL, UK.

\[3\] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

\[4\] University of Bristol, Queen's Building, University Walk, Clifton BS8 1TR, UK.

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEPmbp (Bayesian Estimation of Epidemic Parameters using Model-Based Proposals) is a general-purpose software tool for simulating and performing inference on epidemiological compartmental models. Inference is the method by which suitable model parameters are chosen from available data. To perform inference BEEPmbp accepts a variety of population-based data (time-series giving the rate of transitions, populations in different compartments and marginal distributions). For example, when analysing COVID-19 disease transmission the following data are used: daily hospitalisations, deaths, populations in hospital as well as overall age distributions for these quantities. 

In BEEPmbp priors on model parameters can be specified from a large range of possibilities. The outputs consist of posterior estimates, a pdf containing numerous plots relating to posterior variation in the model parameters and system state as well as diagnostic information. 

## Epidemiological model features:
*	Specify an arbitrary compartmental epidemiological model.
*	Incorporate spatial and age-structured models.
*	Capture time-variation in reproduction number Rt and external force of infection. 
*	Split population into arbitrary demographic classifications (e.g. age and/or sex). 
*	Incorporate susceptibility variation for different demographic groups.
*	Incorporate area-based covariates to modify the force of infection, either fixed (e.g. population density) or time-varying (e.g. temperature). 
*	Incorporate a user specified age-mixing matrix (along with potential time modification). 
*	Specify a matrix for mixing between different areas (along with potential time modification).
*	Perform prediction, counterfactuals as well as posterior predictive checks.

## Data features:
*	Accepts a variety of different data types (informing transition, populations and marginal).
*	Incorporate splines to relate measured data to system properties.

## Software implementation features:
*	Efficient parallel implementation (written in C++ with MPI).
*	Choose from 7 different inference algorithms.
*	Generates pdf report for easy visualisation of outputs.

## Downloading and running

All information about downloading and running BEEPmbp can be found in the [user manual](BEEPmbp_Manual_v1.0.pdf "User guide").

## Requirements

You need to have MPI installed and the mpicxx compiler to compile and run this code.
