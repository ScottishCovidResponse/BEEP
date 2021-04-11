# BEEPmbp

C. M. Pooley† \[1\], I. Hinder \[2\], R. Bailey \[3\], R. Williams\[4\], S. Catterall \[1\],  A. Doeschl-Wilson \[3\] and Glenn Marion \[1\]

\[1\] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK.

\[2\] The University of Manchester, Oxford Rd, Manchester, M13 9PL, UK.

\[3\] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

\[4\] University of Bristol, Queen's Building, University Walk, Clifton BS8 1TR, UK.

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEPmbp (Bayesian Estimation of Epidemic Parameters using Model-Based Proposals) is a general-purpose software tool for simulating and performing inference on compartmental models. Inference is the method by which suitable model parameters are chosen from available data. To perform inference BEEPmbp accepts a variety of population-based data (time-series giving the rate of transitions, populations in different compartments and marginal distributions). For example, when analysing COVID-19 disease transmission the following data are used: daily hospitalisations, deaths, populations in hospital as well as overall age distributions for these quantities. 

In BEEPmbp priors on model parameters can be specified from a large range of possibilities. The outputs consist of posterior estimates, a pdf containing numerous plots relating to posterior variation in the model parameters and system state as well as diagnostic information. 

## Features:
*	Specify an arbitrary compartmental model.
*	Incorporate spatial and age-structured models.
*	Time-varying splines to capture variation in reproduction number R(t) and the external force of infection. 
*	Split population into arbitrary demographic classifications (age, sex). 
*	Incorporate susceptibility variation for different demographic groups.
*	Add area-based covariates to modify the force of infection (e.g. population density). 
*	Incorporate a user specified age-mixing matrix (along with potential time modification). 
*	Specify matrix for mixing between different areas (along with potential time modification).
*	Perform posterior predictive checks as well as analyse counterfactuals.
*	Choose from 6 different inference algorithms.
*	Efficient parallel implementation.

## Downloading and running

All information about downloading and running BEEPmbp can be found in the [user manual](BEEPmbp_Manual_v1.0.pdf "User guide").

## Requirements

You need to have MPI installed and the mpicxx compiler to compile and run this code.
