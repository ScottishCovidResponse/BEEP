# BEEP

C. M. Pooley† \[1\], I. Hinder \[2\], R. Bailey \[3\], R. Williams\[4\], S. Catterall \[1\],  A. Doeschl-Wilson \[3\] and Glenn Marion \[1\]

\[1\] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK.

\[2\] The University of Manchester, Oxford Rd, Manchester, M13 9PL, UK.

\[3\] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

\[4\] University of Bristol, Queen's Building, University Walk, Clifton BS8 1TR, UK.

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEP (Bayesian Estimation of Epidemiological Parameters) is a general-purpose software tool for simulating and performing inference on compartmental models. 

BEEP incorporates three modes of operation:

**Simulation** – Given a set of model parameters, potential system dynamics can be sampled from the model (note, compartmental models are inherently stochastic, so random differences in disease transmission naturally lead to differences in epidemic outcome).

**Inference** – This is the method by which model parameters are estimated from available data (along with associated uncertainties in these estimates). BEEP accepts a variety of different data types: time series measurements giving transition numbers between selected compartments (e.g. daily cases or weekly deaths), populations in different compartments (e.g. hospitalised population measured each week) and marginal data (e.g. distribution for total number of deaths for different age groups). 

**Prediction** – Based on the results of inference, predictions from the model can be made. These can either estimate future behaviour (scenario analysis), or can be used to look at how things would have turned out differently had the model been altered in some specified way (counterfactual analysis).


## Epidemiological model features:
*	Specify arbitrary compartmental epidemiological models.
*	Capture time-variation in reproduction number Rt and external force of infection.
*	Incorporate spatial stratification.
*	Split population into arbitrary demographic classifications (e.g. age and/or sex). 
*	Incorporate susceptibility variation for different demographic groups.
*	Incorporate area-based covariates that modify the force of infection, either fixed effects (e.g. population density) or time-varying effects (e.g. temperature). 
*	Incorporate a user specified age-mixing matrix (along with potential time modification). 
*	Specify a matrix for mixing between different areas (along with potential time modification).
*	Perform predictions, scenario and counterfactual analysis as well as posterior predictive checks.

## Data features:
*	Accepts a variety of different data types (from transitions, populations and marginal distributions).
*	Incorporate splines to relate measured data to system properties (e.g. to account for the fact that only a fraction of true cases are observed).

## Software implementation:
*	Efficient parallel code (written in C++ with MPI).
*	Choose from 8 different inference algorithms.
*	A web browser visualisation tool for viewing results on maps, graphs, histograms and tables.


## Downloading and running

All information about downloading and running BEEP can be found in the [user manual](BEEP_Manual_v1.0.pdf "User guide").

## Requirements

MPI must be installed to allow for parallelisation. Building and compiling the code requires CMake (minimum version 3.13 required) and mpicxx.
