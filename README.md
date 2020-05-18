
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

To simulate from the model use:

./analysis 1024 0
 
The first number gives the number of groups into which the houses are divided (should be a power of 4) and the second number changes the random seed.

Simulation generates the file "Weekly case data.txt".

Inference can then be performed on this data by using:

./analysis 1024 1000 0

Here the first number gives the number of groups into which the houses are divided (should be a power of 4), the second number gives the number of PMCMC samples and the third number changes the random seed.

Model parameter samples are generated in the file "trace.txt".
