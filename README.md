
# CoronaBICI

C. M. Pooley† [1] and Glenn Marion [1]

[1] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK 

[2] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

CoronaPMCMC is a code for analysing coronavirus regional level data. This divides the area under study (e.g. Scotland or the UK) into a large number of geographical grouping of houses (e.g. ~65k groups). The model captures short range and long range disease transmission (as well as household effects and the inclusion of schools / workplaces / carehomes, which  will be added soon). The data to be analysed with be weekly case number at Healthboard level along with national mortality data.

Parameter inference is performed using particle MCMC (PMCMC). This relies on being able to simulate from the model in a fast way. For this purpose we have developed a fast Gillespie algoritm which relies on a multiscale approach.

## Performing analysis

To complie the code use

g++ analysis.cc -o analysis -O3

To simulate from the model use:

 ./analysis 1024 0
 
The first number gives the number of groups into which the houses are divided (should be a power of 4) and the second number changes the random seed.

This generates the file "Weekly case data.txt".

Inference can then be performed on this data by using:

./analysis 1024 1000 0

Here the first number gives the number of groups into which the houses are divided (should be a power of 4), the second number give the number of PMCMC samples and the third number changes the random seed.

Output samples are generated in the file "trace.txt"
