# BEEPmbp

| Branch        | Test status   | Codacy grade |
| ------------- | ------------- | ------------ |
| master        | [![](https://github.com/ScottishCovidResponse/BEEPmbp/workflows/CI/badge.svg?branch=master)](https://github.com/ScottishCovidResponse/BEEPmbp/actions?query=workflow%3ACI)      |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=master)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=master)           |
| dev           | [![](https://github.com/ScottishCovidResponse/BEEPmbp/workflows/CI/badge.svg?branch=dev)](https://github.com/ScottishCovidResponse/BEEPmbp/actions?query=workflow%3ACI)         |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=dev)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=dev)                 |
| chrispooley   | [![](https://github.com/ScottishCovidResponse/BEEPmbp/workflows/CI/badge.svg?branch=chrispooley)](https://github.com/ScottishCovidResponse/BEEPmbp/actions?query=workflow%3ACI) |[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6b91cb37e62409ab926da36727e6f61?branch=chrispooley)](https://www.codacy.com/gh/ScottishCovidResponse/BEEPmbp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ScottishCovidResponse/BEEPmbp&amp;utm_campaign=Badge_Grade?branch=chrispooley) |

C. M. Pooley† \[1\], I. Hinder \[2\], R. Bailey \[3\], R. Williams\[4\], S. Catterall \[1\],  A. Doeschl-Wilson \[3\] and Glenn Marion \[1\]

\[1\] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK.

\[2\] The University of Manchester, Oxford Rd, Manchester, M13 9PL, UK.

\[3\] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

\[4\] University of Bristol, Queen's Building, University Walk, Clifton BS8 1TR, UK.

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BEEPmbp (Bayesian Estimation of Epidemic Parameters using Model-Based Proposals) is a general-purpose software tool for simulating and performing inference on compartmental models. Inference is the method by which suitable model parameters are chosen from available data. To perform inference BEEPmbp accepts a variety of population-based data (time-series giving the rate of transitions, populations in different compartments and marginal distributions). For example, when analysing COVID-19 disease transmission the following data are used: daily hospitalisations, deaths, populations in hospital as well as overall age distributions for these quantities. 

In BEEPmbp priors on model parameters can be specified from a large range of possibilities. The outputs consist of posterior estimates, a pdf containing numerous plots relating to posterior variation in the model parameters and system state as well as diagnostic information. 

## BEEPmbp features:
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

All information about downloading and running BEEPmbp can be found in the [user manual](BEEPmbp_Manual_v1.0.pdf "User guide")

## Requirements

You need to have MPI installed and the mpicxx compiler to compile and run this code.


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
  [Actions](https://github.com/ScottishCovidResponse/BEEPmbp/actions)
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
    which will store the new results from `regression_test_results` into `tests/*/refdata`.  If you ran the tests
  with an uncommitted version of the code, this will fail; you need to commit the code changes and rerun `make test`
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

-   Running `make` stores intermediate build objects in the `build`
  directory, but the executable (for historical and convenience
  reasons) is written to the current directory.

-   Clean the build with `make clean`
