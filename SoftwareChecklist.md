# SCRC Software checklist

This checklist is part of ongoing work on a model scoresheet for SCRC models. It relates to software implementation, and assumes that other documents cover questions about model validation, data provenance and quality, quality of science, and policy readiness.

## Software Details

### Model / software name

> BEEPmbp

### Date

> 31-Jul-2020

### Version identifier

> ebf714f

## Overall statement

Do we have sufficient confidence in the correctness of the software to trust the results?

This is your overall judgement on the level of confidence based on all the aspects of the checklist. There is no formulaic way to arrive at this overall assessment based on the individual checklist answers but please explain how the measures in place combine to reach this level of confidence and make clear any caveats (eg applies for certain ways of using the software and not others).

> -   [ ] Yes
> -   [ ] Yes, with caveats
> -   [X] No
>
> The model has been undergoing heavy feature development over the past few weeks, and is now starting to stabilise, with improvements to the inference algorithm. Due to the need for urgent and wide-ranging development by the model owner, it has not been possible to attempt large-scale code cleanup, refactoring, and implementation of units tests due to the problem of conflicts, and this will be necessary to gain confidence in the quality of the software. The software is not yet integrated in the data pipeline.

## Checklist

Please use a statement from this list: "Sufficiently addressed", "Some work remaining or caveats", or "Needs to be addressed" to begin each response.

Additionally, for each question please explain the situation and include any relevant links (eg tool dashboards, documentation). The sub bullet points are to make the scope of the question clear and should be covered if relevant but do not have to be answered individually.

### Can a run be repeated and reproduce exactly the same results?

-   How is stochasticity handled?
-   Is sufficient meta-data logged to enable a run to be reproduced: Is the exact code version recorded (and whether the repository was "clean"), including versions of dependent libraries (e.g. an environment.yml file or similar) along with all command line arguments and the content of any configuration files? 
-   Is there up-to-date documentation which explains precisely how to run the code to reproduce existing results? 

> -   [ ] Sufficiently addressed
> -   [X] Some work remaining or caveats
> -   [ ] Needs to be addressed
> 
> Version information is recorded in the log file, and will be recorded in the data pipeline. Input parameters are recorded in the log file.  Dependencies are handled as submodules, so are included in the git hash of the main repository.  Random seeds are specified in the parameter file. However, the decision of how many proposals to run on each process is made based on timings, and hence runs cannot be exactly reproduced. There is an alternative reproducible method implemented (propsmethod=fixed), but it is expected to be very inefficient. We [plan](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/602) to store the number of proposals run at each step so that it can be replayed for exact reproduction of an efficient run, but this is not yet implemented.

### Are there appropriate tests?  (And are they automated?)

-   Are there unit tests? What is covered?
-   System and integration tests?  Automated model validation tests?
-   Regression tests? (Which show whether changes to the code lead to changes in the output. Changes to the model will be expected to change the output, but many other changes, such as refactoring and adding new features, should not. Having these tests gives confidence that the code hasn't developed bugs due to unintentional changes.)
-   Is there CI?
-   Is everything you need to run the tests (including documentation) in the repository (or the data pipeline where appropriate)?

> -   [ ] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [X] Needs to be addressed
> 
> There is a regression test framework which catches inadvertent changes to the results. Coverage of the regression test framework is ~80%. We have a unit test framework. Both types of tests are run [regularly in CI](https://github.com/ScottishCovidResponse/BEEPmbp/actions?query=branch%3Adev).  The [README](README.md) describes how to run the tests.  However, there are minimal [unit tests](codetests), as the important parts of the code do not currently lend themselves to unit testing. This will require some refactoring to fix, and unit tests to be [added](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/556).

### Are the scientific results of runs robust to different ways of running the code?

-   Running on a different machine?
-   With different number of processes?
-   With different compilers and optimisation levels?
-   Running in debug mode?

(We don't require bitwise identical results here, but the broad conclusions after looking at the results of the test case should be the same.) 

> -   [ ] Sufficiently addressed
> -   [X] Some work remaining or caveats
> -   [ ] Needs to be addressed
> 
> We have tested on several compilers in macOS and Linux, and find broadly the same results. Regression tests are passed (bitwise identical results of ASCII output files, with default (6 digit) floating point precision output) on both Clang and GCC on Linux, but [not on macOS](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/628) in either. We have not yet investigated the root cause of the differences on macOS. Further testing of optimisation flags needs to be performed.

### Has any sort of automated code checking been applied?

-   For C++, this might just be the compiler output when run with "all warnings". It could also be more extensive static analysis. For other languages, it could be e.g. pylint, StaticLint.jl, etc.
-   If there are possible issues reported by such a tool, have they all been either fixed or understood to not be important?

> -   [ ] Sufficiently addressed
> -   [X] Some work remaining or caveats
> -   [ ] Needs to be addressed
> 
> We enable all compiler warnings with GCC and use Clang Static Analyzer. There are currently [no warnings](https://github.com/ScottishCovidResponse/BEEPmbp/actions/runs/190447662).  We have integrated CODACY and get a B grade. There are 126 total issues, 107 of which are for code style, 6 are "error prone" code patterns, 13 are for performance.  The "error-prone" warnings need to be [fixed](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/702).

### Is the code clean, generally understandable and readable and written according to good software engineering principles?

-   Is it modular?  Are the internal implementation details of one module hidden from other modules?
-   Commented where necessary?
-   Avoiding red flags such as very long functions, global variables, copy and pasted code, etc.?

> -   [ ] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [X] Needs to be addressed
> 
> Code has been split into classes, but all member variables are public (i.e. poor encapsulation). Functions are generally not too long. There are minimal or no global variables. Addressing the public member variables is a high priority.

### Is there sufficient documentation?

-   Is there a readme?
-   Does the code have user documentation?
-   Does the code have developer documentation?
-   Does the code have algorithm documentation? e.g. something that describes how the model is actually simulated, or inference is performed?
-   Is all the documentation up to date? 

> -   [ ] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [X] Needs to be addressed
> 
> There is a README which explains how to run the code. The example input files contain comments next to each parameter. There is a Doxygen configuration file, but Doxygen is not run automatically in CI. The existing documentation is up-to-date. The code has reasonable Doxygen comments.  There is little other developer documentation beyond describing how to run tests etc.  There is no algorithm documentation. The main outstanding issue is that there is no detailed documentation for how to use the code to do something specific.

### Is there suitable collaboration infrastructure?

-   Is the code in a version-controlled repository?
-   Is there a license?
-   Is an issue tracker used?
-   Are there contribution guidelines?

> -   [X] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [ ] Needs to be addressed
> 
> The code is stored in a GitHub repository, and the issue tracker and pull requests are used as part of development, though not for the main work by the model owner.  There is a licence and contribution guidelines.

### Are software dependencies listed and of appropriate quality?

> -   [X] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [ ] Needs to be addressed
> 
> The only external dependency is MPI, a standard interface, and implementations are usually of high quality. The code does not depend on a specific implementation.

### Is input and output data handled carefully?

-   Does the code use the data pipeline for all inputs and outputs?
-   Is the code appropriately parameterized (i.e. have hard coded parameters been removed)?

> -   [ ] Sufficiently addressed
> -   [ ] Some work remaining or caveats
> -   [X] Needs to be addressed
> 
> The code is not yet using the data pipeline ([#366](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/366) and [#367](https://github.com/ScottishCovidResponse/SCRCIssueTracking/issues/367)). All parameters are specified in the input file, not the source code.
