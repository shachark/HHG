# HHG
Heller-Heller-Gorfine Tests of Independence and Equality of Distributions.

This is a mix of R and C++ code relating to the HHG project at the Heller lab, Department of Statistics and Operations Research, School of Mathematics, Tel Aviv University. This includes the development repository for the HHG package true to the point where I handed it off to another maintainer (June, 2014). For the most current version, see CRAN.

This README describes the structure of the project and explains how interested users can set it up on their respective machines.

### Structure

The source code for the development version of the HHG package has the standard R package internal directory structure, see [Writing R Extentions](http://cran.r-project.org/doc/manuals/R-exts.html).

R code for the main functions of the package is under R/HHG.R. 
Extra functions such as simulated data examples, and existing tests (that we don't consider our own contribution) are in R/HHG_extras.R
Some shared utilities are implemented in R/HHG_utils.R
There is some experimental conditional independence code in R/CI.R

Some (most) of the computational heavy-lifting is implemented in C++. It interfaces with the package's R code via the .Call interface (see Writing R Extensions). At some point it might make sense to use the more modern Rcpp infrastructure which could make it easier to follow, but the current implementation is pretty efficient and flexible. We currently have one entry function (HHG_R_C) that handles all the different tests. This is implemented in Package-dev/src/HHG.h and Package-dev/src/HHG.cpp. These files also define some shared data structures and utility functions used throughout the C++ implementation.

When this function is called, it figures out the type of test requested and parses (and preprocesses) the R structures holding the data, allocates an output structure, and then calls SequentialTest to do further processing and populate the output.

SequentialTest does as it name suggests, but also a degenerate sequential test can be a regular permutation test (without Wald test early stopping) or even merely computing a test statistic (with zero permutations). If requested, it handles the Wald sequential logic of early stopping, and knows to distribute the actual computation of resampled test statistics across multiple cores on a single PC (using the common C library pthreads). Each resampling (optional) and computation of a statistic is actually performed by a class called StatsComputer.

StatsComputer should have been a hierarchy of interfaces and classes, but for historical reasons it is one big class with some switches and function pointers (which is basically an old-style C flavor of polymorphism). Redesigning this according to proper OO design principles would be great, but is a lot of work. In any case, the constructor here first allocates temporary buffers needed for the computation, then assigns the function to take care of resampling (e.g., permutation) and test statistic computation.

The tests and their implementations are classified as UV (univariate) and MV (multivariate), and as wither GOF (1-sample), TS (2-sample), KS (k-sample, for k>=2), or IND (independence) which appear in the names of constants and functions. The extended multivariate tests also have a second "test" (the univariate score) defined for them, and it is just the code for computing a UV test statistic.

### Additional files - inventory

This is an attempt to classify and briefly describe the R files under Extras/
 
Software testing: 
- hhg_package_regression_test this is an R script and an RData file with a regression test (the term in software development, not the statistical modeling problem "regression"). The script has three simple modes, and each mode can operate either on the development directory or on the release directory (which is not part of the repository at the time of writing). The three modes are: (1) sanity check (basically, try a simple example of each mode of operation, see if it doesn't crash), (2) generate (creates the data file with results to use as a reference for comparing later), (3) regress (test results using current code against the last saved reference). NOTE: it might need some reviving now, as I haven't run it for a while. It is important to run it prior to CRAN releases! When new features are added to the package they should also be reflected in this file. Still, obviously, this script is only for very basic "dead or alive" testing. More capable testing can be done by rerunning simulations and comparing to results in the various manuscripts.

Simulation:
- ADP_and_MI.r simulation of mutual information estimation using the ADP univariate independence statistic.
- compute-critical-values-gof this goes over null tables generated by generate-nullsim-xdp-gof and generate-nullsim-competition-gof and computes 0.05 etc. quantiles
- compute-critical-values-ks similar to the gof file, but for univariate k-sample tests
- compute-critical-values again similar but for the univariate test of independence
- conditional-sim-finish part of the conditional testing code that as of time of writing is not yet mature
- conditional-sim same
- explore-gof-asymptotic-null-distribution
- generate-nullsims generate null tables for the univariate independence XDP tests
- generate-nullsims-competition generate null tables for some existing univariate independence tests
- generate-nullsims-competition-gof similar but for the GOF tests
- generate-nullsims-competition-ks similar but for the k-sample tests
- generate-nullsims-gof generate null tables for the univariate GOF XDP tests
- generate-nullsims-k-sample generate null tables for the univariate k-sample XDP tests
- hhg-ks-extended-sim new multivariate k-sample test power simulation
- xdp-gof-power-sim\* power simulation for the XDP GOF tests
- xdp_asymptotic_nullsim
- xdp-ks-power-sim power simulation for the XDP k-sample tests
- xdp-power-sim\* power simulation for the XDP independence tests

Real data analysis:
- analyze-aCGH (somewhat boring gene expression data that we analyzed)
- analyze-hapmap-k-sample some k-sample tests run on similar populations in the HapMap genoytping data (the problems were too easy for our tests to be needed)
- hhg-ks-extended-examples new multivariate analysis examples from standard datasets (none were of particular interest so far)
- rosetta\* application of the XDP independence test to gene expression data

Misc.
- competing-k-sample-tests various games with existing k-sample tests
- critical-values\_\*.RData these are data files each containing an R list with critical values for a specific sample size(s) and nominal level, for the various problems of interest. Field names should tell you which value is for which test, but this can also be seen in generate-nullsim-\* and the relevant power simulation code.
- MINE is actually the implementation of the MIC method (Reshef et al.)
- count-rectangles
- ddp-package-example
- fixup_xdp_nullsims
- debug-extended-hhg
- debug-power
- measure-running-time
- try-dyn-slicing
- two_sample_examples
- udf-paper-stuff
- xdp_asymptotic_nullsim_finish
- xdp-k-sample-debug
- hhg-udf-\* additional power simulation for the XDP tests of independence

### To install diectly from GitHub:
```
install.packages('devtools')
devtools::install_github('shachark/HHG')
```
