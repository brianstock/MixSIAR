MixSIAR 3.1.9
=========================

Nothing changed, just trying to get new Zenodo DOI (3.1.8 failed)

MixSIAR 3.1.8
=========================

## MINOR IMPROVEMENTS

* updated manual
* new example: alligator diet, model comparison table with LOO/WAIC weights
* new example: mantis shrimp diet, informative priors + combining sources
* new function 'compare_models' to perform model selection based on LOO/WAIC weights, via 'loo' package
* new function 'combine_sources' to aggregate sources a posteriori
* new function 'summary_stat' to print summary of posterior after combining sources
* new function 'plot_intervals' to plot CI of posterior after combining sources

MixSIAR 3.1.7.9001
=========================

## MINOR IMPROVEMENTS

* Removed `attach.jags` call (now directly use `jags.1$BUGSoutput$sims.list`)
* Added epsilon (multiplicative error term) estimates to output (summary stats, posterior plots)

## BUG FIX

* Fixed bug #129 by using `MCMCpack::rdirichlet` instead of `compositions::rDirichlet.rcomp`
* Small changes to pass R CMD check

MixSIAR 3.1.7.9000
=========================

## MINOR IMPROVEMENTS

* Improved continuous effect output plots: 1) added 95% CI shading, 2) plot each level of fixed effect
* Fixed bug with 1FE + 1RE: empty levels Heidelberger diagnostic (#100)
* Fixed bug with 1FE + 1RE: posterior density plots saved over themselves
* Added return value for output_JAGS function: now returns p.both if 1FE + 1RE

## VERY MINOR IMPROVEMENTS

* Added check/error message: if sources by factor, sources must match for all factor levels
* Updated CITATION and manual

MixSIAR 3.1.7 (2016-08-15)
=========================

## MINOR IMPROVEMENTS

* Changed resid.prop prior from dchisqr(3) to dunif(0,20)

## VERY MINOR IMPROVEMENTS

* Updated install instructions in README (install.packages)
* Fixed error in Cladocera vignette (process vs. resid error)
* Fixed link to manual to CRAN page (DESCRIPTION file)

MixSIAR 3.1.6 (2016-06-24)
=========================

## VERY MINOR IMPROVEMENTS

* CRAN submission
* added files for CRAN submission (cran-comments.md, NEWS.md)
* fixed broken links found when checking for CRAN submission

MixSIAR 3.1.5 (2016-06-21)
=========================

## VERY MINOR IMPROVEMENTS

* Changed "Isotopes" to "Tracers" in GUI

MixSIAR 3.1.4 (2016-06-21)
=========================

## BUG FIX

* Fixed bug #72 by using `utils::globalVariables()`, created when making changes to pass R CMD check for CRAN submission

MixSIAR 3.1.3 (2016-06-01)
=========================

## MINOR IMPROVEMENTS

* Tests via testthat package (see "tests/testthat" folder)
* Changes to pass R CMD check for submission to CRAN

MixSIAR 3.1.2 (2016-03-16)
=========================

## MINOR IMPROVEMENTS

* download latest release by default (instead of master)
* updated citation info

MixSIAR 3.1.1 (2016-03-11)
=========================

## NEW FEATURES

* Added vignettes for script examples, see browseVignettes("MixSIAR")

## MINOR IMPROVEMENTS

* Tested package install on Windows, Mac, Linux
* Revised install instructions
* Couple minor bug fixes
* Revised manual

MixSIAR 3.1.0 (2016-03-09)
=========================

## NEW FEATURES

* Converted to R package structure, so can now run MixSIAR with:

`library(MixSIAR)`
`mixsiar_gui()`

MixSIAR 3.0.2 (2015-11-19)
=========================

## BUG FIX

* Problem loading mix/consumer data in the GUI version (#48)

MixSIAR 3.0.1 (2015-10-29)
=========================

## BUG FIX

* Fixed an error that applies with source data = means and the "Resid*Process" error structure (affects Wolves, Killer whale, and Isopod examples)

MixSIAR 3.0.0 (2015-10-29)
=========================

## NEW FEATURES

* error structures now multivariate, with 3 options: "Residual only", "Residual * Process", "Process only (N=1)"
* fixed coding of fixed effects (cases with 2FE, 1FE + 1RE)
* new function plot_prior to plot uninformative vs. informative priors (included in GUI)
* new function calc_area to calculate Brett (2014) normalized surface area
* normalized continuous covariate
* MCMC chains given different initial values (#35)
* added new examples with fatty acids
* new manual
