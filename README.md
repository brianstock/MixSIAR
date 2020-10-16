MixSIAR
=============
[![cran version](https://www.r-pkg.org/badges/version/MixSIAR)](https://cran.r-project.org/package=MixSIAR)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/MixSIAR?)](https://github.com/r-hub/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1209993.svg)](https://doi.org/10.5281/zenodo.1209993)

MixSIAR is an R package that helps you create and run Bayesian mixing models to analyze biotracer data (i.e. stable isotopes, fatty acids), following the [MixSIAR model framework](https://peerj.com/articles/5096/). MixSIAR represents a collaborative coding project between the investigators behind MixSIR, SIAR, and IsoSource: Brice Semmens, Brian Stock, Eric Ward, Andrew Parnell, Donald Phillips, and Andrew Jackson.

MixSIAR incorporates several years of advances in Bayesian mixing model theory since MixSIR and SIAR, currently:

- Any number of biotracers (examples with 1 isotope, 2 isotope, 8 fatty acids, and 22 fatty acids)
- Source data fit hierarchically within the model
- Source data by categorical covariate (e.g. sources by Region)
- Categorical covariates (up to 2, choice of modeling as random or fixed effects, either nested or independent)
- Continuous covariate (up to 1)
- Error structure options with covariance (Residual * Process, Residual only)
- Concentration dependence
- Plot and include "uninformative"/generalist or informative priors
- Fit multiple models and compare relative support using LOO/WAIC weights

For details, please see the [MixSIAR paper](https://peerj.com/articles/5096/):

- Full description of equations
- Advice/explanation on 4 common issues (error structures, priors, combining sources, covariates)
- Case study highlighting new functionality (model selection with LOO/WAIC weights)

> Stock BC, Jackson AL, Ward EJ, Parnell AC, Phillips DL, Semmens BX. 2018. Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ 6:e5096 https://doi.org/10.7717/peerj.5096

## Installation

The GUI has been removed from the CRAN version of MixSIAR (if desired, see [MixSIARgui](https://github.com/brianstock/MixSIARgui) on GitHub). Running MixSIAR with scripts is easier to install and better for repeated analysis.

1. Download and install/update [R](https://cran.r-project.org/).

2. Download and install [JAGS](http://mcmc-jags.sourceforge.net/).

3. Open R and run:
```
install.packages("MixSIAR", dependencies=TRUE)
library(MixSIAR)
```

If you want the latest changes and bug fixes not yet on CRAN, you can install the GitHub version:
```
remotes::install_github("brianstock/MixSIAR", dependencies=T)
```

## Tutorial

We suggest walking through the [vignettes](http://brianstock.github.io/MixSIAR/articles/index.html) to familiarize yourself with MixSIAR.

There is also an extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```

Alternatively, you can download the manual from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

Clean, runnable `.R` scripts for each vignette are also available in the `example_scripts` folder of the `MixSIAR` package install:
```
library(MixSIAR)
mixsiar.dir <- find.package("MixSIAR")
file.path(mixsiar.dir, "example_scripts")
```

You can then run the Wolves example script with:
```
setwd("choose/where/to/save/output")
source(file.path(mixsiar.dir, "example_scripts", "mixsiar_script_wolves.R"))
```

### Feedback

This software has been improved by the questions, suggestions, and bug reports of the user community. If you have a comment, please use the [Issues](https://github.com/brianstock/MixSIAR/issues) page.

### Citing MixSIAR:

If you use MixSIAR results in publications, please cite the MixSIAR manual as (similar to how you cite R):

> Stock BC and Semmens BX. 2016. MixSIAR GUI User Manual. Version 3.1. https://github.com/brianstock/MixSIAR. doi:10.5281/zenodo.1209993.

The MixSIAR model framework is described in:

> Stock BC, Jackson AL, Ward EJ, Parnell AC, Phillips DL, Semmens BX. 2018. Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ 6:e5096 https://doi.org/10.7717/peerj.5096

The primary citation for _Bayesian mixing models_ (MixSIR):

> Moore, J. W., & Semmens, B. X. (2008). Incorporating uncertainty and prior information into stable isotope mixing models. Ecology Letters, 11(5), 470-480.

If you are using the _residual error term_ (SIAR):

> Parnell, A. C., Inger, R., Bearhop, S., & Jackson, A. L. (2010). Source partitioning using stable isotopes: coping with too much variation. PLoS One, 5(3), e9672.

If you are using a _hierarchical structure/random effects_:

> Semmens, B. X., Ward, E. J., Moore, J. W., & Darimont, C. T. (2009). Quantifying inter-and intra-population niche variability using hierarchical Bayesian stable isotope mixing models. PLoS One, 4(7), e6187.

If you are using _continuous effects_:

> Francis, T. B., Schindler, D. E., Holtgrieve, G. W., Larson, E. R., Scheuerell, M. D., Semmens, B. X., & Ward, E. J. (2011). Habitat structure determines resource use by zooplankton in temperate lakes. Ecology letters, 14(4), 364-372.

If you are using _source fitting_:

> Ward, E. J., Semmens, B. X., & Schindler, D. E. (2010). Including source uncertainty and prior information in the analysis of stable isotope mixing models. Environmental science & technology, 44(12), 4645-4650.

For a detailed description of the math underlying these models, see:

> Parnell, A. C., Phillips, D. L., Bearhop, S., Semmens, B. X., Ward, E. J., Moore, J. W., Jackson, A. L., Grey, J., Kelley, D. J., & Inger, R. (2013). Bayesian stable isotope mixing models. Environmetrics, 24, 387-399.

For an explanation of the _error structures_ ("Process only" vs. "Resid only" vs. "Process * Resid"), see:

> Stock, B. C., & Semmens, B. X. (2016). Unifying error structures in commonly used biotracer mixing models. Ecology, 97(10), 2562–2569.

