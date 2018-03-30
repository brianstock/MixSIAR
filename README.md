MixSIAR
=============
[![cran version](http://www.r-pkg.org/badges/version/MixSIAR)](https://cran.r-project.org/package=MixSIAR)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MixSIAR?)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1209993.svg)](https://doi.org/10.5281/zenodo.1209993)

MixSIAR is an R package that helps you create and run Bayesian mixing models to analyze biotracer data (i.e. stable isotopes, fatty acids), following the MixSIAR model framework. Both graphical user interface (GUI) and script versions are available. 

MixSIAR represents a collaborative coding project between the investigators behind MixSIR and SIAR: Brice Semmens, Brian Stock, Eric Ward, Andrew Parnell, Donald Phillips, Andrew Jackson, Jon Moore, Stuart Bearhop, and Richard Inger.

MixSIAR incorporates several years of advances in Bayesian mixing model theory since MixSIR and SIAR, currently:

- Any number of biotracers (examples with 1 isotope, 2 isotope, 8 fatty acids, and 22 fatty acids)
- Source data fit hierarchically within the model
- Source data by categorical covariate (e.g. sources by Region)
- Categorical covariates (up to 2, choice of modeling as random or fixed effects, either nested or independent)
- Continuous covariate (up to 1)
- Error structure options with covariance (Residual * Process, Residual only)
- Concentration dependence
- Ability to plot and include “uninformative”/generalist or informative priors

### QUICK INSTALL (no GUI)

The script version is easier to install and better for repeated analysis. If you want the **script version only** (no GUI):

1. Download and install/update [R].

2. Download and install [JAGS].

3. Open R and run:
```
install.packages("MixSIAR")
library(MixSIAR)
```

The vignettes can be accessed via:
```
browseVignettes("MixSIAR")
```

There is a more extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```

The manual is also available from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

### FULL INSTALL (with GUI)

Getting the GUI running is more work, but can be a nice introduction to MixSIAR. The install instructions are platform-specific:

#### Windows

1. Download and install/update [R].
2. Download and install [JAGS].
3. Open R. 
4. Install GTK+ dependent packages:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2"))
    ```

5. Load `RGtk2`. You will be prompted to install GTK+. **Follow the automatic prompts and do not interrupt the GTK+ installation!**:

    ```
    library(RGtk2)
    ```

6. Restart R and run:

    ```
    install.packages("MixSIAR", dependencies=TRUE)
    ```

7. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

There is an extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```
Alternatively, you can download the manual from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

#### Mac OS X

1. Download and install/update [R].
2. Download and install [JAGS].
3. Open R. 
4. Install GTK+ dependent R packages:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2"))
    ```

5. Close R.
6. Download and install the newest [GTK+ framework].
7. Install the latest X11 application, [xQuartz].
8. Open R and run:

    ```
    install.packages("MixSIAR", dependencies=TRUE)
    ```
9. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

There is an extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```
Alternatively, you can download the manual from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

#### Linux

1. Download and install/update [R].
2. Download and install [JAGS]. Or, from the terminal: `sudo apt-get install jags r-cran-rjags`.
3. Download and install [GTK+ framework]. From the terminal: `sudo apt-get install libgtk2.0-dev`.
4. Check if GTK+ is installed correctly. Open R, install and load the `RGtk2` package with:

    ```
    install.packages("RGtk2")
    library(RGtk2)
    ```

5. Install MixSIAR:

    ```
    install.packages("MixSIAR", dependencies=TRUE)
    ```

6. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

There is an extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```
Alternatively, you can download the manual from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

### FEEDBACK PLEASE!

This software has been improved by the questions, suggestions, and bug reports of the user community. If you have a comment, ideally use the [Issues] page. You can also post to the [SIAR facebook group] or shoot me an email (_b1stock@ucsd.edu_).

### ON CITING MixSIAR:

If you use MixSIAR results in publications, please cite the MixSIAR manual as (similar to how you cite R):

>B. C. Stock and B. X. Semmens (2016). MixSIAR GUI User Manual. Version 3.1. https://github.com/brianstock/MixSIAR. doi:10.5281/zenodo.1209993.

The primary citation for _Bayesian mixing models_ (MixSIR):

>Moore, J. W., & Semmens, B. X. (2008). Incorporating uncertainty and prior information into stable isotope mixing models. Ecology Letters, 11(5), 470-480.

If you are using the _residual error term_ (SIAR):

>Parnell, A. C., Inger, R., Bearhop, S., & Jackson, A. L. (2010). Source partitioning using stable isotopes: coping with too much variation. PLoS One, 5(3), e9672.

If you are using a _hierarchical structure/random effects_:

>Semmens, B. X., Ward, E. J., Moore, J. W., & Darimont, C. T. (2009). Quantifying inter-and intra-population niche variability using hierarchical Bayesian stable isotope mixing models. PLoS One, 4(7), e6187.

If you are using _continuous effects_:

>Francis, T. B., Schindler, D. E., Holtgrieve, G. W., Larson, E. R., Scheuerell, M. D., Semmens, B. X., & Ward, E. J. (2011). Habitat structure determines resource use by zooplankton in temperate lakes. Ecology letters, 14(4), 364-372.

If you are using _source fitting_:

>Ward, E. J., Semmens, B. X., & Schindler, D. E. (2010). Including source uncertainty and prior information in the analysis of stable isotope mixing models. Environmental science & technology, 44(12), 4645-4650.

For a detailed description of the math underlying these models, see:

>Parnell, A. C., Phillips, D. L., Bearhop, S., Semmens, B. X., Ward, E. J., Moore, J. W., Jackson, A. L., Grey, J., Kelley, D. J., & Inger, R. (2013). Bayesian stable isotope mixing models. Environmetrics, 24, 387-399.

For an explanation of the _error structures_ ("Process only" vs. "Resid only" vs. "Process * Resid"), see:

>Stock, B. C., & Semmens, B. X. (2016). Unifying error structures in commonly used biotracer mixing ­models. Ecology, 97(10), 2562–2569.

Finally... yes, a paper introducing MixSIAR is in the works and will be forthcoming shortly.

### INSTALLATION (GitHub):

If for some reason you can't install using `install.packages`, the GitHub version is another option.

#### Windows (GitHub)

1. Download and install/update [R].
2. Download and install [JAGS].
3. (Optional) If you want to build the vignettes, install [pandoc] or [R Studio].
4. Open R. 
5. Install GTK+ dependent packages:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2", "devtools"))
    ```

6. Load `RGtk2`. You will be prompted to install GTK+. **Follow the automatic prompts and do not interrupt the GTK+ installation!**:

    ```
    library(RGtk2)
    ```

7. Restart R and run:

    ```
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE) # FALSE if no pandoc/R Studio
    ```

8. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

#### Mac OS X (GitHub)

1. Download and install/update [R].
2. Download and install [JAGS].
3. (Optional) If you want to build the vignettes, install [pandoc] or [R Studio].
4. Open R. 
5. Install GTK+ dependent R packages:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2", "devtools"))
    ```

6. Close R.
7. Download and install the newest [GTK+ framework].
8. Install the latest X11 application, [xQuartz].
9. Open R and run:

    ```
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE) # FALSE if no pandoc/R Studio
    ```
10. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

#### Linux (GitHub)

1. Download and install/update [R].
2. Download and install [JAGS]. Or, from the terminal: `sudo apt-get install jags r-cran-rjags`.
3. Download and install [GTK+ framework]. From the terminal: `sudo apt-get install libgtk2.0-dev`.
4. (Optional) If you want to build the vignettes, install [pandoc] or [R Studio].
5. Check if GTK+ is installed correctly. Open R, install and load the `RGtk2` package with:

    ```
    install.packages("RGtk2")
    library(RGtk2)
    ```

6. Install and load devtools, then install MixSIAR:

    ```
    install.packages("devtools")
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE) # FALSE if no pandoc and pandoc-citeproc
    ```

7. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

[pandoc]:https://github.com/jgm/pandoc/releases/
[GTK+ framework]:http://r.research.att.com/#other
[xQuartz]:http://xquartz.macosforge.org/landing/
[R Studio]:https://www.rstudio.com/products/rstudio/download/
[R]:https://cran.r-project.org/bin/
[JAGS]:http://mcmc-jags.sourceforge.net/
[Issues]:https://github.com/brianstock/MixSIAR/issues
[SIAR Facebook group]:https://www.facebook.com/pages/SIAR-Stable-Isotope-Analysis-in-R/148501811896914

