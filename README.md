MixSIAR
=============

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

### INSTALLATION:

The MixSIAR package can be run as a GUI with `mixsiar_gui()`, or as a sequence of R commands (`.R` script). The GUI depends on the `gWidgetsRGtk2` package in R. MixSIAR GUI has been tested and runs on [Windows](#windows), [Mac OS X](#mac-os-x), and [Linux](#linux).

#### Windows

1. Download and install [JAGS].
2. Open R. 
3. Install GTK+ dependent packages. You will be prompted to install GTK+. **Follow the automatic prompts and do not interrupt the GTK+ installation!**:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2", "devtools"))
    ```

4. Restart R and run:

    ```
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE)
    ```

5. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

#### Mac OS X

1. Download and install [JAGS].
2. Open R. 
3. Install GTK+ dependent R packages:

    ```
    install.packages(c("gWidgets", "RGtk2", "gWidgetsRGtk2", "devtools"))
    ```

4. Close R.
5. Download and install the newest GTK+ framework: http://r.research.att.com/#other.
6. Install the latest X11 application (xQuartz): http://xquartz.macosforge.org/landing/.
7. Open R and run:

    ```
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE)
    ```
8. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

#### Linux

1. Download and install [JAGS]. From the terminal: `sudo apt-get install jags r-cran-rjags`.
2. Download and install GTK+. From the terminal: `sudo apt-get install libgtk2.0-dev`. You can also see http://www.gtk.org/download/index.php.
3. Check if GTK+ is installed correctly. Open R, install and load the `RGtk2` package with:

    ```
    install.packages("RGtk2")
    library(RGtk2)
    ```

4. Install and load devtools, then install MixSIAR:

    ```
    install.packages("devtools")
    library(devtools)
    devtools::install_github("brianstock/MixSIAR",
                             dependencies = TRUE, 
                             build_vignettes = TRUE)
    ```

5. Load MixSIAR and run GUI:

    ```
    library(MixSIAR)
    mixsiar_gui()
    ```

### FEEDBACK PLEASE!

This software has been improved by the questions, suggestions, and bug reports of the user community. If you have a comment, ideally use the [Issues] page. You can also post to the [SIAR facebook group] or shoot me an email (_b1stock@ucsd.edu_).

### ON CITING MixSIAR:

If you use MixSIAR GUI results in publications, please cite the MixSIAR GUI manual as (similar to how you cite R):

>B. C. Stock and B. X. Semmens (2013). MixSIAR GUI User Manual, version ##. https://github.com/brianstock/MixSIAR.

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

Finally... yes, a paper introducing MixSIAR is in the works and will be forthcoming shortly.

[JAGS]:http://mcmc-jags.sourceforge.net/
[Issues]:https://github.com/brianstock/MixSIAR/issues
[SIAR Facebook group]:https://www.facebook.com/pages/SIAR-Stable-Isotope-Analysis-in-R/148501811896914