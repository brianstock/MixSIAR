MixSIAR GUI
=============

MixSIAR GUI is a graphical user interface (GUI) that allows you to analyze stable isotope data using the MixSIAR model framework. MixSIAR is a Bayesian mixing model that uses stable isotope data to estimate the proportions of source (prey) contributions to a mixture (consumer).

MixSIAR represents a collaborative coding project between the investigators behind MixSIR and SIAR: Brice Semmens, Brian Stock, Eric Ward, Andrew Parnell, Donald Phillips, Andrew Jackson, Jon Moore, Stuart Bearhop, and Richard Inger.

The MixSIAR GUI incorporates several years of advances in Bayesian mixing model theory since MixSIR and SIAR, currently:

- Any number of isotope/tracer values (usually 1, 2, or 3, but tested up to 20)
- Source fitting (both for raw source data and means+SD data)
- Hierarchical random/fixed effects (up to 2, either hierarchical/nested or not)
- Continuous covariates (up to 1, theoretically any number)
- Option to include Individual random effects
- Error structure options: residual error (SIAR) or not (MixSIR)
- Concentration dependence

### UPDATES (Aug 2014):

Important changes to the MixSIAR GUI 1.0 code (Oct 2013) include:

1. Several bug/error fixes
2. New scripted version, check out "mixsiar_script.r"
3. Altered the GUI function call, now:
```r
source("mixsiar_gui.r")
mixsiar_gui()
```
4. Added option for fixed effects
5. Separated OUTPUT button from RUN MODEL button
6. Simplified MCMC options

### TO DOWNLOAD AND INSTALL:

On the right hand side of this page, click the **"Download ZIP"** button. See the user manual (MixSIAR GUI User Manual 1.0.pdf) for installation instructions and a walk-through of MixSIAR via four example datasets.

### FEEDBACK PLEASE!

This software has been improved by the questions, suggestions, and bug reports of the user community. If you have a comment, ideally use the _issues_ (exclamation point in a circle) or _wiki_ (looks like a book) tabs on the right hand side of this page. You can also post to the [SIAR facebook group] or shoot me an email (_b1stock@ucsd.edu_).

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

[SIAR Facebook group]:https://www.facebook.com/pages/SIAR-Stable-Isotope-Analysis-in-R/148501811896914