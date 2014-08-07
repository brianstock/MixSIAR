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

### TO DOWNLOAD AND INSTALL:

Just above the blue bar, there is a row with "commits", "branches", "releases", and "contributors". Click on **"releases"** and you should see the latest version. Below the release notes, click on the "Source code (zip)" button to download the release. The current user manual will be in the .zip folder with installation instructions and a walk-through of MixSIAR via several example datasets.

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