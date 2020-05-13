# Brian Stock
# March 8, 2016
# Script file to run isopod example without GUI

################################################################################
# Isopod example (8 fatty acids, "Site" as random effect)

# The Isopod Example is from Galloway et al. 2014 and demonstrates MixSIAR applied
# to an 8-dimensional fatty acid dataset. Here the mixture data are isopod
# polyunsaturated fatty acid (PUFA) profiles, with 5 replicates at each of 6 sites
# in Puget Sound, WA:
#   - 8 biotracers (carbon 16.4w3, 18.2w6, 18.3w3, 18.4w3, 20.4w6, 20.5w3, 22.5w3, 22.6w3)
#   - 1 random effect (Site)
#   - source data as means and SDs

# Here we treat Site as a random effect. This makes sense if we are interested in
# the overall population and think of Site as a nuisance factor. Fitting Site as
# a fixed effect would make more sense if we were interested specifically in the
# diet at each Site, as opposed to the overall population diet and variability
# between Sites. This differs from the analysis in Galloway et al. 2014.

# Fatty acid data greatly increase the number of biotracers beyond the typical 2
# stable isotopes, d13C and d15N, which gives the mixing model power to resolve
# more sources. We caution, however, that using fatty acid data is not a panacea
# for the "underdetermined" problem (# sources > # biotracers + 1). As the number
# of sources increases, the "uninformative" prior alpha = 1 has greater influence,
# illustrated in the Cladocera Example, Figure 22 in the Manual.

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "isopod_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("c16.4w3","c18.2w6","c18.3w3","c18.4w3","c20.4w6","c20.5w3","c22.5w3","c22.6w3"),
					 factors="Site",
					 fac_random=TRUE,
					 fac_nested=FALSE,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "isopod_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="means",
						   mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "isopod_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

# Make isospace plot
plot_data(filename="isospace_plot",
		  plot_save_pdf=TRUE,
		  plot_save_png=FALSE,
		  mix,source,discr)

# Plot your prior
plot_prior(alpha.prior=1,source)

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("test" first, then "normal")
jags.1 <- run_model(run="test", mix, source, discr, model_filename,alpha.prior=1)

## "normal" took my laptop ~60 minutes to run
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename,alpha.prior=1)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)

