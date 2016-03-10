# Brian Stock
# March 8, 2016
# Script file to run cladocera example without GUI

########################################################################################
# Cladocera example (22 fatty acids, run each consumer individually as fixed effect)

# The Cladocera Example is from Galloway et al. (2014) and demonstrates MixSIAR applied
# to a 22-dimensional fatty acid dataset. Here the 14 mixture datapoints are Cladocera (wa-
# ter flea) fatty acid profiles from 6 lakes in Finland over 2 seasons. Besides the high di-
# mensionality, the other difference with this analysis is that we fit each mixture datapoint
# individually, because there is no clear covariate structure (some sites have 2 seasons, some
# have 1, some sites in the same lake). We do this by creating an "id" column and treating
# "id" as a fixed effect.

#   - 22 biotracers (carbon 14.0, 16.0, 16.1w9, 16.1w7, 16.2w4, 16.3w3, 16.4w3, 17.0, 18.0,
#                  18.1w9, 18.1w7, 18.2w6, 18.3w6, 18.3w3, 18.4w3, 18.5w3, 20.0, 22.0, 20.4w6,
#                  20.5w3, 22.6w3, BrFA)
#   - Mix datapoints fit independently ("process error" only, MixSIR)
#   - Source data as means and SDs

# Here we fit the "process error" model of MixSIR, which we MUST do when we only
# have one mix datapoint (or here, one mix datapoint per fixed effect). With only
# one datapoint, there is no information to estimate an additional mixture variance
# term, so we have to assume a fixed variance based on the variance of the sources
# (see Moore and Semmens 2008).

# Here we treat “id” as a fixed effect—this will estimate the diet of each mixture
# data point separately (sample size of 1). This makes sense to do when you think
# there will be clear differences between sites/seasons/etc., but only have 1 or 2
# points from each site/season (i.e, you don’t have enough data to estimate a site/
# season effect). If you are interested in the site/season effect, you need
# replicates within each site/season, and then it is best to fit site/season as
# a fixed or random effect.

# Fatty acid data greatly increase the number of biotracers beyond the typical 2
# stable isotopes, d13C and d15N, which gives the mixing model power to resolve
# more sources. We caution, however, that using fatty acid data is not a panacea
# for the "underdetermined" problem (# sources > # biotracers + 1). As the number
# of sources increases, the "uninformative" prior alpha = 1 has greater influence,
# illustrated in the Cladocera Example, Figure 22 in the Manual.

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "cladocera_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("c14.0","c16.0","c16.1w9","c16.1w7","c16.2w4",
					 	         "c16.3w3","c16.4w3","c17.0","c18.0","c18.1w9",
					 	         "c18.1w7","c18.2w6","c18.3w6","c18.3w3","c18.4w3",
					 	         "c18.5w3","c20.0","c22.0","c20.4w6","c20.5w3",
					 	         "c22.6w3","BrFA"),
					 factors="id",
					 fac_random=FALSE,
					 fac_nested=FALSE,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "cladocera_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="means",
						   mix)

# Load discrimination/TDF data

# Note that Galloway et al. [5] conducted feeding trials to create a “resource
# library”. In the mixing model, the sources are actually consumers fed exclusively
# each of the sources. This allowed them to set the discrimination = 0 (see
# “cladocera_discrimination.csv”).

discr.filename <- system.file("extdata", "cladocera_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

# DON'T make isospace plot - MixSIAR makes plots of each pairwise combination
# of tracers.  Here we have 22, so that would be '22 choose 2' = 231 biplots.

# Plot your prior
plot_prior(alpha.prior=1,source)

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("test" first, then "normal")
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)

## "normal" took my laptop ~30 minutes to run
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
# Note that since we fit “id” as a fixed effect, there is no inference on diet at
# the overall population level (no p.global). You should see posterior plots for
# all 14 mixture samples.
output_JAGS(jags.1, mix, source)



