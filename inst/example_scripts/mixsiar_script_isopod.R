# Brian Stock
# March 8, 2016
# Script file to run isopod example without GUI

########################################################################################
# Cladocera example (22 fatty acids, run each consumer individually as fixed effect)

# The Cladocera Example is from Galloway et al. (2014) and demonstrates MixSIAR applied
# to a 22-dimensional fatty acid dataset. Here the 14 mixture datapoints are Cladocera (wa-
# ter flea) fatty acid profiles from 6 lakes in Finland over 2 seasons. Besides the high di-
# mensionality, the other difference with this analysis is that we fit each mixture datapoint
# individually, because there is no clear covariate structure (some sites have 2 seasons, some
# have 1, some sites in the same lake). We do this by creating an “id” column and treating
# “id” as a fixed effect.

# • 22 biotracers (carbon 14.0, 16.0, 16.1ω9, 16.1ω7, 16.2ω4, 16.3ω3, 16.4ω3, 17.0, 18.0,
# 18.1ω9, 18.1ω7, 18.2ω6, 18.3ω6, 18.3ω3, 18.4ω3, 18.5ω3, 20.0, 22.0, 20.4ω6, 20.5ω3,
# 22.6ω3, BrFA)
# 43• Mix datapoints fit independently
# • Source data as means and SDs
# Fatty acid data greatly increase the number of biotracers beyond the typical 2 stable iso-
# topes, δ 13 C and δ 15 N, which gives the mixing model power to resolve more sources. We
# caution, however, that it is not only a matter of # biotracers > # sources + 1. As the num-
# ber of sources increases, the “uninformative” prior α = 1 has greater influence, illus-
# trated in Figure 22.

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "isopod_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("c16.4w3","c18.2w6","c18.3w3","c18.4w3","c20.4w6","c20.5w3","c22.5w3","c22.6w3"),
					 factors="Site",
					 fac_random=FALSE,
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
jags.1 <- run_model(run="test", mix, source, discr, model_filename,
	                alpha.prior=1, resid_err, process_err)

## "normal" took my laptop ~60 minutes to run
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename,
#	                alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)

