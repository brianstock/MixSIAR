# Brian Stock
# March 8, 2016
# Script file to run storm-petrel example without GUI

##########################################################################################
# Storm-petrel example (movement instead of diet, Region = fixed effect)

# The “Storm-petrel Example” is based on Bicknell et al. (2014) and shows MixSIAR in a non-
# diet application: movement. Here the mixture data are juvenile (non-breeding) birds and
# the source data are adult (breeding) birds from 3 breeding colonies. Discrimination is
# set to zero, making the reasonable assumptions that 1) juvenile and adult birds feeding
# in the same region will look isotopically identical, and 2) discrimination is the same for
# juvenile and adult birds.

# We have:
#   2 biotracers (δ 13 C, δ 15 N)
#   1 fixed effect (Region)
#   Raw source data

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "stormpetrel_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("d13C","d15N"),
					 factors="Region",
					 fac_random=FALSE,
					 fac_nested=FALSE,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "stormpetrel_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="raw",
						   mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "stormpetrel_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

# Make isospace plot
plot_data(filename="isospace_plot",
		  plot_save_pdf=TRUE,
		  plot_save_png=FALSE,
		  mix,source,discr)

# Calculate standardized convex hull area
if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

# Plot your prior
plot_prior(alpha.prior=1,source)

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("test" first, then "normal")
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)

