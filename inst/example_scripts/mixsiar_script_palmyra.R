# Brian Stock
# March 8, 2016
# Script file to run palmyra example without GUI

#########################################################################################
# Palmyra example ("Taxa" = fixed effect, "raw" source data)

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "palmyra_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("d13C","d15N"),
					 factors="Taxa",
					 fac_random=FALSE,
					 fac_nested=FALSE,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "palmyra_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="raw",
						   mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "palmyra_discrimination.csv", package = "MixSIAR")
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

# Run the JAGS model ("test" first, then "short")
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)
#jags.1 <- run_model(run="short", mix, source, discr, model_filename, alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)

