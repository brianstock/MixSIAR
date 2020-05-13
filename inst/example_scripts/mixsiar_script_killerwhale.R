# Brian Stock
# March 8, 2016
# Script file to run killer whale example without GUI

################################################################################
# Killer whale example (informative prior)

# The Killer Whale Example demonstrates the difference informative priors make, and
# illustrates how to construct them.

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "killerwhale_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("d13C","d15N"),
					 factors=NULL,
					 fac_random=NULL,
					 fac_nested=NULL,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "killerwhale_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="means",
						   mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "killerwhale_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

# Make isospace plot
plot_data(filename="isospace_plot",
		  plot_save_pdf=TRUE,
		  plot_save_png=FALSE,
		  mix,source,discr)

# Calculate standardized convex hull area
if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

################################################################################
# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
################################################################################
plot_prior(alpha.prior=1, source, filename = "prior_plot_kw_uninf")

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model_kw_uninf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("very long" took ~5 min)
jags.uninf <- run_model(run="test",mix,source,discr,model_filename,alpha.prior = 1, resid_err, process_err)
# jags.uninf <- run_model(run="very long",mix,source,discr,model_filename,alpha.prior = 1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf, mix, source)

################################################################################
# # INFORMATIVE prior (construct alpha from fecal data)
################################################################################
# Let's say we have 14 fecal diet samples that we use to construct alphas...
#   useful in separating some of the sources.

# Our 14 fecal samples were 10, 1, 0, 0, 3
kw.alpha <- c(10,1,0,0,3)

# Generate alpha hyperparameters scaling sum(alpha)=n.sources
kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)

# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
kw.alpha[which(kw.alpha==0)] <- 0.01

# Plot your informative prior
plot_prior(alpha.prior=kw.alpha,
		   source=source,
		   plot_save_pdf=TRUE,
		   plot_save_png=FALSE,
		   filename="prior_plot_kw_inf")

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model_kw_inf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("very long" took ~5 min)
jags.inf <- run_model(run="test",mix,source,discr,model_filename,alpha.prior=kw.alpha)
# jags.inf <- run_model(run="very long",mix,source,discr,model_filename,alpha.prior=kw.alpha)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.inf, mix, source)
