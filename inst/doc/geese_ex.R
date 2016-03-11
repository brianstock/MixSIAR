## ------------------------------------------------------------------------
library(MixSIAR)
mixsiar.dir <- find.package("MixSIAR")
paste0(mixsiar.dir,"/example_scripts")

## ---- eval=FALSE---------------------------------------------------------
#  source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_geese.R"))

## ------------------------------------------------------------------------
library(MixSIAR)

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
mix.filename <- system.file("extdata", "geese_consumer.csv", package = "MixSIAR")

mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("d13C","d15N"),
					 factors="Group",
					 fac_random=FALSE,
					 fac_nested=FALSE,
					 cont_effects=NULL)

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
source.filename <- system.file("extdata", "geese_sources.csv", package = "MixSIAR")

source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=TRUE,
						   data_type="means",
						   mix)

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
discr.filename <- system.file("extdata", "geese_discrimination.csv", package = "MixSIAR")

discr <- load_discr_data(filename=discr.filename, mix)

## ---- eval=FALSE---------------------------------------------------------
#  # Make an isospace plot
#  plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

## ------------------------------------------------------------------------
# Calculate the convex hull area, standardized by source variance
calc_area(source=source,mix=mix,discr=discr)

## ---- eval=FALSE---------------------------------------------------------
#  # default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
#  plot_prior(alpha.prior=1,source)

## ---- eval=FALSE---------------------------------------------------------
#  # Write the JAGS model file
#  model_filename <- "MixSIAR_model.txt"
#  resid_err <- TRUE
#  process_err <- FALSE
#  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

## ---- eval=FALSE---------------------------------------------------------
#  jags.1 <- run_model(run="test", mix, source, discr, model_filename,
#                      alpha.prior = 1, resid_err, process_err)

## ---- eval=FALSE---------------------------------------------------------
#  jags.1 <- run_model(run="short", mix, source, discr, model_filename,
#                      alpha.prior = 1, resid_err, process_err)

## ---- eval=FALSE---------------------------------------------------------
#  output_JAGS(jags.1, mix, source, output_options)

