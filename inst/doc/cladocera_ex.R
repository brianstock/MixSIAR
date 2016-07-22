## ------------------------------------------------------------------------
library(MixSIAR)
mixsiar.dir <- find.package("MixSIAR")
paste0(mixsiar.dir,"/example_scripts")

## ---- eval=FALSE---------------------------------------------------------
#  source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_cladocera.R"))

## ------------------------------------------------------------------------
library(MixSIAR)

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
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

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
source.filename <- system.file("extdata", "cladocera_sources.csv", package = "MixSIAR")

source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="means",
						   mix)

## ------------------------------------------------------------------------
# Replace the system.file call with the path to your file
discr.filename <- system.file("extdata", "cladocera_discrimination.csv", package = "MixSIAR")

discr <- load_discr_data(filename=discr.filename, mix)

## ---- eval=FALSE---------------------------------------------------------
#  # default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
#  plot_prior(alpha.prior=1,source)

## ---- eval=FALSE---------------------------------------------------------
#  # Write the JAGS model file
#  model_filename <- "MixSIAR_model.txt"
#  resid_err <- FALSE
#  process_err <- TRUE
#  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

## ---- eval=FALSE---------------------------------------------------------
#  jags.1 <- run_model(run="test", mix, source, discr, model_filename,
#                      alpha.prior = 1, resid_err, process_err)

## ---- eval=FALSE---------------------------------------------------------
#  jags.1 <- run_model(run="normal", mix, source, discr, model_filename,
#                      alpha.prior = 1, resid_err, process_err)

## ---- eval=FALSE---------------------------------------------------------
#  output_JAGS(jags.1, mix, source, output_options)

