library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename, iso_names=c("d13C","d15N"), factors=c("Region","Pack"), fac_random=c(TRUE,TRUE), fac_nested=c(FALSE,TRUE), cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "wolves_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename, source_factors="Region", conc_dep=FALSE, data_type="means", mix)

# Load discrimination data
discr.filename <- system.file("extdata", "wolves_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

# Plot data
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

# Calculate convex hull area
calc_area(source=source,mix=mix,discr=discr)

# Plot prior
plot_prior(alpha.prior=1,source)

# Define model structure
model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run model using JAGS
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)
# jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)

# Choose output options
output_options <- list(summary_save = TRUE,                 # Save the summary statistics as a txt file?
                       summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
                       sup_post = FALSE,                       # Suppress posterior density plot output in R?
                       plot_post_save_pdf = TRUE,              # Save posterior density plots as pdfs?
                       plot_post_name = "posterior_density",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                       sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                       plot_pairs_save_pdf = TRUE,             # Save pairs plot as pdf?
                       plot_pairs_name = "pairs_plot",         # If yes, specify the base file name (.pdf/.png will be appended later)
                       sup_xy = TRUE,                         # Suppress xy/trace plot output in R?
                       plot_xy_save_pdf = FALSE,                # Save xy/trace plot as pdf?
                       plot_xy_name = "xy_plot",               # If yes, specify the base file name (.pdf/.png will be appended later)
                       gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                       heidel = FALSE,                          # Calculate Heidelberg-Welch diagnostic test?
                       geweke = TRUE,                          # Calculate Geweke diagnostic test?
                       diag_save = TRUE,                       # Save the diagnostics as a txt file?
                       diag_name = "diagnostics",              # If yes, specify the base file name (.txt will be appended later)
                       indiv_effect = FALSE,                   # Is Individual a random effect in the model? (already specified)
                       plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                       plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                       plot_xy_save_png = FALSE)               # Save xy/trace plot as png?

# Create diagnostics, summary statistics, and posterior plots
output_JAGS(jags.1, mix, source, output_options)

