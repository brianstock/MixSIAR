# Brian Stock
# July 2014

# MixSIAR script
# Use MixSIAR without the GUI, using the 4 example datasets

# First load the package dependencies (note that 'gWidgetsRGtk2' is NOT loaded - we're not making the GUI)
require(ggplot2)
require(R2jags)
require(MASS)
require(RColorBrewer)
require(reshape)

# Next we clean up the workspace
rm(list=ls()) 	# deletes everything previously in the workspace
runif(1)		# generates one random number (else JAGS can complain)

# Load all MixSIAR functions into the workspace
source("load_mix_data.r")
source("load_source_data.r")
source("load_discr_data.r")
source("plot_data.r")
source("write_JAGS_model.r")
source("run_model.r")
source("output_JAGS.r")
source("plot_continuous_var.r")

#####################################################################################
# Load mixture data, i.e. your:
#	Consumer isotope values (trophic ecology / diet)
#	Mixed sediment/water tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with mix/consumer data
# 'iso_names' - column headings of the tracers/isotopes you'd like to use
# 'random_effects' - column headings of any random effects
# 'cont_effects' - column headings of any continuous effects

# NOTE: 'iso_names', 'random_effects', and 'cont_effects' can be a subset of your data columns
#   i.e. have 3 isotopes in file but only want MixSIAR to use 2,
#   have data by Region and Pack but only want MixSIAR to use Region

# Wolves example (hierarchical/nested random effects)
mix <- load_mix_data(filename="wolves_consumer.csv", iso_names=c("d13C","d15N"), random_effects=c("Region","Pack"), cont_effects=NULL, fixed_effects=NULL)

# Lake example (continuous effect)
# mix <- load_mix_data(filename="lake_consumer.csv", iso_names=c("d13C","d15N"), random_effects=NULL, cont_effects="Secchi.Mixed", fixed_effects=NULL)

# Geese example (concentration dependence)
# mix <- load_mix_data(filename="geese_consumer.csv", iso_names=c("d13C","d15N"), random_effects="Group", cont_effects=NULL, fixed_effects=NULL)

# Palmyra example
# mix <- load_mix_data(filename="palmyra_consumer.csv", iso_names=c("d13C","d15N"), random_effects="Taxa", cont_effects=NULL, fixed_effects=NULL)

# Storm-petrel example (fixed effect)
# mix <- load_mix_data(filename="7_mix.csv", iso_names=c("d13C","d15N"), random_effects=NULL, cont_effects=NULL, fixed_effects="Region")

# 1-iso example
# mix <- load_mix_data(filename="13_mix.csv", iso_names="d13C", random_effects=NULL, cont_effects=NULL, fixed_effects=NULL)

#####################################################################################
# Load source data, i.e. your:
#	Source isotope values (trophic ecology / diet)
#	Sediment/water source tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with source data
# 'source_random_effects' - column headings of random effects you have source data by
# 'conc_dep' - TRUE or FALSE, do you have concentration dependence data in the file?
# 'data_type' - "means" or "raw", is your source data in the means+SD format, or do you have raw data

# Wolves example
source <- load_source_data(filename="wolves_sources.csv", source_random_effects="Region", conc_dep=FALSE, data_type="means", mix)    

# Lake example
# source <- load_source_data(filename="lake_sources.csv", source_random_effects=NULL, conc_dep=FALSE, data_type="raw", mix)    

# Geese example
# source <- load_source_data(filename="geese_sources.csv", source_random_effects=NULL, conc_dep=TRUE, data_type="means", mix)    

# Palmyra example
# source <- load_source_data(filename="palmyra_sources.csv", source_random_effects=NULL, conc_dep=FALSE, data_type="raw", mix)    

# Storm-petrel example
# source <- load_source_data(filename="7_sources.csv", source_random_effects=NULL, conc_dep=FALSE, data_type="raw", mix)  

# 1-iso example
# source <- load_source_data(filename="13_sources.csv", source_random_effects=NULL, conc_dep=FALSE, data_type="raw", mix)    

#####################################################################################
# Load discrimination data, i.e. your:
#	Trophic Enrichment Factor (TEF) / fractionation values (trophic ecology / diet)
#	xxxxxxxx (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with discrimination data

# Wolves example
discr <- load_discr_data(filename="wolves_discrimination.csv", mix)

# Lake example
# discr <- load_discr_data(filename="lake_discrimination.csv", mix)

# Geese example
# discr <- load_discr_data(filename="geese_discrimination.csv", mix)

# Palmyra example
# discr <- load_discr_data(filename="palmyra_discrimination.csv", mix)

# Storm-petrel example
# discr <- load_discr_data(filename="7_discrimination.csv", mix) 

# 1-iso example
# discr <- load_discr_data(filename="13_discrimination.csv", mix) 

#####################################################################################
# Make isospace plot
# Are the data loaded correctly?
# Is your mixture data in the source polygon?
# Are one or more of your sources confounded/hidden?

# 'filename' - name you'd like MixSIAR to save the isospace plot as (extension will be added automatically)
# 'plot_save_pdf' - TRUE or FALSE, should MixSIAR save the plot as a .pdf?
# 'plot_save_png' - TRUE or FALSE, should MixSIAR save the plot as a .png?

plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

#####################################################################################
# Write JAGS model file
# Model will be saved as 'model_filename' ("MixSIAR_model.txt" is default, but may want to change if in a loop)

# 'indiv_effect' - TRUE or FALSE, do you want Individual as a random effect in the model?
# 'model_filename' - don't need to change, unless you are creating many different model files
# 'nested' - TRUE or FALSE, if you have 2 random effects, is the 2nd nested within the 1st?
# 'resid_err' - TRUE or FALSE, do you want to include residual error in the model (FALSE = MixSIR, TRUE = SIAR)?

# Wolves example
model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
indiv_effect <- TRUE	               # Include Individual as a random effect in the model?
nested <- TRUE                          # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)

# Lake example
# model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
# indiv_effect <- TRUE                    # Include Individual as a random effect in the model?
# nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
# write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)

# Geese and Palmyra examples
# model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
# indiv_effect <- FALSE                   # Include Individual as a random effect in the model?
# nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
# write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)

#####################################################################################
# Run model
# JAGS output will be saved as 'jags.1'

# MCMC run options:
# run <- "test"       	# list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short" 	# list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short"     	# list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal"      	# list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long"  		# list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long" 	# list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme"    	# list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

# Can also set custom MCMC parameters
# run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

# Wolves, Palmyra, Geese examples
jags.1 <- run_model(run="short", indiv_effect,mix,source,discr,model_filename)
# test
# Lake example
# jags.1 <- run_model(run="normal", indiv_effect,mix,source,discr,model_filename)

# You can use 'test' first to check if 1) the data are loaded correctly and 2) the model is specified correctly
# jags.1 <- run_model(run="test", indiv_effect,mix,source,discr,model_filename)

#####################################################################################
# Process JAGS output

# All examples
output_options <- list(summary_save = TRUE,                 # Save the summary statistics as a txt file?
                    summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
                    sup_post = FALSE,                       # Suppress posterior density plot output in R?
                    plot_post_save_pdf = TRUE,              # Save posterior density plots as pdfs?
                    plot_post_name = "posterior_density",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                    sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                    plot_pairs_save_pdf = TRUE,             # Save pairs plot as pdf?
                    plot_pairs_name = "pairs_plot",         # If yes, specify the base file name (.pdf/.png will be appended later)
                    sup_xy = FALSE,                         # Suppress xy/trace plot output in R?
                    plot_xy_save_pdf = TRUE,                # Save xy/trace plot as pdf?
                    plot_xy_name = "xy_plot",               # If yes, specify the base file name (.pdf/.png will be appended later)
                    gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                    heidel = FALSE,                          # Calculate Heidelberg-Welch diagnostic test?
                    geweke = TRUE,                          # Calculate Geweke diagnostic test?
                    diag_save = TRUE,                       # Save the diagnostics as a txt file?
                    diag_name = "diagnostics",              # If yes, specify the base file name (.txt will be appended later)
                    indiv_effect = indiv_effect,            # Is Individual a random effect in the model? (already specified)
                    plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                    plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                    plot_xy_save_png = FALSE)               # Save xy/trace plot as png?
output_JAGS(jags.1, mix, source, output_options)

