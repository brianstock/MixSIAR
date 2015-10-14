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
require(lattice)
require(compositions) # need rdirichlet function to establish initial values for each chain 

# Next we clean up the workspace
rm(list=ls())  # deletes everything previously in the workspace
runif(1)       # generates one random number (else JAGS can complain)

# Load all MixSIAR functions into the workspace
source("load_mix_data_script.r")
source("load_source_data.r")
source("load_discr_data.r")
source("plot_data.r")
source("write_JAGS_model.r")
source("run_model.r")
source("output_JAGS.r")
source("plot_continuous_var.r")
source("plot_prior.r")
source("calc_area.r")

#####################################################################################
# Load mixture data, i.e. your:
#    Consumer isotope values (trophic ecology / diet)
#    Mixed sediment/water tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with mix/consumer data
# 'iso_names' - column headings of the tracers/isotopes you'd like to use
# 'random_effects' - column headings of any random effects
# 'cont_effects' - column headings of any continuous effects

# NOTE: 'iso_names', 'random_effects', and 'cont_effects' can be a subset of your data columns
#   i.e. have 3 isotopes in file but only want MixSIAR to use 2,
#   have data by Region and Pack but only want MixSIAR to use Region

# Wolves example (hierarchical/nested random effects)
mix <- load_mix_data_script(filename="wolves_consumer.csv", iso_names=c("d13C","d15N"), factors=c("Region","Pack"), fac_random=c(TRUE,TRUE), fac_nested=c(FALSE,TRUE), cont_effects=NULL)

# Lake example (continuous effect)
# mix <- load_mix_data_script(filename="lake_consumer.csv", iso_names=c("d13C","d15N"), factors=NULL, fac_random=NULL, fac_nested=NULL, cont_effects="Secchi.Mixed")

# Geese example (concentration dependence)
# mix <- load_mix_data_script(filename="geese_consumer.csv", iso_names=c("d13C","d15N"), factors="Group", fac_random=FALSE, fac_nested=FALSE, cont_effects=NULL)

# Palmyra example (fixed effect)
# mix <- load_mix_data_script(filename="palmyra_consumer.csv", iso_names=c("d13C","d15N"), factors="Taxa", fac_random=FALSE, fac_nested=FALSE, cont_effects=NULL)

# Storm-petrel example (fixed effect)
# mix <- load_mix_data_script(filename="7_mix.csv", iso_names=c("d13C","d15N"), factors="Region", fac_random=FALSE, fac_nested=FALSE, cont_effects=NULL)

# 1-iso example
# mix <- load_mix_data_script(filename="13_mix.csv", iso_names="d13C", factors=NULL, fac_random=NULL, fac_nested=NULL, cont_effects=NULL)

# killer whale - salmon example
# mix <- load_mix_data_script(filename="killerwhale_consumer.csv", iso_names=c("d13C","d15N"), factors=NULL, fac_random=NULL, fac_nested=NULL, cont_effects=NULL)

# Isopod example (8 fatty acids)
# mix <- load_mix_data_script(filename="isopod_consumer.csv", iso_names=c("c16.4w3","c18.2w6","c18.3w3","c18.4w3","c20.4w6","c20.5w3","c22.5w3","c22.6w3"), factors="Site", fac_random=FALSE, fac_nested=FALSE, cont_effects=NULL)

#####################################################################################
# Load source data, i.e. your:
#    Source isotope values (trophic ecology / diet)
#    Sediment/water source tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with source data
# 'source_factors' - column headings of random/fixed effects you have source data by
# 'conc_dep' - TRUE or FALSE, do you have concentration dependence data in the file?
# 'data_type' - "means" or "raw", is your source data in the means+SD format, or do you have raw data

# Wolves example
source <- load_source_data(filename="wolves_sources.csv", source_factors="Region", conc_dep=FALSE, data_type="means", mix)    

# Lake example
# source <- load_source_data(filename="lake_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    

# Geese example
# source <- load_source_data(filename="geese_sources.csv", source_factors=NULL, conc_dep=TRUE, data_type="means", mix)    

# Palmyra example
# source <- load_source_data(filename="palmyra_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    

# Storm-petrel example
# source <- load_source_data(filename="7_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)  

# 1-iso example
# source <- load_source_data(filename="13_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    

# killer whale - salmon example
# source <- load_source_data(filename="killerwhale_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="means", mix)    

# Isopod example
# source <- load_source_data(filename="isopod_sources.csv", source_factors=NULL, conc_dep=FALSE, data_type="means", mix)    

#####################################################################################
# Load discrimination data, i.e. your:
#    Trophic Enrichment Factor (TEF) / fractionation values (trophic ecology / diet)
#    xxxxxxxx (sediment/hydrology fingerprinting)

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

# killer whale - salmon example
# discr <- load_discr_data(filename="killerwhale_discrimination.csv", mix) 

# Isopod example
# discr <- load_discr_data(filename="isopod_discrimination.csv", mix) 

#####################################################################################
# Make isospace plot
# Are the data loaded correctly?
# Is your mixture data in the source polygon?
# Are one or more of your sources confounded/hidden?

# 'filename' - name you'd like MixSIAR to save the isospace plot as (extension will be added automatically)
# 'plot_save_pdf' - TRUE or FALSE, should MixSIAR save the plot as a .pdf?
# 'plot_save_png' - TRUE or FALSE, should MixSIAR save the plot as a .png?

plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

# If 2 isotopes/tracers, calculate the normalized surface area of the convex hull polygon(s)
#   Note that the discrimination SD is added to the source SD (see calc_area.r for details)
#   If source data are by factor (as in wolf ex), computes area for each polygon (one for each of 3 regions in wolf ex)
if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

#####################################################################################
# Write JAGS model file
# Model will be saved as 'model_filename' ("MixSIAR_model.txt" is default, but may want to change if in a loop)

# 'model_filename' - don't need to change, unless you are creating many different model files
# 'resid_err' - TRUE or FALSE, do you want to include residual error in the model (FALSE = MixSIR, TRUE = SIAR)?

# Wolves, Killer whale, Isopod examples
model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err <- FALSE
write_JAGS_model(model_filename, resid_err, mix,source)

# Geese, Palmyra, Lake, Storm-petrel, 1-iso examples
# model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
# write_JAGS_model(model_filename, resid_err=TRUE, mix,source)

#####################################################################################
# Define your prior, and then plot using "plot_prior"
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
# plot_prior(alpha.prior=1,source)

# Killer whale example with INFORMATIVE prior (construct alpha from fecal data)
#   Let's say we have 14 fecal diet samples that we use to construct alphas...useful in separating some of the sources.
# kw.alpha <- c(10,1,0,0,3)   # Our 14 fecal samples were 10, 1, 0, 0, 3
# kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha) # Generate alpha hyperparameters scaling sum(alpha)=n.sources
# kw.alpha[which(kw.alpha==0)] <- 0.001 # the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .001)
# plot_prior(alpha.prior=kw.alpha,source=source,plot_save_pdf=TRUE, plot_save_png=FALSE,filename="prior_plot")

#####################################################################################
# Run model
# JAGS output will be saved as 'jags.1'

# MCMC run options:
# run <- "test"          # list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short"    # list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short"         # list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal"        # list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long"          # list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long"     # list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme"       # list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

# Can also set custom MCMC parameters
# run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

# Good idea to use 'test' first to check if 1) the data are loaded correctly and 2) the model is specified correctly
jags.1 <- run_model(run="test",mix,source,discr,model_filename,alpha.prior = 1,resid_err)

# Wolves, Palmyra, Geese examples
# jags.1 <- run_model(run="short",mix,source,discr,model_filename,alpha.prior = 1,resid_err)

# Lake example
# jags.1 <- run_model(run="normal",mix,source,discr,model_filename,alpha.prior = 1,resid_err)

# Killer whale example with UNINFORMATIVE / GENERALIST prior (alpha = 1)
# jags.1 <- run_model(run="short",mix,source,discr,model_filename,alpha.prior = 1,resid_err) # alpha = 1 by default

# Killer whale example with INFORMATIVE prior (constructed kw.alpha from fecal data above)
# jags.1 <- run_model(run="test",mix,source,discr,model_filename,alpha.prior = kw.alpha,resid_err) # informative prior

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
                    indiv_effect = FALSE,                   # Is Individual a random effect in the model? (already specified)
                    plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                    plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                    plot_xy_save_png = FALSE)               # Save xy/trace plot as png?
output_JAGS(jags.1, mix, source, output_options)

