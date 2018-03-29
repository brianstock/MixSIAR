# Brian Stock
# March 13, 2018
# Script file to run mantis shrimp example

# Load MixSIAR package
library(MixSIAR)

################################################################################
# Load mixture data, i.e. your:
#    Consumer isotope values (trophic ecology / diet)
#    Mixed sediment/water tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with mix/consumer data
# 'iso_names' - column headings of the tracers/isotopes you'd like to use
# 'random_effects' - column headings of any random effects
# 'cont_effects' - column headings of any continuous effects

# 'iso_names', 'random_effects', and 'cont_effects' can be a subset of your columns
#   i.e. have 3 isotopes in file but only want MixSIAR to use 2,
#   or have data by Region and Pack but only want MixSIAR to use Region

# To run on your data, replace the system.file call with the path to your file
mix.filename <- system.file("extdata", "mantis_consumer.csv", package = "MixSIAR")

# Load mixture data
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C","d15N"), 
                     factors="Habitat", 
                     fac_random=FALSE, 
                     fac_nested=FALSE, 
                     cont_effects=NULL)
# mix <- load_mix_data(filename=mix.filename, 
#                      iso_names=c("d13C","d15N"), 
#                      factors=c("Habitat","Season"),
#                      fac_random=c(FALSE,FALSE), 
#                      fac_nested=c(FALSE,FALSE),
#                      cont_effects=NULL)

################################################################################
# Load source data, i.e. your:
#    Source isotope values (trophic ecology / diet)
#    Sediment/water source tracer values (sediment/hydrology fingerprinting)

# 'filename': name of the CSV file with source data
# 'source_factors': column headings of random/fixed effects you have source data by
# 'conc_dep': TRUE or FALSE, do you have concentration dependence data in the file?
# 'data_type': "means" or "raw", is your source data as means+SD, or do you have raw data

# To run on your data, replace the system.file call with the path to your file
source.filename <- system.file("extdata", "mantis_source.csv", package = "MixSIAR")

# Load source data
source <- load_source_data(filename=source.filename, 
                           source_factors="Habitat", 
                           conc_dep=TRUE, 
                           data_type="means", mix)

################################################################################
# Load discrimination data, i.e. your:
#    Trophic Enrichment Factor (TEF) / fractionation values (trophic ecology/diet)
#    xxxxxxxx (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with discrimination data

# To run on your data, replace the system.file call with the path to your file
discr.filename <- system.file("extdata", "mantis_discrimination.csv", package = "MixSIAR")

# Load discrimination data
discr <- load_discr_data(filename=discr.filename, mix)

#####################################################################################
# Make isospace plot
# Are the data loaded correctly?
# Is your mixture data in the source polygon?
# Are one or more of your sources confounded/hidden?

# 'filename' - name you'd like MixSIAR to save the isospace plot as 
#              (extension will be added automatically)
# 'plot_save_pdf' - TRUE or FALSE, should MixSIAR save the plot as a .pdf?
# 'plot_save_png' - TRUE or FALSE, should MixSIAR save the plot as a .png?

plot_data(filename="isospace_plot", 
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix,source,discr)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
#   *Note 2: If source data are by factor (as in wolf ex), computes area for each polygon
#             (one for each of 3 regions in wolf ex)
# if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

################################################################################
# Define your prior, and then plot using "plot_prior"
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
alpha.unif <- rep(1, source$n.sources)
plot_prior(alpha.prior=alpha.unif, 
           source=source,
           filename="prior_uninf")

# But wait... the GENERALIST prior does not reflect our a priori ecological knowledge!
# Current thought is that 'smasher' mantis shrimp SPECIALIZE on hard-shelled prey,
# because they have specialized hammer-like clubs allowing them to break open shells.
# Let's construct a "specialist" informative prior reflecting the expectation
# that N. bredini consumes primarily hard-shelled prey. This will be a more
# conservative test of the hypothesis that N. bredini specializes on hard-shelled prey.

# 'specialist' informative prior (hard-shelled prey get 4x weight of soft-bodied prey,
# based on dietary observations of N. bredini (Caldwell et al. 1989).
alpha.spec <- c(1,1,4,4,1,4)
alpha.spec <- alpha.spec*length(alpha.spec)/sum(alpha.spec)
plot_prior(alpha.prior=alpha.spec, 
           source=source,
           filename="prior_specialist")

# We also have data on prey abundance in the two habitats (coral, seagrass).
# Construct priors using prey abundance so that the alpha_k sum to n_sources
alpha.grass <- c(0.35,1.61,0.43,(51.65+0.26),5.18,40.5)*6/100
plot_prior(alpha.prior=alpha.grass, 
           source=source,
           filename="prior_seagrass")

alpha.coral <- c((14.31+24.74),0.01,15.48,(13.81+4.71),8.44,18.51)*6/100
plot_prior(alpha.prior=alpha.coral, 
           source=source,
           filename="prior_coral")

################################################################################
# Write JAGS model file (define model structure)
# Model will be saved as 'model_filename' ("MixSIAR_model.txt" is default,
#    but may want to change if in a loop)

# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)

# 'model_filename': don't need to change, unless you create many different models
# 'resid_err': include residual error in the model?
# 'process_err': include process error in the model?

#  *Note: If you have only 1 mix datapoint, you have no information about the 
#         mixture/consumer variability. In this case, we ues the original MixSIR
#         error model (which does not fit a residual error term).
#         This is the same behavior as 'siarsolo' in SIAR.

model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

################################################################################
# Run model
# JAGS output will be saved as 'jags.1'

# MCMC run options:
# run <- "test"       # chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE
# run <- "very short" # chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE
# run <- "short"      # chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE
# run <- "normal"     # chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE
# run <- "long"       # chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE
# run <- "very long"  # chainLength=1000000, burn=500000, thin=500, chains=3, calcDIC=TRUE
# run <- "extreme"    # chainLength=3000000, burn=1500000, thin=500, chains=3, calcDIC=TRUE

# Can also set custom MCMC parameters
# run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

# Good idea to use 'test' first to check if
#   1) the data are loaded correctly, and 
#   2) the model is specified correctly
# jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
#                     alpha.prior = alpha.spec, resid_err, process_err)

# After a test run works, increase the MCMC run to a value that may converge
jags.spec <- run_model(run="normal", mix, source, discr, model_filename, 
                    alpha.prior = alpha.spec, resid_err, process_err)

################################################################################
# Process JAGS output

# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

# Create diagnostics, summary statistics, and posterior plots
output_JAGS(jags.spec, mix, source, output_options)

# Show original 6 source names
source$source_names 

# Combine hard-shelled and soft-bodied prey sources
combined <- combine_sources(jags.spec, mix, source, alpha.spec, 
  groups=list(hard=c("clam","crab","snail"), soft=c("alphworm","brittlestar","fish")))
original <- combine_sources(jags.spec, mix, source, alpha.spec, 
  groups=list(alphworm="alphworm",brittlestar="brittlestar",clam="clam",crab="crab",fish="fish",snail="snail"))

summary_stat(original)
summary_stat(combined)
# summary_stat(combined, savetxt=FALSE)
# summary_stat(combined, meanSD=FALSE)
# summary_stat(combined, quantiles=c(.05,.5,.95))
# summary_stat(combined, toprint="fac1", savetxt=FALSE)
# summary_stat(original, toprint="fac1", savetxt=FALSE)
# summary_stat(combined, toprint="p")
# summary_stat(combined, toprint="global")

dev.new()
plot_intervals(combined,toplot="fac1")
dev.new()
plot_intervals(original,toplot="fac1")
dev.new()
plot_intervals(combined,toplot="fac1",levels=1)
dev.new()
plot_intervals(combined,toplot="fac1",levels=2)
