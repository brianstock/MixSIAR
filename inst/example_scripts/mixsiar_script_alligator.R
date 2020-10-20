# Alligator example based on Nifong et al. (2015)
#   Paper (open-access): https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12306
#   Data: https://datadryad.org//resource/doi:10.5061/dryad.j7049

# This case study highlights the main advantage of MixSIAR over previous mixing model software--
# the ability to include fixed and random effects as covariates explaining variability in mixture proportions 
# and calculate relative support for multiple models via information criteria. 

# Nifong et al. (2015) interested in freshwater vs. marine resource use by alligators
# Specific scientific questions:
#   Q1. What is p.marine vs. p.freshwater?
#   Q2. How does p.marine vary with the covariates Length, Sex, and Individual?
#   Q3. How variable are individuals’ diets relative to group-level variability?

# Nifong et al. (2015) approach:
#   - group consumers into 8 subpopulations (all combos of Sex x Size Class)
#     - 2 sexes {male, female}
#     - 4 size classes {small juvenile, large juvenile, subadult, adult}
#   - run 8 separate mixing models for each using SIAR (Parnell et al. 2010). 
#   - to calculate p.marine estimates for the overall population, had to also run a mixing model with all consumers

# MixSIAR approach:
#   - fit several models with fixed and random effects as covariates
#   - evaluate relative support for each model using information criteria (LOO, "compare_models" function)

########################################################################################
# # First, set your working directory to folder containing this script.
# setwd("~/dataS1/")
# # Then, source the mixsiar_script file...
# source("mixsiar_script_alligator.R")

# ------------------------------------------------------------------------------
library(MixSIAR)
mix.filename <- system.file("extdata", "alligator_consumer.csv", package = "MixSIAR")
source.filename <- system.file("extdata", "alligator_sources_simplemean.csv", package = "MixSIAR")
discr.filename <- system.file("extdata", "alligator_TEF.csv", package = "MixSIAR")

# ------------------------------------------------------------------------------
# make a list that holds mix data for all 8 models
#   (models use the same source and TDF data, only differ in covariates on mixture)

# Model 1: NULL
# Model 2: habitat (3 habitats: fresh, intermediate, marine)
# Model 3: sex (male, female)
# Model 4: sclass (4 size classes: small juv, large juv, sub-adult, adult)
# Model 5: length (continuous effect)
# Model 6: sex + sclass (as in Nifong 2015, both as fixed effects)
# Model 7: sex + length (intercept by sex, same slope)
# Model 8: sex_sclass (create new factor = sex * sclass)
#          this is closest to original analysis by Nifong et al. (2015)

n.mod <- 8
mix <- vector("list", n.mod) 
mix[[1]] <- load_mix_data(filename=mix.filename,
                     iso_names=c("d13C","d15N"),
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)
mix[[2]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="habitat",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[3]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[4]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sclass",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[5]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors=NULL,
                          fac_random=NULL,
                          fac_nested=NULL,
                          cont_effects="Length")
mix[[6]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors=c("sex","sclass"),
                          fac_random=c(FALSE,FALSE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects=NULL)
mix[[7]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects="Length")
mix[[8]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex_sclass",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)

# Run the models
source <- vector("list", n.mod)
discr <- vector("list", n.mod)
jags.mod <- vector("list", n.mod)
for(mod in 1:n.mod){
  # create sub-directory and move into it
  mainDir <- getwd()
  subDir <- paste0("model_", mod)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  # load source data
  source[[mod]] <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="means",
                             mix[[mod]])
  
  # load TEF data
  discr[[mod]] <- load_discr_data(filename=discr.filename, mix[[mod]])
  
  # isospace plot
  plot_data(filename=paste0("isospace_plot_", mod),
            plot_save_pdf=TRUE,
            plot_save_png=FALSE,
            mix[[mod]], source[[mod]], discr[[mod]])
  
  # Define model structure and write JAGS model file
  model_filename <- paste0("MixSIAR_model_", mod, ".txt")
  resid_err <- TRUE
  process_err <- TRUE
  write_JAGS_model(model_filename, resid_err, process_err, mix[[mod]], source[[mod]])
  
  # Run the JAGS model
  # "short" MCMC length is plenty long for all models to converge
  jags.mod[[mod]] <- run_model(run="short", mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1)
  
  # Process diagnostics, summary stats, and posterior plots
  output_options=list(
    summary_save = TRUE,                 # Save the summary statistics as a txt file?
    summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
    sup_post = TRUE,                       # Suppress posterior density plot output in R?
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
    plot_xy_save_png = FALSE,
    diag_save_ggmcmc = FALSE,
    return_obj = FALSE)
  
  output_JAGS(jags.mod[[mod]], mix[[mod]], source[[mod]], output_options)
  graphics.off()
  
  # Move back up to root directory
  setwd(mainDir)
}

# Use 'compare_models' to get table with LOOic weights
names(jags.mod) <- c("Null","Habitat","Sex","Size class","Length","Sex + Size class","Sex + Length","Sex : Length")
comparison.table <- compare_models(jags.mod)

# get multiplicative error term estimates, median(xi.C) and median(xi.N)
xi.C <- xi.N <- rep(NA, n.mod)
for(i in 1:n.mod){
  xi.C[i] <- round(median(jags.mod[[i]]$BUGSoutput$sims.list$resid.prop[,1]),1)
  xi.N[i] <- round(median(jags.mod[[i]]$BUGSoutput$sims.list$resid.prop[,2]),1)
}

# add xi.C and xi.N to comparison table (in correct order)
y <- as.numeric(rownames(comparison.table))
comparison.table <- cbind(comparison.table, xi.C[y], xi.N[y])
comparison.table

# Save all results
save.image("alligator_short.RData")
# save only model 5
for(i in c(1:4,6:8)){ jags.mod[[i]] = NULL; mix[[i]] = NULL; source[[i]] = NULL}
save(list=c("jags.mod","mix","source","discr"), file="alligator_mod5.RData")

# modify proportions vs. length plot from best model
library(ggplot2)
output_options$return_obj = TRUE
g <- plot_continuous_var(jags.1 = jags.mod[[5]], mix = mix[[5]], source = source[[5]], output_options)

# legend covers line - move it outside plot area
g[[1]] +
  theme(legend.position = "right")
