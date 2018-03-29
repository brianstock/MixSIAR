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
  jags.mod[[mod]] <- run_model(run="short", mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1, resid_err, process_err)
  
  # Process diagnostics, summary stats, and posterior plots
  # output_JAGS(jags.mod[[mod]], mix[[mod]], source[[mod]])
  output_JAGS(jags.mod[[mod]], mix[[mod]], source[[mod]], output_options=list(
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
                                                  plot_xy_save_png = FALSE))
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
save.image("alligator_ex_allmodels.RData")

###############################################################################
# Plot proportions vs. length from model 5 
R2jags::attach.jags(jags.mod[[5]])
n.sources <- source[[5]]$n.sources
source_names <- source[[5]]$source_names

get_high <- function(x){return(quantile(x,.975))}
get_low <- function(x){return(quantile(x,.025))}
n.plot = 200
chain.len = dim(p.global)[1]
Cont1.plot <- seq(from=round(min(mix[[5]]$CE[[1]]),1), to=round(max(mix[[5]]$CE[[1]]),1), length.out=n.plot)
ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
for(src in 1:n.sources-1){
  for(i in 1:n.plot){
   ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont1[,src]*Cont1.plot[i]
   ilr.low[i,src] <- get_low(ilr.plot[i,src,])
   ilr.median[i,src] <- median(ilr.plot[i,src,])
   ilr.high[i,src] <- get_high(ilr.plot[i,src,])
  }
}

# Transform regression lines from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
for(i in 1:(n.sources-1)){
   e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
   e[,i] <- e[,i]/sum(e[,i])
}
cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.median <- array(data=NA,dim=c(n.plot, n.sources))
cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.low <- array(data=NA,dim=c(n.plot, n.sources))
cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.high <- array(data=NA,dim=c(n.plot, n.sources))
eps.low <- rep(NA, n.plot)
eps.med <- rep(NA, n.plot)
eps.high <- rep(NA, n.plot)    
for(i in 1:n.plot){
  for(j in 1:(n.sources-1)){
    cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
    cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
  }
  for(src in 1:n.sources){
    tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    tmp.p.low[i,src] <- prod(cross.low[i,src,]);
    tmp.p.high[i,src] <- prod(cross.high[i,src,]);
  }
  for(src in 1:n.sources){
    p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
    p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
  }
  eps.med[i] <- calc_eps(p.median[i,])
  eps.low[i] <- calc_eps(p.low[i,])
  eps.high[i] <- calc_eps(p.high[i,])
}
colnames(p.median) <- source_names

Cont1.plot <- Cont1.plot*mix[[5]]$CE_scale + mix[[5]]$CE_center # transform Cont1.plot (x-axis) back to the original scale
df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
colnames(df) <- c("source","median","x","low","high")

# Plot of Diet vs. Cont effect
png("mod5_diet_length.png", height=7, width=7, units='in', res=500)
print(ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
        ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
        ggplot2::ylab("Diet proportion") +
        ggplot2::xlab("Total length (cm)") +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
          panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
          axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
          axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
          legend.justification=c(0,1), legend.title=ggplot2::element_blank()))
dev.off()

########################################################################
# Plot specialization index vs. length from model 5
calc_eps <- function(f){
  n.sources <- length(f)
  gam <- rep(1/n.sources,n.sources)
  phi <- rep(0,n.sources)
  phi[1] <- 1
  sqrt(sum((f-gam)^2))/sqrt(sum((phi-gam)^2))
} 

df.eps <- data.frame(Length=Cont1.plot, med=eps.med, low=eps.low, high=eps.high)
low.inc <- high.inc <- rep(NA,n.plot)
for(i in 2:n.plot){
  low.inc[i] <- ifelse(eps.low[i] > eps.low[i-1], TRUE, FALSE)
  high.inc[i] <- ifelse(eps.high[i] > eps.high[i-1], TRUE, FALSE)
}
low.change <- which(low.inc)[1]
high.change <- which(high.inc)[1]
df.eps$low[low.change:(high.change-1)] <- 0
df.eps$low[high.change:n.plot] <- eps.high[high.change:n.plot]
high.change2 <- which(eps.low > eps.high)[1]
df.eps$high[high.change2:n.plot] <- eps.low[high.change2:n.plot]

png("mod5_epsilon_length.png", height=7, width=7, units='in', res=500)
print(ggplot2::ggplot(data=df.eps, ggplot2::aes(x=Length,y=med)) +
        ggplot2::geom_line(ggplot2::aes(x=Length,y=med),size=1.5) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high), alpha=0.35) +
        ggplot2::ylab(expression(paste("Specialization index (",epsilon,")",sep=""))) +
        ggplot2::xlab("Total length (cm)") +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
          panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
          axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
          axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
          legend.justification=c(0,1), legend.title=ggplot2::element_blank()))
dev.off()
