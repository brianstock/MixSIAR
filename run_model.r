#####################################################################
# Run the model
#####################################################################

# The run_model function is called when the 'RUN MODEL' button is pressed
# It writes the JAGS model file, calls JAGS, plots the JAGS output, and runs diagnostics

# July 8
# Removed f1 dimension from from fac2 variables (cross.fac2, tmp.p.fac2), deleted factor2_empty (lines 62-65)
# changed output_JAGS function call (removed factor1_levels, factor2_levels, factor1_names - they're global)

# August 19
# Added if statement at bottom to run 'plot_continuous_var' if at least one continuous effect exists

# August 28
# Added conc (Concentration Dependence data) to be passed to JAGS
# Changed write_JAGS_model function call to incorporate 'include_conc'

run_model <- function(){
  resid_err <- svalue(resid_err_box)  # If TRUE, add residual/SIAR error to the model.  If FALSE, don't.
  process_err <- svalue(proc_err_box)  # If TRUE, add process/MixSIR error to the model.  If FALSE, don't.

  if(!resid_err && !process_err){
    stop(paste("Neither residual nor process error selected.
  Choose one (or both) of the error terms in Model Error
  Options (top-right corner), and try again.",sep=""))
  }

  if(n.ce > 0 && include_indiv==F){
    stop(paste("In order to analyze a continuous
  covariate, Individual must be included in the model as
  a random effect. Restart MixSIAR and this time make
  sure to check the 'Include Individual as Random Effect'
  box when loading the mixture data."))
  }

  # function that writes the JAGS model file ("MixSIAR model.txt")
  # input: raw_source_data (T/F), sources_by_factor (T/F), n.re (0/1/2), resid_err (T/F), process_err (T/F), include_indiv (T/F), n.ce(non-neg integer), hierarch (T/F), include_conc (T/F)
  # output: JAGS model text file saved in WD ("MixSIAR_model.txt")
  write_JAGS_model(raw_source_data, sources_by_factor, n.re, resid_err, process_err, include_indiv, n.ce, hierarch, include_conc)

  ##############################################################################
  # Make variables to pass to JAGS
  ##############################################################################
  mcmc.chainLength <- as.integer(svalue(txt_chain_length))
  mcmc.burn <- as.integer(svalue(txt_burnin))
  mcmc.thin <- as.integer(svalue(txt_thin))
  mcmc.chains <- as.integer(svalue(txt_num_chains))
  calc.DIC = TRUE
  alpha = rep(1,n.sources)

  # make 'e', an Aitchision-orthonormal basis on the S^d simplex (equation 18, page 292, Egozcue 2003)
  # 'e' is used in the inverse-ILR transform (we pass 'e' to JAGS, where we do the ILR and inverse-ILR)
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }

  # other variables to give JAGS
  cross <- array(data=NA,dim=c(N,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p <- array(data=NA,dim=c(N,n.sources))              # dummy variable for inverse ILR calculation
  #jags.params <- c("p.global", "ilr.global")
  jags.params <- c("p.global")
  if(include_indiv){
    jags.params <- c(jags.params, "ind.sig", "p.ind")
  }
  
  # Random Effect data
  f.data <- character(0)
  if(n.re > 0){
    cross.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources))              # dummy variable for inverse ILR calculation
    jags.params <- c(jags.params, "p.fac1", "fac1.sig")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "cross.fac1", "tmp.p.fac1")
  }
  if(n.re > 1){
    cross.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources))              # dummy variable for inverse ILR calculation
    jags.params <- c(jags.params, "p.fac2", "fac2.sig")
    f.data <- c(f.data, "factor2_levels", "Factor.2", "cross.fac2", "tmp.p.fac2", "factor1_lookup")
  }

  # Source data
  if(raw_source_data){
    s.data <- c("SOURCE_array", "n.rep") # SOURCE_array contains the source data points, n.rep has the number of replicates by source and factor
  } else { # raw_source_data==FALSE
    s.data <- c("MU_array", "SIG2_array", "n_array")  # MU has the source sample means, SIG2 the source variances, n_array the sample sizes
  }
  if(sources_by_factor==T){       # include source factor level data, if we have it
    s.data <- c(s.data, "source_factor_levels")
  }
  if(include_conc){
    s.data <- c(s.data, "conc")   # include Concentration Dependence data, if we have it
  }

  # Continuous Effect data
  c.data <- rep(NA,n.ce)
  if(n.ce > 0){                               # If we have any continuous effects
    for(ce in 1:n.ce){                        # for each continuous effect
      c.data[ce] <- paste("Cont.",ce,sep="")  # add "Cont.ce" to c.data (e.g. c.data[1] <- Cont.1)
      jags.params <- c(jags.params,"ilr.global",paste("ilr.cont",n.ce,sep=""))   # add "ilr.cont(ce)" to jags.params (e.g. ilr.cont1)
    }
  }

  # Always pass JAGS the following data:
  all.data <- c("X_iso", "N", "n.sources", "n.iso", "alpha", "frac_mu", "frac_sig2", "e", "cross", "tmp.p")
  jags.data <- c(all.data, f.data, s.data, c.data)
  # if(resid_err){
  #   jags.params <- c(jags.params,"var.resid")
  # }
  # if(process_err){
  #   jags.params <- c(jags.params,"mix.var")
  # }

  #############################################################################
  # Call JAGS
  #############################################################################
  model.loc <- "MixSIAR_model.txt"
  jags.1 <<- jags(jags.data, parameters.to.save = jags.params, model.file = model.loc,
                n.chains = mcmc.chains, n.burnin = mcmc.burn, n.thin = mcmc.thin,
                n.iter = mcmc.chainLength, DIC = calc.DIC)

  output_options <- list(as.logical(svalue(summary_save)), svalue(summary_name), as.logical(svalue(sup_post)), as.logical(svalue(plot_post_save_pdf)),
                      svalue(plot_post_name), as.logical(svalue(sup_pairs)), as.logical(svalue(plot_pairs_save_pdf)), svalue(plot_pairs_name),
                      as.logical(svalue(sup_xy)), as.logical(svalue(plot_xy_save_pdf)), svalue(plot_xy_name), as.logical(svalue(gelman)),
                      as.logical(svalue(heidel)), as.logical(svalue(geweke)), as.logical(svalue(diag_save)), svalue(diag_name), include_indiv,
                      as.logical(svalue(plot_post_save_png)), as.logical(svalue(plot_pairs_save_png)), as.logical(svalue(plot_xy_save_png)))

  output_JAGS(jags.1, mcmc.chains, n.re, random_effects, n.sources, source_names, output_options)

  if(n.ce > 0){
    plot_continuous_var(output_options)
  }
} # end run_model function

