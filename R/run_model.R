#' Run the JAGS model
#'
#' \code{run_model} calls JAGS to run the mixing model created by
#' \code{\link{write_JAGS_model}}. This happens when the "RUN MODEL" button is
#' clicked in the GUI.
#'
#' @param run list of MCMC parameters (chainLength, burn, thin, chains, calcDIC).
#'   Alternatively, a user can use a pre-defined parameter set by specifying a
#'   valid string:
#'   \itemize{
#'    \item \code{"test"}: chainLength=1000, burn=500, thin=1, chains=3
#'    \item \code{"very short"}: chainLength=10000, burn=5000, thin=5, chains=3
#'    \item \code{"short"}: chainLength=50000, burn=25000, thin=25, chains=3
#'    \item \code{"normal"}: chainLength=100000, burn=50000, thin=50, chains=3
#'    \item \code{"long"}: chainLength=300000, burn=200000, thin=100, chains=3
#'    \item \code{"very long"}: chainLength=1000000, burn=500000, thin=500, chains=3
#'    \item \code{"extreme"}: chainLength=3000000, burn=1500000, thin=500, chains=3
#'   }
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param discr output from \code{\link{load_discr_data}}
#' @param model_filename name of JAGS model file (usually should match \code{filename}
#'   input to \code{\link{write_JAGS_model}}).
#' @param alpha.prior Dirichlet prior on p.global (default = 1, uninformative)
#' @param resid_err include residual error in the model?
#' @param process_err include process error in the model?
#'
#' @return jags.1, a \code{rjags} model object
#'
run_model <- function(run, mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err){
#  resid_err <- mixsiar$resid_err
#  process_err <- mixsiar$process_err

  if(!process_err && !resid_err){
    stop(paste("Invalid error structure, must choose one of:
    1. Residual * Process (resid_err=TRUE, process_err=TRUE)
    2. Residual only (resid_err=TRUE, process_err=FALSE)
    3. Process only (resid_err=FALSE, process_err=TRUE)",sep=""))
  }
  if(resid_err && !process_err) err <- "resid"
  if(process_err && !resid_err) err <- "process"
  if(resid_err && process_err) err <- "mult"
  if(mix$N==1 && err!="process"){
    stop(paste("Invalid error structure. If N=1 mix datapoint,
    must choose Process only error model (MixSIR).
    Set resid_err=FALSE and process_err=TRUE.",sep=""))}
  if(mix$n.fe==1 && mix$N==mix$FAC[[1]]$levels && err!="process"){
    stop(paste("Invalid error structure. If fitting each individual
    mix datapoint separately, must choose Process only error model (MixSIR).
    Set resid_err=FALSE and process_err=TRUE.",sep=""))}

  # Error checks on prior
  if(length(alpha.prior)==1){
    if(alpha.prior==1) alpha.prior = rep(1,source$n.sources)
  }
  if(!is.numeric(alpha.prior)){
    stop(paste("*** Error: Your prior is not a numeric vector of length(n.sources).
        Try again or choose the uninformative prior option. For example,
        c(1,1,1,1) is a valid (uninformative) prior for 4 sources. ***",sep=""))}
  if(length(alpha.prior) != source$n.sources){
    stop(paste("*** Error: Length of your prior does not match the
        number of sources (",source$n.sources,"). Try again. ***",sep=""))}
  if(length(which(alpha.prior==0))!=0){
    stop(paste("*** Error: You cannot set any alpha = 0.
      Instead, set = 0.01.***",sep=""))}

  if(is.numeric(alpha.prior)==F) alpha.prior = 1 # Error checking for user inputted string/ NA
  if(length(alpha.prior)==1) alpha = rep(alpha.prior,source$n.sources) # All sources have same value
  if(length(alpha.prior) > 1 & length(alpha.prior) != source$n.sources) alpha = rep(1,source$n.sources) # Error checking for user inputted string/ NA
  if(length(alpha.prior) > 1 & length(alpha.prior) == source$n.sources) alpha = alpha.prior # All sources have different value inputted by user

  # Cannot set informative prior on fixed effects (no p.global)
  if(!identical(unique(alpha),1) & mix$n.fe>0){
  stop(paste("Cannot set an informative prior with a fixed effect,
  since there is no global/overall population. You can set an
  informative prior on p.global with a random effect.
  To set a prior on each level of a fixed effect you will have to
  modify 'write_JAGS_model.R'",sep=""))}

  # Set mcmc parameters
  if(is.list(run)){mcmc <- run} else { # if the user has entered custom mcmc parameters, use them
    if(run=="test") mcmc <- list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
    if(run=="very short") mcmc <- list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
    if(run=="short") mcmc <- list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
    if(run=="normal") mcmc <- list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
    if(run=="long") mcmc <- list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
    if(run=="very long") mcmc <- list(chainLength=1000000, burn=500000, thin=500, chains=3, calcDIC=TRUE)
    if(run=="extreme") mcmc <- list(chainLength=3000000, burn=1500000, thin=500, chains=3, calcDIC=TRUE)
  }

  n.sources <- source$n.sources
  N <- mix$N
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

  # Random/Fixed Effect data (original)
  # fere <- ifelse(mix$n.effects==2 & mix$n.re < 2,TRUE,FALSE)
  f.data <- character(0)
  if(mix$n.effects > 0 & !mix$fere){
    factor1_levels <- mix$FAC[[1]]$levels
    Factor.1 <- mix$FAC[[1]]$values
    cross.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources))              # dummy variable for inverse ILR calculation
    if(mix$FAC[[1]]$re) jags.params <- c(jags.params, "p.fac1", "fac1.sig") else jags.params <- c(jags.params, "p.fac1")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "cross.fac1", "tmp.p.fac1")
  }
  if(mix$n.effects > 1 & !mix$fere){
    factor2_levels <- mix$FAC[[2]]$levels
    Factor.2 <- mix$FAC[[2]]$values
    if(mix$fac_nested[1]) {factor2_lookup <- mix$FAC[[1]]$lookup; f.data <- c(f.data, "factor2_lookup");}
    if(mix$fac_nested[2]) {factor1_lookup <- mix$FAC[[2]]$lookup; f.data <- c(f.data, "factor1_lookup");}
    cross.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources))              # dummy variable for inverse ILR calculation
    if(mix$FAC[[2]]$re) jags.params <- c(jags.params, "p.fac2", "fac2.sig") else jags.params <- c(jags.params, "p.fac2")
    f.data <- c(f.data, "factor2_levels", "Factor.2", "cross.fac2", "tmp.p.fac2")
  }

  # 2FE or 1FE + 1RE, don't get p.fac2
  # instead, get ilr.both[f1,f2,src], using fac2_lookup (list, each element is vector of fac 2 levels in f1)
  # but do get p.fac1 if fac1=FE and fac2=RE
  if(mix$fere){
    # set up factor 1 as fixed effect (if 1FE + 1RE, fac 1 is fixed effect)
    factor1_levels <- mix$FAC[[1]]$levels
    Factor.1 <- mix$FAC[[1]]$values
    if(mix$n.re==1){ # have p.fac1 (fixed) for 1 FE + 1 RE
      cross.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources))              # dummy variable for inverse ILR calculation
      jags.params <- c(jags.params, "p.fac1")
      f.data <- c(f.data, "cross.fac1", "tmp.p.fac1")
    }

    # set up factor 2
    factor2_levels <- mix$FAC[[2]]$levels
    Factor.2 <- mix$FAC[[2]]$values
    # factor2_lookup <- list()
    # for(f1 in 1:factor1_levels){
    #   factor2_lookup[[f1]] <- unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values==f1)])
    # }
    # f.data <- c(f.data,"factor2_lookup")

    # cross.both <- array(data=NA,dim=c(factor1_levels,factor2_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    # tmp.p.both <- array(data=NA,dim=c(factor1_levels,factor2_levels,n.sources))              # dummy variable for inverse ILR calculation
    # if(mix$FAC[[2]]$re) jags.params <- c(jags.params, "p.both", "fac2.sig") else jags.params <- c(jags.params, "p.both")
    # f.data <- c(f.data, "factor1_levels", "Factor.1", "factor2_levels", "Factor.2", "cross.both", "tmp.p.both")
    if(mix$FAC[[2]]$re) jags.params <- c(jags.params, "fac2.sig")
    jags.params <- c(jags.params,"ilr.global","ilr.fac1","ilr.fac2")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "factor2_levels", "Factor.2")
  }

  # Source data
  if(source$data_type=="raw"){
    SOURCE_array <- source$SOURCE_array
    n.rep <- source$n.rep
    s.data <- c("SOURCE_array", "n.rep") # SOURCE_array contains the source data points, n.rep has the number of replicates by source and factor
  } else { # source$data_type="means"
    MU_array <- source$MU_array
    SIG2_array <- source$SIG2_array
    n_array <- source$n_array
    s.data <- c("MU_array", "SIG2_array", "n_array")  # MU has the source sample means, SIG2 the source variances, n_array the sample sizes
  }
  if(!is.na(source$by_factor)){       # include source factor level data, if we have it
    source_factor_levels <- source$S_factor_levels
    s.data <- c(s.data, "source_factor_levels")
  }
  if(source$conc_dep){
    conc <- source$conc
    s.data <- c(s.data, "conc")   # include Concentration Dependence data, if we have it
  }

  # Continuous Effect data
  c.data <- rep(NA,mix$n.ce)
  if(mix$n.ce > 0){                               # If we have any continuous effects
    for(ce in 1:mix$n.ce){                        # for each continuous effect
      name <- paste("Cont.",ce,sep="")
      assign(name,as.vector(mix$CE[[ce]]))
      c.data[ce] <- paste("Cont.",ce,sep="")  # add "Cont.ce" to c.data (e.g. c.data[1] <- Cont.1)
      jags.params <- c(jags.params,"ilr.global",paste("ilr.cont",ce,sep=""),"p.ind")   # add "ilr.cont(ce)" to jags.params (e.g. ilr.cont1)
    }
  }

  X_iso <- mix$data_iso
  n.iso <- mix$n.iso
  frac_mu <- discr$mu
  frac_sig2 <- discr$sig2
  # Always pass JAGS the following data:
  # all.data <- c("X_iso", "N", "n.sources", "n.iso", "alpha", "frac_mu", "frac_sig2", "e", "cross", "tmp.p")
  all.data <- c("X_iso", "N", "n.sources", "n.iso", "alpha", "frac_mu", "e", "cross", "tmp.p")
  jags.data <- c(all.data, f.data, s.data, c.data)
  # if(resid_err){
  #   jags.params <- c(jags.params,"var.resid")
  # }
  # if(process_err){
  #   jags.params <- c(jags.params,"mix.var")
  # }

  # Error structure objects
  I <- diag(n.iso)
  if(err=="resid" && mix$n.iso>1) jags.data <- c(jags.data,"I")
  if(err!="resid") jags.data <- c(jags.data,"frac_sig2")
  if(err=="mult") jags.params <- c(jags.params,"resid.prop")

  # Set initial values for p.global different for each chain
  jags.inits <- function(){list(p.global=as.vector(compositions::rDirichlet.rcomp(1,alpha)))}

  #############################################################################
  # Call JAGS
  #############################################################################
  jags.1 <- R2jags::jags(jags.data,
                                  inits=jags.inits,
                                  parameters.to.save = jags.params,
                                  model.file = model_filename,
                                  n.chains = mcmc$chains,
                                  n.burnin = mcmc$burn,
                                  n.thin = mcmc$thin,
                                  n.iter = mcmc$chainLength,
                                  DIC = mcmc$calcDIC)

  
  return(jags.1)
} # end run_model function

