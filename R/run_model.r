# Brian Stock
# January 30, 2014

# Function: run_model
#   sets up JAGS objects and calls JAGS ("running the model")
# Usage: jags.1 <- run_model(run,mix,source,discr)
# Input: run              list of mcmc parameters (chainLength, burn, thin, chains, calcDIC)
#                         alternatively, a user can enter "test", "very short", "short", "normal", "long", or "very long"
#                         these have pre-defined mcmc parameter sets (see lines 14-20)
#        mix              output from 'load_mix_data'
#        source           output from 'load_source_data'
#        discr            output from 'load_discrimination_data'
#        model_filename   name of the JAGS model file (created by 'write_JAGS_model.r')
# Output: jags.1, a rjags model object

run_model <- function(run,mix,source,discr,model_filename,alpha.prior = 1,resid_err=FALSE){
  # Set mcmc parameters
  if(run=="test") mcmc <- list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
  if(run=="very short") mcmc <- list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
  if(run=="short") mcmc <- list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
  if(run=="normal") mcmc <- list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
  if(run=="long") mcmc <- list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
  if(run=="very long") mcmc <- list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
  if(run=="extreme") mcmc <- list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)
  if(!exists("mcmc")) mcmc <- run   # if the user has entered custom mcmc parameters, use them
  
  n.sources <- source$n.sources
  N <- mix$N
  # make 'e', an Aitchision-orthonormal basis on the S^d simplex (equation 18, page 292, Egozcue 2003)
  # 'e' is used in the inverse-ILR transform (we pass 'e' to JAGS, where we do the ILR and inverse-ILR)
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(source$n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }

  # other variables to give JAGS
  cross <- array(data=NA,dim=c(N,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p <- array(data=NA,dim=c(N,n.sources))              # dummy variable for inverse ILR calculation
  #jags.params <- c("p.global", "ilr.global")
  jags.params <- c("p.global")
  
  # Random/Fixed Effect data
  f.data <- character(0)
  if(mix$n.effects > 0){
    factor1_levels <- mix$FAC[[1]]$levels
    Factor.1 <- mix$FAC[[1]]$values
    cross.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac1 <- array(data=NA,dim=c(factor1_levels,n.sources))              # dummy variable for inverse ILR calculation
    if(mix$FAC[[1]]$re) jags.params <- c(jags.params, "p.fac1", "fac1.sig") else jags.params <- c(jags.params, "p.fac1")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "cross.fac1", "tmp.p.fac1")
  }
  if(mix$n.effects > 1){
    factor2_levels <- mix$FAC[[2]]$levels
    Factor.2 <- mix$FAC[[2]]$values
    if(mix$fac_nested[1]) {factor2_lookup <- mix$FAC[[1]]$lookup; f.data <- c(f.data, "factor2_lookup");}
    if(mix$fac_nested[2]) {factor1_lookup <- mix$FAC[[2]]$lookup; f.data <- c(f.data, "factor1_lookup");}
    cross.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
    tmp.p.fac2 <- array(data=NA,dim=c(factor2_levels,n.sources))              # dummy variable for inverse ILR calculation
    if(mix$FAC[[2]]$re) jags.params <- c(jags.params, "p.fac2", "fac2.sig") else jags.params <- c(jags.params, "p.fac2")
    f.data <- c(f.data, "factor2_levels", "Factor.2", "cross.fac2", "tmp.p.fac2")
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
  if(source$by_factor==TRUE){       # include source factor level data, if we have it
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
  
  if(is.numeric(alpha.prior)==F) alpha.prior = 1 # Error checking for user inputted string/ NA
  if(length(alpha.prior)==1) alpha = rep(alpha.prior,source$n.sources) # All sources have same value
  if(length(alpha.prior) > 1 & length(alpha.prior) != source$n.sources) alpha = rep(1,source$n.sources) # Error checking for user inputted string/ NA
  if(length(alpha.prior) > 1 & length(alpha.prior) == source$n.sources) alpha = alpha.prior # All sources have different value inputted by user
  
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
  if(resid_err && mix$n.iso>1) jags.data <- c(jags.data,"I")
  if(!resid_err){
    jags.data <- c(jags.data,"frac_sig2")
    jags.params <- c(jags.params,"resid.prop")
  }

  # Set initial values for p.global different for each chain
  jags.inits <- function(){list(p.global=rdirichlet(alpha))}
    
  # jags.inits <- function(){list(p.global=rdirichlet(alpha),fac1.sig=runif(1,min=0,max=20),fac2.sig=runif(1,min=0,max=20))}
  # ilr.cont1=runif(1,min=-2,max=2)

  #############################################################################
  # Call JAGS
  #############################################################################
  jags.1 <- jags(jags.data, inits=jags.inits, parameters.to.save = jags.params, model.file = model_filename,
                n.chains = mcmc$chains, n.burnin = mcmc$burn, n.thin = mcmc$thin,
                n.iter = mcmc$chainLength, DIC = mcmc$calcDIC)
  return(jags.1)
} # end run_model function

