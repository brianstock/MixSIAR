# Brian Stock
# January 30, 2014

# July 17, 2014
# Changed order of input
# Added 'filename' so the user can use the function in a loop and make different model files

# Function: write_JAGS_model
#     creates "MixSIAR_model.txt", the JAGS model file
# Input: 
#   filename        formerly "MixSIAR_model.txt", now the user can change (i.e. batch processing in a for loop)
#   indiv_effect    T/F - include Individual as a random effect in the model?   include_indiv
#   nested          T/F - if multiple random effects are included, is one nested in the other (e.g. Pack within Region)?    hierarch
#   resid_err       T/F - include residual error in the model?            resid_err
#   mix             output from 'load_mix_data'
#     mix$n.re - number of random effects                                 n.re
#     mix$n.ce - number of continuous effects                             n.ce
#     mix$n.fe - number of fixed effects
#     mix$n.effects - n.re + n.fe
#   source          output from 'loac_source_data'
#     source$data_type - "raw" or "means"                                 raw_source_data
#     source$by_factor - T/F, are the source data by a factor?            sources_by_factor
#     source$conc_dep - T/F, do we have concentration dependence data?    include_conc

# Output: JAGS model text file saved to working directory as 'filename' (default is "MixSIAR_model.txt")

write_JAGS_model <- function(filename, resid_err, mix,source){
process_err <- TRUE
resid_err_mult <- FALSE
  if(!resid_err && !process_err){
    stop(paste("Neither residual nor process error selected.
  Choose one (or both) of the error terms in Model Error
  Options (top-right corner), and try again.",sep=""))
  }

cat(paste("# source$data_type: ",source$data_type,sep=""), file=filename)
cat("
", file=filename, append=T)
cat(paste("# source$by_factor: ",source$by_factor,sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# random effects: ",mix$n.re,sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# fixed effects: ",mix$n.fe,sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# nested factors: ",paste(mix$fac_nested,collapse=" "),sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# factors: ",paste(mix$factors,collapse=" "),sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# continuous effects: ",mix$n.ce,sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# resid_err: ",resid_err,sep=""), file=filename, append=T)
# cat("
# ", file=filename, append=T)
# cat(paste("# resid_err_mult: ",resid_err_mult,sep=""), file=filename, append=T)
# cat("
# ", file=filename, append=T)
# cat(paste("# process_err: ",process_err,sep=""), file=filename, append=T)
# cat("
# ", file=filename, append=T)
# cat(paste("# indiv_effect: ",indiv_effect,sep=""), file=filename, append=T)
cat("
", file=filename, append=T)
cat(paste("# source$conc_dep: ",source$conc_dep,sep=""), file=filename, append=T)

if(source$data_type=="raw" && source$by_factor==FALSE){
cat("

var rho[n.sources,n.iso,n.iso], src_cov[n.sources,n.iso,n.iso], src_var[n.sources,n.iso,n.iso], src_Sigma[n.sources,n.iso,n.iso], Sigma.ind[N,n.iso,n.iso], mix.cov[N,n.iso,n.iso];", file=filename, append=T)  
}
if(source$data_type=="raw" && source$by_factor==TRUE){
cat("

var rho[n.sources,source_factor_levels,n.iso,n.iso], src_cov[n.sources,source_factor_levels,n.iso,n.iso], src_var[n.sources,source_factor_levels,n.iso,n.iso], src_Sigma[n.sources,source_factor_levels,n.iso,n.iso], Sigma.ind[N,n.iso,n.iso], mix.cov[N,n.iso,n.iso];", file=filename, append=T)  
}
cat("

model{", file=filename, append=T)

################################################################################
# Setup src_mu and src_tau arrays
# Fit source means (if source$data_type=="raw")
################################################################################

# if(source$data_type=="raw" && source$by_factor==TRUE){ # fit the source means and precisions (by factor 1)
# cat("
#   # uninformed priors on source means and precisions
#   for(src in 1:n.sources){
#     for(f1 in 1:source_factor_levels){
#       for(iso in 1:n.iso){
#         src_mu[src,iso,f1] ~ dnorm(0,.001);
#       }

#     }
#   }
#   # each source data point is distributed normally according to the source means and precisions
#   for(src in 1:n.sources){
#     for(f1 in 1:source_factor_levels){
#       for(r in 1:n.rep[src,f1]){
#         SOURCE_array[src,,f1,r] ~ dmnorm(src_mu[src,,f1],src_Sigma.inv[src,f1,,]);
#       } 
#     }
#   }
# ", file=filename, append=T)
# }
if(source$data_type=="raw" && source$by_factor==TRUE){  # fit the source means and precisions (not by factor 1)
cat("
  # fit source data (big for loop over sources)
  for(src in 1:n.sources){
    for(f1 in 1:source_factor_levels){ # fit source data from each source level separately
      # uninformative priors on source means (src_mu vector)
      for(iso in 1:n.iso){
        src_mu[src,iso,f1] ~ dnorm(0,.001);
      }

      # uninformative priors on source variances (src_tau matrix)
      for(i in 2:n.iso){
        for(j in 1:(i-1)){
          src_tau[src,f1,i,j] <- 0;
          src_tau[src,f1,j,i] <- 0;
        }
      }
      for(i in 1:n.iso){
        src_tau[src,f1,i,i] ~ dgamma(.001,.001);
      }

      # uninformative priors on source correlations (rho matrix)
      for(i in 2:n.iso){
        for(j in 1:(i-1)){
          rho[src,f1,i,j] ~ dunif(-1,1);
          rho[src,f1,j,i] <- rho[src,f1,i,j];
        }
      }
      for(i in 1:n.iso){
        rho[src,f1,i,i] <- 1;
      }

      # Construct source precision matrix (src_Sigma)
      src_var[src,f1,,] <- inverse(src_tau[src,f1,,]);
      src_cov[src,f1,,] <- src_var[src,f1,,] %*% rho[src,f1,,] %*% src_var[src,f1,,];
      src_Sigma[src,f1,,] <- inverse(src_cov[src,f1,,]);

      # each source data point is distributed normally according to the source means and precisions
      for(r in 1:n.rep[src]){
        SOURCE_array[src,,f1,r] ~ dmnorm(src_mu[src,,f1],src_Sigma[src,f1,,]);
      } 
    } # end loop over f1
  } # end source data fitting loop
", file=filename, append=T)  
}

if(source$data_type=="raw" && source$by_factor==FALSE){  # fit the source means and precisions (not by factor 1)
cat("
  # fit source data (big for loop over sources)
  for(src in 1:n.sources){
    # uninformative priors on source means (src_mu vector)
    for(iso in 1:n.iso){
      src_mu[src,iso] ~ dnorm(0,.001);
    }

    # uninformative priors on source variances (src_tau matrix)
    for(i in 2:n.iso){
      for(j in 1:(i-1)){
        src_tau[src,i,j] <- 0;
        src_tau[src,j,i] <- 0;
      }
    }
    for(i in 1:n.iso){
      src_tau[src,i,i] ~ dgamma(.001,.001);
    }

    # uninformative priors on source correlations (rho matrix)
    for(i in 2:n.iso){
      for(j in 1:(i-1)){
        rho[src,i,j] ~ dunif(-1,1);
        rho[src,j,i] <- rho[src,i,j];
      }
    }
    for(i in 1:n.iso){
      rho[src,i,i] <- 1;
    }

    # Construct source precision matrix (src_Sigma)
    src_var[src,,] <- inverse(src_tau[src,,]);
    src_cov[src,,] <- src_var[src,,] %*% rho[src,,] %*% src_var[src,,];
    src_Sigma[src,,] <- inverse(src_cov[src,,]);

    # each source data point is distributed normally according to the source means and precisions
    for(r in 1:n.rep[src]){
      SOURCE_array[src,,r] ~ dmnorm(src_mu[src,],src_Sigma[src,,]);
    } 
  } # end source data fitting loop
", file=filename, append=T)  
}

# Here we fit the source means and variances according to Gelman, p.79-80:
#   u | sig2,y ~ Normal(m, sig2/n), (Eqn 3.8)
#   sig2 | y ~ Inv-X2(n, n-1/n * s2), (Eqn 3.9)
# where u and sig2 are the true source mean and variance,
# m and s2 are the sample mean and variance, and n is the sample size
# We are using an uninformative prior, so k_0 = v_0 = 0
# To simulate sig2, draw tmp.X from chisqr(n), then src_tau = tmp.X/(n-1)*s2  (Gelman p.580)

# If we have only residual error (SIAR model), then we only care about src_mu (src_tau and frac_sig2 go unused)
# If we have process error (MixSIR model), then we use src_mu and src_tau (and frac_sig2, which isn't fit to data yet)

# src_mu and src_tau are used in the X[i,iso] ~ dnorm call at the very bottom of this function
# dim(MU/SIG2_array) is either (n.src,n.iso) or (n.src,n.iso,n.fac), if source$by_factor is FALSE or TRUE, respectively.
if(source$data_type=="means" && source$by_factor==FALSE){ 
cat("
  for(src in 1:n.sources){
    for(iso in 1:n.iso){
      src_mu[src,iso] ~ dnorm(MU_array[src,iso], n_array[src]/SIG2_array[src,iso]);  # Eqn 3.8 but with precision instead of variance
      tmp.X[src,iso] ~ dchisqr(n_array[src]);
      src_tau[src,iso] <- tmp.X[src,iso]/(SIG2_array[src,iso]*(n_array[src] - 1));   # Eqn 3.9, following the simulation on p.580
    }
  }
", file=filename, append=T)
}

if(source$data_type=="means" && source$by_factor==TRUE){ 
cat("
  for(src in 1:n.sources){
    for(f1 in 1:source_factor_levels){
      for(iso in 1:n.iso){
        src_mu[src,iso,f1] ~ dnorm(MU_array[src,iso,f1], n_array[src,f1]/SIG2_array[src,iso,f1]);  # Eqn 3.8 but with precision instead of variance
        tmp.X[src,iso,f1] ~ dchisqr(n_array[src,f1]);
        src_tau[src,iso,f1] <- tmp.X[src,iso,f1]/(SIG2_array[src,iso,f1]*(n_array[src,f1] - 1));   # Eqn 3.9, following the simulation on p.580
      }
    }
  }
", file=filename, append=T) 
} 

###########################################################################
# ILR transform and factors
###########################################################################

cat("
    # draw p.global (global proportion means) from an uninformative Dirichlet,
    # then ilr.global is the ILR-transform of p.global
    p.global[1:n.sources] ~ ddirch(alpha[1:n.sources]);
    for(src in 1:(n.sources-1)){
      gmean[src] <- prod(p.global[1:src])^(1/src);
      ilr.global[src] <- sqrt(src/(src+1))*log(gmean[src]/p.global[src+1]); # page 296, Egozcue 2003
    }
", file=filename, append=T)

if(mix$n.effects > 0 && mix$FAC[[1]]$re){ # factor 1 is random effect
cat("
   fac1.sig ~ dunif(0,20);
   fac1.invSig2 <- 1/(fac1.sig*fac1.sig);
   # draw the fac1 (region) specific ILR terms (random effect)
   for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
          ilr.fac1[f1,src] ~ dnorm(0,fac1.invSig2);
      }
   }   
", file=filename, append=T)}

if(mix$n.effects > 0 && !mix$FAC[[1]]$re){ # factor 1 is fixed effect
cat("
  # draw the fac1 specific ILR terms (fixed effect)
  for(src in 1:(n.sources-1)){
    ilr.fac1[1,src] <- 0;
    for(f1 in 2:factor1_levels){
      ilr.fac1[f1,src] ~ dnorm(0,1);
    }
  }
", file=filename, append=T)}

if(mix$n.effects > 1 && mix$FAC[[2]]$re){ # factor 2 is random effect
cat("
   fac2.sig ~ dunif(0,20);
   fac2.invSig2 <- 1/(fac2.sig*fac2.sig);
   # draw the fac2-specific ILR terms (random effect)
   for(f2 in 1:factor2_levels){
      for(src in 1:(n.sources-1)){
          ilr.fac2[f2,src] ~ dnorm(0,fac2.invSig2);
      }
   }
", file=filename, append=T)}

if(mix$n.effects > 1 && !mix$FAC[[2]]$re){ # factor 2 is fixed effect
cat("
  # draw the fac2 specific ILR terms (fixed effect)
  for(src in 1:(n.sources-1)){
    ilr.fac2[1,src] <- 0;
    for(f2 in 2:factor2_levels){
      ilr.fac2[f2,src] ~ dnorm(0,1);
    }
  }
", file=filename, append=T)}

# Continuous Effects section
ilr.cont.string <- ""
if(mix$n.ce > 0){
cat("
   # Priors on all n.ce continuous effects (slopes for a linear regression in ilr-space)", file=filename, append=T)
for(ce in 1:mix$n.ce){
cat("
   for(src in 1:(n.sources-1)){
      ilr.cont",ce,"[src] ~ dnorm(0,.001)
   }
", file=filename, append=T, sep="")
ilr.cont.string <- paste(ilr.cont.string," + ilr.cont",ce,"[src]*Cont.",ce,"[i]",sep="")}
}

# if(indiv_effect){  # if user has chosen to include Individual as a random effect
# cat("
#    ind.sig ~ dunif(0,20); 
#    ind.invSig2 <- 1/(ind.sig*ind.sig);
#    # generate individual deviates from the global/region/pack mean
#    for(i in 1:N) {
#       for(src in 1:(n.sources-1)) {
#          ilr.ind[i,src] ~ dnorm(0, ind.invSig2);", file=filename, append=T)
# } else { # user has chosen to NOT include Individual as a random effect - remove ind.sig and set ilr.ind=0
cat("
   # DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] <- 0;", file=filename, append=T)
# } # end Individual effect block

# Add ilr.tot line, adding ilr.global, ilr.fac's, ilr.ind, and ilr.cont's
#   This will be different according to n.effects, but *MUST* be directly after the individual effect section!!
#   Can condense in the future by creating a "ilr.rand.string" object in a loop so this line is the same for all n.re
if(mix$n.effects==2){
cat("
         ilr.tot[i,src] <- ilr.global[src] + ilr.fac1[Factor.1[i],src] + ilr.fac2[Factor.2[i],src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file=filename, append=T, sep="")}
if(mix$n.effects==1){
cat("
         ilr.tot[i,src] <- ilr.global[src] + ilr.fac1[Factor.1[i],src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file=filename, append=T, sep="")}
if(mix$n.effects==0){
cat("
         ilr.tot[i,src] <- ilr.global[src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file=filename, append=T, sep="")}

#####################################################################
# Inverse ILR section
#####################################################################
cat("   
   # Inverse ILR math (equation 24, page 294, Egozcue 2003)
   for(i in 1:N){
      for(j in 1:(n.sources-1)){
        cross[i,,j] <- (e[,j]^ilr.tot[i,j])/sum(e[,j]^ilr.tot[i,j]);
      }
      for(src in 1:n.sources){
        tmp.p[i,src] <- prod(cross[i,src,]);
      }
      for(src in 1:n.sources){
        p.ind[i,src] <- tmp.p[i,src]/sum(tmp.p[i,]);
      }
   }
      
   for(src in 1:n.sources) {
      for(i in 1:N){
         # these are weights for variances
         p2[i,src] <- p.ind[i,src]*p.ind[i,src];
      }
   }
", file=filename, append=T)

if(mix$n.effects > 0){
  if(mix$fac_nested[1]){
  cat("
     # Transform ilr.fac1 into p.fac1 (fac1 nested within fac2)
     for(f1 in 1:factor1_levels) {
        for(src in 1:(n.sources-1)) {
            ilr.fac1.tot[f1,src] <- ilr.global[src] + ilr.fac2[factor2_lookup[f1],src] + ilr.fac1[f1,src];
            cross.fac1[f1,,src] <- (e[,src]^ilr.fac1.tot[f1,src])/sum(e[,src]^ilr.fac1.tot[f1,src]);
        }
        for(src in 1:n.sources) {
          tmp.p.fac1[f1,src] <- prod(cross.fac1[f1,src,]);
        }
        for(src in 1:n.sources){
          p.fac1[f1,src] <- tmp.p.fac1[f1,src]/sum(tmp.p.fac1[f1,]);
        }
      }
  ", file=filename, append=T)
  } else {
  cat("   
   # Transform ilr.fac1 into p.fac1 (fac1 not nested within fac2)
   for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
        ilr.fac1.tot[f1,src] <- ilr.global[src] + ilr.fac1[f1,src];
        cross.fac1[f1,,src] <- (e[,src]^ilr.fac1.tot[f1,src])/sum(e[,src]^ilr.fac1.tot[f1,src]);
      }
      for(src in 1:n.sources) {
        tmp.p.fac1[f1,src] <- prod(cross.fac1[f1,src,]);
      }
      for(src in 1:n.sources){
        p.fac1[f1,src] <- tmp.p.fac1[f1,src]/sum(tmp.p.fac1[f1,]);
      }
   }
", file=filename, append=T)
  } # end nested ilr.fac1 section
}

if(mix$n.effects > 1){
  if(mix$fac_nested[2]){
  cat("
     # Transform ilr.fac2 into p.fac2 (fac2 nested within fac1)
      for(f2 in 1:factor2_levels){
        for(src in 1:(n.sources-1)){
            ilr.fac2.tot[f2,src] <- ilr.global[src] + ilr.fac1[factor1_lookup[f2],src] + ilr.fac2[f2,src];
            cross.fac2[f2,,src] <- (e[,src]^ilr.fac2.tot[f2,src])/sum(e[,src]^ilr.fac2.tot[f2,src]);
        }
        for(src in 1:n.sources) {
          tmp.p.fac2[f2,src] <- prod(cross.fac2[f2,src,]);
        }
        for(src in 1:n.sources){
          p.fac2[f2,src] <- tmp.p.fac2[f2,src]/sum(tmp.p.fac2[f2,]);
        }
      }
  ", file=filename, append=T)
  } else {
  cat("
     # Transform ilr.fac2 into p.fac2 (fac2 not nested within fac1)
      for(f2 in 1:factor2_levels){
        for(src in 1:(n.sources-1)){
            ilr.fac2.tot[f2,src] <- ilr.global[src] + ilr.fac2[f2,src];
            cross.fac2[f2,,src] <- (e[,src]^ilr.fac2.tot[f2,src])/sum(e[,src]^ilr.fac2.tot[f2,src]);
        }
        for(src in 1:n.sources) {
          tmp.p.fac2[f2,src] <- prod(cross.fac2[f2,src,]);
        }
        for(src in 1:n.sources){
          p.fac2[f2,src] <- tmp.p.fac2[f2,src]/sum(tmp.p.fac2[f2,]);
        }
      }
  ", file=filename, append=T)
  } # end nested ilr.fac2 section
} # end n.re=2 section

###############################################################################
# mix.mu section
###############################################################################
cat("
    
   # for each isotope and population, calculate the predicted mixtures
   for(iso in 1:n.iso) {
      for(i in 1:N) {
", file=filename, append=T)

if(source$by_factor==T && source$conc_dep==T){
  cat("
         mix.mu[iso,i] <- (inprod(src_mu[,iso,Factor.1[i]],(p.ind[i,]*conc[,iso])) + inprod(frac_mu[,iso],(p.ind[i,]*conc[,iso]))) / inprod(p.ind[i,],conc[,iso]);", file=filename, append=T)
} else if(source$by_factor==T && source$conc_dep==F){
  cat("
         mix.mu[iso,i] <- inprod(src_mu[,iso,Factor.1[i]],p.ind[i,]) + inprod(frac_mu[,iso],p.ind[i,]);", file=filename, append=T)
} else if(source$by_factor==F && source$conc_dep==T){
  cat("  
         mix.mu[iso,i] <- (inprod(src_mu[,iso],(p.ind[i,]*conc[,iso])) + inprod(frac_mu[,iso],(p.ind[i,]*conc[,iso]))) / inprod(p.ind[i,],conc[,iso]);", file=filename, append=T)
} else if(source$by_factor==F && source$conc_dep==F){
  cat("  
         mix.mu[iso,i] <- inprod(src_mu[,iso],p.ind[i,]) + inprod(frac_mu[,iso],p.ind[i,]);", file=filename, append=T)
}
cat("
      }
   }
", file=filename, append=T)

################################################################################
# mix.prcsn section (error terms)
################################################################################

# # Add process/MixSIR error (define mix.var or set equal to 0)
# if(process_err){ # if process_err = TRUE, there are two varieties of mix.var (source$by_factor = T/F)
#   if(source$by_factor){ # source$by_factor = TRUE, include process error as:
# cat("
#          process.var[iso,i] <- inprod(1/src_tau[,iso,Factor.1[i]],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file=filename, append=T)
#   } else {  # source$by_factor = FALSE, include process error as:
# cat("
#          process.var[iso,i] <- inprod(1/src_tau[,iso],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file=filename, append=T)  
#   } # end if/else(source$by_factor)
# } else {  # process_err = FALSE, set process.var[iso,i] = 0:
# cat("
#          process.var[iso,i] <- 0;", file=filename, append=T)
# } # end if/else(process_err)

# # Calculate total variance (mix.var)
# if(resid_err){ # resid_err==TRUE
# cat("
#          mix.var[iso,i] <- process.var[iso,i] + resid.var[iso];
#          mix.prcsn[iso,i] <- 1/mix.var[iso,i];
#       }
#    }
# ", file=filename, append=T)
# } else { # resid_err==FALSE
#   if(resid_err_mult){
# cat("
#          mix.var[iso,i] <- process.var[iso,i] * resid.prop[iso];
#          mix.prcsn[iso,i] <- 1/mix.var[iso,i];
#       }
#    }
# ", file=filename, append=T)
#   } else {  # resid_err==FALSE AND resid_err_mult==FALSE (process_err only = MixSIR)
# cat("
#          mix.var[iso,i] <- process.var[iso,i];
#          mix.prcsn[iso,i] <- 1/mix.var[iso,i];
#       }
#    }
# ", file=filename, append=T)
#   }
# } # end resid_err==FALSE

###############################################################################
# Error structure section
###############################################################################
if(!resid_err && (mix$N > 1)){
  cat("
    # Multiplicative residual error
    for(iso in 1:n.iso){
      resid.prop[iso] ~ dchisqr(3);
    }
", file=filename, append=T)

  if(source$data_type=="means"){
cat("
    
   # Calculate process variance for each isotope and population
   for(iso in 1:n.iso) {
      for(i in 1:N) {
", file=filename, append=T)
    if(source$by_factor){ # source$by_factor = TRUE, include process error as:
cat("
         process.var[iso,i] <- inprod(1/src_tau[,iso,Factor.1[i]],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file=filename, append=T)
    } else {  # source$by_factor = FALSE, include process error as:
cat("
         process.var[iso,i] <- inprod(1/src_tau[,iso],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file=filename, append=T)  
    }
cat("
      }
   }

   # Construct Sigma, the mixture precision matrix
   for(ind in 1:N){
      for(i in 1:n.iso){
        for(j in 1:n.iso){
          Sigma.ind[ind,i,j] <- equals(i,j)*resid.prop[i]/process.var[i,ind];
        }
      }
   }

   # Likelihood
   for(i in 1:N) {
     X_iso[i,] ~ dmnorm(mix.mu[,i], Sigma.ind[i,,]);
   }
} # end model

", file=filename, append=T)
  } # end source$data_type=="means"

  if(source$data_type=="raw"){
cat("

  # resid.prop matrix
  for(i in 1:n.iso){
    for(j in 1:n.iso){
      resid.prop.mat[i,j] <- sqrt(resid.prop[i]*resid.prop[j]);
    }
  }

  # Construct mix covariance
  for(ind in 1:N){
    for(i in 1:n.iso){
      for(j in 1:n.iso){
", file=filename, append=T)
    if(source$by_factor){ # source$by_factor = TRUE, include process error as:
cat("
         mix.cov[ind,i,j] <- equals(i,j)*resid.prop[i]*(inprod(src_cov[,Factor.1[ind],i,j],p2[ind,]) + inprod(frac_sig2[,i],p2[ind,])) + (1-equals(i,j))*inprod(src_cov[,Factor.1[ind],i,j],p2[ind,])*resid.prop.mat[i,j];", file=filename, append=T)  
    } else {  # source$by_factor = FALSE, include process error as:
cat("
         mix.cov[ind,i,j] <- equals(i,j)*resid.prop[i]*(inprod(src_cov[,i,j],p2[ind,]) + inprod(frac_sig2[,i],p2[ind,])) + (1-equals(i,j))*inprod(src_cov[,i,j],p2[ind,])*resid.prop.mat[i,j];", file=filename, append=T)  
    }
cat("
      }
    }
    Sigma.ind[ind,,] <- inverse(mix.cov[ind,,]);
  }

   # Likelihood
  for(i in 1:N){
    X_iso[i,] ~ dmnorm(mix.mu[,i], Sigma.ind[i,,]);
  }
} # end model

", file=filename, append=T)
  } # end source$data_type=="raw"
} # end multiplicative residual error section


if(resid_err && (mix$N > 1)){
cat("
   # Mixture covariance prior (residual error only model)
   Sigma ~ dwish(I,n.iso+1);

   # Likelihood
   for(i in 1:N) {
     X_iso[i,] ~ dmnorm(mix.mu[,i], Sigma);
   }
} # end model

", file=filename, append=T)
}

# if(mix$N == 1){ # If only one datapoint, must fit MixSIR model (no residual error terms)

# }

################################################################################
# Likelihood
################################################################################
# cat("
#    # This section does the likelihood / posterior, N data points
#    for(i in 1:N) {
#       for(iso in 1:n.iso) {
#          X_iso[i,iso] ~ dnorm(mix.mu[iso,i], mix.prcsn[iso,i]);
#       }
#    }
# }


# ", file=filename, append=T)
} # end function write_JAGS_model