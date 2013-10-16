# Function: write_JAGS_model
# Input: raw_source_data (T/F), sources_by_factor (T/F), n.re (0,1,2), resid_err (T/F), process_err (T/F), include_indiv (T/F), n.ce (non-neg integer)
# Output: JAGS model text file saved to WD as "MixSIAR_model.txt"

# May 28
# Added continuous effects section (modified Include Indiv section)
# Changed 'num_factors' to 'n.re' (# random effects) to match the new 'n.ce' variable (# continuous effects)
# Added source mean and var fitting for raw_source_data = FALSE (Ward 2010 paper, "Including Source Uncertainty...")

# July 8
# Removed nesting - removed f1 dimension from from ilr.fac2 and derived objects (ilr.fac2.tot, cross.fac2, tmp.p.fac2, p.fac2)
# Added 'factor1_lookup' to be able to find the fac1 level indexed by fac2, e.g. [1 1 1 2 2 2 2 3]
# Added 'hierarch' (is Factor 2 within Factor 1?)
#   hierarch=T --> ilr.fac2.tot = ilr.global + ilr.fac1 + ilr.fac2
#   hierarch=F --> ilr.fac2.tot = ilr.global + ilr.fac2

# August 28
# Added concentration dependence (modified function input and mix.mu lines)

write_JAGS_model <- function(raw_source_data, sources_by_factor, n.re, resid_err, process_err, include_indiv, n.ce, hierarch, include_conc){
cat(paste("# raw_source_data: ",raw_source_data,sep=""), file="MixSIAR_model.txt")
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# sources_by_factor: ",sources_by_factor,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# # random effects: ",n.re,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# nested random effects: ",hierarch,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# # continuous effects: ",n.ce,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# resid_err: ",resid_err,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# process_err: ",process_err,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# include_indiv: ",include_indiv,sep=""), file="MixSIAR_model.txt", append=T)
cat("
", file="MixSIAR_model.txt", append=T)
cat(paste("# include_conc: ",include_conc,sep=""), file="MixSIAR_model.txt", append=T)
cat("

model{", file="MixSIAR_model.txt", append=T)

################################################################################
# Setup src_mu and src_tau arrays
# Fit source means (if raw_source_data==TRUE)
################################################################################

if(raw_source_data==TRUE && sources_by_factor==TRUE){ # fit the source means and precisions (by factor 1)
cat("
  # uninformed priors on source means and precisions
  for(src in 1:n.sources){
    for(f1 in 1:source_factor_levels){
      for(iso in 1:n.iso){
        src_mu[src,iso,f1] ~ dnorm(0,.001)
        src_tau[src,iso,f1] ~ dgamma(.001,.001)
      }
    }
  }
  # each source data point is distributed normally according to the source means and precisions
  for(src in 1:n.sources){
    for(f1 in 1:source_factor_levels){
      for(r in 1:n.rep[src,f1]){
        for(iso in 1:n.iso){
          SOURCE_array[src,iso,f1,r] ~ dnorm(src_mu[src,iso,f1],src_tau[src,iso,f1])
        }
      } 
    }
  }
", file="MixSIAR_model.txt", append=T)
}

if(raw_source_data==TRUE && sources_by_factor==FALSE){  # fit the source means and precisions (not by factor 1)
cat("
  # uninformed priors on source means and precisions
  for(src in 1:n.sources){
      for(iso in 1:n.iso){
        src_mu[src,iso] ~ dnorm(0,.001)
        src_tau[src,iso] ~ dgamma(.001,.001)
      }
  }
  # each source data point is distributed normally according to the source means and precisions
  for(src in 1:n.sources){
      for(r in 1:n.rep[src]){
        for(iso in 1:n.iso){
          SOURCE_array[src,iso,r] ~ dnorm(src_mu[src,iso],src_tau[src,iso])
        }
      } 
  }
", file="MixSIAR_model.txt", append=T)  
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
# dim(MU/SIG2_array) is either (n.src,n.iso) or (n.src,n.iso,n.fac), if sources_by_factor is FALSE or TRUE, respectively.
if(raw_source_data==FALSE && sources_by_factor==FALSE){ 
cat("
  for(src in 1:n.sources){
    for(iso in 1:n.iso){
      src_mu[src,iso] ~ dnorm(MU_array[src,iso], n_array[src]/SIG2_array[src,iso]);  # Eqn 3.8 but with precision instead of variance
      tmp.X[src,iso] ~ dchisqr(n_array[src]);
      src_tau[src,iso] <- tmp.X[src,iso]/(SIG2_array[src,iso]*(n_array[src] - 1));   # Eqn 3.9, following the simulation on p.580
    }
  }
", file="MixSIAR_model.txt", append=T)
}

if(raw_source_data==FALSE && sources_by_factor==TRUE){ 
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
", file="MixSIAR_model.txt", append=T) 
} 

################################################################################
# ILR transform and factors
################################################################################
# Include 1 of the 3 following sections: n.re = 2, 1, or 0

if(n.re==2){
cat("
    # draw p.global (global proportion means) from an uninformative Dirichlet,
    # then ilr.global is the ILR-transform of p.global
    p.global[1:n.sources] ~ ddirch(alpha[1:n.sources]);
    for(src in 1:(n.sources-1)){
      gmean[src] <- prod(p.global[1:src])^(1/src);
      ilr.global[src] <- sqrt(src/(src+1))*log(gmean[src]/p.global[src+1]); # page 296, Egozcue 2003
    }

   fac1.sig ~ dunif(0,20);
   fac1.invSig2 <- 1/(fac1.sig*fac1.sig);
   # draw the fac1 (region) specific ILR terms
   for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
          ilr.fac1[f1,src] ~ dnorm(0,fac1.invSig2);
      }
   }
   
   fac2.sig ~ dunif(0,20);
   fac2.invSig2 <- 1/(fac2.sig*fac2.sig);
   # draw the fac2-specific ILR terms
   for(f2 in 1:factor2_levels){
      for(src in 1:(n.sources-1)){
          ilr.fac2[f2,src] ~ dnorm(0,fac2.invSig2);
      }
   }
", file="MixSIAR_model.txt", append=T)

ilr.cont.string <- ""
if(n.ce > 0){
cat("
   # Priors on all n.ce continuous effects (slopes for a linear regression in ilr-space)", file="MixSIAR_model.txt", append=T)
for(ce in 1:n.ce){
cat("
   for(src in 1:(n.sources-1)){
      ilr.cont",ce,"[src] ~ dnorm(0,.001)
   }
", file="MixSIAR_model.txt", append=T, sep="")
ilr.cont.string <- paste(ilr.cont.string," + ilr.cont",ce,"[src]*Cont.",ce,"[i]",sep="")
}
}

if(include_indiv){  # if user has chosen to include Individual as a random effect
cat("
   ind.sig ~ dunif(0,20); 
   ind.invSig2 <- 1/(ind.sig*ind.sig);
   # generate individual deviates from the global/region/pack mean
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] ~ dnorm(0, ind.invSig2);", file="MixSIAR_model.txt", append=T)
} else { # user has chosen to NOT include Individual as a random effect - remove ind.sig and set ilr.ind=0
cat("
   # DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] <- 0;", file="MixSIAR_model.txt", append=T)
} # end Individual effect block

# Create the ilr.tot line, adding ilr.global, ilr.fac's, ilr.ind, and ilr.cont's
cat("
         ilr.tot[i,src] <- ilr.global[src] + ilr.fac1[Factor.1[i],src] + ilr.fac2[Factor.2[i],src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file="MixSIAR_model.txt", append=T, sep="")

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
      for(i in 1:N) {
         # these are weights for variances
         p2[i,src] <- p.ind[i,src]*p.ind[i,src]; 
      }
   }
   
   # Transform ilr.fac1 into p.fac1
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
", file="MixSIAR_model.txt", append=T)

if(hierarch){
cat("
   # Transform ilr.fac2 into p.fac2 (hierarchical)
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
", file="MixSIAR_model.txt", append=T)
} else {
cat("
   # Transform ilr.fac2 into p.fac2 (not hierarchical)
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
", file="MixSIAR_model.txt", append=T)
} # end hierarch section  
} # end n.re==2

################################################################################
# n.re==1

if(n.re==1){
cat("
    # draw p.global (global proportion means) from an uninformative Dirichlet,
    # then ilr.global is the ILR-transform of p.global
    p.global[1:n.sources] ~ ddirch(alpha[1:n.sources]);
    for(src in 1:(n.sources-1)){
      gmean[src] <- prod(p.global[1:src])^(1/src);
      ilr.global[src] <- sqrt(src/(src+1))*log(gmean[src]/p.global[src+1]); # page 296, Egozcue 2003
    }

   fac1.sig ~ dunif(0,20);
   fac1.invSig2 <- 1/(fac1.sig*fac1.sig);
   # draw the fac1 (region) specific phi terms
   for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
          ilr.fac1[f1,src] ~ dnorm(0,fac1.invSig2);
      }
   }   
", file="MixSIAR_model.txt", append=T)

ilr.cont.string <- ""
if(n.ce > 0){
cat("
   # Priors on all n.ce continuous effects (slopes for a linear regression in ilr-space)", file="MixSIAR_model.txt", append=T)
for(ce in 1:n.ce){
cat("
   for(src in 1:(n.sources-1)){
      ilr.cont",ce,"[src] ~ dnorm(0,.001)
   }
", file="MixSIAR_model.txt", append=T, sep="")
ilr.cont.string <- paste(ilr.cont.string," + ilr.cont",ce,"[src]*Cont.",ce,"[i]",sep="")
}
}

if(include_indiv){  # if user has chosen to include Individual as a random effect
cat("
   ind.sig ~ dunif(0,20); 
   ind.invSig2 <- 1/(ind.sig*ind.sig);
   # generate individual deviates from the global/region/pack mean
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] ~ dnorm(0, ind.invSig2);", file="MixSIAR_model.txt", append=T)
} else { # user has chosen to NOT include Individual as a random effect - remove ind.sig and set ilr.ind=0
cat("
   # DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] <- 0;", file="MixSIAR_model.txt", append=T)
} # end Individual effect block

# Create the ilr.tot line, adding ilr.global, ilr.fac's, ilr.ind, and ilr.cont's
cat("
         ilr.tot[i,src] <- ilr.global[src] + ilr.fac1[Factor.1[i],src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file="MixSIAR_model.txt", append=T, sep="")

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
      for(i in 1:N) {
         # these are weights for variances
         p2[i,src] <- p.ind[i,src]*p.ind[i,src]; 
      }
   }
   
   # Transform ilr.fac1 into p.fac1
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
", file="MixSIAR_model.txt", append=T)
} # end n.re==1

################################################################################
# n.re==0

if(n.re==0){   # sources_by_factor must be FALSE
cat("
    # draw p.global (global proportion means) from an uninformative Dirichlet
    # then ilr.global is the ILR-transform of p.global
    p.global[1:n.sources] ~ ddirch(alpha[1:n.sources]);
    for(src in 1:(n.sources-1)){
      gmean[src] <- prod(p.global[1:src])^(1/src);
      ilr.global[src] <- sqrt(src/(src+1))*log(gmean[src]/p.global[src+1]); # page 296, Egozcue 2003
    }
", file="MixSIAR_model.txt", append=T)

ilr.cont.string <- ""
if(n.ce > 0){
cat("
   # Priors on all n.ce continuous effects (slopes for a linear regression in ilr-space)", file="MixSIAR_model.txt", append=T)
for(ce in 1:n.ce){
cat("
   for(src in 1:(n.sources-1)){
      ilr.cont",ce,"[src] ~ dnorm(0,.001)
   }
", file="MixSIAR_model.txt", append=T, sep="")
ilr.cont.string <- paste(ilr.cont.string," + ilr.cont",ce,"[src]*Cont.",ce,"[i]",sep="")
}
}

if(include_indiv){  # if user has chosen to include Individual as a random effect
cat("
   ind.sig ~ dunif(0,20); 
   ind.invSig2 <- 1/(ind.sig*ind.sig);
   # generate individual deviates from the global/region/pack mean
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] ~ dnorm(0, ind.invSig2);", file="MixSIAR_model.txt", append=T)
} else { # user has chosen to NOT include Individual as a random effect - remove ind.sig and set ilr.ind=0
cat("
   # DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
   for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
         ilr.ind[i,src] <- 0;", file="MixSIAR_model.txt", append=T)
} # end Individual effect block

# Create the ilr.tot line, adding ilr.global, ilr.fac's, ilr.ind, and ilr.cont's
cat("
         ilr.tot[i,src] <- ilr.global[src]",ilr.cont.string," + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
   }
", file="MixSIAR_model.txt", append=T, sep="")

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
", file="MixSIAR_model.txt", append=T)
} # end n.re==0

###############################################################################
# mix.mu section
###############################################################################
cat("
    
   # for each isotope and population, calculate the predicted mixtures
   for(iso in 1:n.iso) {
      for(i in 1:N) {
", file="MixSIAR_model.txt", append=T)

if(sources_by_factor==T && include_conc==T){
  cat("
         mix.mu[iso,i] <- (inprod(src_mu[,iso,Factor.1[i]],(p.ind[i,]*conc[,iso])) + inprod(frac_mu[,iso],(p.ind[i,]*conc[,iso]))) / inprod(p.ind[i,],conc[,iso]);", file="MixSIAR_model.txt", append=T)
} else if(sources_by_factor==T && include_conc==F){
  cat("
         mix.mu[iso,i] <- inprod(src_mu[,iso,Factor.1[i]],p.ind[i,]) + inprod(frac_mu[,iso],p.ind[i,]);", file="MixSIAR_model.txt", append=T)
} else if(sources_by_factor==F && include_conc==T){
  cat("  
         mix.mu[iso,i] <- (inprod(src_mu[,iso],(p.ind[i,]*conc[,iso])) + inprod(frac_mu[,iso],(p.ind[i,]*conc[,iso]))) / inprod(p.ind[i,],conc[,iso]);", file="MixSIAR_model.txt", append=T)
} else if(sources_by_factor==F && include_conc==F){
  cat("  
         mix.mu[iso,i] <- inprod(src_mu[,iso],p.ind[i,]) + inprod(frac_mu[,iso],p.ind[i,]);", file="MixSIAR_model.txt", append=T)
}

################################################################################
# Error terms
################################################################################

# Add process/MixSIR error (define mix.var or set equal to 0)
if(process_err){ # if process_err = TRUE, there are two varieties of mix.var (sources_by_factor = T/F)
  if(sources_by_factor){ # sources_by_factor = TRUE, include process error as:
cat("
         mix.var[iso,i] <- inprod(1/src_tau[,iso,Factor.1[i]],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file="MixSIAR_model.txt", append=T)
  } else {  # sources_by_factor = FALSE, include process error as:
cat("
         mix.var[iso,i] <- inprod(1/src_tau[,iso],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);", file="MixSIAR_model.txt", append=T)  
  } # end if/else(sources_by_factor)
} else {  # process_err = FALSE, set mix.var[iso,i] = 0:
cat("
         mix.var[iso,i] <- 0;", file="MixSIAR_model.txt", append=T)
} # end if/else(process_err)

# Add residual/SIAR error (define var.resid or set equal to 0)
if(resid_err){ # resid_err==TRUE
cat("
         prec.resid[iso,i] ~ dgamma(.001,.001);
         var.resid[iso,i] <- 1/prec.resid[iso,i];
         mix.totalVar[iso,i] <- mix.var[iso,i] + var.resid[iso,i];
         mix.prcsn[iso,i] <- 1/(mix.totalVar[iso,i]);
      }
   }
", file="MixSIAR_model.txt", append=T)
} else { # resid_err==FALSE
cat("
         var.resid[iso,i] <- 0;
         mix.totalVar[iso,i] <- mix.var[iso,i] + var.resid[iso,i];
         mix.prcsn[iso,i] <- 1/(mix.totalVar[iso,i]);
      }
   }
", file="MixSIAR_model.txt", append=T)
}

################################################################################
# Likelihood
################################################################################
cat("
   # This section does the likelihood / posterior, N data points
   for(i in 1:N) {
      for(iso in 1:n.iso) {
         X_iso[i,iso] ~ dnorm(mix.mu[iso,i], mix.prcsn[iso,i]);
      }
   }
}


", file="MixSIAR_model.txt", append=T)
} # end function write_JAGS_model
