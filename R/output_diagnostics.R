#' Output diagnostics
#'
#' \code{output_diagnostics} returns diagnostics for a fit MixSIAR model
#'
#' @param jags.1 rjags model object, output from \code{\link{run_model}} function
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param output_options list containing options for plots and saving:
#'   \itemize{
#'    \item \code{summary_save}: Save the summary statistics as a txt file? Default = \code{TRUE}
#'    \item \code{summary_name}: Summary statistics file name (.txt will be appended). Default = \code{"summary_statistics"}
#'    \item \code{sup_post}: Suppress posterior density plot output in R? Default = \code{FALSE}
#'    \item \code{plot_post_save_pdf}: Save posterior density plots as pdfs? Default = \code{TRUE}
#'    \item \code{plot_post_name}: Posterior plot file name(s) (.pdf/.png will be appended) Default = \code{"posterior_density"}
#'    \item \code{sup_pairs}: Suppress pairs plot output in R? Default = \code{FALSE}
#'    \item \code{plot_pairs_save_pdf}: Save pairs plot as pdf? Default = \code{TRUE}
#'    \item \code{plot_pairs_name}: Pairs plot file name (.pdf/.png will be appended) Default = \code{"pairs_plot"}
#'    \item \code{sup_xy}: Suppress xy/trace plot output in R? Default = \code{TRUE}
#'    \item \code{plot_xy_save_pdf}: Save xy/trace plot as pdf? Default = \code{FALSE}
#'    \item \code{plot_xy_name}: XY/trace plot file name (.pdf/.png will be appended) Default = \code{"xy_plot"}
#'    \item \code{gelman}: Calculate Gelman-Rubin diagnostic test? Default = \code{TRUE}
#'    \item \code{heidel}: Calculate Heidelberg-Welch diagnostic test? Default = \code{FALSE}
#'    \item \code{geweke}: Calculate Geweke diagnostic test? Default = \code{TRUE}
#'    \item \code{diag_save}: Save the diagnostics as a .txt file? Default = \code{TRUE}
#'    \item \code{diag_name}: Diagnostics file name (.txt will be appended) Default = \code{"diagnostics"}
#'    \item \code{indiv_effect}: artifact, set to FALSE 
#'    \item \code{plot_post_save_png}: Save posterior density plots as pngs? Default = \code{FALSE}
#'    \item \code{plot_pairs_save_png}: Save pairs plot as png? Default = \code{FALSE}
#'    \item \code{plot_xy_save_png}: Save xy/trace plot as png? Default = \code{FALSE}
#'    \item \code{diag_save_ggmcmc}: Save ggmcmc diagnostics as pdf? Default = \code{TRUE}
#'    \item \code{return_obj} Return ggplot objects for later modification? Default = \code{FALSE}
#'   }
#'   
#' @return named list of three data frames (one each for Gelman, Heidelberg-Welch, and Geweke), but only if \code{return_obj = TRUE}
#' 
#' @export
#'   
output_diagnostics <- function(jags.1, mix, source, output_options=list(
                                                  summary_save = TRUE,                 # Save the summary statistics as a txt file?
                                                  summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
                                                  sup_post = FALSE,                       # Suppress posterior density plot output in R?
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
                                                  diag_save_ggmcmc = TRUE,
                                                  return_obj = FALSE)){             # Save ggmcmc diagnostics as pdf?
mcmc.chains <- jags.1$BUGSoutput$n.chains
N <- mix$N
n.re <- mix$n.re
n.effects <- mix$n.effects
if(n.re==1){
  random_effects <- ifelse(mix$FAC[[1]]$re,mix$FAC[[1]]$name,mix$FAC[[2]]$name)
}
if(n.re==2){
  random_effects <- mix$factors
}
n.sources <- source$n.sources
source_names <- source$source_names
# p.global <- ilr.global <- ilr.fac1 <- ilr.fac2 <- fac1.sig <- fac2.sig <- NULL
# ind.sig <- ..scaled.. <- p.fac1 <- p.fac2 <- p.ind <- sources <- NULL
# R2jags::attach.jags(jags.1)
# jags1.mcmc <- coda::as.mcmc(jags.1)
as.mcmc.rjags <- function(x){
  n.chains <- x$n.chains
  sims <- x$sims.array
  n.thin <- x$n.thin
  if(n.chains==1) return(coda::mcmc(sims[, 1, ], thin=n.thin))
  out <- vector("list", length=n.chains)
  for (i in seq(n.chains)) out[[i]] <- coda::mcmc(sims[, i, ], thin=n.thin)
  out <- coda::mcmc.list(out)
  coda::varnames(out) <- dimnames(sims)[[3]]
  return(out)
}
jags1.mcmc <- as.mcmc.rjags(jags.1$BUGSoutput)
n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])

################################################################################
# Calulate diagnostics
################################################################################
# Get number of variables in the model
n.var <- coda::nvar(jags1.mcmc)
# Gelman-Rubin diagnostic
if(output_options[[12]]){  # if Gelman is checked
  if(mcmc.chains == 1){
    gelman <- "*** Error: Gelman diagnostic requires more than one chain ***"
  }
  if(mcmc.chains > 1){    # Gelman diagnostic requires more than one chain
    # Gelman diagnostic, for when the multivariate Gelman fails (matrix not positive definite)
    # Remove the test results for dummy/empty variables
    gelman <- matrix(NA, nrow=n.var, ncol=2)
    for (v in 1:coda::nvar(jags1.mcmc)) {
      gelman[v,] <- coda::gelman.diag(jags1.mcmc[,v])$psrf
    }
    #gelman <- gelman[ind,]
    colnames(gelman) <- c("Point est.","Upper C.I.")
    rownames(gelman) <- coda::varnames(jags1.mcmc)
    #rownames(gelman) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
    gelman.all <- gelman[which(!is.nan(gelman[,1])),] # Remove dummy variables (show up as NA)
    gelman_short <- gelman[order(gelman[,1],decreasing=T),]
    if(n.var>10) gelman_short <- gelman_short[1:10,]
    gelman_fail <- c(length(which(gelman[,1]>1.01)), length(which(gelman[,1]>1.05)), length(which(gelman[,1]>1.1)))
  }
}

# Heidelberger and Welch's diagnostic
# Remove the test results for dummy/empty variables
if(output_options[[13]]){   # if Heidel is checked
  heidel <- coda::heidel.diag(jags1.mcmc)
  w <- which(!is.na(heidel[[1]][,"pvalue"]))  # find all the non-dummy variables
  heidel.all <- data.frame(matrix(NA,nrow=length(w),ncol=3*mcmc.chains))  # create empty data frame
  colstring <- rep(NA,mcmc.chains*3)  # vector of column names
  for(i in 1:mcmc.chains){
    heidel.tmp <- as.data.frame(heidel[[i]][w,c("stest","pvalue","htest")]) # stest, pvalue, and htest are the relevant statistics - get them
    heidel.all[,(3*i-2):(3*i)] <- heidel.tmp
    colstring[(3*i-2):(3*i)] <- c(paste("stest.",i,sep=""), paste("pval.",i,sep=""), paste("hwtest.",i,sep="")) # create the appropriate column names
  }
  #heidel.all <- heidel.all[ind,]
  #rownames(heidel.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
  rownames(heidel.all) <- coda::varnames(jags1.mcmc)[w]
  colnames(heidel.all) <- colstring
  heidel.all <- round(heidel.all,3)
  heidel.all <- replace(heidel.all,heidel.all==0,"fail")  # A normal call to 'heidel.diag' prints "fail" and "pass", for some reason they turn to 0's and 1's
  heidel.all <- replace(heidel.all,heidel.all==1,"pass")  # when you access the statistics directly.  Here we turn the 0's and 1's back into "fail" and "pass"
  # When the stationarity test fails, hwtest returns <NA>...change these NAs to 'fail'
  heidel.all <- replace(heidel.all,is.na(heidel.all),"fail")
  # Count the number of failures (2 tests per chain - 'stationarity' and 'half-width')
  stest_fail <- rep(NA,mcmc.chains); hwtest_fail <- rep(NA,mcmc.chains)
  for(i in 1:mcmc.chains){
    stest_fail[i] <- sum(heidel.all[,3*i-2]=="fail")
    hwtest_fail[i] <- sum(heidel.all[,3*i]=="fail")
  }
  heidel_fail <- rbind(stest_fail,hwtest_fail)
  rownames(heidel_fail) <- c("Stationarity","Half-width")
  colnames(heidel_fail) <- paste("Chain",1:mcmc.chains)
}

# Geweke diagnostic
# Remove the test results for dummy/empty variables
if(output_options[[14]]){ # if Geweke is checked
  geweke <- coda::geweke.diag(jags1.mcmc)
  geweke.all <- data.frame(matrix(NA,nrow=n.var,ncol=mcmc.chains))    # create empty data frame
  colstring <- rep(NA,mcmc.chains)    # vector of column names
  for(i in 1:mcmc.chains){
    geweke.tmp <- as.data.frame(geweke[[i]]$z) # get the relevant geweke statistics
    geweke.all[,i] <- geweke.tmp
    colstring[i] <- c(paste("chain",i,sep=""))  # create the column names "chain1", "chain2", etc.
  }
  #geweke.all <- geweke.all[ind,]
  #rownames(geweke.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
  rownames(geweke.all) <- coda::varnames(jags1.mcmc)
  colnames(geweke.all) <- colstring
  geweke.all <- round(geweke.all,3)
  w <- which(!is.nan(geweke[[1]]$z))  # find all the non-dummy variables
  geweke.all <- geweke.all[w,]
  geweke_fail <- matrix(NA,nrow=1,ncol=mcmc.chains)
  for(i in 1:mcmc.chains){
    geweke_fail[1,i] <- sum(abs(geweke.all[,i])>1.96)
  }
  colnames(geweke_fail) <- paste("Chain",1:mcmc.chains)
  rownames(geweke_fail) <- "Geweke"
}

################################################################################
# Print diagnostics
################################################################################
diags <- list()
if(output_options[[12]]){  # svalue(gelman)
cat("
################################################################################
# Gelman-Rubin Diagnostic
################################################################################

Generally the Gelman diagnostic should be < 1.05

",paste("Out of ",n.var," variables: ",gelman_fail[1]," > 1.01",sep=""),"
                      ",paste(gelman_fail[2]," > 1.05",sep=""),"
                      ",paste(gelman_fail[3]," > 1.1",sep=""),"

The worst variables are:
",sep="")
print(gelman_short)
diags$gelman <- gelman.all

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(getwd(),paste0(output_options[[16]],".txt"))  # svalue(diag_name)
  out <- capture.output(gelman)
  out2 <- capture.output(gelman_short)
  cat("
################################################################################
# Gelman-Rubin Diagnostic
################################################################################

Generally the Gelman diagnostic should be < 1.05

",paste("Out of ",n.var," variables: ",gelman_fail[1]," > 1.01",sep=""),"
                      ",paste(gelman_fail[2]," > 1.05",sep=""),"
                      ",paste(gelman_fail[3]," > 1.1",sep=""),"

The worst variables are:
",out2,"

And here are the Gelman diagnostics for all variables:
",out,sep="\n", file=mypath, append=FALSE)
} # end save Gelman
} # end Gelman printout

if(output_options[[13]]){  # svalue(heidel)
cat("
################################################################################
# Heidelberger and Welch Diagnostic
################################################################################

A few failures is normal and acceptable...
Number of failures in each chain (out of ",n.var," variables):

",sep="")
print(heidel_fail)
#print(heidel.all)
diags$heidel <- heidel.all

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(getwd(),paste0(output_options[[16]],".txt"))  # svalue(diag_name)
  out <- capture.output(heidel.all)
  out2 <- capture.output(heidel_fail)
  cat("
################################################################################
# Heidelberger and Welch Diagnostic
################################################################################

A few failures is normal and acceptable...
Number of failures in each chain (out of ",n.var," variables):

",out2,"

And here are the Heidelberger-Welch diagnostics for all variables:
",out,sep="\n", file=mypath, append=output_options[[12]]) # svalue(gelman)
} # end save Heidel
} # end Heidel printout

if(output_options[[14]]){ # svalue(geweke)
cat("
################################################################################
# Geweke Diagnostic
################################################################################

The Geweke diagnostic is a standard z-score, so we'd expect 5% to be outside +/-1.96
Number of variables outside +/-1.96 in each chain (out of ",n.var,"):

",sep="")
print(geweke_fail)
#print(geweke.all)
diags$geweke <- geweke.all

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(getwd(), paste0(output_options[[16]],".txt"))  # svalue(diag_name)
  out <- capture.output(geweke.all)
  out2 <- capture.output(geweke_fail)
  cat("
################################################################################
# Geweke Diagnostic
################################################################################

The Geweke diagnostic is a standard z-score, so we'd expect 5% to be outside +/-1.96
Number of variables outside +/-1.96 in each chain (out of ",n.var,"):

",out2,"

And here are the Geweke diagnostics for all variables:
",out,sep="\n", file=mypath, append=output_options[[12]]||output_options[[13]]) # svalue(gelman) || svalue(heidel)
} # end Geweke save
} # end Geweke printout

# Use ggmcmc package to create diagnostic plots
if(!is.null(output_options$diag_save_ggmcmc)) if(output_options$diag_save_ggmcmc){
  diag_filename <- file.path(getwd(), paste0(output_options$diag_name,".pdf"))
  ggmcmc::ggmcmc(ggmcmc::ggs(jags1.mcmc), file=diag_filename, plot=c("Rhat","geweke","density","traceplot","running","autocorrelation","crosscorrelation"))
}

if(!is.null(output_options$return_obj)) if(output_options$return_obj) return(diags) else return(NULL)
} # end function output_diagnostics
