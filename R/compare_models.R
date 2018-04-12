#' Compare the predictive accuracy of 2 or more MixSIAR models
#'
#' \code{compare_models} uses the \href{https://CRAN.R-project.org/package=loo}{'loo' package}
#' to compute LOO (leave-one-out cross-validation) or WAIC (widely applicable information criterion)
#' for 2 of more fit MixSIAR models.
#' 
#' LOO and WAIC are "methods for estimating pointwise out-of-sample prediction accuracy
#' from a fitted Bayesian model using the log-likelihood evaluated at the posterior
#' simulations of the parameter values". See \href{https://link.springer.com/article/10.1007/s11222-016-9696-4}{Vehtari, Gelman, & Gabry (2017)}.
#' In brief:
#' 
#' \itemize{
#'  \item LOO and WAIC are preferred over AIC or DIC
#'  \item LOO is more robust than WAIC
#'  \item \code{'loo'} estimates standard errors for the difference in LOO/WAIC between two models
#'  \item We can calculate the relative support for each model using LOO/WAIC weights
#' }
#'
#' @param x list of two or more \code{rjags} model objects (output from \code{\link{run_model}} function)
#' @param loo \code{TRUE/FALSE}: compute LOO if \code{TRUE} (preferred), compute WAIC if \code{FALSE}
#'   
#' @return Data frame with the following columns:
#' 
#' \itemize{
#'  \item \code{Model}: names of \code{x} (input list)
#'  \item \code{LOOic} / \code{WAIC}: LOO information criterion or WAIC
#'  \item \code{se_LOOic} / \code{se_WAIC}: standard error of LOOic / WAIC
#'  \item \code{dLOOic} / \code{dWAIC}: difference between each model and the model with lowest LOOic/WAIC. Best model has dLOOic = 0.
#'  \item \code{se_dLOOic} / \code{se_dWAIC}: standard error of the difference between each model and the model with lowest LOOic/WAIC
#'  \item \code{weight}: relative support for each model, calculated as Akaike weights (p.75 Burnham & Anderson 2002). Interpretation: "an estimate of the probability that the model will make the best predictions on new data, conditional on the set of models considered" (McElreath 2015).
#' }
#'
#' @seealso \href{https://CRAN.R-project.org/package=loo}{'loo' package}
#' @seealso \href{https://link.springer.com/article/10.1007/s11222-016-9696-4}{Vehtari, A, A Gelman, and J Gabry. 2017}. Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing.
#' @seealso Pages 75-88 in \href{http://www.springer.com/us/book/9780387953649}{Burnham, KP and DR Anderson. 2002}. Model selection and multimodel inference: a practical information-theoretic approach. Springer Science & Business Media.
#' @seealso Pages 188-201 in \href{https://www.crcpress.com/Statistical-Rethinking-A-Bayesian-Course-with-Examples-in-R-and-Stan/McElreath/p/book/9781482253443}{McElreath, R. 2016}. Statistical rethinking: a Bayesian course with examples in R and Stan. CRC Press.
#'
#' @examples 
#' \dontrun{
#' # Model 1 = wolf diet by Region + Pack
#' mix.1 <- load_mix_data(filename=mix.filename, 
#'                      iso_names=c("d13C","d15N"), 
#'                      factors=c("Region","Pack"), 
#'                      fac_random=c(TRUE,TRUE), 
#'                      fac_nested=c(FALSE,TRUE), 
#'                      cont_effects=NULL)
#' source.1 <- load_source_data(filename=source.filename, source_factors="Region",
#'                              conc_dep=FALSE, data_type="means", mix.1)
#' discr.1 <- load_discr_data(filename=discr.filename, mix.1)
#' 
#' # Run Model 1
#' jags.1 <- run_model(run="test", mix.1, source.1, discr.1, model_filename, 
#'                     alpha.prior = 1, resid_err=T, process_err=T)
#'                     
#' # Model 2 = wolf diet by Region (no Pack)
#' mix.2 <- load_mix_data(filename=mix.filename, 
#'                      iso_names=c("d13C","d15N"), 
#'                      factors=c("Region"), 
#'                      fac_random=c(TRUE), 
#'                      fac_nested=c(FALSE), 
#'                      cont_effects=NULL)
#' source.2 <- load_source_data(filename=source.filename, source_factors="Region",
#'                              conc_dep=FALSE, data_type="means", mix.2)
#' discr.2 <- load_discr_data(filename=discr.filename, mix.2)
#' 
#' # Run Model 2
#' jags.2 <- run_model(run="test", mix.2, source.2, discr.2, model_filename, 
#'                     alpha.prior = 1, resid_err=T, process_err=T)
#'                     
#' # Compare models 1 and 2 using LOO
#' compare_models(x=list(jags.1, jags.2), loo=TRUE)
#' 
#' # Compare models 1 and 2 using WAIC
#' compare_models(x=list(jags.1, jags.2), loo=FALSE)
#' 
#' # Get WAIC for model 1
#' compare_models(x=list(jags.1), loo=FALSE)
#' 
#' # Create named list of rjags objects to get model names in summary
#' x <- list(jags.1, jags.2)
#' names(x) <- c("Region + Pack", "Region")
#' compare_models(x)
#' }
#' 
compare_models <- function(x, loo=TRUE){
  if(class(x)!="list"){
    stop(paste("x is not a 'list'.
  Retry using e.g. 'compare_models(x=list(jags.1, jags.2))'.",sep=""))    
  }
  if(class(x[[1]])!="rjags"){
    stop(paste("x is not a list of 'rjags' models.
  Retry using e.g. 'compare_models(x=list(jags.1, jags.2))'.",sep=""))    
  }  

  n.mods <- length(x)
  mod.names <- names(x)
  if(n.mods==1){
    if(loo) y <- loo::loo(x[[1]]$BUGSoutput$sims.list$loglik)
    if(!loo) y <- loo::waic(x[[1]]$BUGSoutput$sims.list$loglik)
    cat("Only 1 model selected, returning LOO/WAIC:
If you meant to compare 2+ models, retry using e.g.
'compare_models(x=list(jags.1, jags.2))'.\n\n")    
    return(y)
  }
  
  loo.list <- list()
  nobs <- rep(NA, n.mods)
  for(m in 1:n.mods){
    nobs[m] <- dim(x[[m]]$BUGSoutput$sims.list$loglik)[2]
    if(loo) loo.list[[m]] <- loo::loo(x[[m]]$BUGSoutput$sims.list$loglik)
    if(!loo) loo.list[[m]] <- loo::waic(x[[m]]$BUGSoutput$sims.list$loglik)
  }

  if(length(unique(nobs))!=1){
    stop(paste("Number of mix data points not identical for all models.
  You can only use LOO/WAIC to compare models with the same data.",sep=""))    
  }  

  # Get LOOic, SE, dLOOic, and weights
  if(loo){ # LOO
    LOOic <- sapply(loo.list, function(x) round(x$estimates["looic", "Estimate"],1))
    se_LOOic <- sapply(loo.list, function(x) round(x$estimates["looic", "SE"],1))
  } else { # WAIC
    LOOic <- sapply(loo.list, function(x) round(x$estimates["waic", "Estimate"],1))
    se_LOOic <- sapply(loo.list, function(x) round(x$estimates["waic", "SE"],1))    
  }
  dLOOic <- LOOic - min(LOOic)
  weight <- round(exp(-0.5*dLOOic) / sum(exp(-0.5*dLOOic)), 3)

  # Get SE of pairwise differences in LOOic/WAIC
  dSE.mat <- matrix(NA, nrow=n.mods, ncol=n.mods)
  for(i in 1:(n.mods-1)){
    for(j in (i+1):n.mods){
        dSE.mat[i,j] <- 2*loo::compare(x=list(loo.list[[i]], loo.list[[j]]))['se'] # x2 to put on IC/deviance scale
        dSE.mat[j,i] <- dSE.mat[i,j]
    }#j
  }#i  
  dSE <- round(dSE.mat[, which(dLOOic==0)], 1) # find the model with dLOOic = 0, use that col of dSE.mat
  
  # Produce data frame ordered by LOOic/WAIC
  result <- data.frame(Model=mod.names, LOOic=LOOic, se_LOOic=se_LOOic, dLOOic=dLOOic, se_dLOOic=dSE, weight=weight)
  result <- result[order(result$LOOic),]
  if(!loo) names(result) <- c("Model", "WAIC", "se_WAIC", "dWAIC", "se_dWAIC", "weight")

  return(result)
}
