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
#'  \item 'loo' estimates standard errors for the difference in LOO between two models
#' }
#'
#' @param x list of two or more \code{rjags} model objects (output from \code{\link{run_model}} function)
#' @param loo \code{TRUE/FALSE}: compute LOO if \code{TRUE} (preferred), compute WAIC if \code{FALSE}
#'   
#' @return Output from \code{loo::compare}, a vector or matrix with class 'compare.loo'. 
#' If exactly two MixSIAR models are provided (i.e. if \code{length(x)==2}), then 
#' the difference in expected predictive accuracy and the standard error of the 
#' difference are returned. \emph{A POSITIVE difference in LOO/WAIC means the expected 
#' predictive accuracy for the SECOND model is higher.} If more than two objects 
#' are provided then a matrix of summary information is returned.
#' 
#' @seealso \href{https://CRAN.R-project.org/package=loo}{'loo' package} and \href{https://link.springer.com/article/10.1007/s11222-016-9696-4}{Vehtari, Gelman, & Gabry (2017)}.
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
  if(n.mods==1){
    if(loo) y <- loo::loo(x[[1]]$BUGSoutput$sims.list$loglik)
    if(!loo) y <- loo::waic(x[[1]]$BUGSoutput$sims.list$loglik)
    cat("Only 1 model selected, returning LOO/WAIC: \n
If you meant to compare 2 models, retry using e.g. \n
'compare_models(x=list(jags.1, jags.2))'.")    
    return(y)
  }
  
  loo.list <- list()
  for(m in 1:n.mods){
    if(loo) loo.list[[m]] <- loo::loo(x[[m]]$BUGSoutput$sims.list$loglik)
    if(!loo) loo.list[[m]] <- loo::waic(x[[m]]$BUGSoutput$sims.list$loglik)
  }
  return(loo::compare(x=loo.list))
}
