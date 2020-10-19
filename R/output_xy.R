#' Output XY/trace plots
#'
#' \code{output_xy} makes trace/XY plots for a fit MixSIAR model
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
#' @return NULL
#' 
#' @export
#'   
output_xy <- function(jags.1, mix, source, output_options=list(
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
                                                  return_obj = FALSE)){
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
jags1.mcmc <- coda::as.mcmc(jags.1)
n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])

###########################################################################################
# XY/Trace Plots
###########################################################################################

# XY plots for p.global and factor SD's
if(!output_options[[9]]){  # if 'suppress XY plot' is NOT checked
  # XY plot for p.global
  dev.new()
  print(lattice::xyplot(coda::as.mcmc(jags.1$BUGSoutput$sims.list$p.global),strip=lattice::strip.custom(factor.levels=source_names)))

  # Save the xy p.global plot to file
  if(output_options[[10]]){ # svalue(plot_xy_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_diet_p.pdf",sep=""))  # svalue(plot_xy_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[20]]){ # svalue(plot_xy_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_diet_p.png",sep=""))  # svalue(plot_xy_name)
    dev.copy(png,mypath)
  }

  # XY plot for the factor SDs
  if(output_options[[17]]){ # include_indiv ('ind.sig' is in the model)
    dev.new()
    traceplot_labels <- rep("",length(random_effects)+1)  # +1 because we need to add "Individual SD"
    if(n.re > 0){
      for(i in 1:length(random_effects)){
        traceplot_labels[i] <- paste(random_effects[i]," SD",sep="")
      }
    }
    traceplot_labels[length(random_effects)+1] <- "Individual SD"
    if(n.re==2) print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
    if(n.re==1){
      if(mix$FAC[[1]]$re){
        print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
      } else { # FAC 2 is the 1 random effect
        print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac2.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
      }
    }
    if(n.re==0) print(lattice::xyplot(coda::as.mcmc(jags.1$BUGSoutput$sims.list$ind.sig),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
  } else { # Individual SD is not in the model (no 'ind.sig')
    if(n.re > 0){
      dev.new()
      traceplot_labels <- rep("",length(random_effects))
      for(i in 1:length(random_effects)) { traceplot_labels[i] <- paste(random_effects[i]," SD",sep="") }
      if(n.re==2) print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
      if(n.re==1){
        if(mix$FAC[[1]]$re){
          print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
        } else { # FAC 2 is the 1 random effect
          print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac2.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
        }
      }
    }
  }
  # Save the xy factor SD plot to file
  if(output_options[[10]]){ # svalue(plot_xy_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_SD.pdf",sep=""))  # svalue(plot_xy_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[20]]){ # svalue(plot_xy_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_SD.png",sep=""))  # svalue(plot_xy_name)
    dev.copy(png,mypath)
  }
}

} # end function output_JAGS
