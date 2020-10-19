#' Output pairs plots
#'
#' \code{output_pairs} makes pairs plots for a fit MixSIAR model (only p.global)
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
output_pairs <- function(jags.1, mix, source, output_options=list(
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

# Post-processing for 2 FE or 1FE + 1RE
#   calculate p.both = ilr.global + ilr.fac1 + ilr.fac2
if(mix$fere){
  fac2_lookup <- list()
  for(f1 in 1:mix$FAC[[1]]$levels){
    fac2_lookup[[f1]] <- unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values==f1)])
  }
  ilr.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources-1))
  p.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources))
  cross.both <- array(data=NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels,n.sources,n.sources-1))
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }
  for(i in 1:n.draws){
    for(f1 in 1:mix$FAC[[1]]$levels) {
      for(f2 in fac2_lookup[[f1]]){
        for(src in 1:(n.sources-1)) {
          ilr.both[i,f1,f2,src] <- jags.1$BUGSoutput$sims.list$ilr.global[i,src] + jags.1$BUGSoutput$sims.list$ilr.fac1[i,f1,src] + jags.1$BUGSoutput$sims.list$ilr.fac2[i,f2,src];
          cross.both[i,f1,f2,,src] <- (e[,src]^ilr.both[i,f1,f2,src])/sum(e[,src]^ilr.both[i,f1,f2,src]);
          # ilr.both[,f1,f2,src] <- ilr.global[,src] + ilr.fac1[,f1,src] + ilr.fac2[,f2,src];
        }
        for(src in 1:n.sources) {
          p.both[i,f1,f2,src] <- prod(cross.both[i,f1,f2,src,]);
        }
        p.both[i,f1,f2,] <- p.both[i,f1,f2,]/sum(p.both[i,f1,f2,]);
      } # f2
    } # f1
  }
} # end fere

# Fancy pairs plot of p.global
# Contour plots in the upper right, histograms on the diagonal, correlation coefficients in the lower left
if(!output_options[[6]]){   # if 'suppress pairs plot' is NOT checked
  dev.new()
  # Function: panel.hist (from ?pairs)
  # Purpose: creates histograms on the diagonal of the pairs plot matrix
  panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr), add=TRUE)
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col='blue', xlim=c(0,1),...)
  }
  # Function: panel.cor (from http://personality-project.org/r/r.graphics.html)
  # Purpose: prints correlation coefficients in the lower panel,
  #          scales text sizes to the correlation coefficient magnitudes
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr), add=TRUE)
    par(usr = c(0, 1, 0, 1))
    r = (cor(x, y,use="pairwise"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * abs(r))
  }
  # Function: panel.contour (inspired by http://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay)
  # Purpose: replaces scatterplots with colored contour plots
  panel.contour <- function(x,y){
    n.lines <- 4  # number of contour lines
    my.cols <- rev(RColorBrewer::brewer.pal(n.lines, "RdYlBu"))   # gets some pretty colors
    z <- MASS::kde2d(x,y)   # calculates the 2D kernel density that the contour function needs
    contour(z, drawlabels=FALSE, nlevels=n.lines, col=my.cols, add=TRUE)
  }
  pairs(jags.1$BUGSoutput$sims.list$p.global, labels=source_names, diag.panel=panel.hist, lower.panel=panel.cor, upper.panel=panel.contour)

  # Save the plot to file
  if(output_options[[7]]){ # svalue(plot_pairs_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[8]],".pdf",sep=""))  # svalue(plot_pairs_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[19]]){ # svalue(plot_pairs_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[8]],".png",sep=""))  # svalue(plot_pairs_name)
    dev.copy(png,mypath)
  }
}
} # end function output_pairs
