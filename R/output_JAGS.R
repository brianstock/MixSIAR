#' Process mixing model output from JAGS
#'
#' \code{output_JAGS} processes the mixing model output, prints and saves (in the
#' working directory):
#' \itemize{
#'  \item diagnostics
#'  \item summary statistics
#'  \item posterior density plots
#'  \item pairs plot
#'  \item trace/XY plots
#' }
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
#'   }
#'   
#' @return \code{p.both} -- only if 2 fixed effects OR 1 fixed + 1 random, otherwise \code{NULL}).
#' 
#' \code{p.both} holds the MCMC chains for the estimated proportions at the different factor levels. Dimensions = [n.draws, f1.levels, f2.levels, n.sources].
#' 
#' Calculated by combining the ilr offsets from global intercept:
#'   ilr.both[,f1,f2,src] = ilr.global[,src] + ilr.fac1[,f1,src] + ilr.fac2[,f2,src]
#' And then transforming from ilr- to proportion-space.
#' @export
#'   
output_JAGS <- function(jags.1, mix, source, output_options=list(
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
                                                  diag_save_ggmcmc = TRUE)){             # Save ggmcmc diagnostics as pdf?
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

######################################################################
# Posterior density plots
######################################################################
if(!output_options[[3]]){   # if 'suppress posterior plots' is NOT checked
  n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])   # number of posterior draws
  if(mix$n.fe == 0){ # only if there are no fixed effects, otherwise p.global is meaningless
    # Posterior density plot for p.global
    dev.new()
    df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
    for(i in 1:n.sources){
      df$x[seq(1+n.draws*(i-1),i*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.global[,i]) # fill in the p.global[i] values
      df$sources[seq(1+n.draws*(i-1),i*n.draws)] <- rep(source_names[i],n.draws)  # fill in the source names
    }
    my.title <- "Overall Population"
    print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
            ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
            ggplot2::theme_bw() +
            ggplot2::xlab("Proportion of Diet") +
            ggplot2::ylab("Scaled Posterior Density") +
            ggplot2::xlim(0,1) +
            ggplot2::labs(title = my.title) +
            ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.pdf",sep=""))  # svalue(plot_post_name)
      dev.copy2pdf(file=mypath)
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.png",sep=""))  # svalue(plot_post_name)
      dev.copy(png,mypath)
    }
  }

  if(n.effects >= 1 & mix$n.fe != 2){
    # Posterior density plots for p.fac1's
    for(f1 in 1:mix$FAC[[1]]$levels){    # formerly factor1_levels
      dev.new()
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac1[,f1,src]) # fill in the p.fac1[f1] values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      my.title <- mix$FAC[[1]]$labels[f1]  # formerly factor1_names
      print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
#        geom_density(alpha=.3) +
  ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
  ggplot2::xlim(0,1) +
  ggplot2::theme_bw() +
  ggplot2::xlab("Proportion of Diet") +
  ggplot2::ylab("Scaled Posterior Density") +
  ggplot2::labs(title = my.title) +
  ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],".pdf",sep=""))  # svalue(plot_post_name), factor1_names
        dev.copy2pdf(file=mypath)
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],".png",sep=""))  # svalue(plot_post_name), factor1_names
        dev.copy(png,mypath)
      }
    } # end p.fac1 posterior plots

    if(n.re==2){
      # Posterior density plots for p.fac2's
        for(f2 in 1:mix$FAC[[2]]$levels){  # formerly factor2_levels
          dev.new()
          df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
          for(src in 1:n.sources){
            df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac2[,f2,src]) # fill in the p.fac2 values
            df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
          }
          #my.title <- paste(factor1_names[f1],", ", random_effects[2]," ",f2,sep="") # plot title (ex. "Region 1, Pack 3")
          my.title <- mix$FAC[[2]]$labels[f2] # formerly factor2_names
          print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
#             geom_density(alpha=.3) +
  ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
  ggplot2::theme_bw() +
  ggplot2::xlim(0,1) +
  ggplot2::xlab("Proportion of Diet") +
  ggplot2::ylab("Scaled Posterior Density") +
  ggplot2::labs(title = my.title) +
  ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

          # Save the plot as a pdf file
          if(output_options[[4]]){ # svalue(plot_post_save_pdf)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[2]]$labels[f2],".pdf",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy2pdf(file=mypath)
          }
          if(output_options[[18]]){  # svalue(plot_post_save_png)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[2]]$labels[f2],".png",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy(png,mypath)
          }
        }# end p.fac2 posterior plots
    } # end if(n.re==2)
  } # end if(n.effects >=1 & n.fe != 2)

  # Posterior density plots for p.both (when 2 FE or 1FE + 1RE)
  if(mix$fere){
    for(f1 in 1:mix$FAC[[1]]$levels) {
      for(f2 in fac2_lookup[[f1]]){
        dev.new()
        df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
        for(src in 1:n.sources){
          df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.both[,f1,f2,src]) # fill in the p.both values
          df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
        }
        my.title <- paste(mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],sep=" ") # formerly factor2_names
        print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                ggplot2::theme_bw() +
                ggplot2::xlim(0,1) +
                ggplot2::xlab("Proportion of Diet") +
                ggplot2::ylab("Scaled Posterior Density") +
                ggplot2::labs(title = my.title) +
                ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

        # Save the plot as a pdf file
        if(output_options[[4]]){ # svalue(plot_post_save_pdf)
          mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".pdf",sep="")) #  svalue(plot_post_name), factor2_names
          dev.copy2pdf(file=mypath)
        }
        if(output_options[[18]]){  # svalue(plot_post_save_png)
          mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".png",sep="")) #  svalue(plot_post_name), factor2_names
          dev.copy(png,mypath)
        }
      } # f2
    } # f1
  }

  # Posterior density plot for fac1.sig, fac2.sig, and ind.sig
  if(n.re > 0 || output_options[[17]]){ # only have an SD posterior plot if we have Individual, Factor1, or Factor2 random effects)
    dev.new()
    n.re_ind <- n.re + as.numeric(output_options[[17]]) # this*n.draws will be the length of the plot data frame
    level <- c()
    x <- c()
    if(output_options[[17]]){ # if Individual is in the model, add ind.sig to the SD plot
      level <- c(level,rep("Individual SD",n.draws))
      x <- c(x,jags.1$BUGSoutput$sims.list$ind.sig)
    }
    if(n.re==1){ # if Factor.1 is in the model, add fac1.sig to the SD plot
      if(mix$FAC[[1]]$re){
        level <- c(level,rep(paste(mix$FAC[[1]]$name," SD",sep=""),n.draws))
        x <- c(x,jags.1$BUGSoutput$sims.list$fac1.sig)
      } else { # FAC 2 is the random effect
        level <- c(level,rep(paste(mix$FAC[[2]]$name," SD",sep=""),n.draws))
        x <- c(x,jags.1$BUGSoutput$sims.list$fac2.sig)
      }
    }
    if(n.re==2){ # if Factor.2 is in the model, add fac1.sig and fac2.sig to the SD plot
      level <- c(level,rep(paste(random_effects[1]," SD",sep=""),n.draws), rep(paste(random_effects[2]," SD",sep=""),n.draws))
      x <- c(x,jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig)
    }
    df2 <- data.frame(level=level, x=x) # create the SD plot data frame

    print(ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
#        geom_density(alpha=.3) +
  ggplot2::geom_density(alpha=.3) +
  ggplot2::theme_bw() +
  ggplot2::xlab(expression(sigma)) +
  ggplot2::ylab("Posterior Density") +
  ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))   # + xlim(0,2)

    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_SD.pdf",sep=""))  # svalue(plot_post_name)
      dev.copy2pdf(file=mypath)
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_SD.png",sep=""))  # svalue(plot_post_name)
      dev.copy(png,mypath)
    }
  }
}

# Calculate the summary statistics for the variables we're interested in (p.global's and factor SD's, maybe p.ind's)
# We print them out later, at the very bottom
sig_labels <- NULL; ind_labels <- NULL; fac1_labels <- NULL; fac2_labels <- NULL; sig_stats <- NULL;
getQuant <- function(x) quantile(x,probs=c(.025,.05,.25,.5,.75,.95,.975))
getMeanSD <- function(x) cbind(round(apply(x,2,mean),3),round(apply(x,2,sd),3))

stats <- NULL
sig_stats <- NULL
sig_labels <- NULL
eps_stats <- NULL
eps_labels <- NULL
# print(mix)
# print(mix$n.fe)
if(mix$n.fe == 0){
  global_quants <- t(round(apply(jags.1$BUGSoutput$sims.list$p.global,2,getQuant),3))
  global_means <- getMeanSD(jags.1$BUGSoutput$sims.list$p.global)
  stats <- cbind(global_means, global_quants)
  global_labels <- rep(NA,n.sources)
  for(src in 1:n.sources){
    global_labels[src] <- paste("p.global.",source_names[src],sep="")
  }
  rownames(stats) <- global_labels
}
if(n.effects > 0 & mix$n.fe != 2){
  fac1_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
  fac1_quants <- t(apply(fac1_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  fac1_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),sd),3))$value)
  fac1_stats <- cbind(fac1_means,fac1_quants)
  fac1_labels <- rep(NA,mix$FAC[[1]]$levels*n.sources)
  for(src in 1:n.sources){
    for(f1 in 1:mix$FAC[[1]]$levels){
      fac1_labels[mix$FAC[[1]]$levels*(src-1)+f1] <- paste("p.",mix$FAC[[1]]$labels[f1],".",source_names[src],sep="")
    }
  }
  rownames(fac1_stats) <- fac1_labels
  stats <- rbind(stats,fac1_stats)
  if(mix$FAC[[1]]$re){
    sig_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac1.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac1.sig,2,getQuant),3)))
    sig_labels <- paste(mix$FAC[[1]]$name,".SD",sep="")
  }
}
if(n.re==2){
  fac2_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
  fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  fac2_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),sd),3))$value)
  fac2_stats <- cbind(fac2_means,fac2_quants)
  fac2_labels <- rep(NA,mix$FAC[[2]]$levels*n.sources)
  for(src in 1:n.sources){
    for(f2 in 1:mix$FAC[[2]]$levels){
      fac2_labels[mix$FAC[[2]]$levels*(src-1)+f2] <- paste("p.",mix$FAC[[2]]$labels[f2],".",source_names[src],sep="")
    }
  }
  rownames(fac2_stats) <- fac2_labels
  stats <- rbind(stats,fac2_stats)
  if(mix$FAC[[2]]$re){
    sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
    sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
  }
}
if(mix$fere){
  fac2_quants <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=7)
  fac2_means <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=2)
  fac2_labels <- rep(NA,n.sources*length(unlist(fac2_lookup)))
  i <- 1
  for(f1 in 1:mix$FAC[[1]]$levels) {
    for(f2 in fac2_lookup[[f1]]){
      for(src in 1:n.sources){
        fac2_quants[i,] <- getQuant(p.both[,f1,f2,src])
        fac2_means[i,] <- c(mean(p.both[,f1,f2,src]),sd(p.both[,f1,f2,src]))
        fac2_labels[i] <- paste("p",mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],source_names[src],sep=".")
        i <- i+1
      }
    }
  }
  # fac2_quants <- as.matrix(cast(melt(round(apply(p.both,c(2,3,4),getQuant,na.rm=TRUE),3)),X4+X3+X2~X1)[,-c(1,2)])
  # fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  # fac2_means <- cbind(melt(round(apply(p.fac2,c(2,3),mean),3))$value, melt(round(apply(p.fac2,c(2,3),sd),3))$value)
  fac2_stats <- round(cbind(fac2_means,fac2_quants),3)
  rownames(fac2_stats) <- fac2_labels
  stats <- rbind(stats,fac2_stats)
  if(mix$FAC[[2]]$re){
    sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
    sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
  }
}

if(output_options[[17]]){ # include_indiv (if Individual is in the model)
  ind_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(p.ind,c(2,3),getQuant),3)),X3+X2~X1)[,-c(1,2)])
  ind_quants <- t(apply(ind_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  ind_means <- cbind(reshape2::melt(round(apply(p.ind,c(2,3),mean),3))$value, reshape2::melt(round(apply(p.ind,c(2,3),sd),3))$value)
  ind_stats <- cbind(ind_means,ind_quants)
  ind_labels <- rep(NA,N*n.sources)
  for(src in 1:n.sources){
    for(j in 1:N){
      ind_labels[N*(src-1)+j] <- paste("p.Ind ",j,".",source_names[src],sep="")
    }
  }
  sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$ind.sig),t(round(apply(jags.1$BUGSoutput$sims.list$ind.sig,2,getQuant),3))))
  sig_labels <- c(sig_labels,"Individual.SD")
  rownames(ind_stats) <- ind_labels
  stats <- rbind(stats, ind_stats)
}

# Add SD stats to the top of the summary
rownames(sig_stats) <- sig_labels
stats <- rbind(sig_stats,stats)

# Add epsilon (multiplicative error term) to stat summary
# Also plot posterior density
epsTF <- "resid.prop" %in% names(jags.1$BUGSoutput$sims.list)
if(epsTF){
  eps_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$resid.prop),t(round(apply(jags.1$BUGSoutput$sims.list$resid.prop,2,getQuant),3)))
  eps_labels <- paste0("Epsilon.", 1:mix$n.iso)
  rownames(eps_stats) <- eps_labels
  stats <- rbind(eps_stats,stats)

  # posterior plot
  level <- c()
  x <- c()
  for(j in 1:mix$n.iso){
    level <- c(level,rep(eps_labels[j], n.draws))
    x <- c(x, jags.1$BUGSoutput$sims.list$resid.prop[,j])    
  }
  df2 <- data.frame(level=level, x=x) 
  
  dev.new()
  print(ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
  ggplot2::geom_density(alpha=.3) +
  ggplot2::theme_bw() +
  ggplot2::xlab(expression(epsilon)) +
  ggplot2::ylab("Posterior Density") +
  ggplot2::theme(legend.position=c(.95,.95), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))   # + xlim(0,2)

}
colnames(stats) <- c("Mean","SD","2.5%","5%","25%","50%","75%","95%","97.5%")

# Pack 1 stats only
#stats[grep("Pack 1",rownames(stats)),]

# Region stats only
#stats[grep("Region",rownames(stats)),]

# Region stats, by Region
# byVec <- function(x){ind <- NULL; for(i in 1:length(x)){ ind <- c(ind,grep(x[i],rownames(stats)))}; return(ind)}
# stats[byVec(mix$RE[[1]]$labels),]

# All means
# stats[,"Mean"]

# Region means only
# stats[byVec(mix$RE[[1]]$labels),"Mean"]

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

#print(gelman)

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
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

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
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

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
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

DIC <- jags.1$BUGSoutput$DIC
cat("
################################################################################
# Summary Statistics
################################################################################

DIC = ",DIC,sep="")
out1 <- capture.output(stats)
cat("
",out1,sep="\n")

if(output_options[[1]]){  # svalue(summary_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[2]],".txt",sep=""))  # svalue(summary_name)
  cat("
#################################################################
# Summary Statistics
#################################################################

DIC = ",DIC,sep="", file=mypath, append=FALSE)
cat("
",out1,sep="\n", file=mypath, append=TRUE)
}

# Plot any continuous effects
if(mix$n.ce > 0){
  plot_continuous_var(jags.1,mix,source,output_options)
}

# Use ggmcmc package to create diagnostic plots
if(!is.null(output_options$diag_save_ggmcmc)) if(output_options$diag_save_ggmcmc){
  diag_filename <- paste(getwd(),"/",output_options$diag_name,".pdf",sep="")
  ggmcmc::ggmcmc(ggmcmc::ggs(jags1.mcmc), file=diag_filename, plot=c("Rhat","geweke","density","traceplot","running","autocorrelation","crosscorrelation"))
}

# Return p.both if 2 FE or 1FE + 1RE
if(mix$fere){
  return(p.both)
} else return(NULL) # otherwise return nothing

} # end function output_JAGS
