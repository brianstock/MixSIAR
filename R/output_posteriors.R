#' Output posterior density plots
#'
#' \code{output_posteriors} makes posterior density plots for a fit MixSIAR model
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
#' @return list of ggplot objects (if \code{return_obj = TRUE})
#' 
#' @export
#'   
output_posteriors <- function(jags.1, mix, source, output_options=list(
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

######################################################################
# Posterior density plots
######################################################################
g <- list()

n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])   # number of posterior draws
if(mix$n.fe == 0){ # only if there are no fixed effects, otherwise p.global is meaningless
  # Posterior density plot for p.global
  df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
  for(i in 1:n.sources){
    df$x[seq(1+n.draws*(i-1),i*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.global[,i]) # fill in the p.global[i] values
    df$sources[seq(1+n.draws*(i-1),i*n.draws)] <- rep(source_names[i],n.draws)  # fill in the source names
  }
  my.title <- "Overall Population"
  g$global <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
          ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
          ggplot2::theme_bw() +
          ggplot2::xlab("Proportion") +
          ggplot2::ylab("Scaled Posterior Density") +
          ggplot2::xlim(0,1) +
          ggplot2::labs(title = my.title) +
          ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank())
  if(!output_options[[3]]){ # if NOT suppressing plots
    dev.new()
    print(g$global)
  }

  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_global.pdf"))  # svalue(plot_post_name)
    ggplot2::ggsave(mypath, plot=g$global, width=7, height=5, units='in')
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_global.png"))  # svalue(plot_post_name)
    ggplot2::ggsave(mypath, plot=g$global, width=7, height=5, units='in', dpi=300)
  }
}

if(n.effects >= 1 & mix$n.fe != 2){
  # Posterior density plots for p.fac1's
  g$fac1 <- vector("list", length = mix$FAC[[1]]$levels)
  for(f1 in 1:mix$FAC[[1]]$levels){    # formerly factor1_levels
    df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
    for(src in 1:n.sources){
      df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac1[,f1,src]) # fill in the p.fac1[f1] values
      df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
    }
    my.title <- mix$FAC[[1]]$labels[f1]  # formerly factor1_names
    g$fac1[[f1]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                      ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                      ggplot2::xlim(0,1) +
                      ggplot2::theme_bw() +
                      ggplot2::xlab("Proportion") +
                      ggplot2::ylab("Scaled Posterior Density") +
                      ggplot2::labs(title = my.title) +
                      ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank())
    if(!output_options[[3]]){ # if NOT suppressing plots
      dev.new()
      print(g$fac1[[f1]])
    }
    
    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[1]]$labels[f1],".pdf"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g$fac1[[f1]], width=7, height=5, units='in')
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[1]]$labels[f1],".png"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g$fac1[[f1]], width=7, height=5, units='in', dpi=300)
    }
  } # end p.fac1 posterior plots

  if(n.re==2){
      # Posterior density plots for p.fac2's
      g$fac2 <- vector("list", length = mix$FAC[[2]]$levels)
      for(f2 in 1:mix$FAC[[2]]$levels){  # formerly factor2_levels
        df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
        for(src in 1:n.sources){
          df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac2[,f2,src]) # fill in the p.fac2 values
          df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
        }
        my.title <- mix$FAC[[2]]$labels[f2] # formerly factor2_names
        g$fac2[[f2]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                          ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                          ggplot2::theme_bw() +
                          ggplot2::xlim(0,1) +
                          ggplot2::xlab("Proportion") +
                          ggplot2::ylab("Scaled Posterior Density") +
                          ggplot2::labs(title = my.title) +
                          ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank())
        if(!output_options[[3]]){ # if NOT suppressing plots
          dev.new()
          print(g$fac2[[f2]])
        }
        
        # Save the plot to file
        if(output_options[[4]]){ # svalue(plot_post_save_pdf)
          mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[2]]$labels[f2],".pdf"))  # svalue(plot_post_name), factor1_names
          ggplot2::ggsave(mypath, plot=g$fac2[[f2]], width=7, height=5, units='in')
        }
        if(output_options[[18]]){ # svalue(plot_post_save_png)
          mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[2]]$labels[f2],".png"))  # svalue(plot_post_name), factor1_names
          ggplot2::ggsave(mypath, plot=g$fac2[[f2]], width=7, height=5, units='in', dpi=300)
        }
      }# end p.fac2 posterior plots
  } # end if(n.re==2)
} # end if(n.effects >=1 & n.fe != 2)

# Posterior density plots for p.both (when 2 FE or 1FE + 1RE)
if(mix$fere){
  g$both <- vector("list", length = mix$FAC[[1]]$levels)
  for(f1 in 1:mix$FAC[[1]]$levels) {
    g$both[[f1]] <- vector("list", length = mix$FAC[[2]]$levels)
    for(f2 in fac2_lookup[[f1]]){
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.both[,f1,f2,src]) # fill in the p.both values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      my.title <- paste(mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],sep=" ") # formerly factor2_names
      g$both[[f1]][[f2]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
              ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
              ggplot2::theme_bw() +
              ggplot2::xlim(0,1) +
              ggplot2::xlab("Proportion of Diet") +
              ggplot2::ylab("Scaled Posterior Density") +
              ggplot2::labs(title = my.title) +
              ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank())
      if(!output_options[[3]]){ # if NOT suppressing plots
        dev.new()
        print(g$both[[f1]][[f2]])
      }
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".pdf"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g$both[[f1]][[f2]], width=7, height=5, units='in')
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".png"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g$both[[f1]][[f2]], width=7, height=5, units='in', dpi=300)
      }
    } # f2
  } # f1
} # end p.both

# Posterior density plot for random effect variance terms (fac1.sig, fac2.sig, and ind.sig)
if(n.re > 0 || output_options[[17]]){ # only have an SD posterior plot if we have Individual, Factor1, or Factor2 random effects)
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

  g$sig <- ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
              ggplot2::geom_density(alpha=.3) +
              ggplot2::theme_bw() +
              ggplot2::xlab(expression(sigma)) +
              ggplot2::ylab("Posterior Density") +
              ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank())
  if(!output_options[[3]]){ # if NOT suppressing plots
    dev.new()
    print(g$sig)
  }
  
  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_SD.pdf"))
    ggplot2::ggsave(mypath, plot=g$sig, width=7, height=5, units='in')
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_SD.png"))
    ggplot2::ggsave(mypath, plot=g$sig, width=7, height=5, units='in', dpi=300)
  }
}

# epsilon (multiplicative error term)
epsTF <- "resid.prop" %in% names(jags.1$BUGSoutput$sims.list)
if(epsTF){
  eps_labels <- paste0("Epsilon.", 1:mix$n.iso)
  level <- c()
  x <- c()
  for(j in 1:mix$n.iso){
    level <- c(level,rep(eps_labels[j], n.draws))
    x <- c(x, jags.1$BUGSoutput$sims.list$resid.prop[,j])    
  }
  df2 <- data.frame(level=level, x=x) 
  
  g$epsilon <- ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
                ggplot2::geom_density(alpha=.3) +
                ggplot2::theme_bw() +
                ggplot2::xlab(expression(epsilon)) +
                ggplot2::ylab("Posterior Density") +
                ggplot2::theme(legend.position=c(.95,.95), legend.justification=c(1,1), legend.title=ggplot2::element_blank())

  if(!output_options[[3]]){ # if NOT suppressing plots
    dev.new()
    print(g$epsilon)
  }
  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_epsilon.pdf"))
    ggplot2::ggsave(mypath, plot=g$epsilon, width=7, height=5, units='in')
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(getwd(),paste0(output_options[[5]],"_epsilon.png"))
    ggplot2::ggsave(mypath, plot=g$epsilon, width=7, height=5, units='in', dpi=300)
  }
}

# Plot any continuous effects
if(mix$n.ce > 0){
  g$cont <- plot_continuous_var(jags.1, mix, source, output_options)
}

if(!is.null(output_options$return_obj)) if(output_options$return_obj) return(g) else return(NULL)
} # end function
