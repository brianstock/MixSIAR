#' Plot proportions by a continuous covariate
#'
#' \code{plot_continuous_var} creates a plot of how the mixture proportions
#' change according to a continuous covariate, as well as plots of the mixture
#' proportions for the individuals with minimum, median, and maximum covariate
#' values. Called by \code{\link{output_JAGS}} if any continuous effects are in
#' the model.
#'
#' MixSIAR fits a continuous covariate as a linear regression in ILR/transform-space.
#' Two terms are fit for the proportion of each source: an intercept and a slope.
#' The plotted line uses the posterior median estimates of the intercept and slope, and
#' the lines are curved because of the ILR-transform back into proportion-space. The 
#' 95\% credible intervals are shaded.
#' 
#' If the model contains both a continuous AND a categorical (factor) covariate, MixSIAR
#' fits a different intercept term for each factor level and all levels share the
#' same slope term.
#'
#' @param jags.1 output from \code{\link{run_model}}
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param output_options list containing options for plots and saving, passed from \code{\link{output_JAGS}} or \code{\link{output_posteriors}}
#' @param alphaCI alpha level for credible intervals (default = 0.05, 95\% CI)
#' @param exclude_sources_below don't plot sources with median proportion below this level for entire range of continuous effect variable (default = 0.1)
#'
#' @seealso Francis et al. 2011
#' @export
plot_continuous_var <- function(jags.1, mix, source, output_options, alphaCI=0.05, exclude_sources_below=0.1){
# added only to pass R CMD check
# ilr.global <- x <- p.global <- p.ind <- sources <- ..scaled.. <- NULL
R2jags::attach.jags(jags.1)
n.sources <- source$n.sources
source_names <- source$source_names

for(ce in 1:mix$n.ce){
  if(mix$n.effects == 1){ # if there is a FE/RE in addition to continuous effect
    g <- vector("list", length = mix$FAC[[1]]$levels)
    for(f1 in 1:mix$FAC[[1]]$levels){
      g[[f1]] <- vector("list", 4)
      fac.lab <- mix$FAC[[1]]$labels[f1]
      label <- mix$cont_effects[ce]
      cont <- mix$CE[[ce]]
      ilr.cont <- get(paste("ilr.cont",ce,sep=""))

      n.plot = 200
      chain.len = dim(p.global)[1]
      Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
      ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
      for(src in 1:n.sources-1){
        for(i in 1:n.plot){
         ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i] + ilr.fac1[,f1,src]
        }
      }

      # Transform every draw from ILR-space to p-space
      e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
      for(i in 1:(n.sources-1)){
         e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
         e[,i] <- e[,i]/sum(e[,i])
      }
      # dummy variables for inverse ILR calculation
      cross <- array(data=NA,dim=c(n.plot, chain.len, n.sources, n.sources-1))  
      tmp <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
      p.plot <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
      for(i in 1:n.plot){
        for(d in 1:chain.len){
          for(j in 1:(n.sources-1)){
            cross[i,d,,j] <- (e[,j]^ilr.plot[i,j,d])/sum(e[,j]^ilr.plot[i,j,d]);
          }
          for(src in 1:n.sources){
            tmp[i,d,src] <- prod(cross[i,d,src,]);
          }
          for(src in 1:n.sources){
            p.plot[i,d,src] <- tmp[i,d,src]/sum(tmp[i,d,]);
          }
        }
      }
      # now take quantiles, after ILR transform of every draw
      get_high <- function(x){return(quantile(x, 1-alphaCI/2))}
      get_low <- function(x){return(quantile(x, alphaCI/2))}
      p.low <- apply(p.plot, c(1,3), get_low)
      p.high <- apply(p.plot, c(1,3), get_high)
      p.median <- apply(p.plot, c(1,3), median)
      colnames(p.median) <- source_names
      
      Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
      df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
      colnames(df) <- c("source","median","x","low","high")
      df$source <- factor(df$source, levels=source_names)
      
      # remove sources from plot with very low proportions
      rm.srcs <- apply(p.median, 2, function(x) all(x < exclude_sources_below))
      df <- subset(df, source %in% source_names[!rm.srcs])

      g[[f1]][[1]] <- ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
                        ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
                        ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
                        ggplot2::labs(title = fac.lab) +
                        ggplot2::ylab("Proportion") +
                        ggplot2::xlab(label) +
                        ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                          panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                          axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                          axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                          legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      if(!output_options[[3]]){ # if NOT suppressing plots
        dev.new()
        print(g[[f1]][[1]])
      }
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_",label,"_cont.pdf"))
        ggplot2::ggsave(mypath, plot=g[[f1]][[1]], width=7, height=5, units='in')
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_",label,"_cont.png"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g[[f1]][[1]], width=7, height=5, units='in', dpi=300)
      }      
      
      # Posterior plot for min(Cont1), median(Cont1), and max(Cont1)
      # Page 370 in Francis et al (2011)
      n.draws <- length(p.global[,1])
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      f1.cont <- cont[mix$FAC[[1]]$values == f1,1]
      min_ind <- which(cont==min(f1.cont) & mix$FAC[[1]]$values == f1)[1]   # find the index of min(Cont1)
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,min_ind,src]) # fill in the p.ind values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      cont.lab = format(cont[min_ind]*mix$CE_scale + mix$CE_center, digits=3)
      my.title <- paste0(fac.lab," ",label," = ",cont.lab)
      g[[f1]][[2]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
              ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
              ggplot2::theme_bw() +
              ggplot2::xlab("Proportion") +
              ggplot2::ylab("Scaled Posterior Density") +
              ggplot2::labs(title = my.title) +
              ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                             panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                             axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                             axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                             legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      if(!output_options[[3]]){ # if NOT suppressing plots
        dev.new()
        print(g[[f1]][[2]])
      }
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_min_",label,".pdf"))
        ggplot2::ggsave(mypath, plot=g[[f1]][[2]], width=7, height=5, units='in')
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_min_",label,".png"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g[[f1]][[2]], width=7, height=5, units='in', dpi=300)
      }   

      # median(Cont1)
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      is.odd <- length(f1.cont) %% 2                # mod 2 division - odd length will return 1, even length returns 0
      if(is.odd==1) {med_ind <- which(cont==median(f1.cont) & mix$FAC[[1]]$values == f1)}   # find the index of median(Cont1)
      if(is.odd==0){ # If Cont.1 is even, this finds the index of the value just below the median. Here, median(Cont.1) has no corresponding index.
        if(sum(cont==median(f1.cont)) < length(f1.cont)/2) med_ind <- which(cont==max(sort(f1.cont[which(f1.cont<median(f1.cont))])) & mix$FAC[[1]]$values == f1) else {
          med_ind <- which(cont==median(f1.cont))[1]
        }
      }   
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,med_ind,src]) # fill in the p.ind values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      cont.lab = format(cont[med_ind]*mix$CE_scale + mix$CE_center, digits=3)
      my.title <- paste0(fac.lab," ",label," = ",cont.lab)
      g[[f1]][[3]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
              ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
              ggplot2::theme_bw() +
              ggplot2::xlab("Proportion") +
              ggplot2::ylab("Scaled Posterior Density") +
              ggplot2::labs(title = my.title) +
              ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                             panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                             axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                             axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                             legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      if(!output_options[[3]]){ # if NOT suppressing plots
        dev.new()
        print(g[[f1]][[3]])
      }
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_median_",label,".pdf"))
        ggplot2::ggsave(mypath, plot=g[[f1]][[3]], width=7, height=5, units='in')
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_median_",label,".png"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g[[f1]][[3]], width=7, height=5, units='in', dpi=300)
      }   
      
      # max(Cont1)
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      max_ind <- which(cont==max(f1.cont) & mix$FAC[[1]]$values == f1)[1]   # find the index of max(Cont1)
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,max_ind,src])    # fill in the p.ind values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      cont.lab = format(cont[max_ind]*mix$CE_scale + mix$CE_center, digits=3)
      my.title <- paste0(fac.lab," ",label," = ",cont.lab)
      g[[f1]][[4]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
              ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
              ggplot2::theme_bw() +
              ggplot2::xlab("Proportion") +
              ggplot2::ylab("Scaled Posterior Density") +
              ggplot2::labs(title = my.title) +
              ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                             panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                             axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                             axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                             legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      if(!output_options[[3]]){ # if NOT suppressing plots
        dev.new()
        print(g[[f1]][[4]])
      }
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_max_",label,".pdf"))
        ggplot2::ggsave(mypath, plot=g[[f1]][[4]], width=7, height=5, units='in')
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(getwd(),paste0(output_options[[5]],"_",fac.lab,"_max_",label,".png"))  # svalue(plot_post_name), factor1_names
        ggplot2::ggsave(mypath, plot=g[[f1]][[4]], width=7, height=5, units='in', dpi=300)
      }  
    } # end loop over factor1 levels
  } # end if YES FE/RE

  if(mix$n.effects == 0){
    g <- vector("list", 4)
    label <- mix$cont_effects[ce]
    cont <- mix$CE[[ce]]
    ilr.cont <- get(paste("ilr.cont",ce,sep=""))

    n.plot = 200
    chain.len = dim(p.global)[1]
    Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
    ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
    for(src in 1:n.sources-1){
      for(i in 1:n.plot){
       ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i]
      }
    }

    # Transform regression lines from ILR-space to p-space
    e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
    for(i in 1:(n.sources-1)){
       e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
       e[,i] <- e[,i]/sum(e[,i])
    }
    # dummy variables for inverse ILR calculation
    cross <- array(data=NA,dim=c(n.plot, chain.len, n.sources, n.sources-1))  
    tmp <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
    p.plot <- array(data=NA,dim=c(n.plot, chain.len, n.sources))  
    for(i in 1:n.plot){
      for(d in 1:chain.len){
        for(j in 1:(n.sources-1)){
          cross[i,d,,j] <- (e[,j]^ilr.plot[i,j,d])/sum(e[,j]^ilr.plot[i,j,d]);
        }
        for(src in 1:n.sources){
          tmp[i,d,src] <- prod(cross[i,d,src,]);
        }
        for(src in 1:n.sources){
          p.plot[i,d,src] <- tmp[i,d,src]/sum(tmp[i,d,]);
        }
      }
    }
    # now take quantiles, after ILR transform of every draw
    get_high <- function(x){return(quantile(x, 1-alphaCI/2))}
    get_low <- function(x){return(quantile(x, alphaCI/2))}
    p.low <- apply(p.plot, c(1,3), get_low)
    p.high <- apply(p.plot, c(1,3), get_high)
    p.median <- apply(p.plot, c(1,3), median)
    colnames(p.median) <- source_names
    
    Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
    df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
    colnames(df) <- c("source","median","x","low","high")
    df$source <- factor(df$source, levels=source_names)
    
    # remove sources from plot with very low proportions
    rm.srcs <- apply(p.median, 2, function(x) all(x < exclude_sources_below))
    df <- subset(df, source %in% source_names[!rm.srcs])

    # Plot of Diet vs. Cont effect
    # Page 370 in Francis et al (2011)
    g[[1]] <- ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
            ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
            ggplot2::ylab("Proportion") +
            ggplot2::xlab(label) +
            ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
              panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
              axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
              axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
              legend.justification=c(0,1), legend.title=ggplot2::element_blank())

    if(!output_options[[3]]){ # if NOT suppressing plots
      dev.new()
      print(g[[1]])
    }
    
    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_cont.pdf"))
      ggplot2::ggsave(mypath, plot=g[[1]], width=7, height=5, units='in')
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_cont.png"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g[[1]], width=7, height=5, units='in', dpi=300)
    }    
    
    # Posterior plot for min(Cont1), median(Cont1), and max(Cont1)
    # Page 370 in Francis et al (2011)
    n.draws <- length(p.global[,1])
    df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
    min_ind <- which(cont==min(cont))[1]   # find the index of min(Cont1)
    for(src in 1:n.sources){
      df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,min_ind,src]) # fill in the p.ind values
      df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
    }
    my.title <- paste("Min(",label,")",sep="")
    g[[2]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
            ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
            ggplot2::theme_bw() +
            ggplot2::xlab("Proportion") +
            ggplot2::ylab("Scaled Posterior Density") +
            ggplot2::labs(title = my.title) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                           panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                           axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                           axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                           legend.justification=c(0,1), legend.title=ggplot2::element_blank())
    if(!output_options[[3]]){ # if NOT suppressing plots
      dev.new()
      print(g[[2]])
    }
    
    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_min_",label,".pdf"))
      ggplot2::ggsave(mypath, plot=g[[2]], width=7, height=5, units='in')
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_min_",label,".png"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g[[2]], width=7, height=5, units='in', dpi=300)
    }   

    # median(Cont1)
    df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
    is.odd <- length(cont) %% 2                # mod 2 division - odd length will return 1, even length returns 0
    if(is.odd==1) {med_ind <- which(cont==median(cont))}   # find the index of median(Cont1)
    if(is.odd==0){ # If Cont.1 is even, this finds the index of the value just below the median. Here, median(Cont.1) has no corresponding index.
      if(sum(cont==median(cont)) < length(cont)/2) med_ind <- which(cont==max(sort(cont[which(cont<median(cont))]))) else {
        med_ind <- which(cont==median(cont))[1]
      }
    }   
    for(src in 1:n.sources){
      df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,med_ind,src]) # fill in the p.ind values
      df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
    }
    my.title <- paste("Median(",label,")",sep="")
    g[[3]] <- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
            ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
            ggplot2::theme_bw() +
            ggplot2::xlab("Proportion") +
            ggplot2::ylab("Scaled Posterior Density") +
            ggplot2::labs(title = my.title) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                           panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                           axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                           axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                           legend.justification=c(0,1), legend.title=ggplot2::element_blank())
    if(!output_options[[3]]){ # if NOT suppressing plots
      dev.new()
      print(g[[3]])
    }
    
    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_median_",label,".pdf"))
      ggplot2::ggsave(mypath, plot=g[[3]], width=7, height=5, units='in')
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_median_",label,".png"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g[[3]], width=7, height=5, units='in', dpi=300)
    }   
    
    # max(Cont1)
    df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
    max_ind <- which(cont==max(cont))[1]   # find the index of max(Cont1)
    for(src in 1:n.sources){
      df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,max_ind,src])    # fill in the p.ind values
      df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
    }
    my.title <- paste("Max(",label,")",sep="")
    g[[4]]<- ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
            ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
            ggplot2::theme_bw() +
            ggplot2::xlab("Proportion") +
            ggplot2::ylab("Scaled Posterior Density") +
            ggplot2::labs(title = my.title) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels=c("0", "0.25","0.5","0.75","1")) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                           panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                           axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                           axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                           legend.justification=c(0,1), legend.title=ggplot2::element_blank())
    if(!output_options[[3]]){ # if NOT suppressing plots
      dev.new()
      print(g[[4]])
    }
    
    # Save the plot to file
    if(output_options[[4]]){ # svalue(plot_post_save_pdf)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_max_",label,".pdf"))
      ggplot2::ggsave(mypath, plot=g[[4]], width=7, height=5, units='in')
    }
    if(output_options[[18]]){ # svalue(plot_post_save_png)
      mypath <- file.path(getwd(),paste0(output_options[[5]],"_max_",label,".png"))  # svalue(plot_post_name), factor1_names
      ggplot2::ggsave(mypath, plot=g[[4]], width=7, height=5, units='in', dpi=300)
    } 
  } # end if NO FE/RE
} # end loop over ce

if(!is.null(output_options$return_obj)) if(output_options$return_obj) return(g) else return(NULL)
} #end plot_continuous_var function
