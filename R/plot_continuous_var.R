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
#' The plot uses the posterior median estimates of the intercept and slope, and
#' the lines are curved because of the ILR-transform back into proportion-space.
#'
#' @param jags.1 output from \code{\link{run_model}}
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param output_options list containing options for plots and saving,
#' passed from \code{\link{output_JAGS}}
#'
#' @seealso Francis et al. 2011
plot_continuous_var <- function(jags.1, mix, source, output_options){
# added only to pass R CMD check
# ilr.global <- x <- p.global <- p.ind <- sources <- ..scaled.. <- NULL
R2jags::attach.jags(jags.1)
n.sources <- source$n.sources
source_names <- source$source_names

for(ce in 1:mix$n.ce){
  label <- mix$cont_effects[ce]
  cont <- mix$CE[[ce]]
  ilr.cont <- get(paste("ilr.cont",ce,sep=""))

  get_high <- function(x){return(quantile(x,.975))}
  get_low <- function(x){return(quantile(x,.025))}
  Cont1.plot <- seq(from=round(min(cont),1),to=round(max(cont),1),by=0.1)
  ilr.median <- array(NA,dim=c(length(Cont1.plot),n.sources-1))
  for(src in 1:n.sources-1){
     ilr.median[,src] <- median(ilr.global[,src]) + median(ilr.cont[,src])*Cont1.plot
  }
  N_plot <- length(Cont1.plot)

  # Transform regression lines from ILR-space to p-space
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
     e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
     e[,i] <- e[,i]/sum(e[,i])
  }
  cross.med <- array(data=NA,dim=c(N_plot,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p.med <- array(data=NA,dim=c(N_plot,n.sources))              # dummy variable for inverse ILR calculation
  p.median <- array(data=NA,dim=c(N_plot,n.sources))
  for(i in 1:N_plot){
    for(j in 1:(n.sources-1)){
      cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    }
    for(src in 1:n.sources){
      tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    }
    for(src in 1:n.sources){
      p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    }
  }
  colnames(p.median) <- source_names

  Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
  df <- data.frame(reshape2::melt(p.median)[,2:3],rep(Cont1.plot,n.sources))
  colnames(df) <- c("source","median","x")

  # medians <- data.frame(cont,apply(p.ind,c(2,3),median))
  # colnames(medians) <- c("cont",source_names)
  # medians <- melt(medians,id="cont")

  # Plot of Diet vs. Cont effect
  # Page 370 in Francis et al (2011)
  dev.new()
  print(ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
          ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
           # geom_point(data=medians,aes(x=cont,y=value,colour=variable), guide=F) +
          ggplot2::ylab("Proportion of Diet") +
          ggplot2::xlab(label) +
          ggplot2::ylim(0,1) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank()))

  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",label,".pdf",sep=""))  # svalue(plot_post_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",label,".png",sep=""))  # svalue(plot_post_name)
    dev.copy(png,mypath)
  }

  # Posterior plot for min(Cont1), median(Cont1), and max(Cont1)
  # Page 370 in Francis et al (2011)
  n.draws <- length(p.global[,1])
  # min(Cont1)
  dev.new()
  df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
  min_ind <- which(cont==min(cont))   # find the index of min(Cont1)
  for(src in 1:n.sources){
    df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,min_ind,src]) # fill in the p.ind values
    df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
  }
  my.title <- paste("Diet of min(",label,") individual",sep="")
  print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
          ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
          ggplot2::theme_bw() +
          ggplot2::xlab("Proportion of Diet") +
          ggplot2::ylab("Scaled Posterior Density") +
          ggplot2::labs(title = my.title) +
          ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_min_",label,".pdf",sep=""))  # svalue(plot_post_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_min_",label,".png",sep=""))  # svalue(plot_post_name)
    dev.copy(png,mypath)
  }

  # median(Cont1)
  dev.new()
  df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
  is.odd <- length(cont) %% 2                # mod 2 division - odd length will return 1, even length returns 0
  if(is.odd==1) {med_ind <- which(cont==median(cont))}   # find the index of median(Cont1)
  if(is.odd==0) {med_ind <- which(cont==max(sort(cont[which(cont<median(cont))])))}   # If Cont.1 is even, this finds the index of the value just below the median. Here, median(Cont.1) has no corresponding index.
  for(src in 1:n.sources){
    df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,med_ind,src]) # fill in the p.ind values
    df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
  }
  my.title <- paste("Diet of median(",label,") individual",sep="")
  print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
          ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
          ggplot2::theme_bw() +
          ggplot2::xlab("Proportion of Diet") +
          ggplot2::ylab("Scaled Posterior Density") +
          ggplot2::labs(title = my.title) +
          ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_median_",label,".pdf",sep=""))  # svalue(plot_post_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_median_",label,".png",sep=""))  # svalue(plot_post_name)
    dev.copy(png,mypath)
  }

  # max(Cont1)
  dev.new()
  df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
  max_ind <- which(cont==max(cont))   # find the index of max(Cont1)
  for(src in 1:n.sources){
    df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,max_ind,src])    # fill in the p.ind values
    df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
  }
  my.title <- paste("Diet of max(",label,") individual",sep="")
  print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
          ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
          ggplot2::theme_bw() +
          ggplot2::xlab("Proportion of Diet") +
          ggplot2::ylab("Scaled Posterior Density") +
          ggplot2::labs(title = my.title) +
          ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))

  # Save the plot to file
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_max_",label,".pdf",sep=""))  # svalue(plot_post_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_max_",label,".png",sep=""))  # svalue(plot_post_name)
    dev.copy(png,mypath)
  }

} # end loop over ce
} #end plot_continuous_var function
