#' Plot biotracer data (1-D)
#'
#' \code{plot_data_one_iso} creates a 1-D plot of mix and source tracer data and
#' saves the plot to a file in the working directory
#'
#' An important detail is that \code{plot_data_one_iso} plots the raw mix data
#' and \emph{adds the TDF to the source data}, since this is the polygon that the
#' mixing model uses to determine proportions. The plotted source means are:
#' \deqn{\mu_source + \mu_discr}
#' The source error bars are +/- 1 standard deviation, \emph{calculated as a
#' combination of source and TDF variances:}
#' \deqn{\sqrt{\sigma^2_source + \sigma^2_discr}}
#'
#' \code{plot_data_one_iso} looks for 'C', 'N', 'S', and 'O' in the biotracer column
#' headers and assumes they are stable isotopes, labeling the axes with, e.g.,
#' expression(paste(delta^13, "C (u2030)",sep="")).
#'
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param discr output from \code{\link{load_discr_data}}
#' @param filename name of the plot file(s) to save (e.g. "isospace_plot")
#' @param plot_save_pdf T/F, save the plot(s) as a pdf?
#' @param plot_save_png T/F, save the plot(s) as a png?
#' @param return_obj T/F, whether or not to return ggplot object for further modification, defaults to F
#'
#' @seealso \code{\link{plot_data}}
plot_data_one_iso <- function(mix,source,discr,filename,plot_save_pdf,plot_save_png,return_obj=FALSE){
  # added only to pass R CMD check
  # x <- position_jitter <- scolour <- xmin <- xmax <- label <- NULL

  # Look in the isotope column headers for 'C', 'N', 'S', and 'O'
  # Make the x and y labels for the isospace plot
  if(length(grep("C",mix$iso_names))==1) x_label <- expression(paste(delta^13, "C (\u2030)",sep=""))
  if(length(grep("N",mix$iso_names))==1) x_label <- expression(paste(delta^15, "N (\u2030)",sep=""))
  if(length(grep("S",mix$iso_names))==1) x_label <- expression(paste(delta^34, "S (\u2030)",sep=""))
  if(length(grep("O",mix$iso_names))==1) x_label <- expression(paste(delta^18, "O (\u2030)",sep=""))
  if(length(grep("SP",mix$iso_names))==1) x_label <- expression(paste(delta^15, "N-SP (\u2030)",sep=""))
  if(!exists("x_label")) x_label <- mix$iso_names

  y_data <- 0.5
  y <- rep(y_data,mix$N)
  df <- data.frame(x = mix$data_iso, y = y)
  spacing <- 0.1

  if(!is.na(source$by_factor)){
    source_linetype <- sort(rep(1:source$n.sources,source$S_factor_levels))    # each source gets a different linetype (assumes S_MU is sorted by source and then factor, which it is)
    source_color <- factor(as.numeric(source$S_factor1))  # color sources by factor 1 (ex: region)
    index <- seq(from=1,to=1+(source$n.sources-1)*source$S_factor_levels,by=source$S_factor_levels)  # "index" gets the row in S_MU of the first instance of each source (for making the source labels)
    discr_mu_plot <- array(NA,dim=c(length(source$S_MU[,1]),mix$n.iso))     # Since discr$mu is not by factor, it needs to be expanded out by 'source_factor_levels' to match the dimensions of S_MU.  I.e. if n.sources=10, n.iso=2, and source_factor_levels=3 (condor data), frac_mu is 10x2 and S_MU is 30x2.  This makes frac_mu_plot, a 30x2 matrix.
    discr_sig2_plot <- array(NA,dim=c(length(source$S_MU[,1]),mix$n.iso))   # Same for discr$sig2
    for(i in 1:source$n.sources){
      discr_mu_plot[index[i]:(index[i]+source$S_factor_levels-1),] <- matrix(rep(discr$mu[i],source$S_factor_levels),nrow=source$S_factor_levels,ncol=mix$n.iso,byrow=T)
      discr_sig2_plot[index[i]:(index[i]+source$S_factor_levels-1),] <- matrix(rep(discr$sig2[i],source$S_factor_levels),nrow=source$S_factor_levels,ncol=mix$n.iso,byrow=T)
    }
    y_sources <- seq(y_data+0.2,(source$S_factor_levels*source$n.sources*spacing)-spacing+y_data+0.2,by=spacing)
    MU_plot <- source$S_MU[,mix$iso_names] + discr_mu_plot    # add discrimination mean to the source mean values
    SIG_plot <- sqrt(source$S_SIG[,mix$iso_names]^2 + discr_sig2_plot)  # add discrimination SD to the source SD values
  } else {  # source$by_factor = NA
    source_linetype <- 1:source$n.sources    # each source gets a different linetype
    source_color <- factor(rep("black",source$n.sources))  # this doesn't work...solution was to make separate ggplot calls for by_factor and not_by_factor
    index <- 1:source$n.sources              # "index" gets the row in S_MU of the first instance of each source (since not by factor, only one instance of each source)
    discr_mu_plot <- discr$mu
    discr_sig2_plot <- discr$sig2
    y_sources <- seq(y_data+0.2,(source$n.sources*spacing)-spacing+y_data+0.2,by=spacing)
    source$S_factor_levels <- 0.5   # just for correctly spacing the source labels (y in source.labels, line 61)
    MU_plot <- source$S_MU + discr_mu_plot    # add discrimination mean to the source mean values
    SIG_plot <- sqrt(source$S_SIG^2 + discr_sig2_plot)  # add discrimination SD to the source SD values
  }

  MU_plot <- as.vector(MU_plot)
  SIG_plot <- as.vector(SIG_plot)

  df_sources <- data.frame(x=MU_plot,
                           xmin = MU_plot - SIG_plot,
                           xmax = MU_plot + SIG_plot,
                           y = y_sources,
                           linetype = source_linetype,
                           scolour = source_color)

  source.labels <- data.frame(
    x = MU_plot[index],    # label sources just left
    #y = MU_plot[index,2] + rep(1.5,n.sources),    # and up from their means
    y = y_sources[index] + spacing*source$S_factor_levels,
    label = source$source_names
  )
  .e <- environment()
  dev.new()

  if(mix$n.effects==2){
    # ggplot2 will only make 6 different shapes, so force it to use enough for Factor.2
    shapes <- c(16,17,15,3,7,8,1,6,35,36,37,4,18,14,11,9,13)
    shapes <- shapes[1:mix$FAC[[2]]$levels]  # 1:factor2_levels
    if(!is.na(source$by_factor)){ # sources by factor, want to color the sources by factor1
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                       shape = factor(mix$FAC[[2]]$values)), position=position_jitter(width=.2,height=.1), show.legend=T) +    # Factor.2
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +   # factor1_names
        ggplot2::scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
        ggplot2::geom_point(data=df_sources,
                            ggplot2::aes(x=x, y=y,colour=scolour),
                   size=2,
                   show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                                ggplot2::aes(xmin=xmin,xmax=xmax,colour=scolour),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::ylab("") +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    } else { # sources not by factor (make the sources black)
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                       shape = factor(mix$FAC[[2]]$values)), position=position_jitter(width=.2,height=.1), show.legend=T) +   # Factor.2
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +  # factor1_names
        ggplot2::scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
        ggplot2::geom_point(data=df_sources,
                            ggplot2::aes(x=x, y=y),
                   size=2,
                   show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                       ggplot2::aes(xmin=xmin,xmax=xmax),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::ylab("") +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    }
  } # end n.effects==2
  if(mix$n.effects==1){
    if(!is.na(source$by_factor)){ # sources by factor, want to color the sources by factor1
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), position=position_jitter(width=.2,height=.1), show.legend=T) +  # Factor.1
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +  # factor1_names
        ggplot2::geom_point(data=df_sources,
                            ggplot2::aes(x=x, y=y, colour=scolour),
                   size=2,
                   show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                                ggplot2::aes(xmin=xmin,xmax=xmax,colour=scolour),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::ylab("") +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    } else { # sources not by factor (make the sources black)
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), position=position_jitter(width=.2,height=.1), show.legend=T) +  # Factor.1
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +     # factor1_names
        ggplot2::geom_point(data=df_sources,
                            ggplot2::aes(x = x, y = y),
                   size=2,
                   show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                                ggplot2::aes(xmin=xmin,xmax=xmax),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::ylab("") +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    }
  } # end n.effects==1
  if(mix$n.effects==0){
    g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.2,height=.1)) +
      ggplot2::geom_point(data=df_sources,
                          ggplot2::aes(x=x,y=y),
                 size=2,
                 show.legend=F) +
      ggplot2::geom_errorbarh(data=df_sources,
                              ggplot2::aes(xmin=xmin,xmax=xmax),
                     size=1,
                     height=0,
                     linetype=source_linetype,
                     show.legend=F) +
      ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::ylab("") +
      ggplot2::xlab(x_label) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
    print(g)
  }
  if(plot_save_pdf==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,".pdf",sep=""))
    cairo_pdf(filename=mypath, width=7, height=7)
    print(g)
    dev.off()
  }
  if(plot_save_png==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,".png",sep=""))
    png(filename=mypath)
    print(g)
    dev.off()
  }
  if(return_obj==TRUE) return(g)
}
