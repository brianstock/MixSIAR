#' Plot biotracer data (2-D)
#'
#' \code{plot_data_two_iso} creates a 2-D plot of mix and source tracer data and
#' saves the plot to a file in the working directory
#'
#' An important detail is that \code{plot_data_two_iso} plots the raw mix data
#' and \emph{adds the TDF to the source data}, since this is the polygon that the
#' mixing model uses to determine proportions. The plotted source means are:
#' \deqn{\mu_source + \mu_discr}
#' The source error bars are +/- 1 standard deviation, \emph{calculated as a
#' combination of source and TDF variances:}
#' \deqn{\sqrt{\sigma^2_source + \sigma^2_discr}}
#'
#' \code{plot_data_two_iso} looks for 'C', 'N', 'S', and 'O' in the biotracer column
#' headers and assumes they are stable isotopes, labeling the axes with, e.g.,
#' expression(paste(delta^13, "C (u2030)",sep="")).
#'
#' @param isotopes 2-vector of biotracer indices to plot (e.g. c(1,2) or c(2,3))
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param discr output from \code{\link{load_discr_data}}
#' @param filename name of the plot file(s) to save (e.g. "isospace_plot")
#' @param plot_save_pdf T/F, save the plot(s) as a pdf?
#' @param plot_save_png T/F, save the plot(s) as a png?
#' @param return_obj T/F, whether or not to return ggplot object for further modification, defaults to F
#' @export
#' @seealso \code{\link{plot_data}}
plot_data_two_iso <- function(isotopes,mix,source,discr,filename,plot_save_pdf,plot_save_png,return_obj=FALSE){
  # added only to pass R CMD check
  x <- y <- ymin <- ymax <- scolour <- xmin <- xmax <- label <- NULL

  # Plot the 2 input isotopes (iso1 on x-axis, iso2 on y-axis)
  df <- data.frame(x = mix$data_iso[,isotopes[1]], y = mix$data_iso[,isotopes[2]])
  # Look in the isotope column headers for 'C', 'N', 'S', and 'O'
  # Make the x and y labels for the isospace plot
  if(length(grep("C",mix$iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^13, "C (\u2030)",sep=""))
  if(length(grep("N",mix$iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^15, "N (\u2030)",sep=""))
  if(length(grep("S",mix$iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^34, "S (\u2030)",sep=""))
  if(length(grep("O",mix$iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^18, "O (\u2030)",sep=""))
  if(length(grep("SP",mix$iso_names[isotopes[1]]))==1) y_label <- expression(paste(delta^15, "N-SP (\u2030)",sep=""))
  if(length(grep("C",mix$iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^13, "C (\u2030)",sep=""))
  if(length(grep("N",mix$iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^15, "N (\u2030)",sep=""))
  if(length(grep("S",mix$iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^34, "S (\u2030)",sep=""))
  if(length(grep("O",mix$iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^18, "O (\u2030)",sep=""))
  if(length(grep("SP",mix$iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^15, "N-SP (\u2030)",sep=""))
  if(!exists("x_label")) x_label <- mix$iso_names[isotopes[1]]
  if(!exists("y_label")) y_label <- mix$iso_names[isotopes[2]]

  if(!is.na(source$by_factor)){
    source_linetype <- sort(rep(1:source$n.sources,source$S_factor_levels))    # each source gets a different linetype (assumes source$S_MU is sorted by source and then factor, which it is)
    source_color <- factor(as.numeric(source$S_factor1))  # color sources by factor 1 (ex: region)
    index <- seq(from=1,to=1+(source$n.sources-1)*source$S_factor_levels,by=source$S_factor_levels)  # "index" gets the row in source$S_MU of the first instance of each source (for making the source labels)
    discr_mu_plot <- array(NA,dim=c(length(source$S_MU[,1]),mix$n.iso))     # Since discr$mu is not by factor, it needs to be expanded out by 'source$S_factor_levels' to match the dimensions of source$S_MU.  I.e. if source$n.sources=10, n.iso=2, and source$S_factor_levels=3 (condor data), frac_mu is 10x2 and source$S_MU is 30x2.  This makes frac_mu_plot, a 30x2 matrix.
    discr_sig2_plot <- array(NA,dim=c(length(source$S_MU[,1]),mix$n.iso))   # Same for discr$sig2
    for(i in 1:source$n.sources){
      discr_mu_plot[index[i]:(index[i]+source$S_factor_levels-1),] <- matrix(rep(discr$mu[i,],source$S_factor_levels),nrow=source$S_factor_levels,ncol=mix$n.iso,byrow=T)
      discr_sig2_plot[index[i]:(index[i]+source$S_factor_levels-1),] <- matrix(rep(discr$sig2[i,],source$S_factor_levels),nrow=source$S_factor_levels,ncol=mix$n.iso,byrow=T)
    }
  } else {  # source$by_factor==FALSE
    source_linetype <- 1:source$n.sources    # each source gets a different linetype
    source_color <- factor(rep("black",source$n.sources))  # this doesn't work...solution was to make separate ggplot calls for by_factor and not_by_factor
    index <- 1:source$n.sources              # "index" gets the row in S_MU of the first instance of each source (since not by factor, only one instance of each source)
    discr_mu_plot <- discr$mu
    discr_sig2_plot <- discr$sig2
  }

  MU_plot <- array(NA,dim=c(length(source$S_MU[,1]),2))    # MU_plot will hold the source means adjusted for fractionation/enrichment
  SIG_plot <- array(NA,dim=c(length(source$S_SIG[,1]),2))  # SIG_plot will hold the source sds adjusted for fractionation/enrichment
  #for(src in 1:source$n.sources){
  for(iso in 1:2){
    MU_plot[,iso] <- source$S_MU[,isotopes[iso]] + discr_mu_plot[,isotopes[iso]]    # add fractionation mean to the source mean values
    SIG_plot[,iso] <- sqrt(source$S_SIG[,isotopes[iso]]^2 + discr_sig2_plot[,isotopes[iso]])  # add fractionation sd to the source sd values
  }
  #}

  df_sources <- data.frame(x=MU_plot[,1], y=MU_plot[,2],
                           ymin = MU_plot[,2] - SIG_plot[,2],
                           ymax = MU_plot[,2] + SIG_plot[,2],
                           xmin = MU_plot[,1] - SIG_plot[,1],
                           xmax = MU_plot[,1] + SIG_plot[,1],
                           linetype = source_linetype,
                           scolour = source_color)

  source.labels <- data.frame(
    x = MU_plot[index,1] - rep(0.6,source$n.sources),    # label sources just left
    y = MU_plot[index,2] + rep(0.6,source$n.sources),    # and up from their means
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
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values), # Factor.1
                       shape = factor(mix$FAC[[2]]$values)), size=2.5, show.legend=T) +   # Factor.2
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +  # factor1_names
        ggplot2::scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
        ggplot2::geom_pointrange(data=df_sources,
                        ggplot2::aes(ymin=ymin,ymax=ymax,colour=scolour),
                        size=1,
                        linetype=source_linetype,
                        show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                       ggplot2::aes(xmin=xmin,xmax=xmax,colour=scolour),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::ylab(y_label) +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    } else { # sources not by factor (make the sources black)
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                       shape = factor(mix$FAC[[2]]$values)), size=2.5, show.legend=T) +  # Factor.2
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +    # factor1_names
        ggplot2::scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
        ggplot2::geom_pointrange(data=df_sources,
                        ggplot2::aes(ymin=ymin,ymax=ymax),
                        size=1,
                        linetype=source_linetype,
                        show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                       ggplot2::aes(xmin=xmin,xmax=xmax),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::ylab(y_label) +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    }

  } # end n.effects==2
  if(mix$n.effects==1){
    if(!is.na(source$by_factor)){ # sources by factor, want to color the sources by factor1
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), show.legend=T) +  # Factor.1
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),    # Factor.1
                              labels = mix$FAC[[1]]$labels) +  # factor1_names
        ggplot2::geom_pointrange(data=df_sources,
                        ggplot2::aes(ymin=ymin,ymax=ymax,colour=scolour),
                        size=1,
                        linetype=source_linetype,
                        show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                       ggplot2::aes(xmin=xmin,xmax=xmax,colour=scolour),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::ylab(y_label) +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    } else { # sources not by factor (make the sources black)
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), show.legend=T) +  # Factor.1
        ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                              labels = mix$FAC[[1]]$labels) +    # factor1_names
        ggplot2::geom_pointrange(data=df_sources,
                        ggplot2::aes(ymin=ymin,ymax=ymax),
                        size=1,
                        linetype=source_linetype,
                        show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                       ggplot2::aes(xmin=xmin,xmax=xmax),
                       size=1,
                       height=0,
                       linetype=source_linetype,
                       show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::ylab(y_label) +
        ggplot2::xlab(x_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    }
  } # end n.effects==1
  if(mix$n.effects==0){
    g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_point() +
      ggplot2::geom_pointrange(data=df_sources,
                      ggplot2::aes(ymin=ymin,ymax=ymax),
                      size=1,
                      linetype=source_linetype,
                      show.legend=F) +
      ggplot2::geom_errorbarh(data=df_sources,
                     ggplot2::aes(xmin=xmin,xmax=xmax),
                     size=1,
                     height=0,
                     linetype=source_linetype,
                     show.legend=F) +
      ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
      ggplot2::ylab(y_label) +
      ggplot2::xlab(x_label) +
      ggplot2::theme_bw()
    print(g)
  }
  if(plot_save_pdf==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,"_",isotopes[1],"_",isotopes[2],".pdf",sep=""))
    # dev.copy2pdf(file=mypath)
    cairo_pdf(filename=mypath, width=7, height=7)
    print(g)
    dev.off()
  }
  if(plot_save_png==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,"_",isotopes[1],"_",isotopes[2],".png",sep=""))
    png(filename=mypath)
    print(g)
    dev.off()
  }
  if(return_obj==TRUE) return(g)
} # End plot_data_two_iso function
