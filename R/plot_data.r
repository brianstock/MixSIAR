# Brian Stock
# January 30, 2014

# Jan 30, 2014
# - Created new function 'plot_data' that calls 'plot_data_two_iso' and 'plot_data_one_iso' to actually make the plots
#     'plot_data_two_iso' and 'plot_data_one_iso' are below, in this file

# June 6, 2013
# - Made 'plot_data' able to be called on any combination of 2 isotopes
#   Added input: 'isotopes', a 2-vector of isotope column indices (for X_iso, S_MU, S_SIG, frac_mu, frac_sig2)
#     Now when the "Plot Data" button is presed, plot_data is called in a loop for each combo of isotopes,
#     e.g. for n.iso=2: plot_data(1,2).  For n.iso=3: plot_data(1,2), plot_data(1,3), and plot_data(2,3).
#   Instead of assuming C is on the x-axis and N on the y-axis, plot_data will look for C, N, S, and O in the
#     column headers and label the axes appropriately.

# July 17, 2014
# Changed the order of input to be consistent with other functions (options first, then previously created objects, i.e. 'mix')

# Function: plot_data
#   makes an isospace plot of the consumer and source data, using each combination of isotopes (if more than 2)
#   if n.iso = 1    calls plot_data_one_iso,
#   if n.iso = 2    calls plot_data_two_iso,
#   if n.iso > 2    calls plot_data_two_iso for each pairwise combo
# Usage: plot_data(mix,source,discr,filename,plot_save_pdf,plot_save_png)
# Input: mix            output from 'load_mix_data'
#        source         output from 'load_source_data'
#        discr          output from 'load_discr_data'
#        filename  string, name of the plot file(s) to save (e.g. "isospace_plot")
#        plot_save_pdf  T/F, save the plot(s) as a pdf?
#        plot_save_png  T/F, save the plot(s) as a png?
# Output: produces an isospace plot(s) in R, saves the plot(s) to a file(s) in the working directory

plot_data <- function(filename,plot_save_pdf,plot_save_png,mix,source,discr){
  if(mix$n.iso==1){
    plot_data_one_iso(mix,source,discr,filename,plot_save_pdf,plot_save_png)
  } else {
    for(iso1 in 1:(mix$n.iso-1)){
      for(iso2 in (iso1+1):mix$n.iso){
        plot_data_two_iso(c(iso1,iso2),mix,source,discr,filename,plot_save_pdf,plot_save_png)
      }
    }
  }
}

# Function: plot_data_two_iso
#   makes an isospace plot for two selected isotopes
# Usage: plot_data_two_iso(isotopes,mix,source,discr,filename,plot_save_pdf,plot_save_png)
# Input: isotopes       two-vector of isotope indices to plot (e.g. c(1,2) or c(2,3))
#        mix            output from 'load_mix_data'
#        source         output from 'load_source_data'
#        discr          output from 'load_discr_data'
#        filename  string, name of the plot file to save (e.g. "isospace_plot")
#        plot_save_pdf  T/F, save the plot as a pdf?
#        plot_save_png  T/F, save the plot as a png?
# Output: produces an isospace plot in R, saves the plot to a file in the working directory

plot_data_two_iso <- function(isotopes,mix,source,discr,filename,plot_save_pdf,plot_save_png){
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
  
  if(source$by_factor){
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
    x = MU_plot[index,1] - rep(1,source$n.sources),    # label sources just left
    y = MU_plot[index,2] + rep(0.75,source$n.sources),    # and up from their means
    label = source$source_names
  )
  .e <- environment()
  dev.new()

  if(mix$n.effects==2){
    # ggplot2 will only make 6 different shapes, so force it to use enough for Factor.2
    shapes <- c(16,17,15,3,7,8,1,6,35,36,37,4,18,14,11,9,13)
    shapes <- shapes[1:mix$FAC[[2]]$levels]  # 1:factor2_levels
    if(source$by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
          geom_point(aes(colour = factor(mix$FAC[[1]]$values), # Factor.1
                        shape = factor(mix$FAC[[2]]$values)), size=2.5, show_guide=T) +   # Factor.2
          scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +  # factor1_names
          scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
          geom_pointrange(data=df_sources,
                          aes(ymin=ymin,ymax=ymax,colour=scolour),
                          size=1,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax,colour=scolour),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          ylab(y_label) +
          xlab(x_label) +
          theme_bw() +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
   } else { # sources not by factor (make the sources black)      
    print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
          geom_point(aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                        shape = factor(mix$FAC[[2]]$values)), size=2.5, show_guide=T) +  # Factor.2
          scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +    # factor1_names
          scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
          geom_pointrange(data=df_sources,
                          aes(ymin=ymin,ymax=ymax),
                          size=1,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          ylab(y_label) +
          xlab(x_label) +
          theme_bw() +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    }

  } # end n.effects==2
  if(mix$n.effects==1){
    if(source$by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y),,environment=.e) +
            geom_point(aes(colour = factor(mix$FAC[[1]]$values)), show_guide=T) +  # Factor.1
            scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),    # Factor.1
                            labels = mix$FAC[[1]]$labels) +  # factor1_names
            geom_pointrange(data=df_sources,
                            aes(ymin=ymin,ymax=ymax,colour=scolour),
                            size=1,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_errorbarh(data=df_sources,
                            aes(xmin=xmin,xmax=xmax,colour=scolour),
                            size=1,
                            height=0,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
            ylab(y_label) +
            xlab(x_label) +
            theme_bw() +
            theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    } else { # sources not by factor (make the sources black)
      print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
            geom_point(aes(colour = factor(mix$FAC[[1]]$values)), show_guide=T) +  # Factor.1
            scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +    # factor1_names
            geom_pointrange(data=df_sources,
                            aes(ymin=ymin,ymax=ymax),
                            size=1,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_errorbarh(data=df_sources,
                            aes(xmin=xmin,xmax=xmax),
                            size=1,
                            height=0,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
            ylab(y_label) +
            xlab(x_label) +
            theme_bw() +
            theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    }
  } # end n.effects==1
  if(mix$n.effects==0){
    print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point() +
          geom_pointrange(data=df_sources,
                          aes(ymin=ymin,ymax=ymax),
                          size=1,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          ylab(y_label) +
          xlab(x_label) +
          theme_bw())
  }
  if(plot_save_pdf==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,"_",isotopes[1],"_",isotopes[2],".pdf",sep=""))
    dev.copy2pdf(file=mypath)
  }
  if(plot_save_png==TRUE){
    mypath <- file.path(paste(getwd(),"/",filename,"_",isotopes[1],"_",isotopes[2],".png",sep=""))
    dev.copy(png,mypath)
  }
} # End plot_data_two_iso function

# Function: plot_data_one_iso
#   makes an isospace plot when there is only one isotope
# Usage: plot_data_one_iso(mix,source,discr,filename,plot_save_pdf,plot_save_png)
# Input: mix            output from 'load_mix_data'
#        source         output from 'load_source_data'
#        discr          output from 'load_discr_data'
#        filename  string, name of the plot file to save (e.g. "isospace_plot")
#        plot_save_pdf  T/F, save the plot as a pdf?
#        plot_save_png  T/F, save the plot as a png?
# Output: produces an isospace plot in R, saves the plot to a file in the working directory

plot_data_one_iso <- function(mix,source,discr,plot_filename,plot_save_pdf,plot_save_png){
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

  if(source$by_factor){
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
  } else {  # source$by_factor==FALSE
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
    if(source$by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
          geom_point(aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                        shape = factor(mix$FAC[[2]]$values)), position=position_jitter(width=.2,height=.1), show_guide=T) +    # Factor.2
          scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +   # factor1_names
          scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
          geom_point(data=df_sources,
                          aes(x=x, y=y,colour=scolour),
                          size=2,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax,colour=scolour),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          scale_y_continuous(breaks = NULL) + 
          ylab("") +
          xlab(x_label) +
          theme_bw() +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
   } else { # sources not by factor (make the sources black)      
    print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
          geom_point(aes(colour = factor(mix$FAC[[1]]$values),   # Factor.1
                        shape = factor(mix$FAC[[2]]$values)), position=position_jitter(width=.2,height=.1), show_guide=T) +   # Factor.2
          scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +  # factor1_names
          scale_shape_manual(values=shapes, labels=mix$FAC[[2]]$labels) +  # factor2_names
          geom_point(data=df_sources,
                          aes(x=x, y=y),
                          size=2,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          scale_y_continuous(breaks = NULL) + 
          ylab("") +
          xlab(x_label) +
          theme_bw() +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    }
  } # end n.effects==2
  if(mix$n.effects==1){
    if(source$by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
            geom_point(aes(colour = factor(mix$FAC[[1]]$values)), position=position_jitter(width=.2,height=.1), show_guide=T) +  # Factor.1
            scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +  # factor1_names
            geom_point(data=df_sources,
                            aes(x=x, y=y, colour=scolour),
                            size=2,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_errorbarh(data=df_sources,
                            aes(xmin=xmin,xmax=xmax,colour=scolour),
                            size=1,
                            height=0,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
            scale_y_continuous(breaks = NULL) + 
            ylab("") +
            xlab(x_label) +
            theme_bw() +
            theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    } else { # sources not by factor (make the sources black)
      print(ggplot(data = df,aes(x = x,y = y),environment=.e) +
            geom_point(aes(colour = factor(mix$FAC[[1]]$values)), position=position_jitter(width=.2,height=.1), show_guide=T) +  # Factor.1
            scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
                            labels = mix$FAC[[1]]$labels) +     # factor1_names
            geom_point(data=df_sources,
                            aes(x = x, y = y),
                            size=2,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_errorbarh(data=df_sources,
                            aes(xmin=xmin,xmax=xmax),
                            size=1,
                            height=0,
                            linetype=source_linetype,
                            show_guide=F) +
            geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
            scale_y_continuous(breaks = NULL) + 
            ylab("") +
            xlab(x_label) +
            theme_bw() +
            theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
    }
  } # end n.effects==1
  if(mix$n.effects==0){
    print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point(position=position_jitter(width=.2,height=.1)) +
          geom_point(data=df_sources,
                          aes(x=x,y=y),
                          size=2,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_errorbarh(data=df_sources,
                          aes(xmin=xmin,xmax=xmax),
                          size=1,
                          height=0,
                          linetype=source_linetype,
                          show_guide=F) +
          geom_text(data=source.labels, aes(x=x,y=y,label=label), show_guide=F) +
          scale_y_continuous(breaks = NULL) + 
          ylab("") +
          xlab(x_label) +
          theme_bw() +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))
  }
  if(plot_save_pdf==TRUE){
    mypath <- file.path(paste(getwd(),"/",plot_filename,".pdf",sep=""))
    dev.copy2pdf(file=mypath)
  }
  if(plot_save_png==TRUE){
    mypath <- file.path(paste(getwd(),"/",plot_filename,".png",sep=""))
    dev.copy(png,mypath)
  }
}

