####################################################################
# Isospace Plot
####################################################################

# The function plot_data is called when the 'Plot data' button is pressed,
# probably before the model is run, but doesn't matter.
# plot_data makes an isospace plot of the consumer and source data, using each combination of isotopes (if more than 2)

# June 6
# - Made 'plot_data' able to be called on any combination of 2 isotopes
#   Added input: 'isotopes', a 2-vector of isotope column indices (for X_iso, S_MU, S_SIG, frac_mu, frac_sig2)
#     Now when the "Plot Data" button is presed, plot_data is called in a loop for each combo of isotopes,
#     e.g. for n.iso=2: plot_data(1,2).  For n.iso=3: plot_data(1,2), plot_data(1,3), and plot_data(2,3).
#   Instead of assuming C is on the x-axis and N on the y-axis, plot_data will look for C, N, S, and O in the
#     column headers and label the axes appropriately.

plot_data <- function(isotopes){
  # Plot the 2 input isotopes (iso1 on x-axis, iso2 on y-axis)
  df <- data.frame(x = X_iso[,isotopes[1]], y = X_iso[,isotopes[2]])
  # Look in the isotope column headers for 'C', 'N', 'S', and 'O'
  # Make the x and y labels for the isospace plot 
  if(length(grep("C",iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^13, "C (%)",sep=""))
  if(length(grep("N",iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^15, "N (%)",sep=""))
  if(length(grep("S",iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^34, "S (%)",sep=""))
  if(length(grep("O",iso_names[isotopes[1]]))==1) x_label <- expression(paste(delta^18, "O (%)",sep=""))
  if(length(grep("C",iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^13, "C (%)",sep=""))
  if(length(grep("N",iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^15, "N (%)",sep=""))
  if(length(grep("S",iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^34, "S (%)",sep=""))
  if(length(grep("O",iso_names[isotopes[2]]))==1) y_label <- expression(paste(delta^18, "O (%)",sep=""))
  if(!exists("x_label")) x_label <- iso_names[isotopes[1]]
  if(!exists("y_label")) y_label <- iso_names[isotopes[2]]
  
  if(sources_by_factor){
    source_linetype <- sort(rep(1:n.sources,source_factor_levels))    # each source gets a different linetype (assumes S_MU is sorted by source and then factor, which it is)
    source_color <- factor(as.numeric(S_factor1))  # color sources by factor 1 (ex: region)
    index <- seq(from=1,to=1+(n.sources-1)*source_factor_levels,by=source_factor_levels)  # "index" gets the row in S_MU of the first instance of each source (for making the source labels)
    frac_mu_plot <- array(NA,dim=c(length(S_MU[,1]),n.iso))     # Since frac_mu is not by factor, it needs to be expanded out by 'source_factor_levels' to match the dimensions of S_MU.  I.e. if n.sources=10, n.iso=2, and source_factor_levels=3 (condor data), frac_mu is 10x2 and S_MU is 30x2.  This makes frac_mu_plot, a 30x2 matrix.
    frac_sig2_plot <- array(NA,dim=c(length(S_MU[,1]),n.iso))   # Same for frac_sig2
    for(i in 1:n.sources){
      frac_mu_plot[index[i]:(index[i]+source_factor_levels-1),] <- matrix(rep(frac_mu[i,],source_factor_levels),nrow=source_factor_levels,ncol=n.iso,byrow=T)
      frac_sig2_plot[index[i]:(index[i]+source_factor_levels-1),] <- matrix(rep(frac_sig2[i,],source_factor_levels),nrow=source_factor_levels,ncol=n.iso,byrow=T)
    }
  } else {  # sources_by_factor==FALSE
    source_linetype <- 1:n.sources    # each source gets a different linetype
    source_color <- factor(rep("black",n.sources))  # this doesn't work...solution was to make separate ggplot calls for by_factor and not_by_factor
    index <- 1:n.sources              # "index" gets the row in S_MU of the first instance of each source (since not by factor, only one instance of each source)
    frac_mu_plot <- frac_mu
    frac_sig2_plot <- frac_sig2
  }

  MU_plot <- array(NA,dim=c(length(S_MU[,1]),2))    # MU_plot will hold the source means adjusted for fractionation/enrichment
  SIG_plot <- array(NA,dim=c(length(S_SIG[,1]),2))  # SIG_plot will hold the source sds adjusted for fractionation/enrichment
  #for(src in 1:n.sources){
    for(iso in 1:2){
      MU_plot[,iso] <- S_MU[,isotopes[iso]] + frac_mu_plot[,isotopes[iso]]    # add fractionation mean to the source mean values
      SIG_plot[,iso] <- sqrt(S_SIG[,isotopes[iso]]^2 + frac_sig2_plot[,isotopes[iso]])  # add fractionation sd to the source sd values
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
    x = MU_plot[index,1] - rep(1,n.sources),    # label sources just left
    y = MU_plot[index,2] + rep(0.75,n.sources),    # and up from their means
    label = source_names
  )

  dev.new()

  if(n.re==2){
    # ggplot2 will only make 6 different shapes, so force it to use enough for Factor.2
    shapes <- c(16,17,15,3,7,8,1,6,35,36,37,4,18,14,11,9,13)
    shapes <- shapes[1:factor2_levels]
    if(sources_by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point(aes(colour = factor(Factor.1),
                        shape = factor(Factor.2)), size=2.5, show_guide=T) +                    # X[,random_effects[2]]
          scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +
          scale_shape_manual(values=shapes, labels=factor2_names) + 
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
    print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point(aes(colour = factor(Factor.1),
                        shape = factor(Factor.2)), size=2.5, show_guide=T) +                    # X[,random_effects[2]]
          scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +
          scale_shape_manual(values=shapes, labels=factor2_names) + 
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

  } # end n.re==2
  if(n.re==1){
    if(sources_by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y)) +
            geom_point(aes(colour = factor(Factor.1)), show_guide=T) +
            scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +
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
      print(ggplot(data = df,aes(x = x,y = y)) +
            geom_point(aes(colour = factor(Factor.1)), show_guide=T) +
            scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +            
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
  } # end n.re==1
  if(n.re==0){
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
  if(svalue(plot_save_pdf)==TRUE){
    mypath <- file.path(paste(getwd(),"/",svalue(plot_filename),"_",isotopes[1],"_",isotopes[2],".pdf",sep=""))
    dev.copy2pdf(file=mypath)
  }
  if(svalue(plot_save_png)==TRUE){
    mypath <- file.path(paste(getwd(),"/",svalue(plot_filename),"_",isotopes[1],"_",isotopes[2],".png",sep=""))
    dev.copy(png,mypath)
  }
} # End plot_data function

