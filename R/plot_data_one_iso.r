####################################################################
# Isospace Plot - one isotope
####################################################################

# July 3

# The function plot_data_one_iso is called when the 'Plot data' button is pressed,
# probably before the model is run, but doesn't matter.
# plot_data_one_iso makes an isospace plot of the consumer and source data **when there is only one isotope selected**

plot_data_one_iso <- function(){
	# Look in the isotope column headers for 'C', 'N', 'S', and 'O'
	# Make the x and y labels for the isospace plot 
	if(length(grep("C",iso_names))==1) x_label <- expression(paste(delta^13, "C (%)",sep=""))
	if(length(grep("N",iso_names))==1) x_label <- expression(paste(delta^15, "N (%)",sep=""))
	if(length(grep("S",iso_names))==1) x_label <- expression(paste(delta^34, "S (%)",sep=""))
	if(length(grep("O",iso_names))==1) x_label <- expression(paste(delta^18, "O (%)",sep=""))
	if(!exists("x_label")) x_label <- iso_names

	y_data <- 0.5
	y <- rep(y_data,N)
	df <- data.frame(x = X_iso, y = y)
	spacing <- 0.1

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
		y_sources <- seq(y_data+0.2,(source_factor_levels*n.sources*spacing)-spacing+y_data+0.2,by=spacing)
		MU_plot <- S_MU[,iso_names] + frac_mu_plot    # add fractionation mean to the source mean values
		SIG_plot <- sqrt(S_SIG[,iso_names]^2 + frac_sig2_plot)  # add fractionation sd to the source sd values
	} else {  # sources_by_factor==FALSE
		source_linetype <- 1:n.sources    # each source gets a different linetype
		source_color <- factor(rep("black",n.sources))  # this doesn't work...solution was to make separate ggplot calls for by_factor and not_by_factor
		index <- 1:n.sources              # "index" gets the row in S_MU of the first instance of each source (since not by factor, only one instance of each source)
		frac_mu_plot <- frac_mu
		frac_sig2_plot <- frac_sig2
		y_sources <- seq(y_data+0.2,(n.sources*spacing)-spacing+y_data+0.2,by=spacing)
		source_factor_levels <- 0.5		# just for correctly spacing the source labels (y in source.labels, line 61)
		MU_plot <- S_MU + frac_mu_plot    # add fractionation mean to the source mean values
		SIG_plot <- sqrt(S_SIG^2 + frac_sig2_plot)  # add fractionation sd to the source sd values
	}

  df_sources <- data.frame(x=MU_plot,
                         xmin = MU_plot - SIG_plot,
                         xmax = MU_plot + SIG_plot,
                         y = y_sources,
                         linetype = source_linetype,
                         scolour = source_color)

  source.labels <- data.frame(
    x = MU_plot[index],    # label sources just left
    #y = MU_plot[index,2] + rep(1.5,n.sources),    # and up from their means
    y = y_sources[index] + spacing*source_factor_levels,
    label = source_names
  )

  dev.new()

  if(n.re==2){
    if(sources_by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point(aes(colour = factor(Factor.1),
                        shape = factor(Factor.2)), position=position_jitter(width=.2,height=.1), show_guide=T) +                    # X[,random_effects[2]]
          scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +          
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
    print(ggplot(data = df,aes(x = x,y = y)) +
          geom_point(aes(colour = factor(Factor.1),
                        shape = factor(Factor.2)), position=position_jitter(width=.2,height=.1), show_guide=T) +                    # X[,random_effects[2]]
          scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +
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
  } # end n.re==2
  if(n.re==1){
    if(sources_by_factor){ # sources by factor, want to color the sources by factor1
      print(ggplot(data = df,aes(x = x,y = y)) +
            geom_point(aes(colour = factor(Factor.1)), position=position_jitter(width=.2,height=.1), show_guide=T) +
            scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +
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
      print(ggplot(data = df,aes(x = x,y = y)) +
            geom_point(aes(colour = factor(Factor.1)), position=position_jitter(width=.2,height=.1), show_guide=T) +
            scale_colour_discrete(breaks = levels(factor(Factor.1)),
                            labels = factor1_names) +            
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
  } # end n.re==1
  if(n.re==0){
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
  if(svalue(plot_save_pdf)==TRUE){
    mypath <- file.path(paste(getwd(),"/",svalue(plot_filename),".pdf",sep=""))
    dev.copy2pdf(file=mypath)
  }
  if(svalue(plot_save_png)==TRUE){
    mypath <- file.path(paste(getwd(),"/",svalue(plot_filename),".png",sep=""))
    dev.copy(png,mypath)
  }
}
