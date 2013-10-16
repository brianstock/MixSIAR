# Function: output_JAGS
# Input: code differs based on n.re (0,1,2) and mcmc.chains (1 or more),
#        but also takes in jags.1, random_effects, n.sources, and source_names 
# Output: prints Gelman, Heidelberger, and Geweke diagnostics;
#         summary statistics for p.global, p.fac1, p.fac2, and factor SDs;
#         trace/XY plots for p.global's and factor SD's;
#         pairs plot for p.global;
#         posterior density plots for p.global's and factor SD's

# May 28
# Changed 'num_factors' to 'n.re' (# random effects) to match the new 'n.ce' variable (# continuous effects)

# July 8
# Adjusted factor2 posterior plots (removed f1 dimension, so same as factor1)
# Removed factor1_levels, factor2_levels, and factor1_names from the function call (global vars, don't need to pass btwn functions)
# Changed plot saving to .png instead of .pdf (so you can flip through them quickly using a photo viewer)

# July 24
# Added option to save plots as pdf and/or png (changed plot_post_save to plot_post_save_pdf and plot_post_save_png)
# Changed 'Open plot in new window' to 'Suppress plot output', since most users will want to see the plots (changed 'plot_post_w' to 'sup_post')

# Oct 4
# Added DIC printout to Summary Statistics file
# Changed 'p.fac1[1,1]' labels to be more informative (e.g. p.Region 1.Salmon)

# Oct 13
# Scaled posterior density plots so you can always see each distribution clearly

#output_options is a LIST:
#1) svalue(summary_save), 
#2) svalue(summary_name), 
#3) svalue(sup_post), 
#4) svalue(plot_post_save_pdf),
#5) svalue(plot_post_name), 
#6) svalue(sup_pairs), 
#7) svalue(plot_pairs_save_pdf), 
#8) svalue(plot_pairs_name),
#9) svalue(sup_xy), 
#10) svalue(plot_xy_save_pdf), 
#11) svalue(plot_xy_name), 
#12) svalue(gelman), 
#13) svalue(heidel), 
#14) svalue(geweke),
#15) svalue(diag_save), 
#16) svalue(diag_name),
#17) include_indiv (include_indiv = TRUE means we have 'ind.sig' and 'p.ind' in the model)
#18) svalue(plot_post_save_png)
#19) svalue(plot_pairs_save_png)
#20) svalue(plot_xy_save_png)

output_JAGS <- function(jags.1, mcmc.chains, n.re, random_effects, n.sources, source_names, output_options){
attach.jags(jags.1)
jags1.mcmc<<-as.mcmc(jags.1)

###########################################################################################
# XY/Trace Plots
###########################################################################################

# XY plots for p.global and factor SD's
if(!output_options[[9]]){  # if 'suppress XY plot' is NOT checked
  # XY plot for p.global
  dev.new()
  print(xyplot(as.mcmc(p.global),strip=strip.custom(factor.levels=source_names)))
  
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
    if(n.re==2) print(xyplot(as.mcmc(cbind(fac1.sig,fac2.sig,ind.sig)),strip=strip.custom(factor.levels=traceplot_labels)))
    if(n.re==1) print(xyplot(as.mcmc(cbind(fac1.sig,ind.sig)),strip=strip.custom(factor.levels=traceplot_labels)))
    if(n.re==0) print(xyplot(as.mcmc(ind.sig),strip=strip.custom(factor.levels=traceplot_labels)))
  } else { # Individual SD is not in the model (no 'ind.sig')
    if(n.re > 0){
      dev.new()
      traceplot_labels <- rep("",length(random_effects))  
      for(i in 1:length(random_effects)) { traceplot_labels[i] <- paste(random_effects[i]," SD",sep="") }
      if(n.re==2) print(xyplot(as.mcmc(cbind(fac1.sig,fac2.sig)),strip=strip.custom(factor.levels=traceplot_labels)))
      if(n.re==1) print(xyplot(as.mcmc(cbind(fac1.sig)),strip=strip.custom(factor.levels=traceplot_labels)))
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
    usr <- par("usr"); on.exit(par(usr))
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
    usr <- par("usr"); on.exit(par(usr))
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
    my.cols <- rev(brewer.pal(n.lines, "RdYlBu"))   # gets some pretty colors
    z <- kde2d(x,y)   # calculates the 2D kernel density that the contour function needs
    contour(z, drawlabels=FALSE, nlevels=n.lines, col=my.cols, add=TRUE)
  }
  pairs(p.global, labels=source_names, diag.panel=panel.hist, lower.panel=panel.cor, upper.panel=panel.contour)
  
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
  # Posterior density plot for p.global
  dev.new()
  n.draws <- length(p.global[,1])   # number of posterior draws
  df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
  for(i in 1:n.sources){
    df$x[seq(1+n.draws*(i-1),i*n.draws)] <- as.matrix(p.global[,i]) # fill in the p.global[i] values
    df$sources[seq(1+n.draws*(i-1),i*n.draws)] <- rep(source_names[i],n.draws)  # fill in the source names
  }
  my.title <- "Overall Population"
  print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
    geom_density(alpha=.3) +
    theme_bw() +
    xlab("Proportion of Diet") +
    ylab("Posterior Density") +
    xlim(0,1) +
    labs(title = my.title) +
    theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))
  
  # Save the plot to file  
  if(output_options[[4]]){ # svalue(plot_post_save_pdf)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.pdf",sep=""))  # svalue(plot_post_name)
    dev.copy2pdf(file=mypath)
  }
  if(output_options[[18]]){ # svalue(plot_post_save_png)
    mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.png",sep=""))  # svalue(plot_post_name)
    dev.copy(png,mypath)
  }
  
  if(n.re >= 1){
    # Posterior density plots for p.fac1's
    for(f1 in 1:factor1_levels){ 
      dev.new()
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      for(src in 1:n.sources){
        df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.fac1[,f1,src]) # fill in the p.fac1[f1] values
        df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
      }
      my.title <- factor1_names[f1]
      print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
#        geom_density(alpha=.3) +
        geom_density(alpha=.3, aes(y=..scaled..)) +
        xlim(0,1) +
        theme_bw() +
        xlab("Proportion of Diet") +
        ylab("Scaled Posterior Density") +
        labs(title = my.title) +
        theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))
        
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",factor1_names[f1],".pdf",sep=""))  # svalue(plot_post_name)
        dev.copy2pdf(file=mypath)
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",factor1_names[f1],".png",sep=""))  # svalue(plot_post_name)
        dev.copy(png,mypath)
      }
    } # end p.fac1 posterior plots
    
    if(n.re==2){
      # Posterior density plots for p.fac2's
        for(f2 in 1:factor2_levels){
          dev.new()
          df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
          for(src in 1:n.sources){
            df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.fac2[,f2,src]) # fill in the p.fac2 values
            df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
          }
          #my.title <- paste(factor1_names[f1],", ", random_effects[2]," ",f2,sep="") # plot title (ex. "Region 1, Pack 3")
          my.title <- factor2_names[f2]
          print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
#             geom_density(alpha=.3) +
            geom_density(alpha=.3, aes(y=..scaled..)) +
            theme_bw() +
            xlim(0,1) +
            xlab("Proportion of Diet") +
            ylab("Scaled Posterior Density") +
            labs(title = my.title) +
            theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))
            
          # Save the plot as a pdf file  
          if(output_options[[4]]){ # svalue(plot_post_save_pdf)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",factor2_names[f2],".pdf",sep="")) #  svalue(plot_post_name)
            dev.copy2pdf(file=mypath)
          }
          if(output_options[[18]]){  # svalue(plot_post_save_png)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",factor2_names[f2],".png",sep="")) #  svalue(plot_post_name)
            dev.copy(png,mypath)
          }
        }# end p.fac2 posterior plots
    } # end if(n.re==2)
  } # end if(n.re>=1)
    
  # Posterior density plot for fac1.sig, fac2.sig, and ind.sig
  if(n.re > 0 || output_options[[17]]){ # only have an SD posterior plot if we have Individual, Factor1, or Factor2 random effects)
    dev.new()
    n.re_ind <- n.re + as.numeric(output_options[[17]]) # this*n.draws will be the length of the plot data frame
    level <- c()
    x <- c()
    if(output_options[[17]]){ # if Individual is in the model, add ind.sig to the SD plot
      level <- c(level,rep("Individual SD",n.draws))
      x <- c(x,ind.sig)
    }
    if(n.re==1){ # if Factor.1 is in the model, add fac1.sig to the SD plot
      level <- c(level,rep(paste(random_effects[1]," SD",sep=""),n.draws))
      x <- c(x,fac1.sig)
    }
    if(n.re==2){ # if Factor.2 is in the model, add fac1.sig and fac2.sig to the SD plot
      level <- c(level,rep(paste(random_effects[1]," SD",sep=""),n.draws), rep(paste(random_effects[2]," SD",sep=""),n.draws))
      x <- c(x,fac1.sig,fac2.sig)
    }
    df2 <- data.frame(level=level, x=x) # create the SD plot data frame
    
    print(ggplot(df2, aes(x=x, fill=level, colour=level)) +
#        geom_density(alpha=.3) +
      geom_density(alpha=.3) +
      theme_bw() +
      xlab(expression(sigma)) +
      ylab("Posterior Density") +
      theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))   # + xlim(0,2)
    
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
s <- summary(jags1.mcmc)
stats <- s$stat
quantiles <- s$quantile

# Re-index summary stats using 'ind' to only print out the statistics we're interested in, in a more sensical order
fac2_ind <- grep("^p.fac2",rownames(stats))     # find the rows for p.fac2's
fac1_ind <- grep("^p.fac1",rownames(stats))     # find the rows for p.fac1's
global_ind <- grep("^p.global",rownames(stats)) # find the rows for p.global's
ind_ind <- grep("^p.ind",rownames(stats))       # find the rows for p.ind's
sig_ind <- grep("sig",rownames(stats))          # find the rows for SDs
ind <- c(sig_ind,global_ind,fac1_ind,fac2_ind,ind_ind)
stats <- stats[ind,]
quantiles <- quantiles[ind,]

# Create new labels for SD, fac1, fac2, and ind terms
# Now instead of displaying 'p.fac1[2,3]', it will display 'p.Region 2.Salmon'
sig_labels <- NULL; ind_labels <- NULL; fac1_labels <- NULL; fac2_labels <- NULL;
global_labels <- rep(NA,n.sources)
for(src in 1:n.sources){
  global_labels[src] <- paste("p.global.",source_names[src],sep="")
}
if(n.re > 0){
  sig_labels <- paste(random_effects[1],".SD",sep="")
  fac1_labels <- rep(NA,factor1_levels*n.sources)
  for(src in 1:n.sources){
    for(f1 in 1:factor1_levels){
      fac1_labels[factor1_levels*(src-1)+f1] <- paste("p.",factor1_names[f1],".",source_names[src],sep="")
    }
  }
}
if(n.re > 1){
  sig_labels <- c(sig_labels, paste(random_effects[2],".SD",sep=""))
  fac2_labels <- rep(NA,factor2_levels*n.sources)
  for(src in 1:n.sources){
    for(f2 in 1:factor2_levels){
      fac2_labels[factor2_levels*(src-1)+f2] <- paste("p.",factor2_names[f2],".",source_names[src],sep="")
    }
  }
}
if(output_options[[17]]){ # include_indiv (if Individual is in the model)
  sig_labels <- c(sig_labels,"Individual.SD")
  ind_labels <- rep(NA,N*n.sources)
  for(src in 1:n.sources){
    for(j in 1:N){
      ind_labels[N*(src-1)+j] <- paste("p.Ind ",j,".",source_names[src],sep="")
    }
  }
} 
rownames(stats) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
rownames(quantiles) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)

################################################################################
# Calulate diagnostics
################################################################################
# Gelman-Rubin diagnostic
if(output_options[[12]]){  # if Gelman is checked
  if(mcmc.chains == 1){
    gelman <- "*** Error: Gelman diagnostic requires more than one chain ***"
  } 
  if(mcmc.chains > 1){    # Gelman diagnostic requires more than one chain
    # Gelman diagnostic, for when the multivariate Gelman fails (matrix not positive definite)
    # Remove the test results for dummy/empty variables
    gelman <- matrix(NA, nrow=nvar(jags1.mcmc), ncol=2)
    for (v in 1:nvar(jags1.mcmc)) {
      gelman[v,] <- gelman.diag(jags1.mcmc[,v])$psrf
    }
    gelman <- gelman[ind,]
    colnames(gelman) <- c("Point est.","Upper C.I.")
    #rownames(gelman) <- varnames(jags1.mcmc)
    rownames(gelman) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
    gelman <- gelman[which(!is.nan(gelman[,1])),] # Remove dummy variables (show up as NA)    
  }
}

# Heidelberger and Welch's diagnostic
# Remove the test results for dummy/empty variables
if(output_options[[13]]){   # if Heidel is checked
  heidel <- heidel.diag(jags1.mcmc)
  w <- which(!is.na(heidel[[1]][,"pvalue"]))  # find all the non-dummy variables    
  heidel.all <- data.frame(matrix(NA,nrow=length(w),ncol=3*mcmc.chains))  # create empty data frame
  colstring <- rep(NA,mcmc.chains*3)  # vector of column names
  for(i in 1:mcmc.chains){
    heidel.tmp <- as.data.frame(heidel[[i]][w,c("stest","pvalue","htest")]) # stest, pvalue, and htest are the relevant statistics - get them
    heidel.all[,(3*i-2):(3*i)] <- heidel.tmp
    colstring[(3*i-2):(3*i)] <- c(paste("stest.chain",i,sep=""), paste("p.val.chain",i,sep=""), paste("hwtest.chain",i,sep="")) # create the appropriate column names
  }
  heidel.all <- heidel.all[ind,]
  rownames(heidel.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
  colnames(heidel.all) <- colstring
  heidel.all <- replace(heidel.all,heidel.all==0,"failed")  # A normal call to 'heidel.diag' prints "failed" and "passed", for some reason they turn to 0's and 1's
  heidel.all <- replace(heidel.all,heidel.all==1,"passed")  # when you access the statistics directly.  Here we turn the 0's and 1's back into "failed" and "passed"
}

# Geweke diagnostic
# Remove the test results for dummy/empty variables
if(output_options[[14]]){ # if Geweke is checked
  geweke <- geweke.diag(jags1.mcmc)
  w <- which(!is.nan(geweke[[1]]$z))  # find all the non-dummy variables
  geweke.all <- data.frame(matrix(NA,nrow=length(w),ncol=mcmc.chains))    # create empty data frame
  colstring <- rep(NA,mcmc.chains)    # vector of column names
  for(i in 1:mcmc.chains){
    geweke.tmp <- as.data.frame(geweke[[i]]$z[w]) # get the relevant geweke statistics
    geweke.all[,i] <- geweke.tmp
    colstring[i] <- c(paste("chain",i,sep=""))  # create the column names "chain1", "chain2", etc.
  }
  geweke.all <- geweke.all[ind,]
  rownames(geweke.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
  colnames(geweke.all) <- colstring
}

################################################################################
# Print diagnostics
################################################################################

if(output_options[[12]]){  # svalue(gelman)
cat("
################################################################################
# Gelman-Rubin Diagnostic
################################################################################

")
print(gelman)

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
  out <- capture.output(gelman)
  cat("
################################################################################
# Gelman-Rubin Diagnostic
################################################################################
      ",out,sep="\n", file=mypath, append=FALSE)
} # end save Gelman
} # end Gelman printout

if(output_options[[13]]){  # svalue(heidel)
cat("
################################################################################
# Heidelberger and Welch Diagnostic
################################################################################

")
print(heidel.all)

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
  out <- capture.output(heidel.all)
  cat("
################################################################################
# Heidelberger and Welch Diagnostic
################################################################################
      ",out,sep="\n", file=mypath, append=output_options[[12]]) # svalue(gelman)
} # end save Heidel
} # end Heidel printout

if(output_options[[14]]){ # svalue(geweke)
cat("
################################################################################
# Geweke Diagnostic
################################################################################

")
print(geweke.all)

if(output_options[[15]]){  # svalue(diag_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
  out <- capture.output(geweke.all)
  cat("
#################################################################
# Geweke Diagnostic
#################################################################
  ",out,sep="\n", file=mypath, append=output_options[[12]]||output_options[[13]]) # svalue(gelman) || svalue(heidel)
} # end Geweke save
} # end Geweke printout

DIC <<- jags.1$BUGSoutput$DIC
cat("
################################################################################
# Summary Statistics
################################################################################

DIC = ",DIC,sep="")
out1 <- capture.output(stats)
out2 <- capture.output(quantiles)
cat("
",out1,out2,sep="\n")

if(output_options[[1]]){  # svalue(summary_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[2]],".txt",sep=""))  # svalue(summary_name)
  cat("
#################################################################
# Summary Statistics
#################################################################
  
DIC = ",DIC,sep="", file=mypath, append=FALSE)
cat("
",out1,out2,sep="\n", file=mypath, append=TRUE)
}

} # end function output_JAGS