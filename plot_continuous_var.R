# Function: plot_continuous_var
# Input: output_options (for plot saving and naming)
# Output: prints posterior plots of diet vs. Cont.1,
#           diet of min(Cont.1), median(Cont.1), and max(Cont.1)

# Oct 13
# Fixed even Cont.1 bug - previously gave error if length(Cont.1) was even,
# because there was no index that matches median(Cont.1) in that case.

plot_continuous_var <- function(output_options){

get_high <- function(x){return(quantile(x,.975))}
get_low <- function(x){return(quantile(x,.025))}
Cont1.plot <- seq(from=round(min(Cont.1),1),to=round(max(Cont.1),1),by=0.1)
ilr.median <- array(NA,dim=c(length(Cont1.plot),n.sources-1))
# ilr.high <- array(NA,dim=c(length(Cont1.plot),n.sources-1))
# ilr.low <- array(NA,dim=c(length(Cont1.plot),n.sources-1))
for(src in 1:n.sources-1){
   ilr.median[,src] <- median(ilr.global[,src]) + median(ilr.cont1[,src])*Cont1.plot
   # ilr.high[,src] <- get_high(ilr.global[,src]) + get_high(ilr.cont1[,src])*Cont1.plot
   # ilr.low[,src] <- get_low(ilr.global[,src]) + get_low(ilr.cont1[,src])*Cont1.plot
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
# cross.high <- array(data=NA,dim=c(N_plot,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
# tmp.p.high <- array(data=NA,dim=c(N_plot,n.sources))              # dummy variable for inverse ILR calculation
# p.high <- array(data=NA,dim=c(N_plot,n.sources))
# cross.low <- array(data=NA,dim=c(N_plot,n.sources,n.sources-1))  # dummy variable for inverse ILR calculation
# tmp.p.low <- array(data=NA,dim=c(N_plot,n.sources))              # dummy variable for inverse ILR calculation
# p.low <- array(data=NA,dim=c(N_plot,n.sources))
for(i in 1:N_plot){
  for(j in 1:(n.sources-1)){
    cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    # cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
    # cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
  }
  for(src in 1:n.sources){
    tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    # tmp.p.high[i,src] <- prod(cross.high[i,src,]);
    # tmp.p.low[i,src] <- prod(cross.low[i,src,]);
  }
  for(src in 1:n.sources){
    p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    # p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
    # p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
  }
}

colnames(p.median) <- source_names
# colnames(p.high) <- source_names
# colnames(p.low) <- source_names

# df <- data.frame(melt(p.median)[,2:3],melt(p.high)[,3],melt(p.low)[,3],rep(Cont1.plot,n.sources))
df <- data.frame(melt(p.median)[,2:3],rep(Cont1.plot,n.sources))
colnames(df) <- c("source","median","x")

medians <- data.frame(Cont.1,apply(p.ind,c(2,3),median))
colnames(medians) <- c("Cont.1",source_names)
medians <- melt(medians,id="Cont.1")

# Plot of Diet vs. Cont.1
# Page 370 in Francis et al (2011)
dev.new()
print(ggplot(data=df,aes(x=x,y=median)) +
         #geom_point(aes(colour=variable), guide=F) + 
         #geom_smooth(method="glm", family="binomial", fill=NA, size=1.1, show_guide=F) +
         #geom_ribbon(aes(ymax=high,ymin=low,group=source,fill=source), alpha=0.3) +        ###### working, but leave off for now
         #geom_smooth(method="loess", fill=NA) + 
         geom_line(aes(x=x, y=median,group=source,colour=source),size=1.5) +
         geom_point(data=medians,aes(x=Cont.1,y=value,colour=variable), guide=F) + 
         #geom_line(data=df.pred, aes(x=Cont.1, y=pred_2, colour="blue")) +
         #geom_line(data=df.pred, aes(x=Cont.1, y=pred_3, colour="purple")) +
         #geom_ribbon(aes(ymin=low, ymax=high, fill=variable, alpha=0.3)) +
         #scale_colour_discrete(name="", breaks = c("X1","X2","X3"), labels = source_names) +
         ylab("Proportion of Diet") +
         #xlab("Secchi depth : Mixed layer depth") +
         xlab(cont_effects[1]) +
         ylim(0,1) +
         #xlim(min(Cont1.plot)-0.1,max(Cont1.plot)+0.1) + 
         theme_bw() +
         theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank()))

# Save the plot to file  
if(output_options[[4]]){ # svalue(plot_post_save_pdf)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_Cont1.pdf",sep=""))  # svalue(plot_post_name)
  dev.copy2pdf(file=mypath)
}
if(output_options[[18]]){ # svalue(plot_post_save_png)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_Cont1.png",sep=""))  # svalue(plot_post_name)
  dev.copy(png,mypath)
}

# Posterior plot for min(Cont1), median(Cont1), and max(Cont1)
# Page 370 in Francis et al (2011)
n.draws <- length(p.global[,1])
# min(Cont1)
dev.new()
df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
min_ind <- which(Cont.1==min(Cont.1))   # find the index of min(Cont1)
for(src in 1:n.sources){
  df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,min_ind,src]) # fill in the p.ind values
  df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
}
my.title <- paste("Diet of min(",cont_effects[1],") individual",sep="")
print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
  geom_density(alpha=.3, aes(y=..scaled..)) +
  theme_bw() +
  xlab("Proportion of Diet") +
  ylab("Scaled Posterior Density") +
  labs(title = my.title) +
  theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))

# Save the plot to file  
if(output_options[[4]]){ # svalue(plot_post_save_pdf)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_min_Cont1.pdf",sep=""))  # svalue(plot_post_name)
  dev.copy2pdf(file=mypath)
}
if(output_options[[18]]){ # svalue(plot_post_save_png)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_min_Cont1.png",sep=""))  # svalue(plot_post_name)
  dev.copy(png,mypath)
}

# median(Cont1)
dev.new()
df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
is.odd <- length(Cont.1) %% 2                # mod 2 division - odd length will return 1, even length returns 0
if(is.odd==1) {med_ind <- which(Cont.1==median(Cont.1))}   # find the index of median(Cont1)
if(is.odd==0) {med_ind <- which(Cont.1==max(sort(Cont.1[which(Cont.1<median(Cont.1))])))}   # If Cont.1 is even, this finds the index of the value just below the median. Here, median(Cont.1) has no corresponding index.
for(src in 1:n.sources){
  df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,med_ind,src]) # fill in the p.ind values
  df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
}
my.title <- paste("Diet of median(",cont_effects[1],") individual",sep="")
print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
  geom_density(alpha=.3, aes(y=..scaled..)) +
  theme_bw() +
  xlab("Proportion of Diet") +
  ylab("Scaled Posterior Density") +
  labs(title = my.title) +
  theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))

# Save the plot to file  
if(output_options[[4]]){ # svalue(plot_post_save_pdf)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_median_Cont1.pdf",sep=""))  # svalue(plot_post_name)
  dev.copy2pdf(file=mypath)
}
if(output_options[[18]]){ # svalue(plot_post_save_png)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_median_Cont1.png",sep=""))  # svalue(plot_post_name)
  dev.copy(png,mypath)
}

# max(Cont1)
dev.new()
df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
max_ind <- which(Cont.1==max(Cont.1))   # find the index of max(Cont1)
for(src in 1:n.sources){
  df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.ind[,max_ind,src])    # fill in the p.ind values
  df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
}
my.title <- paste("Diet of max(",cont_effects[1],") individual",sep="")
print(ggplot(df, aes(x=x, fill=sources, colour=sources)) +
  geom_density(alpha=.3, aes(y=..scaled..)) +
  theme_bw() +
  xlab("Proportion of Diet") +
  ylab("Scaled Posterior Density") +
  labs(title = my.title) +
  theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title=element_blank()))

# Save the plot to file  
if(output_options[[4]]){ # svalue(plot_post_save_pdf)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_max_Cont1.pdf",sep=""))  # svalue(plot_post_name)
  dev.copy2pdf(file=mypath)
}
if(output_options[[18]]){ # svalue(plot_post_save_png)
  mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_max_Cont1.png",sep=""))  # svalue(plot_post_name)
  dev.copy(png,mypath)
}

} #end plot_continuous_var function
