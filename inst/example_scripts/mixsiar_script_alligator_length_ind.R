# Brian Stock
# April 21, 2017
# Script file to run alligator example without GUI
# Recreates Figures 6 and 8 in PeerJ paper:
#   https://peerj.com/articles/5096/#fig-6

# ------------------------------------------------------------------------------
# Alligator example
#   continuous effect of Length
#   random effect of Individual

# choose where to save output
# setwd("path/to/save")

library(MixSIAR)
library(tidyr)
library(ggplot2)
mix.filename <- system.file("extdata", "alligator_consumer.csv", package = "MixSIAR")
source.filename <- system.file("extdata", "alligator_sources_simplemean.csv", package = "MixSIAR")
discr.filename <- system.file("extdata", "alligator_TEF.csv", package = "MixSIAR")

# load mix data (cont effect of Length, random effect of Individual)
mix <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="ID",
                          fac_random=TRUE,
                          fac_nested=NULL,
                          cont_effects="Length")

# load source data
source <- load_source_data(filename=source.filename,
                         source_factors=NULL,
                         conc_dep=FALSE,
                         data_type="means",
                         mix=mix)

# load TEF data
discr <- load_discr_data(filename=discr.filename, mix=mix)

# isospace plot
plot_data(filename="isospace_plot",
        plot_save_pdf=TRUE,
        plot_save_png=FALSE,
        mix, source, discr)

# Define model structure and write JAGS model file
model_filename <- paste0("MixSIAR_model_cont_ind.txt")
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model
jags.mod <- run_model(run="short", mix, source, discr, model_filename, alpha.prior=1)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.mod, mix, source)
graphics.off()

load("/home/brian/Documents/Isotopes/MixSIAR_paper_PeerJ/alligator example/mixsiar_analysis_loo/alligator_allmodels_short.RData")
jags.mod <- jags.mod[[11]]
mix <- mix[[11]]
source <- source[[11]]

# ----------------------------------------------------------------
# Process posteriors to make Figure 6
# Have to 
R2jags::attach.jags(jags.mod)
n.sources <- source$n.sources
source_names <- source$source_names

calc_eps <- function(f){
  n.sources <- length(f)
  gam <- rep(1/n.sources,n.sources)
  phi <- rep(0,n.sources)
  phi[1] <- 1
  sqrt(sum((f-gam)^2))/sqrt(sum((phi-gam)^2))
} 

ce=1
label <- mix$cont_effects[ce]
cont <- mix$CE[[ce]]
ilr.cont <- get(paste("ilr.cont",ce,sep=""))

# get_high <- function(x){return(quantile(x,.975))} # 95% CI 
# get_low <- function(x){return(quantile(x,.025))}
get_high <- function(x){return(quantile(x,.95))} # 90% CI 
get_low <- function(x){return(quantile(x,.05))}    
# get_high <- function(x){return(quantile(x,.75))} # 50% CI 
# get_low <- function(x){return(quantile(x,.25))}    
n.plot = 200
chain.len = dim(p.global)[1]
Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
for(src in 1:n.sources-1){
  for(i in 1:n.plot){
   ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i]
   ilr.low[i,src] <- get_low(ilr.plot[i,src,])
   ilr.median[i,src] <- median(ilr.plot[i,src,])
   ilr.high[i,src] <- get_high(ilr.plot[i,src,])
  }
}

# Transform regression lines from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
for(i in 1:(n.sources-1)){
   e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
   e[,i] <- e[,i]/sum(e[,i])
}
cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.median <- array(data=NA,dim=c(n.plot, n.sources))
cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.low <- array(data=NA,dim=c(n.plot, n.sources))
cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
p.high <- array(data=NA,dim=c(n.plot, n.sources))
eps.low <- rep(NA, n.plot)
eps.med <- rep(NA, n.plot)
eps.high <- rep(NA, n.plot)    
for(i in 1:n.plot){
  for(j in 1:(n.sources-1)){
    cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
    cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
  }
  for(src in 1:n.sources){
    tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    tmp.p.low[i,src] <- prod(cross.low[i,src,]);
    tmp.p.high[i,src] <- prod(cross.high[i,src,]);
  }
  for(src in 1:n.sources){
    p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
    p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
  }
  eps.med[i] <- calc_eps(p.median[i,])
  eps.low[i] <- calc_eps(p.low[i,])
  eps.high[i] <- calc_eps(p.high[i,])
}
colnames(p.median) <- source_names

Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
colnames(df) <- c("source","median","x","low","high")

original <- combine_sources(jags.mod, mix, source, alpha.prior=1, 
  groups=list(Freshwater="Freshwater",Marine="Marine"))
col.ind.marine <- grep("p.ind\\[.*,2]", colnames(original$post))
p.ind.marine <- as.data.frame(original$post[,col.ind.marine])
lengths <- mix$CE_orig[[1]]
colnames(p.ind.marine) <- paste0("ind",1:length(lengths))
df.ind <- p.ind.marine %>% gather(Ind, p)
df.ind$Length <- rep(lengths, each=dim(original$post)[1])
df.ind$Ind <- factor(df.ind$Ind)

cols <- RColorBrewer::brewer.pal(9,"Blues")
png("fig6_diet_length_ind.png", width=7, height=4.5, units='in', res=300)
print(ggplot(df.ind) +
  geom_ribbon(data=df[df$source=="Marine",], mapping=aes(x=x, ymin=low, ymax=high), alpha=0.35, fill=cols[6]) +
  geom_line(data=df[df$source=="Marine",], mapping=aes(x=x, y=median), size=1.5, color=cols[9]) +  
  geom_linerange(mapping=aes(x=Length, y=p, group=Ind), stat = "summary", color=cols[6], size=1, 
               fun.ymin = function(z) {quantile(z,0.05)},
               fun.ymax = function(z) {quantile(z,0.95)}) +  
  geom_point(mapping=aes(x=Length, y=p, group=Ind), stat = "summary", shape = 21, color=cols[9], fill=cols[9], size=3,
               fun.y = function(z) {quantile(z,0.5)}) +  
  scale_y_continuous(limits=c(0,1), expand=c(0.01,0.01)) +
  ylab(expression(italic(p)[marine])) +
  xlab("Total length (cm)") +
  theme_bw() +
  theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
    axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
    axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
    legend.justification=c(0,1), legend.title=ggplot2::element_blank()))
dev.off()

# -----------------------------------------------------------
# Figure 8 (histogram of specialization index)
med.p.ind <- apply(p.ind.marine,2,median)
med.q.ind <- 1-med.p.ind
df.fig8 <- data.frame(pmarine=med.p.ind, pfresh=med.q.ind)
df.fig8$eps <- apply(df.fig8, 1, calc_eps)

png("fig8_eps_hist.png", width=7, height=4.5, units='in', res=500)
print(ggplot(df.fig8, aes(x=eps)) + 
  geom_histogram(bins=20) + 
  xlab(expression(paste("Specialization index (",epsilon[ind],")",sep=""))) +
  ylab("Count") +
  labs(title = "") +
  scale_y_continuous(expand = c(0,0), limits=c(0,125)) +  
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = "black"), axis.title=element_text(size=16), 
    axis.text=element_text(size=14), legend.text=element_text(size=14), legend.position=c(.02,1), 
    legend.justification=c(0,1), legend.title=element_blank()))
dev.off()


