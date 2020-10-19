#' Output summary statistics
#'
#' \code{output_stats} returns summary statistics from a fit MixSIAR model
#'
#' @param jags.1 rjags model object, output from \code{\link{run_model}} function
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param output_options list containing options for plots and saving:
#'   \itemize{
#'    \item \code{summary_save}: Save the summary statistics as a txt file? Default = \code{TRUE}
#'    \item \code{summary_name}: Summary statistics file name (.txt will be appended). Default = \code{"summary_statistics"}
#'    \item \code{sup_post}: Suppress posterior density plot output in R? Default = \code{FALSE}
#'    \item \code{plot_post_save_pdf}: Save posterior density plots as pdfs? Default = \code{TRUE}
#'    \item \code{plot_post_name}: Posterior plot file name(s) (.pdf/.png will be appended) Default = \code{"posterior_density"}
#'    \item \code{sup_pairs}: Suppress pairs plot output in R? Default = \code{FALSE}
#'    \item \code{plot_pairs_save_pdf}: Save pairs plot as pdf? Default = \code{TRUE}
#'    \item \code{plot_pairs_name}: Pairs plot file name (.pdf/.png will be appended) Default = \code{"pairs_plot"}
#'    \item \code{sup_xy}: Suppress xy/trace plot output in R? Default = \code{TRUE}
#'    \item \code{plot_xy_save_pdf}: Save xy/trace plot as pdf? Default = \code{FALSE}
#'    \item \code{plot_xy_name}: XY/trace plot file name (.pdf/.png will be appended) Default = \code{"xy_plot"}
#'    \item \code{gelman}: Calculate Gelman-Rubin diagnostic test? Default = \code{TRUE}
#'    \item \code{heidel}: Calculate Heidelberg-Welch diagnostic test? Default = \code{FALSE}
#'    \item \code{geweke}: Calculate Geweke diagnostic test? Default = \code{TRUE}
#'    \item \code{diag_save}: Save the diagnostics as a .txt file? Default = \code{TRUE}
#'    \item \code{diag_name}: Diagnostics file name (.txt will be appended) Default = \code{"diagnostics"}
#'    \item \code{indiv_effect}: artifact, set to FALSE 
#'    \item \code{plot_post_save_png}: Save posterior density plots as pngs? Default = \code{FALSE}
#'    \item \code{plot_pairs_save_png}: Save pairs plot as png? Default = \code{FALSE}
#'    \item \code{plot_xy_save_png}: Save xy/trace plot as png? Default = \code{FALSE}
#'    \item \code{diag_save_ggmcmc}: Save ggmcmc diagnostics as pdf? Default = \code{TRUE}
#'    \item \code{return_obj} Return ggplot objects for later modification? Default = \code{FALSE}
#'   }
#'   
#' @return data frame of summary statistics (if \code{return_obj = TRUE})
#' 
#' @export
#'   
output_stats <- function(jags.1, mix, source, output_options=list(
                                                  summary_save = TRUE,                 # Save the summary statistics as a txt file?
                                                  summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
                                                  sup_post = FALSE,                       # Suppress posterior density plot output in R?
                                                  plot_post_save_pdf = TRUE,              # Save posterior density plots as pdfs?
                                                  plot_post_name = "posterior_density",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                                  sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                                                  plot_pairs_save_pdf = TRUE,             # Save pairs plot as pdf?
                                                  plot_pairs_name = "pairs_plot",         # If yes, specify the base file name (.pdf/.png will be appended later)
                                                  sup_xy = TRUE,                         # Suppress xy/trace plot output in R?
                                                  plot_xy_save_pdf = FALSE,                # Save xy/trace plot as pdf?
                                                  plot_xy_name = "xy_plot",               # If yes, specify the base file name (.pdf/.png will be appended later)
                                                  gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                                                  heidel = FALSE,                          # Calculate Heidelberg-Welch diagnostic test?
                                                  geweke = TRUE,                          # Calculate Geweke diagnostic test?
                                                  diag_save = TRUE,                       # Save the diagnostics as a txt file?
                                                  diag_name = "diagnostics",              # If yes, specify the base file name (.txt will be appended later)
                                                  indiv_effect = FALSE,                   # Is Individual a random effect in the model? (already specified)
                                                  plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                                                  plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                                                  plot_xy_save_png = FALSE,
                                                  diag_save_ggmcmc = TRUE,
                                                  return_obj = FALSE)){
mcmc.chains <- jags.1$BUGSoutput$n.chains
N <- mix$N
n.re <- mix$n.re
n.effects <- mix$n.effects
if(n.re==1){
  random_effects <- ifelse(mix$FAC[[1]]$re,mix$FAC[[1]]$name,mix$FAC[[2]]$name)
}
if(n.re==2){
  random_effects <- mix$factors
}
n.sources <- source$n.sources
source_names <- source$source_names
# p.global <- ilr.global <- ilr.fac1 <- ilr.fac2 <- fac1.sig <- fac2.sig <- NULL
# ind.sig <- ..scaled.. <- p.fac1 <- p.fac2 <- p.ind <- sources <- NULL
# R2jags::attach.jags(jags.1)
jags1.mcmc <- coda::as.mcmc(jags.1)
n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])

# Post-processing for 2 FE or 1FE + 1RE
#   calculate p.both = ilr.global + ilr.fac1 + ilr.fac2
if(mix$fere){
  fac2_lookup <- list()
  for(f1 in 1:mix$FAC[[1]]$levels){
    fac2_lookup[[f1]] <- unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values==f1)])
  }
  ilr.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources-1))
  p.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources))
  cross.both <- array(data=NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels,n.sources,n.sources-1))
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }
  for(i in 1:n.draws){
    for(f1 in 1:mix$FAC[[1]]$levels) {
      for(f2 in fac2_lookup[[f1]]){
        for(src in 1:(n.sources-1)) {
          ilr.both[i,f1,f2,src] <- jags.1$BUGSoutput$sims.list$ilr.global[i,src] + jags.1$BUGSoutput$sims.list$ilr.fac1[i,f1,src] + jags.1$BUGSoutput$sims.list$ilr.fac2[i,f2,src];
          cross.both[i,f1,f2,,src] <- (e[,src]^ilr.both[i,f1,f2,src])/sum(e[,src]^ilr.both[i,f1,f2,src]);
          # ilr.both[,f1,f2,src] <- ilr.global[,src] + ilr.fac1[,f1,src] + ilr.fac2[,f2,src];
        }
        for(src in 1:n.sources) {
          p.both[i,f1,f2,src] <- prod(cross.both[i,f1,f2,src,]);
        }
        p.both[i,f1,f2,] <- p.both[i,f1,f2,]/sum(p.both[i,f1,f2,]);
      } # f2
    } # f1
  }
} # end fere

# Calculate the summary statistics for the variables we're interested in (p.global's and factor SD's, maybe p.ind's)
# We print them out later, at the very bottom
sig_labels <- NULL; ind_labels <- NULL; fac1_labels <- NULL; fac2_labels <- NULL; sig_stats <- NULL;
getQuant <- function(x) quantile(x,probs=c(.025,.05,.25,.5,.75,.95,.975))
getMeanSD <- function(x) cbind(round(apply(x,2,mean),3),round(apply(x,2,sd),3))

stats <- NULL
sig_stats <- NULL
sig_labels <- NULL
eps_stats <- NULL
eps_labels <- NULL
# print(mix)
# print(mix$n.fe)
if(mix$n.fe == 0){
  global_quants <- t(round(apply(jags.1$BUGSoutput$sims.list$p.global,2,getQuant),3))
  global_means <- getMeanSD(jags.1$BUGSoutput$sims.list$p.global)
  stats <- cbind(global_means, global_quants)
  global_labels <- rep(NA,n.sources)
  for(src in 1:n.sources){
    global_labels[src] <- paste("p.global.",source_names[src],sep="")
  }
  rownames(stats) <- global_labels
}
if(n.effects > 0 & mix$n.fe != 2){
  fac1_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
  fac1_quants <- t(apply(fac1_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  fac1_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),sd),3))$value)
  fac1_stats <- cbind(fac1_means,fac1_quants)
  fac1_labels <- rep(NA,mix$FAC[[1]]$levels*n.sources)
  for(src in 1:n.sources){
    for(f1 in 1:mix$FAC[[1]]$levels){
      fac1_labels[mix$FAC[[1]]$levels*(src-1)+f1] <- paste("p.",mix$FAC[[1]]$labels[f1],".",source_names[src],sep="")
    }
  }
  rownames(fac1_stats) <- fac1_labels
  stats <- rbind(stats,fac1_stats)
  if(mix$FAC[[1]]$re){
    sig_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac1.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac1.sig,2,getQuant),3)))
    sig_labels <- paste(mix$FAC[[1]]$name,".SD",sep="")
  }
}
if(n.re==2){
  fac2_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
  fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  fac2_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),sd),3))$value)
  fac2_stats <- cbind(fac2_means,fac2_quants)
  fac2_labels <- rep(NA,mix$FAC[[2]]$levels*n.sources)
  for(src in 1:n.sources){
    for(f2 in 1:mix$FAC[[2]]$levels){
      fac2_labels[mix$FAC[[2]]$levels*(src-1)+f2] <- paste("p.",mix$FAC[[2]]$labels[f2],".",source_names[src],sep="")
    }
  }
  rownames(fac2_stats) <- fac2_labels
  stats <- rbind(stats,fac2_stats)
  if(mix$FAC[[2]]$re){
    sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
    sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
  }
}
if(mix$fere){
  fac2_quants <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=7)
  fac2_means <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=2)
  fac2_labels <- rep(NA,n.sources*length(unlist(fac2_lookup)))
  i <- 1
  for(f1 in 1:mix$FAC[[1]]$levels) {
    for(f2 in fac2_lookup[[f1]]){
      for(src in 1:n.sources){
        fac2_quants[i,] <- getQuant(p.both[,f1,f2,src])
        fac2_means[i,] <- c(mean(p.both[,f1,f2,src]),sd(p.both[,f1,f2,src]))
        fac2_labels[i] <- paste("p",mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],source_names[src],sep=".")
        i <- i+1
      }
    }
  }
  # fac2_quants <- as.matrix(cast(melt(round(apply(p.both,c(2,3,4),getQuant,na.rm=TRUE),3)),X4+X3+X2~X1)[,-c(1,2)])
  # fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  # fac2_means <- cbind(melt(round(apply(p.fac2,c(2,3),mean),3))$value, melt(round(apply(p.fac2,c(2,3),sd),3))$value)
  fac2_stats <- round(cbind(fac2_means,fac2_quants),3)
  rownames(fac2_stats) <- fac2_labels
  stats <- rbind(stats,fac2_stats)
  if(mix$FAC[[2]]$re){
    sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
    sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
  }
}

if(output_options[[17]]){ # include_indiv (if Individual is in the model)
  ind_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(p.ind,c(2,3),getQuant),3)),X3+X2~X1)[,-c(1,2)])
  ind_quants <- t(apply(ind_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
  ind_means <- cbind(reshape2::melt(round(apply(p.ind,c(2,3),mean),3))$value, reshape2::melt(round(apply(p.ind,c(2,3),sd),3))$value)
  ind_stats <- cbind(ind_means,ind_quants)
  ind_labels <- rep(NA,N*n.sources)
  for(src in 1:n.sources){
    for(j in 1:N){
      ind_labels[N*(src-1)+j] <- paste("p.Ind ",j,".",source_names[src],sep="")
    }
  }
  sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$ind.sig),t(round(apply(jags.1$BUGSoutput$sims.list$ind.sig,2,getQuant),3))))
  sig_labels <- c(sig_labels,"Individual.SD")
  rownames(ind_stats) <- ind_labels
  stats <- rbind(stats, ind_stats)
}

# Add SD stats to the top of the summary
rownames(sig_stats) <- sig_labels
stats <- rbind(sig_stats,stats)

# Add epsilon (multiplicative error term) to stat summary
epsTF <- "resid.prop" %in% names(jags.1$BUGSoutput$sims.list)
if(epsTF){
  eps_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$resid.prop),t(round(apply(jags.1$BUGSoutput$sims.list$resid.prop,2,getQuant),3)))
  eps_labels <- paste0("Epsilon.", 1:mix$n.iso)
  rownames(eps_stats) <- eps_labels
  stats <- rbind(eps_stats,stats)
}
colnames(stats) <- c("Mean","SD","2.5%","5%","25%","50%","75%","95%","97.5%")

# Pack 1 stats only
#stats[grep("Pack 1",rownames(stats)),]

# Region stats only
#stats[grep("Region",rownames(stats)),]

# Region stats, by Region
# byVec <- function(x){ind <- NULL; for(i in 1:length(x)){ ind <- c(ind,grep(x[i],rownames(stats)))}; return(ind)}
# stats[byVec(mix$RE[[1]]$labels),]

# All means
# stats[,"Mean"]

# Region means only
# stats[byVec(mix$RE[[1]]$labels),"Mean"]

# ------------------------------------------------------
# print summary stats
DIC <- jags.1$BUGSoutput$DIC
cat("
################################################################################
# Summary Statistics
################################################################################

DIC = ",DIC,sep="")
out1 <- capture.output(stats)
cat("
",out1,sep="\n")

if(output_options[[1]]){  # svalue(summary_save)
  mypath <- file.path(paste(getwd(),"/",output_options[[2]],".txt",sep=""))  # svalue(summary_name)
  cat("
#################################################################
# Summary Statistics
#################################################################

DIC = ",DIC,sep="", file=mypath, append=FALSE)
cat("
",out1,sep="\n", file=mypath, append=TRUE)
}

if(!is.null(output_options$return_obj)) if(output_options$return_obj) return(stats) else return(NULL)
} # end function output_stats
