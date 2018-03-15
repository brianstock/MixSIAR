#' Summary statistics from posterior of MixSIAR model
#'
#' \code{summary_stat} prints and saves summary statistics
#'
#' @param combined list, output from \code{\link{combine_sources}} function
#' @param toprint vector, which parameters to print? Options are: \code{"p"} to print stats for proportions only (all factors),
#'        \code{"global"} to only print overall proportions, \code{"fac1"} to only print factor 1 proportions, \code{"fac2"} to print
#'        factor 2 proportions. Set = \code{"epsilon"} to print multiplicative error term only.
#'        Default = \code{"all"}, prints stats for all model parameters.
#' @param groupby character, group stats by "factor" or "source"? I.e. in wolves example, group proportions by Region 1, Region 2, Region 3
#'        (\code{groupby="factor"}) vs. Deer, Marine Mammals, Salmon (\code{groupby="source"}). Currently only \code{"factor"} is implemented.
#' @param meanSD \code{TRUE/FALSE}, print mean and SD for the parameters?
#' @param quantiles vector, which quantiles to print. Default = \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}.
#' @param savetxt \code{TRUE/FALSE}, save results as .txt file (in working directory)?
#' @param filename character, file name to save results as (\code{.txt} will be appended automatically)
#'  
#' @return NULL
#'  
#' @seealso \code{\link{combine_sources}} and \code{\link{plot_intervals}}

#' @examples 
#' \dontrun{
#' # first run mantis shrimp example

#' # combine 6 sources into 2 groups of interest (hard-shelled vs. soft-bodied)
#' #   'hard' = 'clam' + 'crab' + 'snail'           # group 1 = hard-shelled prey
#' #   'soft' = 'alphworm' + 'brittlestar' + 'fish' # group 2 = soft-bodied prey
#' combined <- combine_sources(jags.1, mix, source, alpha.prior=alpha, 
#'                 groups=list(hard=c("clam","crab","snail"), soft=c("alphworm","brittlestar","fish")))
#'	
#' summary_stat(combined)
#' summary_stat(combined, savetxt=FALSE)
#' summary_stat(combined, meanSD=FALSE)
#' summary_stat(combined, quantiles=c(.05,.5,.95))
#' summary_stat(combined, toprint="fac1")
#' summary_stat(combined, toprint="p")
#' summary_stat(combined, toprint="global")
#' }
summary_stat <- function(combined, toprint="all", groupby="factor", meanSD=TRUE, quantiles=c(.025,.25,.5,.75,.975), savetxt = TRUE, filename="summary_statistics"){
	# Error check for groupby
	if(!(groupby %in% c("factor","source"))){
      stop(paste0("groupby must be 'factor' or 'source'.
See ?summary_stat"))
  	}
	if(groupby=="source") warning("My apologies... groupby='source' not implemented yet.")  	

	# Error check for toprint
	if(!(toprint %in% c("all","p","global","fac1","fac2","epsilon"))){
	  stop(paste0("toprint options are 'all', 'p', 'global', 'fac1', 'fac2', 'epsilon'.
See ?summary_stat"))
	}
	if(combined$mix$n.fe != 0 & toprint=="global"){
	  stop(paste0("toprint='global' but no global proportions in model.
You have fixed effects in your model, so intercept terms are meaningless.
Set toprint='fac1' to print factor 1 proportions,
or toprint='fac2' to print factor 2 proportions."))		
	}

	post <- combined$post
	n.sources <- combined$source.new$n.sources
	source_names <- combined$source.new$source_names

	# easier - just get the proportions from post as new objects
	p.global <- post[,grep("p.global",colnames(post))]
	p.fac1 <- post[,grep("p.fac1",colnames(post))]
	p.fac2 <- post[,grep("p.fac2",colnames(post))]
	p.both <- post[,grep("p.both",colnames(post))]

	# groupby="source" -- need to reorder p.fac1, p.fac2, p.both

	stats <- NULL; eps_stats <- NULL; eps_labels <- NULL;
	sig_labels <- NULL; fac1_labels <- NULL; fac2_labels <- NULL; sig_stats <- NULL;
	getQuant <- function(x) quantile(x, probs=quantiles)
	getMeanSD <- function(x) cbind(round(apply(x,2,mean),3),round(apply(x,2,sd),3))

	if(combined$mix$n.fe == 0 & (toprint %in% c("all","p","global"))){
	  global_quants <- t(round(apply(p.global,2,getQuant),3))
	  if(meanSD){ global_means <- getMeanSD(p.global)} else {global_means <- NULL}
	  stats <- cbind(global_means, global_quants)
	  global_labels <- paste("p.",source_names,sep="")
	  rownames(stats) <- global_labels
	}
	if(combined$mix$n.effects > 0 & combined$mix$n.fe != 2 & (toprint %in% c("all","p","fac1"))){
	  if(groupby=="factor"){
	  	fac1_quants <- t(round(apply(p.fac1,2,getQuant),3))
	  	if(meanSD){ fac1_means <- getMeanSD(p.fac1)} else {fac1_means <- NULL}
	    fac1_stats <- cbind(fac1_means,fac1_quants)
		fac1_labels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[1]]$labels),sep=".")))
		rownames(fac1_stats) <- fac1_labels
		stats <- rbind(stats,fac1_stats)
	  }
	  if(combined$mix$FAC[[1]]$re & (toprint %in% c("all","fac1"))){
	    sig_stats <- cbind(getMeanSD(combined$jags.1$BUGSoutput$sims.list$fac1.sig),t(round(apply(combined$jags.1$BUGSoutput$sims.list$fac1.sig,2,getQuant),3)))
	    sig_labels <- paste(combined$mix$FAC[[1]]$name,".SD",sep="")
	  }
	}
	if(combined$mix$n.re==2 & (toprint %in% c("all","p","fac2"))){
	  if(groupby=="factor"){
	  	fac2_quants <- t(round(apply(p.fac2,2,getQuant),3))
	  	if(meanSD){ fac2_means <- getMeanSD(p.fac2)} else {fac2_means <- NULL}
	    fac2_stats <- cbind(fac2_means,fac2_quants)
		fac2_labels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[2]]$labels),sep=".")))
		rownames(fac2_stats) <- fac2_labels
		stats <- rbind(stats,fac2_stats)
	  }
	  if(combined$mix$FAC[[2]]$re & (toprint %in% c("all","fac2"))){
	    sig_stats <- rbind(sig_stats,cbind(getMeanSD(combined$jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(combined$jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
	    sig_labels <- c(sig_labels,paste(combined$mix$FAC[[2]]$name,".SD",sep=""))
	  }
	}
	if(combined$mix$fere & (toprint %in% c("all","p","fac1","fac2"))){
		fac2_quants <- t(round(apply(p.both,2,getQuant),3))
		if(meanSD){ fac2_means <- getMeanSD(p.both)} else {fac2_means <- NULL}
		fac2_labels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[2]]$labels,combined$mix$FAC[[1]]$labels),sep=".")))
		fac2_stats <- round(cbind(fac2_means,fac2_quants),3)
		rownames(fac2_stats) <- fac2_labels
		stats <- rbind(stats,fac2_stats)
		if(combined$mix$FAC[[2]]$re & (toprint %in% c("all","fac2","fac1"))){
		  sig_stats <- rbind(sig_stats,cbind(getMeanSD(combined$jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(combined$jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
		  sig_labels <- c(sig_labels,paste(combined$mix$FAC[[2]]$name,".SD",sep=""))
		}
	}

	# Add SD stats to the top of the summary
	rownames(sig_stats) <- sig_labels
	stats <- rbind(sig_stats,stats)

	# Add epsilon (multiplicative error term) to stat summary
	# Also plot posterior density
	epsTF <- "resid.prop" %in% names(combined$jags.1$BUGSoutput$sims.list)
	if(epsTF & (toprint %in% c("all","epsilon"))){
	  if(meanSD) eps_stats <- cbind(getMeanSD(combined$jags.1$BUGSoutput$sims.list$resid.prop),t(round(apply(combined$jags.1$BUGSoutput$sims.list$resid.prop,2,getQuant),3)))
	  if(!meanSD) eps_stats <- t(round(apply(combined$jags.1$BUGSoutput$sims.list$resid.prop,2,getQuant),3))
	  eps_labels <- paste0("Epsilon.", 1:combined$mix$n.iso)
	  rownames(eps_stats) <- eps_labels
	  stats <- rbind(eps_stats,stats)
	}

	if(meanSD) colnames(stats) <- c("Mean","SD", paste(round(100*quantiles, 1), "%", sep=""))
	if(!meanSD) colnames(stats) <- paste(round(100*quantiles, 1), "%", sep="")

out1 <- capture.output(stats)
cat("
################################################################
# Summary Statistics
################################################################
",out1,sep="\n")

if(savetxt){  # svalue(summary_save)
  mypath <- file.path(paste(getwd(),"/",filename,".txt",sep=""))  
  cat("
#################################################################
# Summary Statistics
#################################################################
",out1,sep="\n", file=mypath, append=TRUE)
}

}
