#' Plot posterior uncertainty intervals from a MixSIAR model
#'
#' \code{plot_intervals} plots the posterior interval estimates (quantile-based) from the MCMC draws in a MixSIAR model.
#' Calls \href{http://mc-stan.org/bayesplot/reference/MCMC-intervals.html}{bayesplot::mcmc_intervals}.
#'
#' @param combined list, output from \code{\link{combine_sources}} function
#' @param toplot vector, which parameters to plot? Options are similar to \code{\link{summary_stat}}: 
#' \itemize{
#'  \item \code{"p"}: plots all proportions (default)
#'  \item \code{"global"}: plots overall proportions
#'  \item \code{"fac1"}: plots factor 1 proportions
#'  \item \code{"fac2"}: plots factor 2 proportions
#'  \item \code{"epsilon"}: plots multiplicative error terms
#'  \item \code{"sd"}: plots random effect SD terms
#' }
#' @param levels vector if \code{toplot="fac1"} or \code{toplot="fac2"}, which level(s) to plot? Plots all levels if \code{level=NULL} (default). Specify levels as a vector, e.g. in wolves ex, \code{levels=1} to plot Region 1, \code{levels=c(1,2)} to plot Regions 1 and 2.
#' @param groupby character, group by "factor" or "source"? I.e. in wolves example, group proportions by Region 1, Region 2, Region 3
#'        (\code{groupby="factor"}) vs. Deer, Marine Mammals, Salmon (\code{groupby="source"}). Currently only \code{"factor"} is implemented.
#' @param savepdf \code{TRUE/FALSE}, save plot as .pdf file (in working directory)?
#' @param filename character, file name to save results as (\code{.pdf} will be appended automatically)
#' @param ... additional arguments to pass to \href{http://mc-stan.org/bayesplot/reference/MCMC-intervals.html}{bayesplot::mcmc_intervals}. For example:
#'   \itemize{
#'    \item \code{prob}: sets inner (thick) interval (default = 50\%)
#'    \item \code{prob_outer}: sets outer (thin) interval (default = 90\%)
#'    \item \code{point_est}: what point estimate to use (dot), default = \code{"median"}, can also use \code{"mean"} or \code{"none"}
#'   } 
#' 
#' @return NULL
#' @export
#' @seealso \code{\link{combine_sources}} and \code{\link{summary_stat}}
#' 
#' @examples 
#' \dontrun{
#' # 1. run mantis shrimp example
#' original <- combine_sources(jags.1, mix, source, alpha, 
#'                 groups=list(alphworm="alphworm",brittlestar="brittlestar",clam="clam",
#'	                           crab="crab",fish="fish",snail="snail"))

#' # 2. combine 6 sources into 2 groups of interest (hard-shelled vs. soft-bodied)
#' #   'hard' = 'clam' + 'crab' + 'snail'           # group 1 = hard-shelled prey
#' #   'soft' = 'alphworm' + 'brittlestar' + 'fish' # group 2 = soft-bodied prey
#' combined <- combine_sources(jags.1, mix, source, alpha.prior=alpha, 
#'                 groups=list(hard=c("clam","crab","snail"), soft=c("alphworm","brittlestar","fish")))
#'	
#' plot_intervals(combined,toplot="fac1")
#' plot_intervals(original,toplot="fac1")
#' plot_intervals(combined,toplot="fac1",levels=1)
#' plot_intervals(combined,toplot="fac1",levels=2)
#'
#' }
plot_intervals <- function(combined, toplot="p", levels=NULL, groupby="factor", savepdf=FALSE, filename="post_intervals", ...){
	# Error check for groupby
	if(!(groupby %in% c("factor","source"))){
      stop(paste0("groupby must be 'factor' or 'source'.
See ?plot_intervals"))
  	}
	if(groupby=="source") warning("My apologies... groupby='source' not implemented yet.")  	

	# Error check for toplot
	if(!(toplot %in% c("sd","p","global","fac1","fac2","epsilon"))){
	  stop(paste0("toplot options are 'p', 'global', 'fac1', 'fac2', 'epsilon', 'sd'.
See ?plot_intervals"))
	}
	if(combined$mix$n.fe != 0 & toplot=="global"){
	  stop(paste0("toplot='global' but no global proportions in model.
You have fixed effects in your model, so intercept terms are meaningless.
Set toplot='p' to plot all proportions,
toplot='fac1' to plot factor 1 proportions,
or toplot='fac2' to plot factor 2 proportions."))		
	}
	if(combined$mix$n.re != 0 & toplot=="sd"){
	  stop(paste0("toplot='sd' but no random effects in model."))		
	}	

	n.sources <- combined$source.new$n.sources
	source_names <- combined$source.new$source_names

	if(toplot=="p"){
		print(bayesplot::mcmc_intervals(combined$post, regex_pars="p\\.", ...) + 
		  ggplot2::xlab("Proportion") + 
		  # ggplot2::scale_y_discrete(name="Source",labels=combined$source.new$source_names) + 
		  ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)))
	}
	if(toplot=="global"){
		labs <- grep("global", colnames(combined$post), value=TRUE)
		print(bayesplot::mcmc_intervals(combined$post, regex_pars="global") + 
		  ggplot2::xlab("Proportion") + 
		  ggplot2::scale_y_discrete(name="Source",labels=combined$source.new$source_names) + 
		  ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
		  ggplot2::ggtitle("Global proportions"))
	}
	if(toplot=="fac1" & combined$mix$n.effects > 0){
		if(is.null(levels)) levels <- 1:combined$mix$FAC[[1]]$levels
		if(length(levels)==1){
			yname <- "Source"
			ylabels <- source_names
			titlelab <- paste0(combined$mix$FAC[[1]]$labels[levels], " proportions")
		} else {
			yname <- ""
			ylabels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[1]]$labels[levels]),sep=".")))
			titlelab <- paste0("Proportions by ",combined$mix$FAC[[1]]$name)
		}
		print(bayesplot::mcmc_intervals(combined$post, regex_pars=paste0("p.fac1\\[",levels, collapse="|")) + 
		  ggplot2::xlab("Proportion") + 
		  ggplot2::scale_y_discrete(name=yname, labels=ylabels) + 
		  ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
		  ggplot2::ggtitle(titlelab))
	}
	if(toplot=="fac2" & combined$mix$n.effects > 1 & !combined$mix$fere){
		if(is.null(levels)) levels <- 1:combined$mix$FAC[[2]]$levels
		if(length(levels)==1){
			yname <- "Source"
			ylabels <- source_names
			titlelab <- paste0(combined$mix$FAC[[2]]$labels[levels], " proportions")
		} else {
			yname <- ""
			ylabels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[2]]$labels[levels]),sep=".")))
			titlelab <- paste0("Proportions by ",combined$mix$FAC[[2]]$name)
		}
		print(bayesplot::mcmc_intervals(combined$post, regex_pars=paste0("p.fac2\\[",levels, collapse="|")) + 
		  ggplot2::xlab("Proportion") + 
		  ggplot2::scale_y_discrete(name=yname, labels=ylabels) + 
		  ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
		  ggplot2::ggtitle(titlelab))
	}	
	if((toplot %in% c("fac1","fac2")) & combined$mix$fere){
		yname <- ""
		ylabels <- paste0("p.",do.call(paste, c(expand.grid(source_names,combined$mix$FAC[[2]]$labels,combined$mix$FAC[[1]]$labels),sep=".")))
		titlelab <- paste0("Proportions by ",combined$mix$FAC[[1]]$name," and ",combined$mix$FAC[[2]]$name)
		print(bayesplot::mcmc_intervals(combined$post, regex_pars="p.both") + 
		  ggplot2::xlab("Proportion") + 
		  ggplot2::scale_y_discrete(name=yname, labels=ylabels) + 
		  ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
		  ggplot2::ggtitle(titlelab))		
	}
	if(toplot=="epsilon"){
		ylabels <- paste0("Epsilon ", 1:combined$mix$n.iso)
		xmax <- max(combined$post[,grep("resid.prop", colnames(combined$post))])
		print(bayesplot::mcmc_intervals(combined$post, regex_pars="resid.prop") + 
		  ggplot2::xlab("") + 
		  ggplot2::scale_y_discrete(name="", labels=ylabels) + 
		  ggplot2::scale_x_continuous(limits=c(0,xmax), expand=c(0.01,0.01)) +
		  ggplot2::ggtitle("Epsilon error terms"))		
	}	
	if(toplot=="sd"){
		xmax <- max(combined$post[,grep("\\.sig", colnames(combined$post))])
		print(bayesplot::mcmc_intervals(combined$post, regex_pars="\\.sig") + 
		  ggplot2::xlab("") + 
		  ggplot2::scale_y_discrete(name="", labels=ylabels) + 
		  ggplot2::scale_x_continuous(limits=c(0,xmax), expand=c(0.01,0.01)) +
		  ggplot2::ggtitle("Random effect SD terms"))		
	}
}
