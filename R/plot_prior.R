#' Plot prior
#'
#' \code{plot_prior} plots your prior on the global diet proportions (p.global)
#' and the uninformative prior side-by-side. Your prior is in red, and the
#' "uninformative"/generalist prior (alpha = 1) in dark grey.
#'
#' @param alpha.prior vector of alpha (dirichlet hyperparameters, none can be = 0)
#' @param source output from \code{\link{load_source_data}}
#' @param plot_save_pdf T/F, save the plot as a pdf?
#' @param plot_save_png T/F, save the plot as a png?
#' @param filename name of the file to save (e.g. "prior_plot")
#'
plot_prior <- function(alpha.prior = 1,source,plot_save_pdf=TRUE, plot_save_png=FALSE, filename="prior_plot"){
	# Error check for alpha = 0
	if(length(which(alpha.prior==0))!=0){
      stop(paste("*** Error: You cannot set any alpha = 0.
      Instead, set = 0.01.***",sep=""))
	}

	n.sources <- source$n.sources
	if(is.numeric(alpha.prior)==F) alpha.prior = 1 # Error checking for user inputted string/ NA
	if(length(alpha.prior)==1) alpha = rep(alpha.prior,n.sources) # All sources have same value
	if(length(alpha.prior) > 1 & length(alpha.prior) != n.sources) alpha = rep(1,n.sources) # Error checking for user inputted string/ NA
	if(length(alpha.prior) > 1 & length(alpha.prior) == n.sources) alpha = alpha.prior # All sources have different value inputted by user

	alpha.unif <- rep(1,n.sources)
	alpha.jeff <- rep(1/n.sources,n.sources)
	p = compositions::rDirichlet.rcomp(10000, alpha)
	p.unif = compositions::rDirichlet.rcomp(10000, alpha.unif)
	# p.jeff = compositions::rDirichlet.rcomp(10000, alpha.jeff)

	alpha_lab <- paste0("(",paste0(round(alpha,2),collapse=","),")",sep="")
	alpha.unif_lab <- paste0("(",paste0(round(alpha.unif,2),collapse=","),")",sep="")
	# alpha.jeff_lab <- paste0("(",paste0(round(alpha.jeff,2),collapse=","),")",sep="")

	dev.new()
	layout(matrix(c(seq(1:(2*n.sources)),(2*n.sources)+1,(2*n.sources)+1), ncol=2, byrow=TRUE), heights=c(rep(3,n.sources),2))
	par(mai=rep(0.3,4))
	for(i in 1:n.sources){
		hist(p[,i], breaks = seq(0,1,length.out=40),col="red", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
		hist(p.unif[,i], breaks = seq(0,1,length.out=40),col="darkgrey", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
		# hist(p.jeff[,i], breaks = seq(0,1,length.out=40),col="lightgrey", main = paste0("Source ",i,": ",alpha.jeff_lab),xlab=expression(p[i]),xlim=c(0,1))
	}

	par(mai=c(0,0,0,0))
	plot.new()
	legend(x="center", ncol=2,legend=c(paste0("Your prior: ",alpha_lab),paste0("\"Uninformative\" prior",alpha.unif_lab)),
	       fill=c("red","darkgrey"),bty = "n",cex=1.5)

	if(plot_save_pdf==TRUE){
		mypath <- file.path(paste(getwd(),"/",filename,".pdf",sep=""))
		cairo_pdf(filename=mypath, width=7, height=7)
		layout(matrix(c(seq(1:(2*n.sources)),(2*n.sources)+1,(2*n.sources)+1), ncol=2, byrow=TRUE), heights=c(rep(3,n.sources),2))
		par(mai=rep(0.3,4))
		for(i in 1:n.sources){
			hist(p[,i], breaks = seq(0,1,length.out=40),col="red", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
			hist(p.unif[,i], breaks = seq(0,1,length.out=40),col="darkgrey", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
		}

		par(mai=c(0,0,0,0))
		plot.new()
		legend(x="center", ncol=2,legend=c(paste0("Your prior: ",alpha_lab),paste0("\"Uninformative\" prior",alpha.unif_lab)),
		       fill=c("red","darkgrey"),bty = "n",cex=1.5)
		dev.off()
	}
	if(plot_save_png==TRUE){
		mypath <- file.path(paste(getwd(),"/",filename,".png",sep=""))
		png(filename=mypath)
		layout(matrix(c(seq(1:(2*n.sources)),(2*n.sources)+1,(2*n.sources)+1), ncol=2, byrow=TRUE), heights=c(rep(3,n.sources),2))
		par(mai=rep(0.3,4))
		for(i in 1:n.sources){
			hist(p[,i], breaks = seq(0,1,length.out=40),col="red", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
			hist(p.unif[,i], breaks = seq(0,1,length.out=40),col="darkgrey", main = paste0("Source ",i),xlab=expression(p[i]),xlim=c(0,1))
		}

		par(mai=c(0,0,0,0))
		plot.new()
		legend(x="center", ncol=2,legend=c(paste0("Your prior: ",alpha_lab),paste0("\"Uninformative\" prior",alpha.unif_lab)),
		       fill=c("red","darkgrey"),bty = "n",cex=1.5)
		dev.off()
	}
}
