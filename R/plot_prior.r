# Brian Stock
# October 8, 2015

# Function: plot_prior
#   Plots the prior on diet proportions (p.global)
# Usage: plot_prior(alpha.prior)
# Input: alpha.prior 		vector of alpha (dirichlet hyperparameters)
# 		 source 		 	loaded source data (need n.sources)
# Output: none (displays plots)
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

plot_prior <- function(alpha.prior = 1,source,plot_save_pdf=TRUE, plot_save_png=FALSE, filename="prior_plot"){
	n.sources <- source$n.sources
	if(is.numeric(alpha.prior)==F) alpha.prior = 1 # Error checking for user inputted string/ NA
	if(length(alpha.prior)==1) alpha = rep(alpha.prior,n.sources) # All sources have same value
	if(length(alpha.prior) > 1 & length(alpha.prior) != n.sources) alpha = rep(1,n.sources) # Error checking for user inputted string/ NA
	if(length(alpha.prior) > 1 & length(alpha.prior) == n.sources) alpha = alpha.prior # All sources have different value inputted by user

	alpha.unif <- rep(1,n.sources)
	alpha.jeff <- rep(1/n.sources,n.sources)
	p = rDirichlet.rcomp(10000, alpha)
	p.unif = rDirichlet.rcomp(10000, alpha.unif)
	# p.jeff = rDirichlet.rcomp(10000, alpha.jeff)

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
