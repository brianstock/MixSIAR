# Brian Stock
# October 2015

# Function: calc_area
#   Calculates the normalized surface area of the SOURCE + TDF convex hull (for 2 dimensions)
# Usage: area <- calc_area(source)
# Input: source           output from 'load_source_data'
# Output: scalar (if source$by_factor == FALSE) 	normalized surface area/volume of the convex hull
#         vector (if source$by_factor == TRUE)		for each factor, the normalized surface area/volume of the convex hull

calc_area <- function(source,mix,discr){
	if(mix$n.iso!=2){
      stop(paste("*** Error: MixSIAR can only calculate normalized polygon
      	area if there are 2 dimensions (tracers/isotopes).***",sep=""))}

	library(splancs)
	if(mix$n.iso==2){
		if(!source$by_factor){ # == FALSE
			MU_plot <- source$S_MU + discr$mu   			# source means adjusted for fractionation/enrichment
			SIG_plot <- sqrt(source$S_SIG^2 + discr$sig2) 	# source sds adjusted for fractionation/enrichment
			x <- MU_plot[,1]
			y <- SIG_plot[,2]
			ind <- chull(x,y)
			area <- splancs::areapl(cbind(x[ind],y[ind]))
			bot <- prod(apply(SIG_plot,2,mean))
			val <- area/bot
		}
		if(source$by_factor){ # == TRUE
			MU_plot <- SIG_plot <- source$MU_array
			val <- rep(NA,source$S_factor_levels)
			for(f1 in 1:source$S_factor_levels){
				MU_plot[,,f1] <- source$MU_array[,,f1] + discr$mu 				# source means adjusted for fractionation/enrichment
				SIG_plot[,,f1] <- sqrt(source$SIG2_array[,,f1] + discr$sig2) 	# source sds adjusted for fractionation/enrichment
				x <- MU_plot[,1,f1]
				y <- SIG_plot[,2,f1]
				ind <- chull(x,y)
				area <- splancs::areapl(cbind(x[ind],y[ind]))
				bot <- prod(apply(SIG_plot[,,f1],2,mean))
				val[f1] <- area/bot
			}
		}
		return(val)
	}
}
