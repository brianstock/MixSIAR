#' Calculate the normalized surface area of the source convex hull
#'
#' \code{calc_area()} calculates the normalized surface area of the SOURCE + TDF
#'  convex hull, only if there are exactly 2 biotracers.
#'
#' Important detail is that, unlike in Brett (2014), \code{calc_area} uses the
#' combined SOURCE + TDF variance to normalize the surface area:
#' \deqn{\sqrt{\sigma^2_source + \sigma^2_discr}}
#' This is the variance used in fitting the mixing model.
#'
#' \code{calc_area()} relies on the \code{splancs::areapl()} function from the \code{splancs}
#' package. If \code{splancs} is not installed, a WARNING message will appear.
#'
#' @param source output from \code{\link{load_source_data}}
#' @param mix output from \code{\link{load_mix_data}}
#' @param discr output from \code{\link{load_discr_data}}
#'
#' @return If source$by_factor = FALSE, \code{calc_area} returns a scalar, the
#'   normalized surface area of the SOURCE + TDF convex hull
#' @return If source$by_factor = TRUE, \code{calc_area} returns a vector, where
#'   the entries are the normalized surface areas of the convex hull of each
#'   source factor level (e.g. source data by 3 Regions, returns a 3-vector of
#'   the areas of the Region 1 convex hull, Region 2 convex hull, etc.)
#'
#' @seealso Brett (2014): \url{https://www.researchgate.net/profile/Michael_Brett/publication/269873625_Resource_polygon_geometry_predicts_Bayesian_stable_isotope_mixing_model_bias/links/549884090cf2519f5a1de635.pdf}
#' @export
calc_area <- function(source,mix,discr){
  if (!"splancs" %in% installed.packages()){
    stop(paste("*** Error: 'splancs' package not installed. 'splancs' is
        required to run the 'calc_area' function. Install 'splancs' by
        entering 'install.packages('splancs')' in the R console. ***",sep=""))}

  if(mix$n.iso!=2){
      stop(paste("*** Error: MixSIAR can only calculate normalized polygon
      	area if there are 2 dimensions (tracers/isotopes).***",sep=""))}

	# library(splancs)
	if(mix$n.iso==2){
		if(is.na(source$by_factor)){ # == FALSE
			MU_plot <- source$S_MU + discr$mu   			# source means adjusted for fractionation/enrichment
			SIG_plot <- sqrt(source$S_SIG^2 + discr$sig2) 	# source sds adjusted for fractionation/enrichment
			x <- MU_plot[,1]
			y <- MU_plot[,2]
			ind <- chull(x,y)
			area <- splancs::areapl(cbind(x[ind],y[ind]))
			bot <- prod(apply(SIG_plot,2,mean))
			val <- area/bot
		}
		if(!is.na(source$by_factor)){ # == TRUE
			MU_plot <- SIG_plot <- source$MU_array
			val <- rep(NA,source$S_factor_levels)
			for(f1 in 1:source$S_factor_levels){
				MU_plot[,,f1] <- source$MU_array[,,f1] + discr$mu 				# source means adjusted for fractionation/enrichment
				SIG_plot[,,f1] <- sqrt(source$SIG2_array[,,f1] + discr$sig2) 	# source sds adjusted for fractionation/enrichment
				x <- MU_plot[,1,f1]
				y <- MU_plot[,2,f1]
				ind <- chull(x,y)
				area <- splancs::areapl(cbind(x[ind],y[ind]))
				bot <- prod(apply(SIG_plot[,,f1],2,mean))
				val[f1] <- area/bot
			}
		}
		return(val)
	}
}
