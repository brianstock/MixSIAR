#' Combine sources from a finished MixSIAR model (\emph{a posteriori})
#'
#' \code{combine_sources} aggregates the proportions from multiple sources.
#' Proportions are summed across posterior draws, since the source proportions
#' are correlated.
#'
#' @param jags.1 \code{rjags} model object, output from \code{\link{run_model}}
#' @param mix list, output from \code{\link{load_mix_data}}
#' @param source list, output from \code{\link{load_source_data}}
#' @param alpha.prior vector with length = n.sources, Dirichlet prior on p.global (default = 1, uninformative)
#' @param groups list, which sources to combine, and what names to give the new combined sources. See example.
#'   
#' @return \code{combined}, a list including:
#' \itemize{
#'  \item \code{combined$post}: matrix, posterior draws with new source groupings
#'  \item \code{combined$source.new}: list, original \code{source} list with modified entries for \code{n.sources} and \code{source_names}
#'  \item \code{combined$groups}: (input) list, shows original and combined sources
#'  \item \code{combined$jags.1}: (input) \code{rjags} model object
#'  \item \code{combined$source.old}: (input) list of original source data
#'  \item \code{combined$mix}: (input) list of original mix data
#'  \item \code{combined$alpha.prior}: (input) prior vector
#' }
#'
#' @seealso \code{\link{summary_stat}} and \code{\link{plot_intervals}}

#' @examples 
#' \dontrun{
#' # first run mantis shrimp example

#' # combine 6 sources into 2 groups of interest (hard-shelled vs. soft-bodied)
#' #   'hard' = 'clam' + 'crab' + 'snail'           # group 1 = hard-shelled prey
#' #   'soft' = 'alphworm' + 'brittlestar' + 'fish' # group 2 = soft-bodied prey
#' combined <- combine_sources(jags.1, mix, source, alpha.prior=alpha, 
#'                 groups=list(hard=c("clam","crab","snail"), soft=c("alphworm","brittlestar","fish")))
#'
#' # get posterior medians for new source groupings
#' apply(combined$post, 2, median)
#' summary_stat(combined, meanSD=FALSE, quantiles=c(.025,.5,.975), savetxt=FALSE)
#' }
#' 
combine_sources <- function(jags.1, mix, source, alpha.prior=1, groups){
  old.source.names <- sort(unlist(groups, use.names=FALSE))
  if(!identical(old.source.names, source$source_names)){
    stop(paste("Source names in 'groups' list do not match those in
  'source$source_names'. All previous sources must appear in 'groups'.",sep=""))    
  }
  
  # New source object (only names and # of sources changed, not data) 
  n.sources <- source$n.sources
  source.new <- source
  source.new$n.sources <- length(groups)
  source.new$source_names <- names(groups)

  # Get matrix of posterior draws
  # post <- jags.1$BUGSoutput$sims.array
  # names(dimnames(post)) <- c("iterations", "chains","parameters")
  post.mat <- jags.1$BUGSoutput$sims.matrix
  old.other <- post.mat[,-c(grep("p.global",colnames(post.mat)), grep("p.fac1",colnames(post.mat)), grep("p.fac2",colnames(post.mat)))]
  new.fac1 <- new.fac2 <- new.both <- NULL

  # Combine posterior draws into new source groupings
  old.global <- post.mat[,grep("p.global",colnames(post.mat))]
  n.draws <- dim(old.global)[1]
  new.global <- matrix(NA, nrow=n.draws, ncol=source.new$n.sources)
  colnames(new.global) <- paste0("p.global[",1:source.new$n.sources,"]")
  for(i in 1:source.new$n.sources){
    old <- groups[[i]]
    old.levels <- match(old, source$source_names)
    new.global[,i] <- apply(as.matrix(old.global[,old.levels]), 1, sum)
  }

  if(!mix$fere){
    # combine factor 1
    if(mix$n.effects > 0){
      for(f1 in 1:mix$FAC[[1]]$levels){
        new.fac1.tmp <- matrix(NA, nrow=n.draws, ncol=source.new$n.sources)
        colnames(new.fac1.tmp) <- paste0("p.fac1[",f1,",",1:source.new$n.sources,"]")
        old.fac1.tmp <- post.mat[,grep(paste0("p.fac1\\[",f1), colnames(post.mat))]
        for(i in 1:source.new$n.sources){
          old <- groups[[i]]
          old.levels <- match(old, source$source_names)
          new.fac1.tmp[,i] <- apply(as.matrix(old.fac1.tmp[,old.levels]), 1, sum)
        }
        new.fac1 <- cbind(new.fac1, new.fac1.tmp)
      }
    }
    # combine factor 2
    if(mix$n.effects > 1){
      for(f2 in 1:mix$FAC[[2]]$levels){
        new.fac2.tmp <- matrix(NA, nrow=n.draws, ncol=source.new$n.sources)
        colnames(new.fac2.tmp) <- paste0("p.fac2[",f2,",",1:source.new$n.sources,"]")
        old.fac2.tmp <- post.mat[,grep(paste0("p.fac2\\[",f2), colnames(post.mat))]
        for(i in 1:source.new$n.sources){
          old <- groups[[i]]
          old.levels <- match(old, source$source_names)
          new.fac2.tmp[,i] <- apply(as.matrix(old.fac2.tmp[,old.levels]), 1, sum)
        }
        new.fac2 <- cbind(new.fac2, new.fac2.tmp)
      }
    }
  }

  # Post-processing for 2 FE or 1FE + 1RE, calculate p.both = ilr.global + ilr.fac1 + ilr.fac2
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

    # Now combine sources for p.both
    for(f1 in 1:mix$FAC[[1]]$levels){
      for(f2 in 1:mix$FAC[[2]]$levels){
        new.both.tmp <- matrix(NA, nrow=n.draws, ncol=source.new$n.sources)
        colnames(new.both.tmp) <- paste0("p.both[",f1,",",f2,",",1:source.new$n.sources,"]")
        old.both.tmp <- p.both[,f1,f2,]
        for(i in 1:source.new$n.sources){
          old <- groups[[i]]
          old.levels <- match(old, source$source_names)
          new.both.tmp[,i] <- apply(as.matrix(old.both.tmp[,old.levels]), 1, sum)
        }
        new.both <- cbind(new.both, new.both.tmp)
      }
    }
  } # end fere

  # Combine posterior matrices
  post.new <- cbind(old.other, new.global, new.fac1, new.fac2, new.both)

  return(list(post=post.new, source.new=source.new, groups=groups, jags.1=jags.1, source.old=source, mix=mix, alpha.prior=alpha.prior))
}
