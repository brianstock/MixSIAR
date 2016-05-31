#' Load trophic discrimination factor (TDF) data
#'
#' \code{load_discr_data} loads the trophic discrimination factor (TDF) data.
#' TDF is the amount that a consumer's tissue biotracer values are modified
#' (enriched/depleted) \emph{after} consuming a source. If tracers are conservative,
#' then set TDF = 0 (ex. essential fatty acids, fatty acid profile data,
#' element concentrations).
#'
#' @param filename csv file with the discrimination data
#' @param mix output from \code{\link{load_mix_data}}
#'
#' @return discr, a list including:
#' \itemize{
#'  \item \code{discr$mu}, matrix of discrimination means
#'  \item \code{discr$sig2}, matrix of discrimination variances
#' }
load_discr_data <- function(filename,mix){
  DISCR <- read.csv(filename)
  row.names(DISCR)<-DISCR[,1]     # store the row names of DISCR (sources)
  DISCR <- as.matrix(DISCR[-1])   # remove source names column of DISCR
  DISCR <- DISCR[order(rownames(DISCR)),]  # rearrange DISCR so sources are in alphabetical order

  # Make sure the iso columns of discr_mu and discr_sig2 are in the same order as S_MU and S_SIG
  # check that MU_names and SIG_names are in colnames(DISCR)
  if(sum(is.na(match(mix$MU_names,colnames(DISCR)))) > 0){
    stop(paste("*** Error: Discrimination mean column names mislabeled.
    Should be 'Mean' + iso_names from mix data file, e.g. 'Meand13C' if
    mix$iso_names = 'd13C'. Please ensure headings in discr data file match
    this format and try again.",sep=""))}
  if(sum(is.na(match(mix$SIG_names,colnames(DISCR)))) > 0){
    stop(paste("*** Error: Discrimination SD column names mislabeled.
    Should be 'SD' + iso_names from mix data file, e.g. 'SDd13C' if
    mix$iso_names = 'd13C'. Please ensure headings in discr data file match
    this format and try again.",sep=""))}

  discr_mu_cols <- match(mix$MU_names,colnames(DISCR))   # get the column numbers of DISCR that correspond to the means
  discr_sig_cols <- match(mix$SIG_names,colnames(DISCR))   # get the column numbers of DISCR that correspond to the SDs
  discr_mu <- as.matrix(DISCR[,discr_mu_cols])                            # DISCR means
  discr_sig2 <- as.matrix(DISCR[,discr_sig_cols]*DISCR[,discr_sig_cols])    # DISCR variances

  return(list(
    mu = discr_mu,
    sig2 = discr_sig2))
}

