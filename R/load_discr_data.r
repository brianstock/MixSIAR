# Brian Stock
# January 29, 2014

# Function: load_discr_data
# Usage: discr <- load_discr_data(filename,mix)
# Input: filename   csv file with the discrimination data
#        mix              output from 'load_mix_data'
# Output: discr, a list including:
#        discr$mu     matrix of discrimination means
#        discr$sig2   matrix of discrimination variances

load_discr_data <- function(filename,mix){
  DISCR <- read.csv(filename)
  row.names(DISCR)<-DISCR[,1]     # store the row names of DISCR (sources)
  DISCR <- as.matrix(DISCR[-1])   # remove source names column of DISCR
  DISCR <- DISCR[order(rownames(DISCR)),]  # rearrange DISCR so sources are in alphabetical order

  # Make sure the iso columns of discr_mu and discr_sig2 are in the same order as S_MU and S_SIG
  discr_mu_cols <- match(mix$MU_names,colnames(DISCR))   # get the column numbers of DISCR that correspond to the means
  discr_sig_cols <- match(mix$SIG_names,colnames(DISCR))   # get the column numbers of DISCR that correspond to the SDs
  discr_mu <- as.matrix(DISCR[,discr_mu_cols])                            # DISCR means
  discr_sig2 <- as.matrix(DISCR[,discr_sig_cols]*DISCR[,discr_sig_cols])    # DISCR variances

  return(list(
    mu = discr_mu,
    sig2 = discr_sig2))
}

