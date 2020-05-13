#' Load source data
#'
#' \code{load_source_data} specifies the source data structure (factors,
#' concentration dependence, data type) and loads the source data file. \emph{Sources
#' are sorted alphabetically.}
#'
#' WARNING messages check for:
#' \itemize{
#'  \item More than one source factor selected
#'  \item Source factor not in mixture data
#'  \item Source sample sizes missing or entered incorrectly
#'  \item Source SD = 0
#' }
#'
#' @param filename character, csv file with the source data.
#' @param source_factors character, column heading in 'filename' that matches
#'   a Fixed or Random Effect from the mixture data (\code{mix$factors}).
#'   Only used if you have source data by a factor (e.g. "Region"), otherwise \code{NULL}.
#' @param conc_dep T/F, \code{TRUE} indicates you have concentration dependence
#'   data in 'filename'.
#' @param data_type \code{"raw"} or \code{"means"}. "Raw" source data are repeated
#'  source biotracer measurements, "means" data are source biotracer values as
#'  means, SDs, and sample size. See manual for formatting.
#' @param mix list, output from \code{\link{load_mix_data}}.
#'
#' @return \code{source}, a list including:
#' \itemize{
#'   \item \code{source$n.sources}: integer, number of sources
#'   \item \code{source$source_names}: vector, source names/labels
#'   \item \code{source$S_MU}: matrix, source means used for plotting - NOT
#'   passed to JAGS. If sources are by factor, then the third column of S_MU
#'   will be the factor values (e.g. for 4 sources and 3 Regions:
#'   1 2 3 1 2 3 1 2 3 1 2 3)
#'   \item \code{source$S_SIG}: matrix, source SDs used for plotting - NOT passed
#'    to JAGS. Same structure as S_MU.
#'   \item \code{source$S_factor1}: factor or NULL, factor values if sources are
#'    by factor.
#'   \item \code{source$S_factor_levels}: scalar or NULL, number of \code{S_factor1}
#'    levels if sources are by factor.
#'   \item \code{source$conc}: matrix or NULL, concentration dependence values
#'   for each isotope
#'   \item \code{source$MU_array}: array of source means, dim(src,iso,f1) or
#'   dim(src,iso) if data_type="means", NULL if data_type="raw".
#'   \item \code{source$SIG2_array}: array of source variances, dim(src,iso,f1)
#'   or dim(src,iso) if data_type="means", NULL if data_type="raw".
#'   \item \code{source$n_array}: vector/matrix of source sample sizes,
#'   dim(src,f1) or dim(src) if data_type="means", NULL if data_type="raw".
#'   \item \code{source$SOURCE_array}: array of source data, dim(src,iso,f1,replicate)
#'    or dim(src,iso,replicate) if data_type="raw", NULL if data_type="means".
#'   \item \code{source$n.rep}: vector/matrix of source sample sizes, dim(src,f1)
#'    or dim(src) if data_type="raw", NULL if data_type="means".
#'   \item \code{source$by_factor}: NA or factor number, are the source data by a Fixed or Random Effect?
#'   \item \code{source$data_type}: \code{"raw"} or \code{"means"}, same as input.
#'   \item \code{source$conc_dep}: T/F, same as input.
#' }
#'
#' @seealso \code{\link{load_mix_data}} and \code{\link{load_discr_data}}
#' @export
#'
load_source_data <- function(filename,source_factors=NULL,conc_dep,data_type,mix){
  SOURCE <- read.csv(filename)
  source.fac <- length(source_factors)
  if(source.fac > 1){
    stop(paste("More than one source factor.
    MixSIAR can only fit source data by up to ONE factor.
    Please specify 0 or 1 source factor and try again.",sep=""))
  }
  # source factor must be in colnames(SOURCE)
  if(sum(is.na(match(source_factors,colnames(SOURCE)))) > 0){
    stop(paste("Your 'source_factors' do not match column names in your
        source data file (case sensitive). Please check your source .csv data
        file and load_source_data line, then try again. ***",sep=""))
  }
  test_fac <- match(source_factors,mix$factors)
  if(source.fac==1 && (length(test_fac)==0 || is.na(test_fac))){
    stop(paste("Source factor not in mix$factors.
    You cannot model a source random effect that is not included
    as a random/fixed effect for the mixture/consumer. Either
     1) remove the source factor (reload source data), or
     2) include the random/fixed effect in the mixture (reload mix data).
    Could be a mismatch between column headings in the mix and source data files.",sep=""))
  }
  if(source.fac==0) by_factor <- NA else by_factor <- match(source_factors, mix$factors)

  # turn source names into numbers
  src <- as.factor(SOURCE[,1])
  source_names <- levels(src)   # first save the source names in source_names
  n.sources <- length(source_names)    # n.sources is the number of sources, which is the length of the source names vector
  levels(src) <- 1:n.sources    # convert the source names in SOURCE into numbers
  SOURCE[,1] <- as.numeric(src)
  # sorts SOURCE by source name, then by the source factors (if present)
  source_factor_cols <- match(source_factors,colnames(SOURCE)) # the column number(s) of the user-selected random effects
  S_factor_levels <- rep(0,length(source_factor_cols))                # the number of levels of each source random effect
  if(!is.na(by_factor)){
    SOURCE <- SOURCE[order(SOURCE[,1],SOURCE[,source_factor_cols]),]
  } else {  # by_factor = NA
    SOURCE <- SOURCE[order(SOURCE[,1]),]
  }
  # Concentration Dependence section (must be after SOURCE is sorted)
  if(conc_dep){   # if we have concentration dependence data
    CONC_names <- paste("Conc",mix$iso_names,sep="")
    # check that CONC_names are in colnames(SOURCE)
    if(sum(is.na(match(CONC_names,colnames(SOURCE)))) > 0){
    stop(paste("Concentration dependence column names mislabeled.
    Should be 'Conc' + iso_names, e.g. 'Concd13C' if iso_names = 'd13C'.
    Please ensure Conc headings in source data file match iso_names
    in mix data file and try again. Alternatively, you may have set conc_dep=T
    by mistake.",sep=""))}

    CONC_iso_cols <- match(CONC_names,colnames(SOURCE))
    conc <- do.call(rbind,lapply(split(SOURCE[,CONC_iso_cols],list(SOURCE[,1])),colMeans))  # calculate conc means for each [src,iso]
  } else conc <- NULL

  # colSds function
  # calculates the sd of each column of a matrix or data frame (x)
  # formerly from the 'matrixStats' package... but that stopped working in R 3.1
  col_sd <- function(x){
    sds <- apply(x,2,sd)
    return(sds)
  }

  if(data_type=="raw"){
    if(sum(is.na(match(mix$iso_names,colnames(SOURCE)))) > 0){
      stop(paste("With raw source data, the iso_names in mix data
        file must be in the column headings of source data file. Please check
        your source and mix .csv data files and try again. ***",sep=""))}
    S_iso_cols <- match(mix$iso_names,colnames(SOURCE))   # find the column numbers of the user-selected isotopes
    # Create S_MU and S_SIG - the source means and sds by isotopes and fac1 (if source data are by factor)
    if(!is.na(by_factor)){  # if we have raw source data BY FACTOR
      for(fac in 1:length(source_factor_cols)){
        S_factor_levels[fac] <- length(levels(SOURCE[,source_factor_cols[fac]]))
      }

      # if by_factor, test that the sources are IDENTICAL for each level
      list.sources.bylev <- vector("list", S_factor_levels)
      test.sources.identical <- rep(TRUE, S_factor_levels)
      S_factor_levels_names <- levels(as.factor(SOURCE[,source_factor_cols[fac]]))
      for(lev in 1:S_factor_levels){
        list.sources.bylev[[lev]] <- SOURCE[SOURCE[,source_factors] == S_factor_levels_names[lev], 1]
      }
      if(length(unique(list.sources.bylev))!=1){
        stop(paste("Sources for each level of *",source_factors,"* do not match.
        If you have different sources (or # of sources) for levels of *",source_factors,"*,
        you must fit separate MixSIAR models for each level.",sep=""))        
      }

      if(mix$n.iso > 1) mu <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),colMeans))  # calculate the means for each iso/source/fac1 combination
      if(mix$n.iso==1) mu <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),mean))
      S_MU <- cbind(mu,rep(1:S_factor_levels,n.sources))   # S_MU has source means by (columns): isotopes and fac1.  Sorted by source name and then factor 1
      if(mix$n.iso > 1) sig <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),col_sd))   # calculate the sds for each iso/source/fac1 combination
      if(mix$n.iso==1) sig <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),sd))
      S_SIG <- cbind(sig,rep(1:S_factor_levels,n.sources))  # S_SIG has source sds by (columns): isotopes and fac1.  Sorted by source name and then factor 1
      colnames(S_MU) <- c(mix$iso_names, source_factors)
      S_MU_factor_col <- match(source_factors,colnames(S_MU))
      S_factor1 <- S_MU[,S_MU_factor_col]
    }
    if(is.na(by_factor)){ # if we have raw source data NOT by factor
      if(mix$n.iso > 1) S_MU <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),colMeans))   # calculate the means for each iso/source combination
      if(mix$n.iso==1) S_MU <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),mean))
      if(mix$n.iso > 1) S_SIG <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),col_sd))   # calculate the sds for each iso/source combination
      if(mix$n.iso==1) S_SIG <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),sd))
      colnames(S_MU) <- c(mix$iso_names)
      S_factor1 <- NULL; S_factor_levels <- NULL
    }

    # make 'n.rep': the array of source/factor replicates
    n.rep <- array(0,dim=c(n.sources,S_factor_levels))
    if(!is.na(by_factor)){
      for(src in 1:n.sources){
        for(f1 in 1:S_factor_levels){
          n.rep[src,f1] <- table(SOURCE[which(SOURCE[,1]==src),source_factor_cols])[f1]
        }
      }
    } else {  # by_factor = NA
        for(src in 1:n.sources){
        n.rep[src] <- length(which(SOURCE[,1]==src))
      }
    }
    max.rep <- max(n.rep)

    # Make 'SOURCE_array' to pass to JAGS: the SOURCE data indexed by source, isotope, factor1, and replicate
    # For this to work SOURCE needs to be sorted by source and then factor1 (done above)
    SOURCE_array <- array(NA,dim=c(n.sources,mix$n.iso,S_factor_levels,max.rep))
    count <- 1
    if(!is.na(by_factor)){
      for(src in 1:n.sources){
        for(f1 in 1:S_factor_levels){
          for(r in 1:n.rep[src,f1]){
            for(iso in 1:mix$n.iso){
              SOURCE_array[src,iso,f1,r] <- SOURCE[count,mix$iso_names[iso]]  # columns of SOURCE are: source number, iso1, iso2, fac1 (we defined SOURCE above, so we know the column indicies are in this format)
            }
            count <- count+1
          }
        }
      }
    } else {  # by_factor = NA
        for(src in 1:n.sources){
          for(r in 1:n.rep[src]){
            for(iso in 1:mix$n.iso){
              SOURCE_array[src,iso,r] <- SOURCE[count,mix$iso_names[iso]]   # columns of SOURCE are: source number, iso1, iso2 (we defined SOURCE above, so we know the column indicies are in this format)
            }
            count <- count+1
          }
        }
    }

    # Create NULL objects that would be used in "means" source loading
    MU_array <- NULL; SIG2_array <- NULL; n_array <- NULL
  } # end RAW data loading

  if(data_type=="means"){
    # check that MU_names and SIG_names are in colnames(SOURCE)
    if(sum(is.na(match(mix$MU_names,colnames(SOURCE)))) > 0){
      stop(paste("Source mean column names mislabeled.
    Should be 'Mean' + iso_names from mix data file, e.g. 'Meand13C' if
    mix$iso_names = 'd13C'. Please ensure headings in source data file match
    this format and try again. Alternatively, if you have raw source data,
    you should set data_type='raw'.",sep=""))}
    if(sum(is.na(match(mix$SIG_names,colnames(SOURCE)))) > 0){
      stop(paste("Source SD column names mislabeled.
    Should be 'SD' + iso_names from mix data file, e.g. 'SDd13C' if
    mix$iso_names = 'd13C'. Please ensure headings in source data file match
    this format and try again. Alternatively, if you have raw source data,
    you should set data_type='raw'.",sep=""))}

    S_MU_iso_cols <- match(mix$MU_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes
    S_SIG_iso_cols <- match(mix$SIG_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes

    sample_size_col <- match("n",colnames(SOURCE))    # Find the column titled "n" be in the source means file, where n is the sample size for each source
    if(is.na(sample_size_col)){
      stop(paste("Source sample sizes missing or entered incorrectly.
        Check your sources .csv data file to be sure you have a column
        titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
    }
    S_sample_size <- SOURCE[,sample_size_col]         # Get the sample sizes

    if(is.na(by_factor)){
      S_MU <- as.matrix(SOURCE[,c(S_MU_iso_cols[],source_factor_cols)]) # rearrange S_MU columns to be the selected isotopes (in order) and source random effects
      S_SIG <- as.matrix(SOURCE[,c(S_SIG_iso_cols[],source_factor_cols)]) # rearrange S_SIG columns to be the selected isotopes (in order) and source random effects
      MU_array <- S_MU                      # MU_array is passed to JAGS
      SIG2_array <- S_SIG*S_SIG             # SIG2_array is passed to JAGS
      n_array <- S_sample_size              # n_array is passed to JAGS
      S_factor1 <- NULL; S_factor_levels <- NULL
    } else {  # if by_factor = 1 or 2
      for(fac in 1:length(source_factor_cols)){
        S_factor_levels[fac] <- length(levels(as.factor(SOURCE[,source_factor_cols[fac]])))
      }

      # if by_factor, test that the sources are IDENTICAL for each level
      list.sources.bylev <- vector("list", S_factor_levels)
      test.sources.identical <- rep(TRUE, S_factor_levels)
      S_factor_levels_names <- levels(as.factor(SOURCE[,source_factor_cols[fac]]))
      for(lev in 1:S_factor_levels){
        list.sources.bylev[[lev]] <- SOURCE[SOURCE[,source_factors] == S_factor_levels_names[lev], 1]
      }
      if(length(unique(list.sources.bylev))!=1){
        stop(paste("Sources for each level of *",source_factors,"* do not match.
        If you have different sources (or # of sources) for levels of *",source_factors,"*,
        you must fit separate MixSIAR models for each level.",sep=""))        
      }

      if(length(S_sample_size) != (n.sources*S_factor_levels)){
        stop(paste("Source sample sizes missing or entered incorrectly.
        Check your source_means.csv data file to be sure you have a column
        titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
      }
      S_factor1 <- factor(SOURCE[,source_factor_cols])               # save the Factor1 labels for the isospace plot
      SOURCE[,source_factor_cols] <- as.numeric(factor(SOURCE[,source_factor_cols]))  # turn the Factor1 labels into numbers (already have them saved as 'factor1_levels' from X)
      S_MU <- as.matrix(SOURCE[,c(S_MU_iso_cols[],source_factor_cols)]) # rearrange S_MU columns to be the selected isotopes (in order) and source random effects
      S_SIG <- as.matrix(SOURCE[,c(S_SIG_iso_cols[],source_factor_cols)]) # rearrange S_SIG columns to be the selected isotopes (in order) and source random effects

      # Create MU_array, SIG2_array, and n_array
      MU_array <- array(NA,dim=c(n.sources,mix$n.iso,S_factor_levels))
      SIG2_array <- array(NA,dim=c(n.sources,mix$n.iso,S_factor_levels))
      n_array <- array(NA,dim=c(n.sources,S_factor_levels))
      count <- 1
      for(src in 1:n.sources){
        for(f1 in 1:S_factor_levels){
          for(iso in 1:mix$n.iso){
            MU_array[src,iso,f1] <- S_MU[count,iso]      # MU_array is the source mean data, but indexed as [src,iso,f1]. Passed to JAGS
            SIG2_array[src,iso,f1] <- S_SIG[count,iso]*S_SIG[count,iso]   # same for SIG2_array, the source sd data
          }
          n_array[src,f1] <- S_sample_size[count]        # n_array is the source sample sizes but indexed as [src,f1].  Passed to JAGS
          count <- count+1
        }
      }
    } # end if(by_factor)==TRUE

    # Create NULL objects that would be used in "raw" source loading
    SOURCE_array <- NULL; n.rep <- NULL
  } # end MEANS data loading

  # Error check for zero SD
  if(length(which(S_SIG==0))>0){
    stop(paste("You have at least one source SD = 0.
    Check your source data file to be sure each SD entry is non-zero.",sep=""))
  }

  return(list(
      n.sources = n.sources,
      source_names = source_names,
      S_MU = S_MU,
      S_SIG = S_SIG,
      S_factor1 = S_factor1,
      S_factor_levels = S_factor_levels,
      conc = conc,
      MU_array = MU_array,
      SIG2_array = SIG2_array,
      n_array = n_array,
      SOURCE_array = SOURCE_array,
      n.rep = n.rep,
      by_factor = by_factor,
      data_type = data_type,
      conc_dep = conc_dep))
}


