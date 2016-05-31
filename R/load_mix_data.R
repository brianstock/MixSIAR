#' Load mixture data
#'
#' \code{load_mix_data} loads the mixture data file and names the biotracers and
#' any Fixed, Random, or Continuous Effects.
#'
#' @param filename csv file with the mixture/consumer data
#' @param iso_names vector of isotope column headings in 'filename'
#' @param factors vector of random/fixed effect column headings in 'filename'.
#'          NULL if no factors.
#' @param fac_random vector of TRUE/FALSE, TRUE if factor is random effect, FALSE
#'          if fixed effect. NULL if no factors.
#' @param fac_nested vector of TRUE/FALSE, TRUE if factor is nested within the
#'          other. Only applies if 2 factors. NULL otherwise.
#' @param cont_effects vector of continuous effect column headings in 'filename'
#'
#' @return \code{mix}, a list including:
#' \itemize{
#'  \item \code{mix$data}: dataframe, raw mix/consumer data (all columns in 'filename'),
#'  \item \code{mix$data_iso}: matrix, mix/consumer biotracer/isotope values (those
#'        specified in 'iso_names'),
#'  \item \code{mix$n.iso}: integer, number of biotracers/isotopes,
#'  \item \code{mix$n.re}: integer, number of random effects,
#'  \item \code{mix$n.ce}: integer, number of continuous effects,
#'  \item \code{mix$FAC}: list of fixed/random effect values, each of which contains:
#'    \itemize{
#'      \item \code{values}: factor, values of the effect for each mix/consumer point
#'      \item \code{levels}: numeric vector, total number of values
#'      \item \code{labels}: character vector, names for each factor level
#'      \item \code{lookup}: numeric vector, if 2 factors and Factor.2 is nested
#'      within Factor.1, stores Factor.1 values for each level of Factor.2 (e.g.
#'      Wolf Ex has 8 Packs in 3 Regions, and \code{mix$FAC[[2]]$lookup =
#'      c(1,1,1,2,2,2,2,3)}, the Regions each Pack belongs to).
#'      \item \code{re}: T/F, is the factor a Random Effect? (FALSE = Fixed Effect)
#'      \item \code{name}: character, name of the factor (e.g. "Region")
#'    }
#'  \item \code{mix$CE}: list of length \code{n.ce}, contains the \code{cont_effects}
#'   values centered (subtract the mean) and scaled (divide by SD)
#'  \item \code{mix$CE_orig}: list of length \code{n.ce}, contains the original
#'   (unscaled) \code{cont_effects} values
#'  \item \code{mix$CE_center}: vector of length \code{n.ce}, means of each \code{cont_effects}
#'  \item \code{mix$CE_scale}: vector of length \code{n.ce}, SD of each \code{cont_effects}
#'  \item \code{mix$cont_effects}: vector of length \code{n.ce}, names of each \code{cont_effects}
#'  \item \code{mix$MU_names}: vector of biotracer/iso MEAN column headings to look for
#'  in the source and discrimination files (e.g. 'd13C' in \code{iso_names},
#'  'Meand13C' here)
#'  \item \code{mix$SIG_names}: vector of biotracer/iso SD column headings to look for
#'  in the source and discrimination files (e.g. 'd13C' in \code{iso_names},
#'  'SDd13C' here)
#'  \item \code{mix$iso_names}: vector of isotope column headings in 'filename' (same
#'  as input)
#'  \item \code{mix$N}: integer, number of mix/consumer data points
#'  \item \code{mix$n.fe}: integer, number of Fixed Effects
#'  \item \code{mix$n.effects}: integer, number of Fixed Effects + Random Effects
#'  \item \code{mix$factors}: vector of length \code{n.effects}, names of the
#'  Fixed and Random Effects
#'  \item \code{mix$fac_random}: T/F vector of length \code{n.effects} indicating
#'  which effects are Random (= TRUE) and Fixed (= FALSE)
#'  \item \code{mix$fac_nested}: T/F vector of length \code{n.effects} indicating
#'  which effects are nested within the other, if any
#'  \item \code{mix$fere}: TRUE if there are 2 Fixed Effects or 1 Fixed Effect and
#'  1 Random Effect, FALSE otherwise. Used by \code{write_JAGS_model}.
#' }
#'
#' If no biotracer/isotope columns are specified, a WARNING prompts the user to
#'  select 2, 1, or 0.
#'
#' If more than 2 Fixed/Random Effects are selected, a WARNING prompts the user
#' to select 2, 1, or 0.
#'
#' If more than 1 Continuous Effect is selected, a WARNING prompts the user to
#' select 1 or 0.
#'
load_mix_data <- function(filename,iso_names,factors,fac_random,fac_nested,cont_effects){
  X <- read.csv(filename)         # raw consumer data
  n.iso <- length(iso_names)      # number of isotopes
  # if n.iso = 0, error
  if(n.iso == 0){
  stop(paste("*** Error: No isotopes/tracers selected. Please select 1 or more
        isotopes/tracers, and then load your consumer/mixture data again. ***",sep=""))}
  if(length(fac_random) != length(factors)){
  stop(paste("*** Error: You have specified factors to include without saying
        if they are random or fixed effects (length of fac_random should
        match length of factors). Please check your load_mix_data line and try again. ***",sep=""))}
  if(length(factors)==2 && length(fac_nested)!=2){
  stop(paste("*** Error: You have specified factors to include without saying
        if they are nested or independent (length of fac_nested should
        match length of factors). Please check your load_mix_data line and try again. ***",sep=""))}
  if(length(factors)==2){
    if(!is.na(fac_nested[1])){
      if(fac_nested[1]==TRUE && fac_nested[2]==TRUE){
        stop(paste("*** Error: Both factors cannot be nested within each other. Please check
              the fac_nested argument in your load_mix_data line and try again. ***",sep=""))}
    }
  }

  n.effects <- length(factors)
  n.re <- sum(fac_random)  # number of random effects
  n.fe <- n.effects-n.re   # number of fixed effects
  fere <- ifelse(n.effects==2 & n.re < 2,TRUE,FALSE) # either 2 FE or 1FE + 1RE
  if(n.effects==1) fac_nested <- FALSE

  # if n.effects > 2, error
  if(n.effects > 2){
  stop(paste("*** Error: More than 2 random/fixed effects selected (MixSIAR can only
        currently handle 0, 1, or 2 random/fixed effects). Please choose 0, 1,
        or 2 random/fixed effects and then load your consumer/mixture data again. ***",sep=""))}

  n.ce <- length(cont_effects)    # number of continuous effects
  # if n.ce > 1, error
  if(n.ce > 1){
  stop(paste("*** Error: More than 1 continuous effect selected (MixSIAR can only
        currently handle 0 or 1 continuous effects). Please choose 0 or 1
        continuous effects and then load your consumer/mixture data again. ***",sep=""))}

  # check that iso_names, factors, and cont_effects are in colnames(X)
  if(sum(is.na(match(iso_names,colnames(X)))) > 0){
  stop(paste("*** Error: Your 'iso_names' do not match column names in your
        mixture data file (case sensitive). Please check your mix .csv data 
        file and load_mix_data line, then try again. ***",sep=""))}
  if(sum(is.na(match(factors,colnames(X)))) > 0){
  stop(paste("*** Error: Your 'factors' do not match column names in your
        mixture data file (case sensitive). Please check your mix .csv data
        file and load_mix_data line, then try again. ***",sep=""))}
  if(sum(is.na(match(cont_effects,colnames(X)))) > 0){
  stop(paste("*** Error: Your 'cont_effects' do not match column names in your
        mixture data file (case sensitive). Please check your mix .csv data 
        file and load_mix_data line, then try again. ***",sep=""))}

  N <- dim(X)[1]                  # number of consumer data points
  X_iso_cols <- match(iso_names,colnames(X))   # find the column indicies of the user-selected isotopes in X
  X_iso <- as.matrix(X[,X_iso_cols[]])         # keep the original X but create 'X_iso' to pass JAGS that only has the selected consumer isotope values (in the order of selection)
  MU_names <- paste("Mean",iso_names,sep="")  # vector of iso_names to look for means in the source and discrimination files (e.g. 'Meand13C')
  SIG_names <- paste("SD",iso_names,sep="")   # vector of iso_names to look for SDs in the source and discrimination files (e.g. 'SDd13C')

  FAC <- replicate(n.effects, NULL)     # FE is a list of length=n.fe that will contain the fixed_effects values, levels, and labels
  if(n.effects > 0){
    for(i in 1:n.effects){
      re <- fac_random[i] # is the ith factor a random effect (TRUE) or fixed effect (FALSE)?
      fac_values <- X[,factors[i]]
      fac_name <- factors[i]
      fac_levels <- length(unique(fac_values))
      if(is.numeric(fac_values)){ # if the factor was input as numbers, add the column name to each number (e.g. "Region 1", "Region 2", "Region 3")
        fac_labels <- paste(rep(factors[i],fac_levels),levels(factor(fac_values)),sep=" ")
      } else {  # if factor was input as text names, save them as is
        fac_labels <- levels(factor(fac_values))
      }
      fac_values <- as.numeric(factor(fac_values))
      FAC[[i]] <- list(values = fac_values,
                      levels = fac_levels,
                      labels = fac_labels,
                      lookup = NULL,
                      re = re,
                      name = fac_name)
    }
    if(n.re==2 & !is.na(fac_nested[1])){ # Only script version here. GUI sets nested = NA because need additional question answered
      if(n.re==2 & fac_nested[2]){
        for(lev in 1:FAC[[2]]$levels){
          FAC[[2]]$lookup[lev] <- FAC[[1]]$values[which(FAC[[2]]$values==lev)][1]
        }
      }
      if(n.re==2 & fac_nested[1]){
        for(lev in 1:FAC[[1]]$levels){
          FAC[[1]]$lookup[lev] <- FAC[[2]]$values[which(FAC[[1]]$values==lev)][1]
        }
      }
    }
    if(n.fe==1 & n.re==1 & fac_random[1]){ # make fac2 the random effect, fac1 the fixed effect
      tmp <- FAC[[1]]
      FAC[[1]] <- FAC[[2]]
      FAC[[2]] <- tmp
      factors <- rev(factors)
      fac_random <- rev(fac_random)
      fac_nested <- rev(fac_nested)
    }

  } # end random/fixed effects loop

  CE_orig <- replicate(n.ce, NULL)
  CE <- replicate(n.ce, NULL)  # CE is a list of length=n.ce that will contain the cont_effects values
  CE_center <- rep(NA,n.ce) # mean of CE
  CE_scale <- rep(NA,n.ce) # sd of CE
  if(n.ce > 0){ # If we have any continuous effects
    for(i in 1:n.ce){ # For each continuous effect selected, put the data from X into CE
      CE_orig[[i]] <- X[,cont_effects[i]]  # Get the values for CE[[1]] from X (raw consumer file)
      CE[[i]] <- scale(X[,cont_effects[i]],center=TRUE,scale=TRUE)
      CE_center[i] <- attributes(CE[[i]])$"scaled:center"
      CE_scale[i] <- attributes(CE[[i]])$"scaled:center"
    }
  } # end continuous effects loop

  # return a list of data, data_iso, n.iso, n.re, n.ce, RE, CE, random_effects, MU_names, SIG_names, iso_names, ...
  return(list(
      data = X,
      data_iso = X_iso,
      n.iso = n.iso,
      n.re = n.re,
      n.ce = n.ce,
      FAC = FAC,
      CE = CE,
      CE_orig = CE_orig,
      CE_center = CE_center,
      CE_scale = CE_scale,
      cont_effects = cont_effects,
      MU_names = MU_names,
      SIG_names = SIG_names,
      iso_names = iso_names,
      N = N,
      n.fe = n.fe,
      n.effects = n.effects,
      factors = factors,
      fac_random = fac_random,
      fac_nested = fac_nested,
      fere = fere))
} # end load_mix_data function


