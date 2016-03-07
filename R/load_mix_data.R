# Brian Stock
# January 27, 2014

# Function: load_mix_data
# Usage: mix <- load_mix_data(filename,iso_names,random_effects,cont_effects)
# Input: filename (csv file with the mixture/consumer data),
#         iso_names (vector of isotope column headings in 'filename'),
#         factors (vector of categorical covariate column headings in 'filename'),
#         fac_random (vector of length(factors), is each factor a random effect?  TRUE = random effect, FALSE = fixed effect)
#         fac_nested (vector of length(factors), is each factor nested in the other?)
#           ex: 2 factors Region and Pack, both random effects, Pack nested within Region
#                 factors=c("Region","Pack")
#                 fac_random=c(TRUE,TRUE)
#                 fac_nested=c(FALSE,TRUE)
#           ex: fixed effect of Region and random effect of Pack, Region and Pack independent
#                 factors=c("Region","Pack")
#                 fac_random=c(FALSE,TRUE)
#                 fac_nested=c(FALSE,FALSE)
#         cont_effects (vector of continuous effect column headings in 'filename'),
#
#
# Output: mix, a list including:
#        mix$data (X)           raw consumer data
#        mix$data_iso (X_iso)   consumer isotope values (those specified in 'iso_names')
#        mix$n.iso              number of isotopes
#        mix$n.re               number of random effects
#        mix$n.ce               number of continuous effects
#        mix$RE                 list of random effect values (Factor.2), levels (factor2_levels), labels (factor2_levels), and lookup (factor1_lookup)
#        mix$random_effects     vector of included random effects (e.g. "Region" "Pack")
#        mix$CE                 list of continuous effect values (Cont.1)
#        mix$cont_effects       vector of included continuous effects (e.g. "Secchi.Mixed")
#        mix$MU_names           vector of mean column headings to look for in the source and discrimination files (e.g. 'Meand13C')
#        mix$SIG_names          vector of SD column headings to look for in the source and discrimination files (e.g. 'SDd13C')
#        mix$iso_names          vector of included isotopes (e.g. "d13C" "d15N")
#        mix$N                  scalar, number of mixture/consumer data points
#        mix$n.fe               number of fixed effects
#        mix$FE                 list of fixed effect values
#        mix$fixed_effects      vector of included fixed effects
#        mix$n.effects          number of random + fixed effects (n.re + n.fe)

# load_mix_data <- function(filename,iso_names,random_effects,fixed_effects,cont_effects){
load_mix_data <- function(filename,iso_names,factors,fac_random,fac_nested,cont_effects){
  X <- read.csv(filename)         # raw consumer data
  n.iso <- length(iso_names)      # number of isotopes
  # if n.iso = 0, error
  if(n.iso == 0){
  stop(paste("*** Error: No isotopes/tracers selected. Please select 1 or more
        isotopes/tracers, and then load your consumer/mixture data again. ***",sep=""))}
  if(length(fac_random) < length(factors)){
  stop(paste("*** Error: You have specified factors to include without saying
        if they are random or fixed effects (length of fac_random should
        match length of factors). Please check your load_mix_data line and try again. ***",sep=""))}
  if(length(factors)==2 && length(fac_nested)!=2){
  stop(paste("*** Error: You have specified factors to include without saying
        if they are nested or independent (length of fac_nested should
        match length of factors). Please check your load_mix_data line and try again. ***",sep=""))}
  if(length(factors)==2 & !is.na(fac_nested[1])){
    if(fac_nested[1]==TRUE && fac_nested[2]==TRUE){
      stop(paste("*** Error: Both factors cannot be nested within each other. Please check
            the fac_nested argument in your load_mix_data line and try again. ***",sep=""))}
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


