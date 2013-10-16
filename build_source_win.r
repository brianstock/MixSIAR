#' Build the window to read-in SOURCE data
#'
#' This function is run when a user clicks the "Load source data" button in the main MixSIAR GUI,
#' and creates a separate GUI window (\code{source_win}) where the user loads the source data. It contains:
#'  1) Questions the user if their source data is by any of their Random Effects (if they have Random Effects)
#'  2) Questions the user if they have Concentration Dependence data in their source file,
#'  3) Individual Effect gcheckbox ("Include 'Individual' as a Random Effect"),
#'  4) "I'm finished" gbutton that closes \code{mix_win} and manipulates \code{X} into all the objects we use later.
#' If more than 2 Random Effects are selected, a separate WARNING gwindow prompts the user to select 2, 1, or 0.
#' If 2 Random Effects are selected, a separate gwindow asks the user if the effects are hierarchical/nested.
#' If more than 1 Continuous Effect is selected, a separate WARNING gwindow prompts the user to select 1 or 0.
#' Finally, the function adds a green check image if the data is successfully loaded, or a red x image if not.
build_source_win <- function(){
  source_win <- gwindow("Read in your SOURCE data")
  source_grp_all <- ggroup(horizontal=F, container=source_win)
  if(n.re > 0){    # if we have factors in the mix data
    re_list <<- vector("list",n.re)
    for(i in 1:n.re){      # for each random effect
      cur <- random_effects[i]
      grp <- ggroup(cont=source_grp_all,horizontal=T)
      lbl <- glabel(
          paste("Does your source data vary by ",cur,"?",sep=""),
          cont=grp)
      re_list[[i]] <<- gradio(c("Yes","No"), cont=grp, horizontal=T)
      svalue(re_list[[i]]) <<- "No"
    }
  }
  grp_conc <- ggroup(cont=source_grp_all,horizontal=T)
  lbl_conc <- glabel("Do you have Concentration Dependence data?",cont=grp_conc)
  conc_rad <- gradio(c("Yes","No"),cont=grp_conc,horizontal=T)
  svalue(conc_rad) <- "No"

  addSpring(source_grp_all)
  lbl_data_type <- glabel("Do you have raw source data, or source means and SDs?",
                        cont=source_grp_all)
  grp_btns <- ggroup(cont=source_grp_all, horizontal=T)
  addSpring(grp_btns)
  grp_raw <- ggroup(cont=grp_btns, horizontal=T)
  raw_data <- gbutton(
    text = "Load raw source data",
    cont = grp_raw,
    expand = T,
    handler = function(h, ...){
      gfile(
        text = "Load raw source data file",
        type = "open",
        action = "read.csv",
        handler = function(h, ...){
          tryCatch(   # reads source data into 'SOURCE'
            {
            data_frame_name <- make.names("SOURCE")
            the_data <- do.call(h$action, list(h$file, header=T))
            assign(data_frame_name, the_data, envir = globalenv())
            add(grp_raw,gimage("check.png"))
            },
            error = function(e){
              svalue(status_bar) <- "Could not load data"
              add(grp_raw,gimage("red x.png"))
            }
          )
        }
      )
      source_random_effects <- character(0)  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
      if(n.re > 0){    # if we have random effects in the mix data
        for(i in 1:n.re){  # for each random effect
          if(svalue(re_list[[i]])=="Yes"){   # if random_effect[i] is selected by the user, add it to source_random_effects
            source_random_effects <- c(source_random_effects,random_effects[i])
          }
        }
      }
      if(length(source_random_effects)==0) sources_by_factor <<- FALSE else sources_by_factor <<- TRUE

      S_iso_cols <- match(iso_names,colnames(SOURCE))   # find the column numbers of the user-selected isotopes
      source_factor_cols <- match(source_random_effects,colnames(SOURCE)) # the column number(s) of the user-selected random effects
      source_factor_levels <- rep(0,length(source_factor_cols))                # the number of levels of each source random effect
      n.sources <- length(levels(SOURCE[,1]))        # n.sources is the number of sources, which is the length of the source names vector
                                                      # the first column of SOURCE MUST be the source names

      # Create S_MU and S_SIG - the source means and sds by isotopes and fac1 (if source data are by factor)
      if(sources_by_factor){  # if we have raw source data BY FACTOR
        for(fac in 1:length(source_factor_cols)){
          source_factor_levels[fac] <- length(levels(SOURCE[,source_factor_cols[fac]]))
        }
        mu <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),colMeans))  # calculate the means for each iso/source/fac1 combination
        S_MU <- cbind(mu,rep(1:source_factor_levels,n.sources))   # S_MU has source means by (columns): isotopes and fac1.  Sorted by source name and then factor 1
        sig <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,source_factor_cols[1]],SOURCE[,1])),colSds))   # calculate the sds for each iso/source/fac1 combination
        S_SIG <- cbind(sig,rep(1:source_factor_levels,n.sources))  # S_SIG has source sds by (columns): isotopes and fac1.  Sorted by source name and then factor 1
        colnames(S_MU) <- c(iso_names, source_random_effects)
        S_MU_factor_col <- match(source_random_effects,colnames(S_MU))
        S_factor1 <- S_MU[,S_MU_factor_col]
        assign("S_factor1", S_factor1, envir = .GlobalEnv)
      }
      if(!sources_by_factor){ # if we have raw source data NOT by factor
        S_MU <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),colMeans))   # calculate the means for each iso/source combination
        S_SIG <- do.call(rbind,lapply(split(SOURCE[,S_iso_cols],list(SOURCE[,1])),colSds))   # calculate the sds for each iso/source combination
        colnames(S_MU) <- c(iso_names)
      }
      assign("S_MU", S_MU, envir = .GlobalEnv)
      assign("S_SIG", S_SIG, envir = .GlobalEnv)

      # turn source names into numbers
      source_names <- levels(SOURCE[,1])   # first save the source names in source_names
      levels(SOURCE[,1]) <- 1:n.sources    # convert the source names in SOURCE into numbers
      SOURCE[,1] <- as.numeric(SOURCE[,1])
      # sorts SOURCE by source name, then by the source factors (if present)
      if(sources_by_factor){
        SOURCE <- SOURCE[order(SOURCE[,1],SOURCE[,source_factor_cols]),]
      } else {  # sources_by_factor==FALSE
        SOURCE <- SOURCE[order(SOURCE[,1]),]
      }

      # make 'n.rep': the array of source/factor replicates
      n.rep <- array(0,dim=c(n.sources,source_factor_levels))
      if(sources_by_factor){
        for(src in 1:n.sources){
          for(f1 in 1:source_factor_levels){
            n.rep[src,f1] <- table(SOURCE[which(SOURCE[,1]==src),source_factor_cols])[f1]
          }
        }
      } else {  # sources_by_factor==FALSE
          for(src in 1:n.sources){
          n.rep[src] <- length(which(SOURCE[,1]==src))
        }
      }
      assign("n.rep", n.rep, envir = .GlobalEnv)
      max_rep <<- max(n.rep)

      # Make 'SOURCE_array' to pass to JAGS: the SOURCE data indexed by source, isotope, factor1, and replicate
      # For this to work SOURCE needs to be sorted by source and then factor1 (done above)
      SOURCE_array <- array(NA,dim=c(n.sources,n.iso,source_factor_levels,max_rep))
      count <- 1
      if(sources_by_factor){
        for(src in 1:n.sources){
          for(f1 in 1:source_factor_levels){
            for(r in 1:n.rep[src,f1]){
              for(iso in 1:n.iso){
                SOURCE_array[src,iso,f1,r] <- SOURCE[count,iso_names[iso]]  # columns of SOURCE are: source number, iso1, iso2, fac1 (we defined SOURCE above, so we know the column indicies are in this format)
              }
              count <- count+1
            }
          }
        }
      } else {  # sources_by_factor==FALSE
          for(src in 1:n.sources){
            for(r in 1:n.rep[src]){
              for(iso in 1:n.iso){
                SOURCE_array[src,iso,r] <- SOURCE[count,iso_names[iso]]   # columns of SOURCE are: source number, iso1, iso2 (we defined SOURCE above, so we know the column indicies are in this format)
              }
              count <- count+1
            }
          }
      }

      # Concentration Dependence section
      if(svalue(conc_rad)=="Yes") include_conc <- TRUE else include_conc <- FALSE   # do we have concentration dependence data?
      if(include_conc){   # if we have concentration dependence data
        CONC_names <- paste("Conc",iso_names,sep="")
        CONC_iso_cols <- match(CONC_names,colnames(SOURCE))
        conc <- do.call(rbind,lapply(split(SOURCE[,CONC_iso_cols],list(SOURCE[,1])),colMeans))  # calculate conc means for each [src,iso]
        assign("conc", conc, envir=.GlobalEnv)
      }
      assign("include_conc", include_conc, envir=.GlobalEnv)

# S_MU was created with columns of 'iso_names' 'source_factor_cols'...      S_iso_cols <- match(iso_names,colnames(S_MU))   # If we have raw source data, S_MU and S_SIG are only used for the isospace plot.  We need to know which columns of S_MU correspond to the isotopes
#      source_factor_cols <- match(source_random_effects,colnames(S_MU)) # We need to know which columns of S_MU (and S_SIG) correspond to the source random effect (if it exists)

      assign("n.sources", n.sources, envir = .GlobalEnv)
      assign("source_names", source_names, envir = .GlobalEnv)
#      assign("S_iso_cols", S_iso_cols, envir = .GlobalEnv)
#      assign("source_factor_cols", source_factor_cols, envir = .GlobalEnv)
      assign("source_factor_levels", source_factor_levels, envir = .GlobalEnv)
      assign("SOURCE_array", SOURCE_array, envir = .GlobalEnv)
      raw_source_data <<- TRUE  # We have raw source data (the user clicked the "I have raw source data" button)
    }
  )

  lbl_or <- glabel(" OR ", cont=grp_btns)

  grp_means_sds <- ggroup(cont=grp_btns, horizontal=T)
  means_sds <- gbutton(
    text = "Load source means and SDs",
    cont = grp_means_sds,
    expand = T,
    handler	= function(h, ...){
      gfile(
        text = "Load source data file",
        type = "open",
        action = "read.csv",
        handler = function(h, ...){
          tryCatch(   # reads source data into 'SOURCE'
            {
            data_frame_name <- make.names("SOURCE")
            the_data <- do.call(h$action, list(h$file, header=T))
            assign(data_frame_name, the_data, envir = globalenv())
            },
            error = function(e){
              svalue(status_bar) <- "Could not load data"
              add(grp_means_sds,gimage("red x.png"))
            }
          )
        }
      )

      # sample_size_col <- match("n",colnames(S_MU))      # Find the column titled "n" be in the source means file, where n is the sample size for each source
      # if(is.na(sample_size_col)){
      #   stop(paste("*** Error: Source sample sizes missing.
      #   Check your source_means.csv data file to be sure you have a column
      #   titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
      # }
      # S_sample_size <- S_MU[,sample_size_col]         # Save the sample sizes
      # if(length(S_sample_size) != (n.sources*source_factor_levels)){
      #   stop(paste("*** Error: Source sample sizes missing or entered incorrectly.
      #   Check your source_means.csv data file to be sure you have a column
      #   titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
      # }

      # Error check: zero SD
      # if(length(which(S_SIG==0))>0){
      #   stop(paste("*** Error: Zero standard deviation.
      #   Check your source_sds.csv data file to be sure each entry is non-zero.",sep=""))
      # }

      source_random_effects <<- character(0)
      if(n.re > 0){    # if we have random effects in the mix data
        for(i in 1:n.re){  # for each random effect
          if(svalue(re_list[[i]])=="Yes"){   # if random_effect[i] is selected, add it to source_random_effects
            source_random_effects <<- c(source_random_effects,random_effects[i])  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
          }
        }
      }
      if(length(source_random_effects)==0) sources_by_factor <<- FALSE else sources_by_factor <<- TRUE  # do we have any source random effects?

      if(!sources_by_factor){
        SOURCE <- SOURCE[order(SOURCE[,1]),] # SOURCE is the source data file - sort it by the source names (alphabetically)
        
        # Concentration Dependence section (must be done after SOURCE is sorted but before SOURCE[,1] is deleted)
        if(svalue(conc_rad)=="Yes") include_conc <- TRUE else include_conc <- FALSE   # do we have concentration dependence data?
        if(include_conc){   # if we have concentration dependence data
            CONC_names <- paste("Conc",iso_names,sep="")
            CONC_iso_cols <- match(CONC_names,colnames(SOURCE))
            conc <- do.call(rbind,lapply(split(SOURCE[,CONC_iso_cols],list(SOURCE[,1])),colMeans))  # calculate conc means for each [src,iso]
            assign("conc", conc, envir=.GlobalEnv)
        }
        assign("include_conc", include_conc, envir=.GlobalEnv)
        
        row.names(SOURCE) <- SOURCE[,1]    # temporarily store the source names
        SOURCE <- as.matrix(SOURCE[-1])    # remove the column of source names from SOURCE

        sample_size_col <- match("n",colnames(SOURCE))    # Find the column titled "n" be in the source means file, where n is the sample size for each source
        if(is.na(sample_size_col)){
          stop(paste("*** Error: Source sample sizes missing or entered incorrectly.
    Check your sources.csv data file to be sure you have a column
    titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
        }
        S_sample_size <- SOURCE[,sample_size_col]         # Get the sample sizes
        n.sources <- length(row.names(SOURCE))            # n.sources is the number of sources, which is the length of the source names vector

        S_MU_iso_cols <- match(MU_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes
        S_SIG_iso_cols <- match(SIG_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes
        source_factor_cols <- match(source_random_effects,colnames(SOURCE)) # get the S_MU column number(s) of the user-selected source random effect(s)
        S_MU <- as.matrix(SOURCE[,c(S_MU_iso_cols[],source_factor_cols)]) # rearrange S_MU columns to be the selected isotopes (in order) and source random effects
        S_SIG <- as.matrix(SOURCE[,c(S_SIG_iso_cols[],source_factor_cols)]) # rearrange S_SIG columns to be the selected isotopes (in order) and source random effects
#        S_iso_cols <- match(iso_names,colnames(S_MU))             # We just rearranged the S_MU columns, so re-calculate where the S_iso_cols are

        source_names <- rownames(SOURCE)   # source_names is a vector of the SORTED source names

        MU_array <- S_MU                      # MU_array is passed to JAGS
        SIG2_array <- S_SIG*S_SIG             # SIG2_array is passed to JAGS
        n_array <- S_sample_size              # n_array is passed to JAGS
      } # end if(sources_by_factor)==FALSE
      if(sources_by_factor){
        source_factor_cols <- match(source_random_effects,colnames(SOURCE))   # get the S_MU column number(s) of the user-selected source random effect(s)
        source_factor_levels <- rep(0,length(source_factor_cols))                # get the number of levels of each source random effect
        for(fac in 1:length(source_factor_cols)){
          source_factor_levels[fac] <- length(levels(as.factor(SOURCE[,source_factor_cols[fac]])))
        }
        assign("source_factor_levels", source_factor_levels, envir = .GlobalEnv)

        SOURCE <- SOURCE[order(SOURCE[,1],SOURCE[,source_factor_cols]),]     # sort SOURCE by the source names and by factor(s) within each source
        
        # Concentration Dependence section (must be done after SOURCE is sorted but before SOURCE[,1] is deleted)
        if(svalue(conc_rad)=="Yes") include_conc <- TRUE else include_conc <- FALSE   # do we have concentration dependence data?
        if(include_conc){   # if we have concentration dependence data
            CONC_names <- paste("Conc",iso_names,sep="")
            CONC_iso_cols <- match(CONC_names,colnames(SOURCE))
            conc <- do.call(rbind,lapply(split(SOURCE[,CONC_iso_cols],list(SOURCE[,1])),colMeans))  # calculate conc means for each [src,iso]
            assign("conc", conc, envir=.GlobalEnv)
        }
        assign("include_conc", include_conc, envir=.GlobalEnv)

        source_names <- levels(unique(SOURCE[,1]))          # source_names is a vector of the SORTED source names (as is, the first column of S_MU and S_SIG MUST contain the source names)
        n.sources <- length(source_names)                   # get the number of sources
        SOURCE <- SOURCE[-1]                                # remove the column of source names from SOURCE
        
        sample_size_col <- match("n",colnames(SOURCE))      # Find the column titled "n" be in the source means file, where n is the sample size for each source
        if(is.na(sample_size_col)){
          stop(paste("*** Error: Source sample sizes missing.
          Check your source_means.csv data file to be sure you have a column
          titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
        }
        S_sample_size <- SOURCE[,sample_size_col]         # Save the sample sizes
        if(length(S_sample_size) != (n.sources*source_factor_levels)){
          stop(paste("*** Error: Source sample sizes missing or entered incorrectly.
          Check your source_means.csv data file to be sure you have a column
          titled \"n\" with the sample sizes for each source isotope estimate.",sep=""))
        }
        S_MU_iso_cols <- match(MU_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes
        S_SIG_iso_cols <- match(SIG_names,colnames(SOURCE))             # get the S_MU column numbers of the user-selected isotopes
        source_factor_cols <- match(source_random_effects,colnames(SOURCE)) # get the S_MU column number(s) of the user-selected source random effect(s)
        S_factor1 <- factor(SOURCE[,source_factor_cols])               # save the Factor1 labels for the isospace plot
        SOURCE[,source_factor_cols] <- as.numeric(factor(SOURCE[,source_factor_cols]))  # turn the Factor1 labels into numbers (already have them saved as 'factor1_levels' from X)
        S_MU <- as.matrix(SOURCE[,c(S_MU_iso_cols[],source_factor_cols)]) # rearrange S_MU columns to be the selected isotopes (in order) and source random effects
        S_SIG <- as.matrix(SOURCE[,c(S_SIG_iso_cols[],source_factor_cols)]) # rearrange S_SIG columns to be the selected isotopes (in order) and source random effects
#        S_iso_cols <- match(iso_names,colnames(S_MU))             # We just rearranged the S_MU columns, so re-calculate where the S_iso_cols are
#        source_factor_cols <- match(source_random_effects,colnames(S_MU))  # same for the source_factor_cols
        
        assign("S_factor1", S_factor1, envir = .GlobalEnv)

        MU_array <- array(NA,dim=c(n.sources,n.iso,source_factor_levels))
        SIG2_array <- array(NA,dim=c(n.sources,n.iso,source_factor_levels))
        n_array <- array(NA,dim=c(n.sources,source_factor_levels))
        count <- 1
        for(src in 1:n.sources){
          for(f1 in 1:source_factor_levels){
            for(iso in 1:n.iso){
              MU_array[src,iso,f1] <- S_MU[count,iso]      # MU_array is the source mean data, but indexed as [src,iso,f1]. Passed to JAGS
              SIG2_array[src,iso,f1] <- S_SIG[count,iso]*S_SIG[count,iso]   # same for SIG2_array, the source sd data
            }
            n_array[src,f1] <- S_sample_size[count]        # n_array is the source sample sizes but indexed as [src,f1].  Passed to JAGS
            count <- count+1
          }
        }
      } # end if(sources_by_factor)==TRUE

      assign("n.sources", n.sources, envir = .GlobalEnv)
      assign("source_names", source_names, envir = .GlobalEnv)
#      assign("S_iso_cols", S_iso_cols, envir = .GlobalEnv)
#      assign("source_factor_cols", source_factor_cols, envir = .GlobalEnv)
      assign("S_MU", S_MU, envir = .GlobalEnv)
      assign("S_SIG", S_SIG, envir = .GlobalEnv)
      assign("MU_array", MU_array, envir = .GlobalEnv)
      assign("SIG2_array", SIG2_array, envir = .GlobalEnv)
      assign("n_array", n_array, envir = .GlobalEnv)
      raw_source_data <<- FALSE
      add(grp_means_sds,gimage("check.png"))
    } # end 'load source means and sds' loop
  )
  addSpring(grp_btns)

  addSpring(source_grp_all)
  btn_close_source <- gbutton(
    text		= "I'm finished",
    container	= source_grp_all,
    expand = F,
    handler	= function(h, ...){
      visible(source_win) <- FALSE

    # need to check that the important objects are correct:
    # S_MU, S_SIG, SOURCE_array (or MU_array and SIG_array)

      if(exists("S_MU") && exists("S_SIG")){
        add(grp_source,gimage("check.png"))
        svalue(status_bar) <- "Source data successfully loaded"
      } else {
          svalue(status_bar) <- "Could not load source data"
          add(grp_source,gimage("red x.png"))
      }
    }
  )
  addSpring(source_grp_all)
  visible(source_win) <- TRUE
}

