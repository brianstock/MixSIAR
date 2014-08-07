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
  mix <- mixsiar$mix
  if(!exists("mix")){
  stop(paste("*** Error: Consumer/mixture data not yet loaded. First load your 
      consumer/mix data, then try to load source data again. ***",sep=""))}
  n.re <- mix$n.re
  n.effects <- mix$n.effects
  random_effects <- mix$random_effects
  fixed_effects <- mix$fixed_effects

  source_win <- gwindow("Read in your SOURCE data")
  source_grp_all <- ggroup(horizontal=F, container=source_win)
  if(n.effects > 0){    # if we have factors in the mix data
    re_list <- vector("list",n.re); assign("re_list",re_list,envir=mixsiar);
    for(i in 1:n.effects){      # for each fixed/random effect
      cur <- mix$FAC[[i]]$name
      grp <- ggroup(cont=source_grp_all,horizontal=T)
      lbl <- glabel(
          paste("Does your source data vary by ",cur,"?",sep=""),
          cont=grp)
      mixsiar$re_list[[i]] <- gradio(c("Yes","No"), cont=grp, horizontal=T)
      svalue(mixsiar$re_list[[i]]) <- "No"
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
            assign(data_frame_name, the_data, envir = mixsiar)
            add(grp_raw,gimage("check.png"))
            source_filename <- h$file; assign("source_filename",source_filename,envir=mixsiar)
            },
            error = function(e){
              svalue(mixsiar$status_bar) <- "Could not load data"
              add(grp_raw,gimage("red x.png"))
            }
          )
        }
      )
      source_random_effects <- character(0)  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
      if(n.effects > 0){    # if we have random/fixed effects in the mix data
        for(i in 1:n.effects){  # for each random/fixed effect
          if(svalue(mixsiar$re_list[[i]])=="Yes"){   # if random/fixed effect[i] is selected by the user, add it to source_random_effects
            source_random_effects <- c(source_random_effects,mix$FAC[[i]]$name)
          }
        }
      }
      data_type <- "raw"
      if(svalue(conc_rad)=="Yes") conc_dep <- TRUE else conc_dep <- FALSE
      source <- load_source_data(mixsiar$source_filename,source_random_effects,conc_dep,data_type,mix)
      assign("source", source, envir=mixsiar)
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
            assign(data_frame_name, the_data, envir = mixsiar)
            add(grp_means_sds,gimage("check.png"))
            source_filename <- h$file; assign("source_filename",source_filename,envir=mixsiar)
            },
            error = function(e){
              svalue(mixsiar$status_bar) <- "Could not load data"
              add(grp_means_sds,gimage("red x.png"))
            }
          )
        }
      )
      source_random_effects <- character(0)  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
      if(n.effects > 0){    # if we have random/fixed effects in the mix data
        for(i in 1:n.effects){  # for each random/fixed effect
          if(svalue(mixsiar$re_list[[i]])=="Yes"){   # if random/fixed effect[i] is selected by the user, add it to source_random_effects
            source_random_effects <- c(source_random_effects,mix$FAC[[i]]$name)
          }
        }
      }
      data_type <- "means"
      if(svalue(conc_rad)=="Yes") conc_dep <- TRUE else conc_dep <- FALSE
      source <- load_source_data(mixsiar$source_filename,source_random_effects,conc_dep,data_type,mix)
      assign("source", source, envir=mixsiar)
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
      test <- get("source",envir=mixsiar)
      if(exists("test")){
        add(mixsiar$grp_source,gimage("check.png"))
        svalue(mixsiar$status_bar) <- "Source data successfully loaded"
      } else {
          svalue(mixsiar$status_bar) <- "Could not load source data"
          add(mixsiar$grp_source,gimage("red x.png"))
      }
    }
  )
  addSpring(source_grp_all)
  visible(source_win) <- TRUE
}
