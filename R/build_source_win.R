#' Build the window to read-in source data
#'
#' This function is run when a user clicks the "Load source data" button in the
#' main MixSIAR GUI, and creates a separate GUI window (\code{source_win}) where
#' the user loads the source data. \code{build_source_win} calls
#' \code{\link{load_source_data}} when the user clicks the "I'm finished" button.
#'
#' The "Load source data" window contains:
#' \enumerate{
#'   \item Questions the user if their source data is "by factor" (if they have Fixed or Random Effects),
#'   \item Questions the user if they have Concentration Dependence data,
#'   \item Buttons to either load "raw" or "means + SD" source data.
#' }
#'
#' If \code{mix} does not yet exist, a WARNING appears.
#'
#' A green check appears if the source data is successfully loaded, or a red_x
#' if not.
#'
#' @seealso \code{\link{load_source_data}}, which is run when the "Load source data"
#'          window is closed by clicking the "I'm finished" button.
build_source_win <- function(){
  mix <- mixsiar$mix
  if(!exists("mix")){
  stop(paste("*** Error: Consumer/mixture data not yet loaded. First load your
      consumer/mix data, then try to load source data again. ***",sep=""))}
  n.re <- mix$n.re
  n.effects <- mix$n.effects
  random_effects <- mix$random_effects
  fixed_effects <- mix$fixed_effects

  source_win <- gWidgets::gwindow("Read in your SOURCE data")
  source_grp_all <- gWidgets::ggroup(horizontal=F, container=source_win)
  if(n.effects > 0){    # if we have factors in the mix data
    re_list <- vector("list",n.re); assign("re_list",re_list,envir=mixsiar);
    for(i in 1:n.effects){      # for each fixed/random effect
      cur <- mix$FAC[[i]]$name
      grp <- gWidgets::ggroup(cont=source_grp_all,horizontal=T)
      lbl <- gWidgets::glabel(
          paste("Does your source data vary by ",cur,"?",sep=""),
          cont=grp)
      mixsiar$re_list[[i]] <- gWidgets::gradio(c("Yes","No"), cont=grp, horizontal=T)
      gWidgets::svalue(mixsiar$re_list[[i]]) <- "No"
    }
  }
  grp_conc <- gWidgets::ggroup(cont=source_grp_all,horizontal=T)
  lbl_conc <- gWidgets::glabel("Do you have Concentration Dependence data?",cont=grp_conc)
  conc_rad <- gWidgets::gradio(c("Yes","No"),cont=grp_conc,horizontal=T)
  gWidgets::svalue(conc_rad) <- "No"

  gWidgets::addSpring(source_grp_all)
  lbl_data_type <- gWidgets::glabel("Do you have raw source data, or source means and SDs?",
                        cont=source_grp_all)
  grp_btns <- gWidgets::ggroup(cont=source_grp_all, horizontal=T)
  gWidgets::addSpring(grp_btns)
  grp_raw <- gWidgets::ggroup(cont=grp_btns, horizontal=T)
  raw_data <- gWidgets::gbutton(
    text = "Load raw source data",
    cont = grp_raw,
    expand = T,
    handler = function(h, ...){
      gWidgets::gfile(
        text = "Load raw source data file",
        type = "open",
        action = "read.csv",
        handler = function(h, ...){
          tryCatch(   # reads source data into 'SOURCE'
            {
            data_frame_name <- make.names("SOURCE")
            the_data <- do.call(h$action, list(h$file, header=T))
            assign(data_frame_name, the_data, envir = mixsiar)
            gWidgets::add(grp_raw,gWidgets::gimage(system.file("extdata", "check.png", package = "MixSIAR")))
            source_filename <- h$file; assign("source_filename",source_filename,envir=mixsiar)
            },
            error = function(e){
              gWidgets::svalue(mixsiar$status_bar) <- "Could not load data"
              gWidgets::add(grp_raw,gWidgets::gimage(system.file("extdata", "red_x.png", package = "MixSIAR")))
            }
          )
        }
      )
      source_random_effects <- character(0)  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
      if(n.effects > 0){    # if we have random/fixed effects in the mix data
        for(i in 1:n.effects){  # for each random/fixed effect
          if(gWidgets::svalue(mixsiar$re_list[[i]])=="Yes"){   # if random/fixed effect[i] is selected by the user, add it to source_random_effects
            source_random_effects <- c(source_random_effects,mix$FAC[[i]]$name)
          }
        }
      }
      data_type <- "raw"
      if(gWidgets::svalue(conc_rad)=="Yes") conc_dep <- TRUE else conc_dep <- FALSE
      source <- load_source_data(mixsiar$source_filename,source_random_effects,conc_dep,data_type,mix)
      assign("source", source, envir=mixsiar)
    }
  )

  lbl_or <- gWidgets::glabel(" OR ", cont=grp_btns)

  grp_means_sds <- gWidgets::ggroup(cont=grp_btns, horizontal=T)
  means_sds <- gWidgets::gbutton(
    text = "Load source means and SDs",
    cont = grp_means_sds,
    expand = T,
    handler	= function(h, ...){
      gWidgets::gfile(
        text = "Load source data file",
        type = "open",
        action = "read.csv",
        handler = function(h, ...){
          tryCatch(   # reads source data into 'SOURCE'
            {
            data_frame_name <- make.names("SOURCE")
            the_data <- do.call(h$action, list(h$file, header=T))
            assign(data_frame_name, the_data, envir = mixsiar)
            gWidgets::add(grp_means_sds,gWidgets::gimage(system.file("extdata", "check.png", package = "MixSIAR")))
            source_filename <- h$file; assign("source_filename",source_filename,envir=mixsiar)
            },
            error = function(e){
              gWidgets::svalue(mixsiar$status_bar) <- "Could not load data"
              gWidgets::add(grp_means_sds,gWidgets::gimage(system.file("extdata", "red_x.png", package = "MixSIAR")))
            }
          )
        }
      )
      source_random_effects <- character(0)  # vector of the names of any source random effects, e.g. "Region" or "character(0)" or ("Region", "Pack")
      if(n.effects > 0){    # if we have random/fixed effects in the mix data
        for(i in 1:n.effects){  # for each random/fixed effect
          if(gWidgets::svalue(mixsiar$re_list[[i]])=="Yes"){   # if random/fixed effect[i] is selected by the user, add it to source_random_effects
            source_random_effects <- c(source_random_effects,mix$FAC[[i]]$name)
          }
        }
      }
      data_type <- "means"
      if(gWidgets::svalue(conc_rad)=="Yes") conc_dep <- TRUE else conc_dep <- FALSE
      source <- load_source_data(mixsiar$source_filename,source_random_effects,conc_dep,data_type,mix)
      assign("source", source, envir=mixsiar)
    } # end 'load source means and sds' loop
  )
  gWidgets::addSpring(grp_btns)

  gWidgets::addSpring(source_grp_all)
  btn_close_source <- gWidgets::gbutton(
    text		= "I'm finished",
    container	= source_grp_all,
    expand = F,
    handler	= function(h, ...){
      gWidgets::visible(source_win) <- FALSE
      test <- get("source",envir=mixsiar)
      if(exists("test")){
        gWidgets::add(mixsiar$grp_source,gWidgets::gimage(system.file("extdata", "check.png", package = "MixSIAR")))
        gWidgets::svalue(mixsiar$status_bar) <- "Source data successfully loaded"
      } else {
        gWidgets::svalue(mixsiar$status_bar) <- "Could not load source data"
        gWidgets::add(mixsiar$grp_source,gWidgets::gimage(system.file("extdata", "red_x.png", package = "MixSIAR")))
      }
    }
  )
  gWidgets::addSpring(source_grp_all)
  gWidgets::visible(source_win) <- TRUE
}
