#' Build the window to read-in mixture data.
#'
#' This function is run when a user clicks the "Load mixture data" button in the
#' main MixSIAR GUI, and creates a separate GUI window (\code{mix_win}) where
#' the user loads the mixture data file. \code{build_mix_win} calls
#' \code{\link{load_mix_data}} when the user clicks the "I'm finished" button.
#'
#' The "Load mix data" window contains:
#' \enumerate{
#'   \item "Load consumer data file" gbutton which loads the .csv file into \code{X},
#'   \item gtables where the user selects which columns of \code{X} are Isotopes
#'         and Fixed/Random/Continuous Effects.
#'   \item (formerly) Individual Effect gcheckbox ("Include 'Individual' as a Random Effect"),
#'   \item "I'm finished" gbutton that closes \code{mix_win} and calls
#'         \code{\link{load_mix_data}}.
#' }
#'
#' If more than 2 Fixed/Random Effects are selected, a separate WARNING gwindow
#' prompts the user to select 2, 1, or 0.
#'
#' If more than 1 Continuous Effect is selected, a separate WARNING gwindow
#' prompts the user to select 1 or 0.
#'
#' If 2 Fixed/Random Effects are selected, a separate gwindow asks the user if
#' the effects are hierarchical/nested.
#'
#' Finally, the function adds a green check image if the data is successfully
#' loaded, or a red x image if not.
#'
#' @seealso \code{\link{load_mix_data}}, which is run when the "Load MIX data"
#'          window is closed by clicking the "I'm finished" button.
build_mix_win <- function(){
mix_win <- gwindow("Read in your MIXTURE data",visible=FALSE) # set visible=F here so that it shows up all at once when the function is finished building the window (visible=T at end of function)
mix_grp_all <- ggroup(horizontal=F, container=mix_win)
mix_grp_cons <- ggroup(horizontal=T, container=mix_grp_all)
mix_btn_cons <- gbutton(
  text		= "Load mixture data file",
  container	= mix_grp_cons,
  expand = T,
  handler	= function(h, ...)
  {
    gfile(
      text	= "Load mixture data file",
      type	= "open",
      action= "read.csv",
      handler= function(h, ...)
      {
        tryCatch(
        {
          data_frame_name <- make.names("X")
          the_data <- do.call(h$action, list(h$file, header=T))
          assign(data_frame_name, the_data, envir = mixsiar)
          # Add a real check for X
          add(mix_grp_cons,gimage("check.png"))
          svalue(mix_status_bar) <- "Consumer data file successfully loaded"
          tbl_data_col[] <- colnames(mixsiar$X)
          mix_filename <- h$file; assign("mix_filename",mix_filename,envir=mixsiar)
        },
        error = function(e){
          svalue(mix_status_bar) <- "Could not load data"
          add(mix_grp_cons,gimage("red x.png"))
        }
        )
      }
    )
  }
)

grp_middle <- ggroup(container=mix_grp_all, horizontal=T, expand=T)
grp_data_col <- gframe(text="Data Columns", container=grp_middle, expand=T)
tbl_data_col <- gtable(character(0), container=grp_data_col, multiple=T, expand=T)

grp_right <- ggroup(container=grp_middle, horizontal=F, expand=T)
grp_iso <- gframe(text="Isotopes", container=grp_right, horizontal=T, expand=T)
grp_iso_btns <- ggroup(container=grp_iso, horizontal=F)
add_iso <- gbutton(
              text = ">>",
              container = grp_iso_btns,
              handler = function(h,...){
                tbl_iso[] <- c(tbl_iso[], svalue(tbl_data_col))
                tbl_data_col[] <- tbl_data_col[which(tbl_data_col[] != svalue(tbl_data_col))]
              }
            )
remove_iso <- gbutton(
                text = "<<",
                container = grp_iso_btns,
                handler = function(h,...){
                  tbl_data_col[] <- c(tbl_data_col[], svalue(tbl_iso))
                  tbl_iso[] <- tbl_iso[which(tbl_iso[] != svalue(tbl_iso))]
                }
              )
tbl_iso <- gtable(character(0), container=grp_iso, multiple=T, expand=T)

grp_re <- gframe(text="Random Effects", container=grp_right, horizontal=T, expand=T)
grp_re_btns <- ggroup(container=grp_re, horizontal=F)
add_re <- gbutton(
                text = ">>",
                container = grp_re_btns,
                handler = function(h,...){
                  tbl_re[] <- c(tbl_re[], svalue(tbl_data_col))
                  tbl_data_col[] <- tbl_data_col[which(tbl_data_col[] != svalue(tbl_data_col))]
                }
              )
remove_re <- gbutton(
                  text = "<<",
                  container = grp_re_btns,
                  handler = function(h,...){
                    tbl_data_col[] <- c(tbl_data_col[], svalue(tbl_re))
                    tbl_re[] <- tbl_re[which(tbl_re[] != svalue(tbl_re))]
                  }
                )
tbl_re <- gtable(character(0), container=grp_re, multiple=T, expand=T)

grp_fe <- gframe(text="Fixed Effects", container=grp_right, horizontal=T, expand=T)
grp_fe_btns <- ggroup(container=grp_fe, horizontal=F)
add_fe <- gbutton(
                text = ">>",
                container = grp_fe_btns,
                handler = function(h,...){
                  tbl_fe[] <- c(tbl_fe[], svalue(tbl_data_col))
                  tbl_data_col[] <- tbl_data_col[which(tbl_data_col[] != svalue(tbl_data_col))]
                }
              )
remove_fe <- gbutton(
                  text = "<<",
                  container = grp_fe_btns,
                  handler = function(h,...){
                    tbl_data_col[] <- c(tbl_data_col[], svalue(tbl_fe))
                    tbl_fe[] <- tbl_fe[which(tbl_fe[] != svalue(tbl_fe))]
                  }
                )
tbl_fe <- gtable(character(0), container=grp_fe, multiple=T, expand=T)

grp_cont <- gframe(text="Continuous Effects", container=grp_right, horizontal=T, expand=T)
grp_cont_btns <- ggroup(container=grp_cont, horizontal=F)
add_cont <- gbutton(
                text = ">>",
                container = grp_cont_btns,
                handler = function(h,...){
                  tbl_cont[] <- c(tbl_cont[], svalue(tbl_data_col))
                  tbl_data_col[] <- tbl_data_col[which(tbl_data_col[] != svalue(tbl_data_col))]
                }
              )
remove_cont <- gbutton(
                  text = "<<",
                  container = grp_cont_btns,
                  handler = function(h,...){
                    tbl_data_col[] <- c(tbl_data_col[], svalue(tbl_cont))
                    tbl_cont[] <- tbl_cont[which(tbl_cont[] != svalue(tbl_cont))]
                  }
                )
tbl_cont <- gtable(character(0), container=grp_cont, multiple=T, expand=T)

grp_bottom <- ggroup(container=mix_grp_all, horizontal=F)
mix_status_bar <- gstatusbar("", progress.bar="gui", container=grp_bottom, expand=T)

grp_close <- ggroup(horizontal=T, container=mix_grp_all)
btn_close <- gbutton(
  text		= "I'm finished",
  container	= grp_close,
  expand = T,
  handler	= function(h, ...){
    visible(mix_win) <- FALSE # hide the "Read in your CONSUMER data" window

    # Read in the selected Isotopes, Random Effects, and Continuous Effects
    iso_names <- tbl_iso[]
    random_effects <- tbl_re[]; n.re <- length(random_effects);
    fixed_effects <- tbl_fe[]; n.fe <- length(fixed_effects);
    cont_effects <- tbl_cont[]
    n.effects <- n.re+n.fe

    mix <- load_mix_data(mixsiar$mix_filename,iso_names,random_effects,fixed_effects,cont_effects)
    assign("mix", mix, envir = mixsiar)

    # Hierarchical question box: if the user has included 2 random/fixed effects, ask if the model should be hierarchical (Factor 2 within Factor 1)
    # Creates 'nested', a T/F variable
    #   hierarch=T --> ilr.fac2.tot = ilr.global + ilr.fac1 + ilr.fac2
    #   hierarch=F --> ilr.fac2.tot = ilr.global + ilr.fac2
    nested <- rep(FALSE,n.effects)
    mixsiar$mix$fac_nested <- nested
    assign("nested", nested, envir = mixsiar)
    if(mixsiar$mix$n.effects==2){
      hierarch_win <- gwindow("QUESTION: Hierarchical/Nested Data?", visible=T)
      hierarch_grp_all <- ggroup(horizontal=F, cont=hierarch_win)
      hierarch_msg1 <- glabel(paste("You have 2 random/fixed effects: ",mix$FAC[[1]]$name," and ",mix$FAC[[2]]$name,sep=""),cont=hierarch_grp_all)
      hierarch_msg2 <- glabel("Should MixSIAR run a hierarchical analysis?",cont=hierarch_grp_all)
      h_yes21 <- paste("Yes (",mix$FAC[[2]]$name," nested within ",mix$FAC[[1]]$name,")",sep="")
      h_yes12 <- paste("Yes (",mix$FAC[[1]]$name," nested within ",mix$FAC[[2]]$name,")",sep="")
      h_no <- paste("No (",mix$FAC[[1]]$name,", ",mix$FAC[[2]]$name," independent)",sep="")
      hierarch_rad <- gradio(c(h_yes21,h_yes12,h_no),cont=hierarch_grp_all,horizontal=F)
      svalue(hierarch_rad) <- h_no
      btn_done_hierarch <- gbutton(
        text    = "I'm finished",
        container = hierarch_grp_all,
        expand = F,
        handler = function(h, ...){
          nested <- mixsiar$nested
          if(svalue(hierarch_rad)==h_yes21) nested <- c(FALSE,TRUE)
          if(svalue(hierarch_rad)==h_yes12) nested <- c(TRUE,FALSE)
          visible(hierarch_win) <- FALSE
          mixsiar$mix$fac_nested <- nested

          if(mixsiar$mix$n.re==2 & mixsiar$mix$fac_nested[2]){
            for(lev in 1:mixsiar$mix$FAC[[2]]$levels){
              mixsiar$mix$FAC[[2]]$lookup[lev] <- mixsiar$mix$FAC[[1]]$values[which(mixsiar$mix$FAC[[2]]$values==lev)][1]
            }
          }
          if(mixsiar$mix$n.re==2 & mixsiar$mix$fac_nested[1]){
            for(lev in 1:mixsiar$mix$FAC[[1]]$levels){
              mixsiar$mix$FAC[[1]]$lookup[lev] <- mixsiar$mix$FAC[[2]]$values[which(mixsiar$mix$FAC[[1]]$values==lev)][1]
            }
          }
          if(mixsiar$mix$n.fe==1 & mixsiar$mix$n.re==1 & mixsiar$mix$fac_random[1]){ # make fac2 the random effect, fac1 the fixed effect
            tmp <- mixsiar$mix$FAC[[1]]
            mixsiar$mix$FAC[[1]] <- mixsiar$mix$FAC[[2]]
            mixsiar$mix$FAC[[2]] <- tmp
            mixsiar$mix$factors <- rev(mixsiar$mix$factors)
            mixsiar$mix$fac_random <- rev(mixsiar$mix$fac_random)
            mixsiar$mix$fac_nested <- rev(mixsiar$mix$fac_nested)
          }

        }
      )
      addSpring(hierarch_grp_all)
    }

    # Need to make this check more robust
    test <- get("mix",envir=mixsiar)
    if(exists("test")){
      add(mixsiar$grp_cons,gimage("check.png"))
      svalue(mixsiar$status_bar) <- "Mixture data successfully loaded"
    } else {
        svalue(mixsiar$status_bar) <- "Could not load mixture data"
        add(mixsiar$grp_cons,gimage("red x.png"))
    }
  } # end "I'm finished" button handler/actions
) # end "I'm finished" button
visible(mix_win) <- TRUE
} # end function build_mix_win

