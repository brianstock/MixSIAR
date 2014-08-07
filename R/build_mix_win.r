#' Build the window to read-in MIX data
#'
#' This function is run when a user clicks the "Load mixture data" button in the main MixSIAR GUI,
#' and creates a separate GUI window (\code{mix_win}) where the user loads the mixture data. It contains:
#'  1) "Load consumer data file" gbutton which loads the .csv file into \code{X},
#'  2) gtables where the user selects which columns of \code{X} are Isotopes, Random Effects, and Continuous Effects,
#'  3) Individual Effect gcheckbox ("Include 'Individual' as a Random Effect"),
#'  4) "I'm finished" gbutton that closes \code{mix_win} and manipulates \code{X} into all the objects we use later.
#' If more than 2 Random Effects are selected, a separate WARNING gwindow prompts the user to select 2, 1, or 0.
#' If 2 Random Effects are selected, a separate gwindow asks the user if the effects are hierarchical/nested.
#' If more than 1 Continuous Effect is selected, a separate WARNING gwindow prompts the user to select 1 or 0.
#' Finally, the function adds a green check image if the data is successfully loaded, or a red x image if not.
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
indiv_box <- gcheckbox("Include 'Individual' as a Random Effect", cont=grp_bottom)
svalue(indiv_box) <- TRUE    # Default is to include 'Individual', but can change that here
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

    # Was individual effect checked?
    indiv_effect <- svalue(indiv_box)
    assign("indiv_effect", indiv_effect, envir = mixsiar)

    # Hierarchical question box: if the user has included 2 random effects, ask if the model should be hierarchical (Factor 2 within Factor 1)
    # Creates 'nested', a T/F variable
    #   hierarch=T --> ilr.fac2.tot = ilr.global + ilr.fac1 + ilr.fac2
    #   hierarch=F --> ilr.fac2.tot = ilr.global + ilr.fac2
    nested <- FALSE
    assign("nested", nested, envir = mixsiar)
    if(n.re==2){ 
      hierarch_win <- gwindow("QUESTION: Hierarchical Data?", visible=T)
      hierarch_grp_all <- ggroup(horizontal=F, cont=hierarch_win)
      hierarch_msg1 <- glabel(paste("You have 2 random effects: ",random_effects[1]," and ",random_effects[2],sep=""),cont=hierarch_grp_all)
      hierarch_msg2 <- glabel("Should MixSIAR run a hierarchical analysis?",cont=hierarch_grp_all)
      h_yes <- paste("Yes (",random_effects[2]," within ",random_effects[1],")",sep="")
      h_no <- paste("No (",random_effects[1],", ",random_effects[2]," independent)",sep="")
      hierarch_rad <- gradio(c(h_yes,h_no),cont=hierarch_grp_all,horizontal=F)
      svalue(hierarch_rad) <- h_yes
      btn_done_hierarch <- gbutton(
        text    = "I'm finished",
        container = hierarch_grp_all,
        expand = F,
        handler = function(h, ...){
          nested <- mixsiar$nested
          if(svalue(hierarch_rad)==h_yes){nested <- TRUE} else {nested <- FALSE}
          visible(hierarch_win) <- FALSE
          mixsiar$nested <- nested
        }
      )
      addSpring(hierarch_grp_all)
    }
    
    mix <- load_mix_data(mixsiar$mix_filename,iso_names,random_effects,cont_effects,fixed_effects)
    assign("mix", mix, envir = mixsiar)

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

