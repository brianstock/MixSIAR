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
          assign(data_frame_name, the_data, envir = globalenv())
          # Add a real check for X
          add(mix_grp_cons,gimage("check.png"))
          svalue(mix_status_bar) <- "Consumer data file successfully loaded"
          tbl_data_col[] <- colnames(X)
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
indiv_effect <- gcheckbox("Include 'Individual' as a Random Effect", cont=grp_bottom)
svalue(indiv_effect) <- TRUE    # Default is to include 'Individual'. If we comment this line, default is to NOT include 'Individual'
mix_status_bar <- gstatusbar("", progress.bar="gui", container=grp_bottom, expand=T)

grp_close <- ggroup(horizontal=T, container=mix_grp_all)
btn_close <- gbutton(
  text		= "I'm finished",
  container	= grp_close,
  expand = T,
  handler	= function(h, ...){
    visible(mix_win) <- FALSE # hide the "Read in your CONSUMER data" window

    # Read in the selected Isotopes
    iso_names <- tbl_iso[]
    n.iso <- length(iso_names)
    MU_names <- paste("Mean",iso_names,sep="")
    SIG_names <- paste("SD",iso_names,sep="")
    assign("iso_names", iso_names, envir = .GlobalEnv)
    assign("n.iso", n.iso, envir = .GlobalEnv)
    assign("MU_names", MU_names, envir = .GlobalEnv)
    assign("SIG_names", SIG_names, envir = .GlobalEnv)

    # Read in the selected Random Effects
    random_effects <- tbl_re[]
    n.re <- length(random_effects)
    assign("random_effects", random_effects, envir = .GlobalEnv)
    assign("n.re", n.re, envir = .GlobalEnv)
    if(n.re > 2){ # if the user has included more than 2 random effects, show a warning message and ask which 2 random effects to use
      re_message_win <- gwindow("WARNING: More than 2 Random Effects selected", visible=T)
      re_grp_all <- ggroup(horizontal=F, cont=re_message_win)
      re_msg1 <- glabel("Note: The MixSIAR model can only analyze up to 2 Random Effects.",cont=re_grp_all)
      re_msg2 <- glabel("Please select either 0, 1, or 2 Random Effects to include in the analysis:",cont=re_grp_all)

      grp_choose_re <- ggroup(horizontal=T,cont=re_grp_all, expand=T)
      re_old <- gtable(random_effects, container=grp_choose_re, multiple=T, expand=T)
      grp_new_re_btns <- ggroup(container=grp_choose_re, horizontal=F)
      add_re_new <- gbutton(
                    text = ">>",
                    container = grp_new_re_btns,
                    handler = function(h,...){
                      re_new[] <- c(re_new[], svalue(re_old))
                      re_old[] <- re_old[which(re_old[] != svalue(re_old))]
                    }
                  )
      remove_re_new <- gbutton(
                      text = "<<",
                      container = grp_new_re_btns,
                      handler = function(h,...){
                        re_old[] <- c(re_old[][], svalue(re_new))
                        re_new[] <- re_new[which(re_new[] != svalue(re_new))]
                      }
                    )
      re_new <- gtable(character(0), container=grp_choose_re, multiple=T, expand=T)
      btn_done_re <- gbutton(
        text		= "I'm finished",
        container	= re_grp_all,
        expand = T,
        handler	= function(h, ...){
          random_effects <- re_new[]
          n.re <- length(random_effects)
          assign("random_effects", random_effects, envir = .GlobalEnv)
          assign("n.re", n.re, envir = .GlobalEnv)
          visible(re_message_win) <- FALSE
          if(n.re > 2){
            visible(re_message_win) <- TRUE # If they selected more than 2 Random Effects AGAIN, show them this window AGAIN
          }
        }
      )
    } # End >2 RE warning box

    # Hierarchical question box: if the user has included 2 random effects, ask if the model should be hierarchical (Factor 2 within Factor 1)
    # Creates 'hierarch', a T/F variable
    #   hierarch=T --> ilr.fac2.tot = ilr.global + ilr.fac1 + ilr.fac2
    #   hierarch=F --> ilr.fac2.tot = ilr.global + ilr.fac2
    hierarch <- FALSE
    assign("hierarch", hierarch, envir = .GlobalEnv)
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
          if(svalue(hierarch_rad)==h_yes){hierarch <- TRUE} else {hierarch <- FALSE}
          visible(hierarch_win) <- FALSE
          assign("hierarch", hierarch, envir = .GlobalEnv)
        }
      )
      addSpring(hierarch_grp_all)
    }
    
    # X is the consumer data
    # Sort X by the random effects (needed now?)
    #if(n.re==1) X <- X[order(X[,random_effects[1]]),]                         # Sort X by Factor.1
    #if(n.re==2) X <- X[order(X[,random_effects[1]],X[,random_effects[2]]),]   # Sort X by Factor.1 and then Factor.2 within Factor.1
    N <- dim(X)[1]       # N is the number of consumer data points                                                    
    X_iso_cols <- match(iso_names,colnames(X))   # find the column indicies of the user-selected isotopes in X
    X_iso <- as.matrix(X[,X_iso_cols[]])        # keep the original X but create 'X_iso' to pass JAGS that only has the selected consumer isotope values (in the order of selection)
    assign("N", N, envir = .GlobalEnv)
    assign("X_iso_cols", X_iso_cols, envir = .GlobalEnv)
    assign("X_iso", X_iso, envir = .GlobalEnv)
    if(n.re > 0){
      Factor.1 <- X[,random_effects[1]]
      factor1_levels <- length(unique(Factor.1))
      # get the factor 1 level names (e.g. 'Mainland', 'Inner Island', and 'Outer Island' wolves)
      if(is.numeric(Factor.1)){ # if factor 1 was input as numbers, add the column name to each number (e.g. "Region 1", "Region 2", "Region 3")
        factor1_names <- paste(rep(random_effects[1],factor1_levels),levels(factor(Factor.1)),sep=" ")
      } else {  # if factor 1 was input as text names, save them as is
        factor1_names <- levels(factor(Factor.1))
      }
      Factor.1 <- as.numeric(factor(Factor.1))
      assign("Factor.1", Factor.1, envir = .GlobalEnv)
      assign("factor1_levels", factor1_levels, envir = .GlobalEnv)
      assign("factor1_names", factor1_names, envir = .GlobalEnv)
      if(n.re > 1){
        Factor.2 <- X[,random_effects[2]]
        factor2_levels <- length(unique(Factor.2))
        #factor2_levels <- rep(0,factor1_levels)
        #for(f1 in 1:factor1_levels){                                                      # for each level of Factor.1 ...
        #  factor2_levels[f1] <- length(unique(X[which(Factor.1==f1),random_effects[2]]))  # get the levels of Factor.2 within that level of Factor.1 (e.g. get the # of Packs within Regions 1,2, and 3)
        #}
        #factor2_names <- vector("list", factor1_levels)   # get the factor 2 level names (e.g. 'Mainland', 'Inner Island', and 'Outer Island' wolves)
        # for(f1 in 1:factor1_levels){
        #   if(is.numeric(Factor.2)){ # if factor 2 was input as numbers, add the column name to each number (e.g. "Pack 1", "Pack 2", "Pack 3")
        #     #changed 6/26...factor2_names[[f1]] <- paste(rep(random_effects[2],factor2_levels[f1]),levels(factor(Factor.2)),sep=" ")
        #     factor2_names[[f1]] <- paste(rep(random_effects[2],factor2_levels[f1]),levels(factor(Factor.2[which(Factor.1==f1)])),sep=" ")
        #     tmp_fac <- factor(Factor.2[which(Factor.1==f1)])
        #     levels(tmp_fac) <- 1:factor2_levels[f1]
        #     Factor.2[which(Factor.1==f1)] <- tmp_fac
        #   } else {  # if factor 2 was input as text names, save them as is and turn Factor.2 into numbers
        #     factor2_names[[f1]] <- levels(factor(Factor.2))
        #     Factor.2 <- as.numeric(factor(Factor.2))
        #   } 
        # }
        if(is.numeric(Factor.2)){ # if factor 1 was input as numbers, add the column name to each number (e.g. "Region 1", "Region 2", "Region 3")
          factor2_names <- paste(rep(random_effects[2],factor2_levels),levels(factor(Factor.2)),sep=" ")
        } else {  # if factor 1 was input as text names, save them as is
          factor2_names <- levels(factor(Factor.2))
        }
        Factor.2 <- as.numeric(factor(Factor.2))
        factor1_lookup <- 1:factor2_levels
        for(f2 in 1:factor2_levels){
          factor1_lookup[f2] <- Factor.1[which(Factor.2==f2)][1]
        }
       # factor2_empty <- rep(max(factor2_levels),length(factor2_levels)) - factor2_levels                # 
       # assign("factor2_empty", factor2_empty, envir = .GlobalEnv)
        assign("factor2_levels", factor2_levels, envir = .GlobalEnv)
        assign("Factor.2", Factor.2, envir = .GlobalEnv)
        assign("factor2_names", factor2_names, envir = .GlobalEnv)
        assign("factor1_lookup", factor1_lookup, envir = .GlobalEnv)
      }
    }
    include_indiv <- svalue(indiv_effect)
    assign("include_indiv", include_indiv, envir = .GlobalEnv)

    # Read in the selected Continuous Effects
    cont_effects <- tbl_cont[] # Vector of Continuous Effects
    n.ce <- length(cont_effects) # number of continuous effects
    assign("cont_effects", cont_effects, envir = .GlobalEnv)
    assign("n.ce", n.ce, envir = .GlobalEnv)
    # More than 1 Continuous Effect warning box
    if(n.ce > 1){ # if the user has included more than 1 continuous effect, show a warning message and ask which 1 to use
      ce_message_win <- gwindow("WARNING: More than 1 Continuous Effect selected", visible=T)
      ce_grp_all <- ggroup(horizontal=F, cont=ce_message_win)
      ce_msg1 <- glabel("Note: The MixSIAR model can only analyze up to 1 Continuous Effect.",cont=ce_grp_all)
      ce_msg2 <- glabel("Please select either 0 or 1 Continuous Effect to include in the analysis:",cont=ce_grp_all)

      grp_choose_ce <- ggroup(horizontal=T,cont=ce_grp_all, expand=T)
      ce_old <- gtable(cont_effects, container=grp_choose_ce, multiple=T, expand=T)
      grp_new_ce_btns <- ggroup(container=grp_choose_ce, horizontal=F)
      add_ce_new <- gbutton(
                    text = ">>",
                    container = grp_new_ce_btns,
                    handler = function(h,...){
                      ce_new[] <- c(ce_new[], svalue(ce_old))
                      ce_old[] <- ce_old[which(ce_old[] != svalue(ce_old))]
                    }
                  )
      remove_ce_new <- gbutton(
                      text = "<<",
                      container = grp_new_ce_btns,
                      handler = function(h,...){
                        ce_old[] <- c(ce_old[][], svalue(ce_new))
                        ce_new[] <- ce_new[which(ce_new[] != svalue(ce_new))]
                      }
                    )
      ce_new <- gtable(character(0), container=grp_choose_ce, multiple=T, expand=T)
      btn_done_ce <- gbutton(
        text    = "I'm finished",
        container = ce_grp_all,
        expand = T,
        handler = function(h, ...){
          cont_effects <- ce_new[]
          n.ce <- length(cont_effects)
          assign("cont_effects", cont_effects, envir = .GlobalEnv)
          assign("n.ce", n.ce, envir = .GlobalEnv)
          visible(ce_message_win) <- FALSE
          if(n.ce > 1){
            visible(ce_message_win) <- TRUE # If they selected more than 2 Random Effects AGAIN, show them this window AGAIN
          }
        }
      )
    } # End >1 CE warning box
    if(n.ce > 0){ # If we have any continuous effects
      for(i in 1:n.ce){ # For each continuous effect selected, get the data from X, label it, and assign it globally
        cont.name <- paste("Cont.",i,sep="")  # Call each of the continuous effects vectors "Cont.1", "Cont.2", etc.
        cont.vec <- X[,cont_effects[i]]       # X is temporarily sorted same as X_iso (by fac1, fac2), so we can plot Cont.1 vs. p.ind. X reverts to originally entered form after function ends.
        assign(cont.name, cont.vec, envir = .GlobalEnv)
      }
    }

    # Need to make this check more robust
    if(exists("X_iso")){
      add(grp_cons,gimage("check.png"))
      svalue(status_bar) <- "Mixture data successfully loaded"
    } else {
        svalue(status_bar) <- "Could not load mixture data"
        add(grp_cons,gimage("red x.png"))
    }
  } # end "I'm finished" button handler/actions
) # end "I'm finished" button
visible(mix_win) <- TRUE
} # end function build_mix_win

