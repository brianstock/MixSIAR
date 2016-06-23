#' Run the GUI version of MixSIAR
#'
#' \code{mixsiar_gui} creates the GUI version of MixSIAR.
#'
#' \emph{Before running this function}, a new user will need to install JAGS
#' and GTK+. See Section 2 of the manual for install instructions:
#' \url{https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_3.1_small.pdf}
#'
#' \code{mixsiar_gui} calls most of the other \code{MixSIAR} functions:
#'  \itemize{
#'    \item \code{\link{load_mix_data}}
#'    \item \code{\link{load_source_data}}
#'    \item \code{\link{load_discr_data}}
#'    \item \code{\link{plot_data}}
#'    \item \code{\link{plot_prior}}
#'    \item \code{\link{write_JAGS_model}}
#'    \item \code{\link{run_model}}
#'    \item \code{\link{output_JAGS}}
#'  }
#'
#' @return \code{mixsiar_gui} itself returns nothing, but using the GUI creates
#' R objects that you can access from the console:
#'  \itemize{
#'    \item \code{mixsiar$mix}: mixture data, output of \code{\link{load_mix_data}})
#'    \item \code{mixsiar$source}: source data, output of \code{\link{load_source_data}})
#'    \item \code{mixsiar$discr}: discrimination (TDF) data, output of \code{\link{load_discr_data}}
#'    \item \code{mixsiar$jags.1}: rjags model object with MCMC chains, output of \code{\link{run_model}}
#'  }
#'
#' @examples
#' mixsiar_gui()
#'
#' @seealso \code{\link{load_mix_data}} loads the mixture data file,
#' @seealso \code{\link{load_source_data}} loads the source data file,
#' @seealso \code{\link{load_discr_data}} loads the TDF data file,
#' @seealso \code{\link{plot_data}} creates an isospace plot,
#' @seealso \code{\link{plot_prior}} plots your prior and the uninformative prior,
#' @seealso \code{\link{write_JAGS_model}} creates a JAGS model file (.txt),
#' @seealso \code{\link{run_model}} sets up JAGS objects and calls JAGS to run
#'          the model,
#' @seealso \code{\link{output_JAGS}} processes the JAGS output (prints/saves
#'          diagnostics, summary statistics, and plots)
mixsiar_gui <- function(){
  # only run if RGtk2 can be loaded (suggested dependency)
  if(requireNamespace("gWidgetsRGtk2", quietly=TRUE)) {
  if(requireNamespace("gWidgets", quitely=TRUE)){
    runif(1)
    old <- options("guiToolkit"="RGtk2")
    on.exit(options(old), add = TRUE)

    win<-gWidgets::gwindow("MixSIAR GUI", visible=FALSE)
    grp_all <- gWidgets::ggroup(cont=win, horizontal=FALSE)
    grp_input <- gWidgets::ggroup(cont=grp_all, horizontal=TRUE)

    ############################################################################
    # Read in Data
    ############################################################################
    grp_readin_data <- gWidgets::gframe(text="Read in data", cont=grp_input, horizontal=F)
    grp_cons <- gWidgets::ggroup(horizontal=TRUE, cont=grp_readin_data); assign("grp_cons",grp_cons,envir=mixsiar);
    btn_cons <- gWidgets::gbutton(
      text = "Load mixture data",
      cont = mixsiar$grp_cons,
      expand = T,
      handler	= function(h, ...){
        build_mix_win()
      }
    )

    grp_source <- gWidgets::ggroup(horizontal=TRUE, cont = grp_readin_data); assign("grp_source",grp_source,envir=mixsiar)
    btn_source <- gWidgets::gbutton(
      text = "Load source data",
      cont = mixsiar$grp_source,
      expand = T,
      handler	= function(h, ...){
        build_source_win()
      }
    )

    grp_frac <- gWidgets::ggroup(horizontal=TRUE, cont = grp_readin_data)
    btn_frac <- gWidgets::gbutton(
      text = "Load discrimination data",
      cont = grp_frac,
      expand = T,
      handler	= function(h, ...){
        gWidgets::gfile(
          text = "Load discrimination data file",
          type = "open",
          action = "read.csv",
          handler = function(h, ...)
          {
            tryCatch(   # reads discrimination/fractionation/enrichment means data into 'FRAC'
              {
                data_frame_name <- make.names("FRAC")
                the_data <- do.call(h$action, list(h$file))
                assign(data_frame_name, the_data, envir = mixsiar)
                gWidgets::addSpring(grp_frac)
                gWidgets::svalue(mixsiar$status_bar) <- "Discrimination data successfully loaded"
                discr_filename <- h$file; assign("discr_filename",discr_filename,envir=mixsiar);
              },
              error = function(e){
                gWidgets::svalue(mixsiar$status_bar) <- "Could not load data"
                gWidgets::add(grp_frac,gWidgets::gimage(system.file("extdata", "red_x.png", package = "MixSIAR")))
              }
            )
          }
        )
        discr <- load_discr_data(mixsiar$discr_filename, mixsiar$mix)
        assign("discr", discr, envir = mixsiar)
        gWidgets::add(grp_frac,gWidgets::gimage(system.file("extdata", "check.png", package = "MixSIAR")))
      }
    )
    #############################################################################
    # User-specified MCMC parameters
    #############################################################################
    grp_set_mcmc <- gWidgets::gframe(text="MCMC run length", cont=grp_input, horizontal=TRUE, expand=TRUE)
    grp_mcmc <- gWidgets::ggroup(cont=grp_set_mcmc, horizontal=FALSE)
    mcmc_run <- gWidgets::gradio(c("test","very short","short","normal","long","very long"), cont=grp_mcmc, horizontal=FALSE)
    assign("mcmc_run",mcmc_run,envir=mixsiar)
    gWidgets::svalue(mixsiar$mcmc_run) <- "test"

    ####################################################################
    # Model Error Options
    ####################################################################
    grp_error_priors <- gWidgets::ggroup(cont=grp_input,horizontal=FALSE)
    grp_error <- gWidgets::gframe(text="Error structure", cont=grp_error_priors, horizontal=FALSE,expand=TRUE)
    error_option <- gWidgets::gradio(c("Resid * Process","Residual only","Process only (N=1)"), cont=grp_error, horizontal=FALSE)
    assign("error_option",error_option,envir=mixsiar)
    gWidgets::svalue(mixsiar$error_option) <- "Resid * Process"

    ####################################################################
    # Specify Prior
    ####################################################################
    grp_prior <- gWidgets::gframe(text="Specify prior", cont=grp_error_priors, horizontal=FALSE)
    prior_option <- gWidgets::gradio(c("\"Uninformative\"/Generalist","Informative"),cont=grp_prior,horizontal=FALSE)
    inf_prior <- gWidgets::gedit("", width=15, cont = grp_prior); assign("inf_prior",inf_prior,envir=mixsiar)
    assign("prior_option",prior_option,envir=mixsiar)
    gWidgets::svalue(mixsiar$prior_option) <- "\"Uninformative\"/Generalist"

    status_bar <- gWidgets::gstatusbar("", progress.bar="gui", cont=grp_all, expand=TRUE)
    assign("status_bar",status_bar,envir=mixsiar)

    ####################################################################
    # Isospace Plot
    ####################################################################

    # The 'Make isospace plot' button calls the plot_data function to make an isospace plot
    grp_plot <- gWidgets::ggroup(cont=grp_all, horizontal=T)
    plot_button <- gWidgets::gbutton(
      text = "Make isospace plot",
      cont = grp_plot,
      expand = TRUE,
      handler = function(h, ...){
        plot_data(gWidgets::svalue(mixsiar$plot_filename), gWidgets::svalue(mixsiar$plot_save_pdf), gWidgets::svalue(mixsiar$plot_save_png), mixsiar$mix,mixsiar$source,mixsiar$discr)
        gWidgets::svalue(mixsiar$status_bar) <- "If isospace plot looks good, you can now click 'RUN MODEL' at bottom"
      }
    )

    grp_plot_name <- gWidgets::ggroup(cont=grp_plot, horizontal=T, expand=T)
    plot_lbl <- gWidgets::glabel("Save plot as:",cont=grp_plot_name)
    plot_filename <- gWidgets::gedit("isospace_plot", width=15, cont = grp_plot_name); assign("plot_filename",plot_filename,envir=mixsiar)
    plot_save_pdf <- gWidgets::gcheckbox("pdf", cont = grp_plot_name); assign("plot_save_pdf",plot_save_pdf,envir=mixsiar)
    plot_save_png <- gWidgets::gcheckbox("png", cont = grp_plot_name); assign("plot_save_png",plot_save_png,envir=mixsiar)
    gWidgets::svalue(mixsiar$plot_save_pdf) <- TRUE

    ####################################################################
    # Plot Prior
    ####################################################################

    # The 'Plot prior' button calls the plot_prior function
    grp_plot_prior <- gWidgets::ggroup(cont=grp_all, horizontal=T)
    gWidgets::addSpring(grp_plot_prior)
    plot_button <- gWidgets::gbutton(
      text = "Plot prior",
      cont = grp_plot_prior,
      expand = TRUE,
      handler = function(h, ...){
        if(gWidgets::svalue(mixsiar$prior_option) == "\"Uninformative\"/Generalist"){
          alpha.prior <- rep(1,mixsiar$source$n.sources)
        } else { # prior_option = "Informative"
          alpha.prior <- eval(parse(text=gWidgets::svalue(mixsiar$inf_prior)))
        }
        if(!is.numeric(alpha.prior)){
          stop(paste("*** Error: Your prior is not a numeric vector of length(n.sources).
                     Try again or choose the uninformative prior option. For example,
                     c(1,1,1,1) is a valid (uninformative) prior for 4 sources. ***",sep=""))}
        if(length(alpha.prior) != mixsiar$source$n.sources){
          stop(paste("*** Error: Length of your prior does not match the
                     number of sources (",mixsiar$source$n.sources,"). Try again. ***",sep=""))}
        gWidgets::svalue(mixsiar$status_bar) <- "Success. Your prior is in RED, uninformative/generalist is DARK GREY"
        plot_prior(alpha.prior,mixsiar$source,gWidgets::svalue(plot_save_pdf_prior),gWidgets::svalue(plot_save_png_prior),gWidgets::svalue(plot_filename_prior))
        }
          )

    grp_plot_name_prior <- gWidgets::ggroup(cont=grp_plot_prior, horizontal=T, expand=T)
    plot_lbl_prior <- gWidgets::glabel("Save plot as:",cont=grp_plot_name_prior)
    plot_filename_prior <- gWidgets::gedit("prior_plot", width=15, cont = grp_plot_name_prior); assign("plot_filename_prior",plot_filename_prior,envir=mixsiar)
    plot_save_pdf_prior <- gWidgets::gcheckbox("pdf", cont = grp_plot_name_prior); assign("plot_save_pdf_prior",plot_save_pdf_prior,envir=mixsiar)
    plot_save_png_prior <- gWidgets::gcheckbox("png", cont = grp_plot_name_prior); assign("plot_save_png_prior",plot_save_png_prior,envir=mixsiar)
    gWidgets::svalue(mixsiar$plot_save_pdf_prior) <- TRUE

    ####################################################################
    # Output options
    ####################################################################
    grp_output <- gWidgets::gframe(text="Output options", cont=grp_all, horizontal=F)

    # Summary Statistics options
    grp_summary <- gWidgets::ggroup(cont=grp_output, horizontal=T)
    lbl_summary <- gWidgets::glabel("Summary Statistics", cont=grp_summary, expand=T)
    grp_summary_right <- gWidgets::ggroup(cont=grp_summary, horizontal=F)
    grp_summary_save <- gWidgets::ggroup(cont=grp_summary_right, horizontal=T)
    summary_save <- gWidgets::gcheckbox("Save summary statistics to file: ", cont=grp_summary_save); assign("summary_save",summary_save,envir=mixsiar);
    summary_name <- gWidgets::gedit("summary_statistics", width=20, cont=grp_summary_save); assign("summary_name",summary_name,envir=mixsiar);
    gWidgets::svalue(mixsiar$summary_save) <- TRUE  # Default is to save the summary statistics

    # Posterior Density Plot options
    grp_posterior <- gWidgets::ggroup(cont=grp_output, horizontal=T)
    lbl_posterior <- gWidgets::glabel("Posterior Density Plot", cont=grp_posterior, expand=T)
    grp_post_opt <- gWidgets::ggroup(cont = grp_posterior, horizontal=F)
    # sup_post <- gcheckbox("Suppress plot output", cont = grp_post_opt); assign("sup_post",sup_post,envir=mixsiar);
    grp_post_name <- gWidgets::ggroup(cont = grp_post_opt, horizontal=T, expand=T)
    plot_post_lbl <- gWidgets::glabel("Save plot as:", cont = grp_post_name)
    plot_post_name <- gWidgets::gedit("posterior_density", width=20, cont = grp_post_name); assign("plot_post_name",plot_post_name,envir=mixsiar);
    plot_post_save_pdf <- gWidgets::gcheckbox("pdf", cont = grp_post_name); assign("plot_post_save_pdf",plot_post_save_pdf,envir=mixsiar);
    gWidgets::svalue(mixsiar$plot_post_save_pdf) <- TRUE
    plot_post_save_png <- gWidgets::gcheckbox("png", cont = grp_post_name); assign("plot_post_save_png",plot_post_save_png,envir=mixsiar);

    # Pairs Plot options
    grp_pairs <- gWidgets::ggroup(cont=grp_output, horizontal=T)
    lbl_pairs <- gWidgets::glabel("Pairs Plot", cont = grp_pairs, expand=T)
    grp_pairs_opt <- gWidgets::ggroup(cont = grp_pairs, horizontal=F)
    # sup_pairs <- gcheckbox("Suppress plot output", cont = grp_pairs_opt); assign("sup_pairs",sup_pairs,envir=mixsiar);
    grp_pairs_name <- gWidgets::ggroup(cont = grp_pairs_opt, horizontal=T, expand=T)
    plot_pairs_lbl <- gWidgets::glabel("Save plot as:", cont = grp_pairs_name)
    plot_pairs_name <- gWidgets::gedit("pairs_plot", width=20, cont = grp_pairs_name); assign("plot_pairs_name",plot_pairs_name,envir=mixsiar);
    plot_pairs_save_pdf <- gWidgets::gcheckbox("pdf", cont = grp_pairs_name); assign("plot_pairs_save_pdf",plot_pairs_save_pdf,envir=mixsiar);
    gWidgets::svalue(mixsiar$plot_pairs_save_pdf) <- TRUE
    plot_pairs_save_png <- gWidgets::gcheckbox("png", cont = grp_pairs_name); assign("plot_pairs_save_png",plot_pairs_save_png,envir=mixsiar);

    # # XY Plot options
    # grp_xy <- ggroup(cont=grp_output, horizontal=T)
    # lbl_xy <- glabel("XY Plot", cont = grp_xy, expand=T)
    # grp_xy_opt <- ggroup(cont = grp_xy, horizontal=F)
    # # sup_xy <- gcheckbox("Suppress plot output", cont = grp_xy_opt); assign("sup_xy",sup_xy,envir=mixsiar);
    # grp_xy_name <- ggroup(cont = grp_xy_opt, horizontal=T, expand=T)
    # plot_xy_lbl <- glabel("Save plot as:", cont = grp_xy_name); assign("plot_xy_lbl",plot_xy_lbl,envir=mixsiar);
    # plot_xy_name <- gedit("xy_plot", width=20, cont = grp_xy_name); assign("plot_xy_name",plot_xy_name,envir=mixsiar);
    # plot_xy_save_pdf <- gcheckbox("pdf", cont = grp_xy_name); assign("plot_xy_save_pdf",plot_xy_save_pdf,envir=mixsiar);
    # svalue(mixsiar$plot_xy_save_pdf) <- TRUE
    # plot_xy_save_png <- gcheckbox("png", cont = grp_xy_name); assign("plot_xy_save_png",plot_xy_save_png,envir=mixsiar);

    # Diagnostics options
    grp_diag <- gWidgets::gframe(text="Diagnostics", cont=grp_output, horizontal=F)
    grp_diag_opts <- gWidgets::ggroup(cont=grp_diag, horizontal=T)
    gelman <- gWidgets::gcheckbox("Gelman-Rubin (must have > 1 chain)", cont=grp_diag_opts); assign("gelman",gelman,envir=mixsiar);
    gWidgets::svalue(mixsiar$gelman) <- TRUE
    geweke <- gWidgets::gcheckbox("Geweke", cont=grp_diag_opts); assign("geweke",geweke,envir=mixsiar);
    gWidgets::svalue(mixsiar$geweke) <- TRUE
    grp_diag_save <- gWidgets::ggroup(cont=grp_diag, horizontal=T)
    diag_save <- gWidgets::gcheckbox("Save diagnostics to file:", cont=grp_diag_save); assign("diag_save",diag_save,envir=mixsiar);
    diag_name <- gWidgets::gedit("diagnostics", width=20, cont=grp_diag_save); assign("diag_name",diag_name,envir=mixsiar);
    lbl_diag <- gWidgets::glabel("Note: diagnostics will print in the R command line if you do not choose to save to file",cont=grp_diag)
    gWidgets::svalue(mixsiar$diag_save) <- TRUE   # Default is to save the diagnostics

    ###########################################################################
    # RUN MODEL and Process Output
    ###########################################################################
    grp_bot <- gWidgets::ggroup(cont=grp_all,horizontal=TRUE)

    grp_run_model <- gWidgets::ggroup(cont=grp_bot, horizontal=TRUE, expand=TRUE)
    go_button <- gWidgets::gbutton(text="RUN MODEL", cont=grp_run_model, expand=TRUE,
                                   handler = function(h, ...){
                                     # Error checks on prior
                                     if(gWidgets::svalue(mixsiar$prior_option) == "\"Uninformative\"/Generalist"){
                                       alpha.prior <- rep(1,mixsiar$source$n.sources)
                                     } else { # prior_option = "Informative"
                                       alpha.prior <- eval(parse(text=gWidgets::svalue(mixsiar$inf_prior)))
                                     }

                                     if(gWidgets::svalue(mixsiar$error_option)=="Resid * Process"){resid_err <- TRUE; process_err <- TRUE;}
                                     if(gWidgets::svalue(mixsiar$error_option)=="Residual only"){resid_err <- TRUE; process_err <- FALSE;}
                                     if(gWidgets::svalue(mixsiar$error_option)=="Process only (N=1)"){resid_err <- FALSE; process_err <- TRUE;}
                                     write_JAGS_model("MixSIAR_model.txt", resid_err, process_err, mixsiar$mix, mixsiar$source)

                                     run <- gWidgets::svalue(mixsiar$mcmc_run)
                                     jags.1 <- run_model(run, mixsiar$mix, mixsiar$source, mixsiar$discr, "MixSIAR_model.txt", alpha.prior, resid_err, process_err)
                                     assign("jags.1",jags.1,envir=mixsiar)

                                     test <- get("jags.1",envir=mixsiar)
                                     if(exists("test")){
                                       gWidgets::add(grp_run_model,gWidgets::gimage(system.file("extdata", "check.png", package = "MixSIAR")))
                                       gWidgets::svalue(mixsiar$status_bar) <- "Model run successful. Click 'Process output' for diagnostics, plots, and stats"
                                     }
                                   })

    grp_output <- gWidgets::ggroup(cont=grp_bot, horizontal=TRUE, expand=TRUE)
    output_button <- gWidgets::gbutton(text="Process output", cont=grp_output, expand=TRUE,
                                       handler = function(h, ...){
                                         output_options <- list(gWidgets::svalue(mixsiar$summary_save),     # Save the summary statistics as a txt file?
                                                                gWidgets::svalue(mixsiar$summary_name),            # If yes, specify the base file name (.txt will be appended later)
                                                                # svalue(mixsiar$sup_post),                # Suppress posterior density plot output in R?
                                                                FALSE,
                                                                gWidgets::svalue(mixsiar$plot_post_save_pdf),      # Save posterior density plots as pdfs?
                                                                gWidgets::svalue(mixsiar$plot_post_name),          # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                                                # svalue(mixsiar$sup_pairs),               # Suppress pairs plot output in R?
                                                                FALSE,
                                                                gWidgets::svalue(mixsiar$plot_pairs_save_pdf),     # Save pairs plot as pdf?
                                                                gWidgets::svalue(mixsiar$plot_pairs_name),         # If yes, specify the base file name (.pdf/.png will be appended later)
                                                                # svalue(mixsiar$sup_xy),                  # Suppress xy/trace plot output in R?
                                                                TRUE,
                                                                # svalue(mixsiar$plot_xy_save_pdf),        # Save xy/trace plot as pdf?
                                                                FALSE,
                                                                # svalue(mixsiar$plot_xy_name),            # If yes, specify the base file name (.pdf/.png will be appended later)
                                                                NULL,
                                                                gWidgets::svalue(mixsiar$gelman),                  # Calculate Gelman-Rubin diagnostic test?
                                                                FALSE,                                   # Calculate Heidelberg-Welch diagnostic test?
                                                                gWidgets::svalue(mixsiar$geweke),                  # Calculate Geweke diagnostic test?
                                                                gWidgets::svalue(mixsiar$diag_save),               # Save the diagnostics as a txt file?
                                                                gWidgets::svalue(mixsiar$diag_name),               # If yes, specify the base file name (.txt will be appended later)
                                                                FALSE,                                   # Is Individual a random effect in the model? (already specified)
                                                                gWidgets::svalue(mixsiar$plot_post_save_png),      # Save posterior density plots as pngs?
                                                                gWidgets::svalue(mixsiar$plot_pairs_save_png),     # Save pairs plot as png?
                                                                # svalue(mixsiar$plot_xy_save_png))
                                                                FALSE)
                                         output_JAGS(mixsiar$jags.1, mixsiar$mix, mixsiar$source, output_options)
                                       })

    # Show the GUI once all the code has run
    gWidgets::visible(win) <- TRUE
  }
  } else {
    stop(paste("*** gWidgetsRGtk2 package not able to be loaded. ***
        If 'library('gWidgetsRGtk2')' does not work, MixSIAR GUI will not run.
        On Windows/Linux, try 'install.packages('gWidgetsRGtk2')'.
        On Mac, close R, download and install GTK+ from:
        http://r.research.att.com/#other. Then install latest X11 application
        (xQuartz) from http://xquartz.macosforge.org/landing/.

        If installing GTK+ continues to be problematic, consider using the
        script version of MixSIAR. See manual for help and examples.",sep=""))
  }
}
