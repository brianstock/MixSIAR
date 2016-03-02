#' Run the GUI version of MixSIAR
#'
#' \code{mixsiar_gui} creates the GUI version of MixSIAR.
#'
#' \emph{Before running this function}, a new user will need to install JAGS
#' and GTK+. See Section 2 of the manual for install instructions:
#' \url{https://github.com/brianstock/MixSIAR/blob/master/mixsiar_manual_3.0.pdf}
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

# mixsiar <- new.env() # moved to new file: 'create_mixsiar_envir.R'
mixsiar_gui <- function(){

runif(1)
old <- options("guiToolkit"="RGtk2")
on.exit(options(old), add = TRUE)

win<-gwindow("MixSIAR GUI", visible=FALSE)
grp_all <- ggroup(cont=win, horizontal=FALSE)
grp_input <- ggroup(cont=grp_all, horizontal=TRUE)

############################################################################
# Read in Data
############################################################################
grp_readin_data <- gframe(text="Read in data", cont=grp_input, horizontal=F)
grp_cons <- ggroup(horizontal=TRUE, cont=grp_readin_data); assign("grp_cons",grp_cons,envir=mixsiar);
btn_cons <- gbutton(
  text = "Load mixture data",
  cont = mixsiar$grp_cons,
  expand = T,
  handler	= function(h, ...){
    build_mix_win()
  }
)

grp_source <- ggroup(horizontal=TRUE, cont = grp_readin_data); assign("grp_source",grp_source,envir=mixsiar)
btn_source <- gbutton(
  text = "Load source data",
  cont = mixsiar$grp_source,
  expand = T,
  handler	= function(h, ...){
    build_source_win()
  }
)

grp_frac <- ggroup(horizontal=TRUE, cont = grp_readin_data)
btn_frac <- gbutton(
  text = "Load discrimination data",
  cont = grp_frac,
  expand = T,
  handler	= function(h, ...){
    gfile(
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
          addSpring(grp_frac)
          svalue(mixsiar$status_bar) <- "Discrimination data successfully loaded"
          discr_filename <- h$file; assign("discr_filename",discr_filename,envir=mixsiar);
          },
          error = function(e){
            svalue(mixsiar$status_bar) <- "Could not load data"
            add(grp_frac,gimage("red x.png"))
          }
        )
      }
    )
    discr <- load_discr_data(mixsiar$discr_filename, mixsiar$mix)
    assign("discr", discr, envir = mixsiar)
    add(grp_frac,gimage("check.png"))
  }
)
#############################################################################
# User-specified MCMC parameters
#############################################################################
grp_set_mcmc <- gframe(text="MCMC run length", cont=grp_input, horizontal=TRUE, expand=TRUE)
grp_mcmc <- ggroup(cont=grp_set_mcmc, horizontal=FALSE)
mcmc_run <- gradio(c("test","very short","short","normal","long","very long"), cont=grp_mcmc, horizontal=FALSE)
assign("mcmc_run",mcmc_run,envir=mixsiar)
svalue(mixsiar$mcmc_run) <- "test"

####################################################################
# Model Error Options
####################################################################
grp_error_priors <- ggroup(cont=grp_input,horizontal=FALSE)
grp_error <- gframe(text="Error structure", cont=grp_error_priors, horizontal=FALSE,expand=TRUE)
error_option <- gradio(c("Resid * Process","Residual only","Process only (N=1)"), cont=grp_error, horizontal=FALSE)
assign("error_option",error_option,envir=mixsiar)
svalue(mixsiar$error_option) <- "Resid * Process"

####################################################################
# Specify Prior
####################################################################
grp_prior <- gframe(text="Specify prior", cont=grp_error_priors, horizontal=FALSE)
prior_option <- gradio(c("\"Uninformative\"/Generalist","Informative"),cont=grp_prior,horizontal=FALSE)
inf_prior <- gedit("", width=15, cont = grp_prior); assign("inf_prior",inf_prior,envir=mixsiar)
assign("prior_option",prior_option,envir=mixsiar)
svalue(mixsiar$prior_option) <- "\"Uninformative\"/Generalist"

status_bar <- gstatusbar("", progress.bar="gui", cont=grp_all, expand=TRUE)
assign("status_bar",status_bar,envir=mixsiar)

####################################################################
# Isospace Plot
####################################################################

# The 'Make isospace plot' button calls the plot_data function to make an isospace plot
grp_plot <- ggroup(cont=grp_all, horizontal=T)
plot_button <- gbutton(
  text = "Make isospace plot",
  cont = grp_plot,
  expand = TRUE,
  handler = function(h, ...){
    plot_data(svalue(mixsiar$plot_filename), svalue(mixsiar$plot_save_pdf), svalue(mixsiar$plot_save_png), mixsiar$mix,mixsiar$source,mixsiar$discr)
    svalue(mixsiar$status_bar) <- "If isospace plot looks good, you can now click 'RUN MODEL' at bottom"
  }
)

grp_plot_name <- ggroup(cont=grp_plot, horizontal=T, expand=T)
plot_lbl <- glabel("Save plot as:",cont=grp_plot_name)
plot_filename <- gedit("isospace_plot", width=15, cont = grp_plot_name); assign("plot_filename",plot_filename,envir=mixsiar)
plot_save_pdf <- gcheckbox("pdf", cont = grp_plot_name); assign("plot_save_pdf",plot_save_pdf,envir=mixsiar)
plot_save_png <- gcheckbox("png", cont = grp_plot_name); assign("plot_save_png",plot_save_png,envir=mixsiar)
svalue(mixsiar$plot_save_pdf) <- TRUE

####################################################################
# Plot Prior
####################################################################

# The 'Plot prior' button calls the plot_prior function
grp_plot_prior <- ggroup(cont=grp_all, horizontal=T)
addSpring(grp_plot_prior)
plot_button <- gbutton(
  text = "Plot prior",
  cont = grp_plot_prior,
  expand = TRUE,
  handler = function(h, ...){
    if(svalue(mixsiar$prior_option) == "\"Uninformative\"/Generalist"){
      alpha.prior <- rep(1,mixsiar$source$n.sources)
    } else { # prior_option = "Informative"
      alpha.prior <- eval(parse(text=svalue(mixsiar$inf_prior)))
    }
    if(!is.numeric(alpha.prior)){
      stop(paste("*** Error: Your prior is not a numeric vector of length(n.sources).
        Try again or choose the uninformative prior option. For example,
        c(1,1,1,1) is a valid (uninformative) prior for 4 sources. ***",sep=""))}
    if(length(alpha.prior) != mixsiar$source$n.sources){
      stop(paste("*** Error: Length of your prior does not match the
        number of sources (",mixsiar$source$n.sources,"). Try again. ***",sep=""))}
    svalue(mixsiar$status_bar) <- "Success. Your prior is in RED, uninformative/generalist is DARK GREY"
    plot_prior(alpha.prior,mixsiar$source,svalue(plot_save_pdf_prior),svalue(plot_save_png_prior),svalue(plot_filename_prior))
  }
)

grp_plot_name_prior <- ggroup(cont=grp_plot_prior, horizontal=T, expand=T)
plot_lbl_prior <- glabel("Save plot as:",cont=grp_plot_name_prior)
plot_filename_prior <- gedit("prior_plot", width=15, cont = grp_plot_name_prior); assign("plot_filename_prior",plot_filename_prior,envir=mixsiar)
plot_save_pdf_prior <- gcheckbox("pdf", cont = grp_plot_name_prior); assign("plot_save_pdf_prior",plot_save_pdf_prior,envir=mixsiar)
plot_save_png_prior <- gcheckbox("png", cont = grp_plot_name_prior); assign("plot_save_png_prior",plot_save_png_prior,envir=mixsiar)
svalue(mixsiar$plot_save_pdf_prior) <- TRUE

####################################################################
# Output options
####################################################################
grp_output <- gframe(text="Output options", cont=grp_all, horizontal=F)

# Summary Statistics options
grp_summary <- ggroup(cont=grp_output, horizontal=T)
lbl_summary <- glabel("Summary Statistics", cont=grp_summary, expand=T)
grp_summary_right <- ggroup(cont=grp_summary, horizontal=F)
grp_summary_save <- ggroup(cont=grp_summary_right, horizontal=T)
summary_save <- gcheckbox("Save summary statistics to file: ", cont=grp_summary_save); assign("summary_save",summary_save,envir=mixsiar);
summary_name <- gedit("summary_statistics", width=20, cont=grp_summary_save); assign("summary_name",summary_name,envir=mixsiar);
svalue(mixsiar$summary_save) <- TRUE  # Default is to save the summary statistics

# Posterior Density Plot options
grp_posterior <- ggroup(cont=grp_output, horizontal=T)
lbl_posterior <- glabel("Posterior Density Plot", cont=grp_posterior, expand=T)
grp_post_opt <- ggroup(cont = grp_posterior, horizontal=F)
# sup_post <- gcheckbox("Suppress plot output", cont = grp_post_opt); assign("sup_post",sup_post,envir=mixsiar);
grp_post_name <- ggroup(cont = grp_post_opt, horizontal=T, expand=T)
plot_post_lbl <- glabel("Save plot as:", cont = grp_post_name)
plot_post_name <- gedit("posterior_density", width=20, cont = grp_post_name); assign("plot_post_name",plot_post_name,envir=mixsiar);
plot_post_save_pdf <- gcheckbox("pdf", cont = grp_post_name); assign("plot_post_save_pdf",plot_post_save_pdf,envir=mixsiar);
svalue(mixsiar$plot_post_save_pdf) <- TRUE
plot_post_save_png <- gcheckbox("png", cont = grp_post_name); assign("plot_post_save_png",plot_post_save_png,envir=mixsiar);

# Pairs Plot options
grp_pairs <- ggroup(cont=grp_output, horizontal=T)
lbl_pairs <- glabel("Pairs Plot", cont = grp_pairs, expand=T)
grp_pairs_opt <- ggroup(cont = grp_pairs, horizontal=F)
# sup_pairs <- gcheckbox("Suppress plot output", cont = grp_pairs_opt); assign("sup_pairs",sup_pairs,envir=mixsiar);
grp_pairs_name <- ggroup(cont = grp_pairs_opt, horizontal=T, expand=T)
plot_pairs_lbl <- glabel("Save plot as:", cont = grp_pairs_name)
plot_pairs_name <- gedit("pairs_plot", width=20, cont = grp_pairs_name); assign("plot_pairs_name",plot_pairs_name,envir=mixsiar);
plot_pairs_save_pdf <- gcheckbox("pdf", cont = grp_pairs_name); assign("plot_pairs_save_pdf",plot_pairs_save_pdf,envir=mixsiar);
svalue(mixsiar$plot_pairs_save_pdf) <- TRUE
plot_pairs_save_png <- gcheckbox("png", cont = grp_pairs_name); assign("plot_pairs_save_png",plot_pairs_save_png,envir=mixsiar);

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
grp_diag <- gframe(text="Diagnostics", cont=grp_output, horizontal=F)
grp_diag_opts <- ggroup(cont=grp_diag, horizontal=T)
gelman <- gcheckbox("Gelman-Rubin (must have > 1 chain)", cont=grp_diag_opts); assign("gelman",gelman,envir=mixsiar);
svalue(mixsiar$gelman) <- TRUE
geweke <- gcheckbox("Geweke", cont=grp_diag_opts); assign("geweke",geweke,envir=mixsiar);
svalue(mixsiar$geweke) <- TRUE
grp_diag_save <- ggroup(cont=grp_diag, horizontal=T)
diag_save <- gcheckbox("Save diagnostics to file:", cont=grp_diag_save); assign("diag_save",diag_save,envir=mixsiar);
diag_name <- gedit("diagnostics", width=20, cont=grp_diag_save); assign("diag_name",diag_name,envir=mixsiar);
lbl_diag <- glabel("Note: diagnostics will print in the R command line if you do not choose to save to file",cont=grp_diag)
svalue(mixsiar$diag_save) <- TRUE   # Default is to save the diagnostics

###########################################################################
# RUN MODEL and Process Output
###########################################################################
grp_bot <- ggroup(cont=grp_all,horizontal=TRUE)

grp_run_model <- ggroup(cont=grp_bot, horizontal=TRUE, expand=TRUE)
go_button <- gbutton(text="RUN MODEL", cont=grp_run_model, expand=TRUE,
  handler = function(h, ...){
    # Error checks on prior
    if(svalue(mixsiar$prior_option) == "\"Uninformative\"/Generalist"){
      alpha.prior <- rep(1,mixsiar$source$n.sources)
    } else { # prior_option = "Informative"
      alpha.prior <- eval(parse(text=svalue(mixsiar$inf_prior)))
    }

    if(svalue(mixsiar$error_option)=="Resid * Process"){resid_err <- TRUE; process_err <- TRUE;}
    if(svalue(mixsiar$error_option)=="Residual only"){resid_err <- TRUE; process_err <- FALSE;}
    if(svalue(mixsiar$error_option)=="Process only (N=1)"){resid_err <- FALSE; process_err <- TRUE;}
    write_JAGS_model("MixSIAR_model.txt", resid_err, process_err, mixsiar$mix, mixsiar$source)

    run <- svalue(mixsiar$mcmc_run)
    jags.1 <- run_model(run, mixsiar$mix, mixsiar$source, mixsiar$discr, "MixSIAR_model.txt", alpha.prior)
    assign("jags.1",jags.1,envir=mixsiar)

    test <- get("jags.1",envir=mixsiar)
    if(exists("test")){
      add(grp_run_model,gimage("check.png"))
      svalue(mixsiar$status_bar) <- "Model run successful. Click 'Process output' for diagnostics, plots, and stats"
    }
  })

grp_output <- ggroup(cont=grp_bot, horizontal=TRUE, expand=TRUE)
output_button <- gbutton(text="Process output", cont=grp_output, expand=TRUE,
  handler = function(h, ...){
    output_options <- list(svalue(mixsiar$summary_save),     # Save the summary statistics as a txt file?
                    svalue(mixsiar$summary_name),            # If yes, specify the base file name (.txt will be appended later)
                    # svalue(mixsiar$sup_post),                # Suppress posterior density plot output in R?
                    FALSE,
                    svalue(mixsiar$plot_post_save_pdf),      # Save posterior density plots as pdfs?
                    svalue(mixsiar$plot_post_name),          # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                    # svalue(mixsiar$sup_pairs),               # Suppress pairs plot output in R?
                    FALSE,
                    svalue(mixsiar$plot_pairs_save_pdf),     # Save pairs plot as pdf?
                    svalue(mixsiar$plot_pairs_name),         # If yes, specify the base file name (.pdf/.png will be appended later)
                    # svalue(mixsiar$sup_xy),                  # Suppress xy/trace plot output in R?
                    TRUE,
                    # svalue(mixsiar$plot_xy_save_pdf),        # Save xy/trace plot as pdf?
                    FALSE,
                    # svalue(mixsiar$plot_xy_name),            # If yes, specify the base file name (.pdf/.png will be appended later)
                    NULL,
                    svalue(mixsiar$gelman),                  # Calculate Gelman-Rubin diagnostic test?
                    FALSE,                                   # Calculate Heidelberg-Welch diagnostic test?
                    svalue(mixsiar$geweke),                  # Calculate Geweke diagnostic test?
                    svalue(mixsiar$diag_save),               # Save the diagnostics as a txt file?
                    svalue(mixsiar$diag_name),               # If yes, specify the base file name (.txt will be appended later)
                    FALSE,                                   # Is Individual a random effect in the model? (already specified)
                    svalue(mixsiar$plot_post_save_png),      # Save posterior density plots as pngs?
                    svalue(mixsiar$plot_pairs_save_png),     # Save pairs plot as png?
                    # svalue(mixsiar$plot_xy_save_png))
                    FALSE)
    output_JAGS(mixsiar$jags.1, mixsiar$mix, mixsiar$source, output_options)
  })

# Show the GUI once all the code has run
visible(win) <- TRUE
}
