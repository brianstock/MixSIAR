# Brian Stock
# Oct 3, 2013
# 
# Previous version:
#   MixSIAR GUI 0.2

# Major changes:
#   Convert MixSIAR into a package-able form: function mixsiar() that builds GUI
#   Made global all gWidgets that are used in other functions...quick fix for now

#   Revised wording:
#     - "Discrimination" instead of "Fractionation"
#     - "Mixture" instead of "Consumer"

# Before running this script, a brand new user will need to install the latest
# versions of R and JAGS.  The install.packages("gWidgetsRGtk2") command
# will prompt the user to install GTK+, which also needs to happen.

# Output file: output_JAGS.r
# Model file: write_JAGS_model.r
# Auxillary files:
#   - build_mix_win.r
#   - build_source_win.r
#   - run_model.r
#   - plot_data.r
#   - plot_data_one_iso.r
#   - plot_continuous_var.r

mixsiar <- function(){

# Check for R Version 3.0 or 2.15
minor <- as.numeric(R.Version()$minor)
major <- as.numeric(R.Version()$major)
if(major<3 && minor<15){
  stop(paste("*** Error: You are running ", R.Version()$version.string, ", which is 
                  out of date. Please update R to 3.0 or 2.15 and try again. ***",sep=""))
}

if (!"ggplot2" %in% installed.packages()) install.packages("ggplot2")
if (!"gWidgetsRGtk2" %in% installed.packages()) install.packages("gWidgetsRGtk2")
if (!"runjags" %in% installed.packages()) install.packages("runjags")
if (!"R2jags" %in% installed.packages()) install.packages("R2jags")
if (!"matrixStats" %in% installed.packages()) install.packages("matrixStats")
if (!"MASS" %in% installed.packages()) install.packages("MASS")
if (!"RColorBrewer" %in% installed.packages()) install.packages("RColorBrewer")
if (!"reshape" %in% installed.packages()) install.packages("reshape")

if (!"gWidgetsRGtk2" %in% installed.packages()) stop("*** Error: GTK+ is not installed ***")
if (!"R2jags" %in% installed.packages()) stop("*** Error: JAGS is not installed ***")

require(ggplot2)
require(gWidgetsRGtk2)
require(runjags)
require(R2jags)
require(matrixStats)
require(MASS)
require(RColorBrewer)
require(reshape)

runif(1)
rm(list=ls())
options("guiToolkit"="RGtk2")
source("output_JAGS.r")
source("write_JAGS_model.r")
source("build_mix_win.r")
source("build_source_win.r")
source("run_model.r")
source("plot_data.r")
source("plot_data_one_iso.r")
source("plot_continuous_var.r")

win<-gwindow("MixSIAR GUI", visible=FALSE)
grp_all <- ggroup(cont=win, horizontal=FALSE)
grp_input <- ggroup(cont=grp_all, horizontal=TRUE)

############################################################################
# Read in Data
############################################################################
grp_readin_data <- gframe(text="Read in data", cont=grp_input, horizontal=F)
grp_cons <<- ggroup(horizontal=TRUE, cont=grp_readin_data)
btn_cons <- gbutton(
  text = "Load mixture data",
  cont = grp_cons,
  expand = T,
  handler	= function(h, ...){
    build_mix_win()
  }
)

grp_source <<- ggroup(horizontal=TRUE, cont = grp_readin_data)
btn_source <- gbutton(
  text = "Load source data",
  cont = grp_source,
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
          assign(data_frame_name, the_data, envir = globalenv())
          addSpring(grp_frac)
          svalue(status_bar) <- "Discrimination data successfully loaded"
          },
          error = function(e){
            svalue(status_bar) <- "Could not load data"
            add(grp_frac,gimage("red x.png"))
          }
        )
      }
    )

    row.names(FRAC)<-FRAC[,1]     # store the row names of FRAC (sources)
    FRAC <- as.matrix(FRAC[-1])   # remove source names column of FRAC
    FRAC <- FRAC[order(rownames(FRAC)),]  # rearrange FRAC so sources are in alphabetical order

    # Make sure the iso columns of frac_mu and frac_sig2 are in the same order as S_MU and S_SIG
    frac_mu_cols <- match(MU_names,colnames(FRAC))   # get the column numbers of FRAC that correspond to the means
    frac_sig_cols <- match(SIG_names,colnames(FRAC))   # get the column numbers of FRAC that correspond to the SDs
    frac_mu <- FRAC[,frac_mu_cols]                            # FRAC means
    frac_sig2 <- FRAC[,frac_sig_cols]*FRAC[,frac_sig_cols]    # FRAC variances
 
    assign("frac_mu", frac_mu, envir = .GlobalEnv)    # Assign 'frac_mu' and 'frac_sig2' to the global environment
    assign("frac_sig2", frac_sig2, envir = .GlobalEnv)
    add(grp_frac,gimage("check.png"))
  }
)
#############################################################################
# User-specified MCMC parameters
#############################################################################
grp_set_parameters <- gframe(text="Specify MCMC parameters", 
                             cont=grp_input, 
                             horizontal=FALSE)
# Number of chains value is 'txt_num_chains'
grp_num_chains <- ggroup(cont=grp_set_parameters, horizontal=TRUE)
lbl_num_chains <- glabel(
  "# of Chains: ",
  cont = grp_num_chains
)
addSpring(grp_num_chains)
txt_num_chains <<- gedit("3", width=10, cont = grp_num_chains) 

# Chain length value is 'txt_chain_length'
grp_chain_length <- ggroup(cont=grp_set_parameters, horizontal=TRUE)
lbl_chain_length <- glabel(
  "Chain Length: ",
  cont = grp_chain_length
)
addSpring(grp_chain_length)
txt_chain_length <<- gedit("10000", width=10, cont = grp_chain_length) 

# Burn-in value is 'txt_burnin'
grp_burnin <- ggroup(cont=grp_set_parameters, horizontal=TRUE)
lbl_burnin <- glabel(
  "Burn-in: ",
  cont = grp_burnin
)
addSpring(grp_burnin)
txt_burnin <<- gedit("5000", width=10, cont = grp_burnin)   

# Thinning value is 'txt_thin'
grp_thin <- ggroup(cont=grp_set_parameters, horizontal=TRUE)
lbl_thin <- glabel(
  "Thin: ",
  cont = grp_thin
)
addSpring(grp_thin)
txt_thin <<- gedit("10", width=10, cont = grp_thin)

####################################################################
# Model Error Options
####################################################################
grp_error <- gframe(text="Model Error Options", cont=grp_input, horizontal=F)
resid_err_box <<- gcheckbox("Include residual error", cont=grp_error) # Does the user want to include residual/SIAR error in the model? resid_var = 1/resid_tau, resid_tau ~ dgamma(.001,.001)
proc_err_box <<- gcheckbox("Include process error", cont=grp_error)    # Does the user want to include process/MixSIR error in the model?  process_var = p2[iso,i]*sig2_source[iso,i] + p2[iso,i]*sig2_frac[iso,i]
svalue(resid_err_box) <- TRUE                                           # Default is residual/SIAR error AND process error
svalue(proc_err_box) <- TRUE

status_bar <<- gstatusbar("", progress.bar="gui", cont=grp_all, expand=TRUE)

####################################################################
# Isospace Plot
####################################################################

# The 'Plot data' button calls the plot_data function to make an isospace plot
grp_plot <- ggroup(cont=grp_all, horizontal=T)
plot_button <- gbutton(
  text = "Make isospace plot",
  cont = grp_plot,
  expand = TRUE,
  handler = function(h, ...){
    if(n.iso==1){
      plot_data_one_iso()
    } else {
      for(iso1 in 1:(n.iso-1)){
        for(iso2 in (iso1+1):n.iso){
          plot_data(c(iso1,iso2))
        }
      }
    }
  }
)

grp_plot_name <- ggroup(cont=grp_plot, horizontal=T, expand=T)
plot_lbl <- glabel("Save plot as:",cont=grp_plot_name)
plot_filename <<- gedit("isospace_plot", width=15, cont = grp_plot_name)
plot_save_pdf <<- gcheckbox("pdf", cont = grp_plot_name)
plot_save_png <<- gcheckbox("png", cont = grp_plot_name)
svalue(plot_save_pdf) <- TRUE

####################################################################
# Output options
####################################################################
grp_output <- gframe(text="Output options", cont=grp_all, horizontal=F)

## Choose working directory
#grp_dir <- ggroup(cont=grp_output, horizontal=T)
#chooseWD <- function(h,...){
#  old_dir <- getwd()
#  setwd(tclvalue(tkchooseDirectory()))
#  new_dir <- getwd()
#  file.copy(from=file.path(paste(old_dir,"/MixSIAR.txt",sep="")),
#              to=file.path(paste(new_dir,"/MixSIAR.txt",sep="")))
#}
#choose_wd <- gbutton(
#  text = "Choose working directory to save output",
#  cont = grp_dir,
#  expand = F,
#  handler = chooseWD
#)
#save_work <- gcheckbox("Save R workspace", cont=grp_dir)

# Summary Statistics options
grp_summary <- ggroup(cont=grp_output, horizontal=T)
lbl_summary <- glabel("Summary Statistics", cont=grp_summary, expand=T)
grp_summary_right <- ggroup(cont=grp_summary, horizontal=F)
grp_summary_save <- ggroup(cont=grp_summary_right, horizontal=T)
summary_save <<- gcheckbox("Save summary statistics to file: ", cont=grp_summary_save)
summary_name <<- gedit("summary_statistics", width=20, cont=grp_summary_save)
svalue(summary_save) <- TRUE  # Default is to save the summary statistics

# Posterior Density Plot options
grp_posterior <- ggroup(cont=grp_output, horizontal=T)
lbl_posterior <- glabel("Posterior Density Plot", cont=grp_posterior, expand=T)
grp_post_opt <- ggroup(cont = grp_posterior, horizontal=F)
sup_post <<- gcheckbox("Suppress plot output", cont = grp_post_opt)
grp_post_name <- ggroup(cont = grp_post_opt, horizontal=T, expand=T)
plot_post_lbl <- glabel("Save plot as:", cont = grp_post_name)
plot_post_name <<- gedit("posterior_density", width=20, cont = grp_post_name)
plot_post_save_pdf <<- gcheckbox("pdf", cont = grp_post_name)
svalue(plot_post_save_pdf) <- TRUE
plot_post_save_png <<- gcheckbox("png", cont = grp_post_name)

# Pairs Plot options
grp_pairs <- ggroup(cont=grp_output, horizontal=T)
lbl_pairs <- glabel("Pairs Plot", cont = grp_pairs, expand=T)
grp_pairs_opt <- ggroup(cont = grp_pairs, horizontal=F)
sup_pairs <<- gcheckbox("Suppress plot output", cont = grp_pairs_opt)
grp_pairs_name <- ggroup(cont = grp_pairs_opt, horizontal=T, expand=T)
plot_pairs_lbl <- glabel("Save plot as:", cont = grp_pairs_name)
plot_pairs_name <<- gedit("pairs_plot", width=20, cont = grp_pairs_name)
plot_pairs_save_pdf <<- gcheckbox("pdf", cont = grp_pairs_name)
svalue(plot_pairs_save_pdf) <- TRUE
plot_pairs_save_png <<- gcheckbox("png", cont = grp_pairs_name)

# XY Plot options
grp_xy <- ggroup(cont=grp_output, horizontal=T)
lbl_xy <- glabel("XY Plot", cont = grp_xy, expand=T)
grp_xy_opt <- ggroup(cont = grp_xy, horizontal=F)
sup_xy <<- gcheckbox("Suppress plot output", cont = grp_xy_opt)
grp_xy_name <- ggroup(cont = grp_xy_opt, horizontal=T, expand=T)
plot_xy_lbl <<- glabel("Save plot as:", cont = grp_xy_name)
plot_xy_name <<- gedit("xy_plot", width=20, cont = grp_xy_name)
plot_xy_save_pdf <<- gcheckbox("pdf", cont = grp_xy_name)
svalue(plot_xy_save_pdf) <- TRUE
plot_xy_save_png <<- gcheckbox("png", cont = grp_xy_name)

# Diagnostics options
grp_diag <- gframe(text="Diagnostics", cont=grp_output, horizontal=F)
grp_diag_opts <- ggroup(cont=grp_diag, horizontal=T)
gelman <<- gcheckbox("Gelman-Rubin (must have > 1 chain)", cont=grp_diag_opts)
svalue(gelman) <- TRUE
heidel <<- gcheckbox("Heidelberg-Welch", cont=grp_diag_opts)
svalue(heidel) <- TRUE
geweke <<- gcheckbox("Geweke", cont=grp_diag_opts)
svalue(geweke) <- TRUE
grp_diag_save <- ggroup(cont=grp_diag, horizontal=T)
diag_save <<- gcheckbox("Save diagnostics to file:", cont=grp_diag_save)
diag_name <<- gedit("diagnostics", width=20, cont=grp_diag_save)
lbl_diag <- glabel("Note: diagnostics will print in the R command line if you do not choose to save to file",cont=grp_diag)
svalue(diag_save) <- TRUE   # Default is to save the diagnostics

# The 'RUN MODEL' button calls the main 'run_model' function, which writes the 
# JAGS model file, calls JAGS, plots the JAGS output, and runs diagnostics
go_button <- gbutton(
  text = "RUN MODEL",
  cont = grp_all,
  expand = TRUE,
  handler = function(h, ...){
    run_model()
  }
)

# Show the GUI once all the code has run
visible(win) <- TRUE
}
