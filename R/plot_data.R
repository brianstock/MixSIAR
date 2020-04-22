#' Plot biotracer data
#'
#' \code{plot_data} creates plot(s) of the biotracer data and saves the plot(s)
#' to file(s) in the working directory. All 3 required data files must have been
#' loaded by \code{\link{load_mix_data}}, \code{\link{load_source_data}},
#' and \code{\link{load_discr_data}}. Behavior depends on the number of tracers:
#' \itemize{
#'  \item 1 tracer: calls \code{\link{plot_data_one_iso}} to create a 1-D plot.
#'  \item 2 tracers: calls \code{\link{plot_data_two_iso}} to create a biplot.
#'  \item >2 tracers: calls \code{\link{plot_data_two_iso}} in a loop to create
#'  biplots for each pairwise combination of biotracers.
#' }
#'
#' An important detail is that \code{plot_data_two_iso} and \code{plot_data_one_iso}
#' plot the raw mix data and \emph{add the TDF to the source data}, since this is
#' the polygon that the mixing model uses to determine proportions. The plotted
#' source means are:
#' \deqn{\mu_source + \mu_discr}
#' The source error bars are +/- 1 standard deviation, \emph{calculated as a
#' combination of source and TDF variances:}
#' \deqn{\sqrt{\sigma^2_source + \sigma^2_discr}}
#'
#' \code{plot_data} looks for 'C', 'N', 'S', and 'O' in the biotracer column
#' headers and assumes they are stable isotopes, labeling the axes with, e.g.,
#' expression(paste(delta^13, "C (u2030)",sep="")).
#'
#' @param filename name of the plot file(s) to save (e.g. "isospace_plot")
#' @param plot_save_pdf T/F, save the plot(s) as a pdf?
#' @param plot_save_png T/F, save the plot(s) as a png?
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param discr output from \code{\link{load_discr_data}}
#' @param return_obj T/F, whether or not to return ggplot object for further modification, defaults to F
#'
#' @seealso \code{\link{plot_data_two_iso}}, \code{\link{plot_data_one_iso}}
#' @export
plot_data <- function(filename,plot_save_pdf,plot_save_png,mix,source,discr, return_obj=FALSE){
  # check that discr rownames match source_names
  if(!identical(rownames(discr$mu),source$source_names)){
    stop(paste("*** Error: Source names do not match in source and discr
    data files. Please check your source and discr data file row names.",sep=""))}
  if(!identical(rownames(discr$sig2),source$source_names)){
    stop(paste("*** Error: Source names do not match in source and discr
    data files. Please check your source and discr data file row names.",sep=""))}

  if(mix$n.iso==1){
    plot_data_one_iso(mix,source,discr,filename,plot_save_pdf,plot_save_png,return_obj)
    if(return_obj == TRUE) {
      g = plot_data_one_iso(mix,source,discr,filename,plot_save_pdf,plot_save_png,return_obj=return_obj)
    }
  } else {
    for(iso1 in 1:(mix$n.iso-1)){
      for(iso2 in (iso1+1):mix$n.iso){
        plot_data_two_iso(isotopes=c(iso1,iso2),mix=mix,source=source,
          discr=discr,
          filename=filename,
          plot_save_pdf=plot_save_pdf,
          plot_save_png=plot_save_png,return_obj=return_obj)
        if(return_obj == TRUE) {
          g = plot_data_two_iso(isotopes=c(iso1,iso2),mix=mix,source=source,
            discr=discr,
            filename=filename,
            plot_save_pdf=plot_save_pdf,
            plot_save_png=plot_save_png,return_obj=return_obj)
        }
      }
    }
  }
  if(return_obj == TRUE) {
    return(g)
  }
} # end plot_data function
