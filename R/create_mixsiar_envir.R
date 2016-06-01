#' mixsiar
#'
#' \code{mixsiar} is an R environment the GUI version of MixSIAR uses to store
#' objects and pass values through functions. Internal---only used
#' in the event a user wants to extract/see MixSIAR objects (e.g.
#' loaded mix data are accessed via \code{mixsiar$mix}, and the finished
#' model output, \code{rjags} object, is \code{mixsiar$jags.1})
#'
mixsiar <- new.env()
