#' mixsiar
#'
#' \code{mixsiar} is an R environment the GUI version of MixSIAR uses to store
#' objects and pass values through functions. Internal---only used
#' in the event a user wants to extract/see MixSIAR objects (e.g.
#' loaded mix data are accessed via \code{mixsiar$mix}, and the finished
#' model output, \code{rjags} object, is \code{mixsiar$jags.1})
#'
mixsiar <- new.env()
# Added only to pass R CMD check with 0 NOTEs
if(getRversion() >= "2.15.1")  utils::globalVariables(c("p.global","ilr.global","ilr.fac1","ilr.fac2","fac1.sig","fac2.sig","ind.sig","..scaled..","p.fac1","p.fac2","p.ind","sources","x","position_jitter","scolour","xmin","xmax","label","y","ymin","ymax","low","high"))
