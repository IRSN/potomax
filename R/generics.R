## ****************************************************************************
##' Return periods associated with chosen quantiles.
##'
##' @title Return Periods associated with Chosen Quantiles
##'
##' @param object An object defining \emph{return periods} associated
##' to some "quantiles" values.
##'
##' @param ... Arguments passed to methods.
##'
##' @return An object inheriting from the \code{"data.frame"} class
##' contining the return periods.
##' 
RP <- function(object, ...) {   
    UseMethod("RP", object)
}
