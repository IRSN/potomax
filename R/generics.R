

## ****************************************************************************
##' Return levels associated with chosen periods.
##'
##' @title Return Levels associated with Chosen Periods
##'
##' @param object An object defining \emph{return levels} associated
##' to some period values.
##'
##' @param ... Arguments passed to methods.
##'
##' @return An object inheriting from the \code{"data.frame"} class
##' containing the return levels.
##' 
RL <- function(object, ...) {   
    UseMethod("RL", object)
}



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
##' containing the return periods.
##' 
RP <- function(object, ...) {   
    UseMethod("RP", object)
}

## ****************************************************************************
##' Maximum-Likelihood estimation of a (usually incomplete) object.
##'
##' @title Maximum-Likelihood Estimation
##'
##' @param object An object defining a statistical model, usually
##' incomplete because the estimation has not yet been performed.
##'
##' @param ... Arguments passed to methods.
##'
##' @return A list with the results of the estimation. This include a
##' \code{estimate} element and usually a \code{cov} element.
##' 
MLE <- function(object, ...) {   
    UseMethod("MLE", object)
}
