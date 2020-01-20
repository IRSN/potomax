

## ****************************************************************************
##' \code{RL} is a generic function which computes the return levels
##' associated with chosen periods using a suitable model or data set.
##'
##' @title Return Levels associated with Chosen Periods
##'
##' @param object An object that can be used to define \emph{return
##' levels} associated to some period values.
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
##' \code{RP} is a generic function which computes the return periods
##' associated with chosen return levels or quantile, using a suitable model or
##' data set.
##'
##' @title Return Periods associated with Chosen Return Levels or
##' Quantiles
##'
##' @param object An object that can be used to define \emph{return
##' periods} associated to some "quantiles" values.
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
##' \code{MLE} is a generic function performing Maximum-Likelihood
##' estimation of a (usually incomplete) object.
##'
##' @title Maximum-Likelihood Estimation
##'
##' @param object An object defining a statistical model, usually
##' incomplete because the estimation has not yet been performed.
##'
##' @param ... Arguments passed to methods.
##'
##' @return A list with the results of the estimation. This include an
##' \code{estimate} element, and usually a \code{cov} element. These
##' should be extracted by the \code{coef} and \code{vcov} methods.
##'
##' @seealso \code{\link[stats]{coef}} and \code{\link[stats]{vcov}}.
MLE <- function(object, ...) {   
    UseMethod("MLE", object)
}
