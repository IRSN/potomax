
## ******************************************************************************

##' List of distibutions
##'
##' The definitions could easily be put into a loop.
##'
##' o The first argument is the classical argument of 'd', 'p', q' and
##' 'r' functions. See e.g., \code{link[stats]{Exponential}} distribution.
##'
##' o The second argument is named \code{theta} and all the
##' parameters in the right order.
##'
##' @section Caution: The position of the parameters in the vector
##'     \code{theta} may matter for some uses for which \code{theta}
##'     can be unnamed. The first parameter is always assumed to be a
##'     scale parameter.
##'
##' @import nieve
##' 
##' @noRd
Excd <- list("GPD2" = list(p = 2,
                           parNames = c("scale", "shape")),
             "exp1" = list(p = 1,
                           parNames = c("scale")))

## ==============================================================================
## GPD2 distribution
## ==============================================================================

Excd[["GPD2"]]$dFun <- function(x, theta, log = FALSE,
                                deriv = FALSE, hessian = FALSE) {
    nieve::dGPD2(x = x, scale = theta[1], shape = theta[2],
                 log = log,
                 deriv = deriv, hessian = hessian)
}

Excd[["GPD2"]]$pFun <- function(q, theta, lower.tail = TRUE, ## log.p = FALSE,
                                deriv = FALSE, hessian = FALSE) {
    nieve::pGPD2(q = q , scale = theta[1], shape = theta[2],
                 lower.tail = lower.tail,
                 ## log = log,
                 deriv = deriv, hessian = hessian)
}

Excd[["GPD2"]]$qFun <- function(p, theta, lower.tail = TRUE, ## log.p = FALSE,
                            log = FALSE, deriv = FALSE, hessian = FALSE) {
    nieve::qGPD2(p = p , scale = theta[1], shape = theta[2],
                 lower.tail = lower.tail,
                 deriv = deriv, hessian = hessian)
}

Excd[["GPD2"]]$rFun <- function(n, dist, theta) {
    nieve::rGPD2(n = n , scale = theta[1], shape = theta[2])
}

## ==============================================================================
## exponential distribution
## ==============================================================================

Excd[["exp1"]]$dFun <- function(x, theta, log = FALSE,
                                deriv = FALSE, hessian = FALSE) {
    nieve::dexp1(x = x, scale = theta[1],
                 log = log,
                 deriv = deriv, hessian = hessian)
}

Excd[["exp1"]]$pFun <- function(q, theta, lower.tail = TRUE, ## log.p = FALSE,
                                deriv = FALSE, hessian = FALSE) {
    nieve::pexp1(q = q , scale = theta[1],
                 lower.tail = lower.tail,
                                        # log = log,
                 deriv = deriv, hessian = hessian)
}

Excd[["exp1"]]$qFun <- function(p, theta, lower.tail = TRUE, ## log.p = FALSE,
                                log = FALSE, deriv = FALSE, hessian = FALSE) {
    nieve::qexp1(p = p , scale = theta[1],
                 lower.tail = lower.tail,
                 deriv = deriv, hessian = hessian)
}

Excd[["exp1"]]$rFun <- function(n, theta) {
    nieve::rexp1(n = n , scale = theta[1])   
}
