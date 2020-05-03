## ****************************************************************************
##' Transform Poisson-GP parameters into Point-Process (PP) parameters. 
##'
##' @details In the POT framework the three parameters are the rate
##' \code{lambda} of the Poisson process in time: \eqn{\lambda}, and
##' the two GP parameters: \code{scale} \eqn{\sigma} and \code{shape}
##' \eqn{\xi}. The vector \code{loc} contains the fixed threshold and
##' \code{w} the fixed block duration. These parameters are converted
##' into the vector of three parameters of the GEV distribution for
##' the maximum of the marks \eqn{Y_i}{Yi} on a time interval with
##' duration \code{w}, the number \eqn{N} of these marks being a
##' r.v. with Poisson distribution. More precisely, the GEV
##' distibution applies when \eqn{N > 0}.
##'
##' @usage
##' poisGP2PP(lambda, loc = 0.0, scale = 1.0, shape = 0.0,
##'           w = 1.0, deriv = FALSE) 
##' 
##' @title Transform Poisson-GP Parameters into Point-Process Parameters
##'
##' @param lambda A numeric vector containing the Poisson rate(s).
##'
##' @param loc A numeric vector containing the Generalised Pareto
##' location, i.e. the threshold in the POT framework.
##'
##' @param scale,shape Numeric vectors containing the Generalised Pareto
##' scale and shape parameters.
##'
##' @param w The block duration. Its physical dimension is time and
##' the product \eqn{\lambda \times w}{\lambda * w} is dimensionless.
##'
##' @param deriv Logical. If \code{TRUE} the derivative (Jacobian) of
##' the transformation is computed and returned as an attribute named
##' \code{"gradient"} of the attribute.
##'
##' @return A matrix with three columns representing the Point-Process
##' parameters \code{loc}, \code{scale} and \code{shape}.
##'
##' @author Yves Deville
##'
##' @note This function is essentially a re-implementation in C of the
##' function \code{\link[Renext]{Ren2gev}} of \bold{Renext}. As a
##' major improvement, this function is "vectorised" w.r.t. the
##' parameters so it can transform efficiently a large number of
##' Poisson-GP parameter vectors as it can be required e.g. in a MCMC
##' Bayesian inference. Note also that this function copes with values
##' near zero for the shape parameter: it suitably computes then both
##' the function value and its derivatives.
##' 
##' @references
##' 
##' Yves Deville (2020). \emph{Renext Computing Details}. Technical
##' Report. (Chap. 3 and Appendix A)
##' 
##' @seealso \code{\link{PP2poisGP}} for the reciprocal
##' transformation.
##' 
poisGP2PP <- function(lambda, loc = 0.0, scale = 1.0, shape = 0.0,
                      w = 1.0,
                      deriv = FALSE) {
    
    if (length(w) != 1) {
        stop("'w' must for now have length one")
    }

    if (any(is.na(lambda)) || any(is.na(loc)) || any(is.na(scale)) ||
        any(is.na(shape)) || is.na(w)) {
        stop("NA are not allowed in 'lambda', 'loc', 'scale', 'shape' or 'w'")
    }

    if (any(lambda <= 0.0) || any(scale <= 0.0)) {
        stop("'lambda' and 'scale' must contain positive values")
    }
    
    res <- .Call(Call_poisGP2PP,
                 as.double(lambda),
                 as.double(loc),
                 as.double(scale),
                 as.double(shape),
                 as.double(w),
                 as.integer(deriv))

    g <- attr(res, "gradient")
    
    n <- length(res) / 3
    nm2 <- c("loc", "scale", "shape")
    res <- array(res, dim = c(n, 3),
                 dimnames = list(NULL, nm2))
    
    if (deriv) {
        nm3 <- c("lambda", "loc", "scale", "shape")
        attr(res, "gradient") <-
            array(g,
                  dim = c(n, 3L, 4L),
                  dimnames = list(NULL, nm2, nm3))
        
        
    }
 
    
    return(res)
    
}
