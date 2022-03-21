
## ***********************************************************************
##' Density, distribution function, quantile function and random
##' generation for the one-parameter Exponential Distribution
##' distribution with scale parameter \code{scale}.
##'
##' @name Exp1
##' @rdname Exp1
##'
##' @aliases pexp1 qexp1 rexp1
##' 
##' @title Density, Distribution Function, Quantile Function and
##' Random Generation for the One-Parameter Exponential Distribution
##'
##' @usage
##' dexp1(x, scale = 1.0, log = FALSE,
##'       deriv = FALSE, hessian = FALSE)
##' 
##' pexp1(q, scale = 1.0, lower.tail = TRUE,
##'       deriv = FALSE, hessian = FALSE) 
##'
##' qexp1(p, scale = 1.0, lower.tail = TRUE,
##'       deriv = FALSE, hessian = FALSE)
##' 
##' rexp1(n, scale = 1.0)
##' 
##' @param scale Scale parameter. Numeric vector with suitable length,
##' see \bold{Details}.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##' returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##' computed value w.r.t. the parameter vector is computed, and
##' returned as a \code{"gradient"} attribute of the result. This is a
##' numeric array with dimension \code{c(n, 1)} where \code{n} is the
##' length of the first argument, i.e. \code{x}, \code{p} or \code{q},
##' depending on the function.
##'
##' @param hessian Logical. If \code{TRUE}, the Hessian of each
##' computed value w.r.t. the parameter vector is computed, and
##' returned as a \code{"hessian"} attribute of the result. This is a
##' numeric array with dimension \code{c(n, 1, 1)} where \code{n} is
##' the length of the first argument, i.e. \code{x}, \code{p} or
##' depending on the function. 
##'
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
##' P[X <= x], otherwise, P[X > x].
##'
##' @return A numeric vector with length equal to the maximum of the
##'     two lengths: that of the first argument and that of the
##'     parameter \code{scale}. When \code{deriv} is \code{TRUE}, the
##'     returned value has an attribute named \code{"gradient"} which
##'     is a matrix with \eqn{n} lines and \eqn{1} column containing
##'     the derivative. A row contains the partial derivative of the
##'     corresponding element w.r.t. the parameter \code{"scale"}.
##'
##' @details The probability functions \code{d}, \code{p} and \code{q}
##'     all allow the parameter \code{scale} to be a vector. Then the
##'     recycling rule is used to get two vectors of the same length,
##'     corresponding to the first argument and to the scale
##'     parameter. This behaviour is the standard one for the
##'     probability functions of the \strong{stats} package but is
##'     unusual in R packages devoted to Extreme Value in which the
##'     parameters must generally have length one. Note that the
##'     provided functions can be used e.g. to evaluate the quantile
##'     with a given probability for a large number of values of the
##'     parameter vector \code{shape}. This is frequently required in
##'     he Bayesian framework with MCMC inference.
##'
##' @note The attributes \code{"gradient"} and \code{"hessian"} have
##'     dimension \code{c(n, 1)} and \code{c(n, 1, 1)}, even when
##'     \code{n} equals \code{1}. Use the \code{drop} method on these
##'     objects to drop the extra dimension if wanted i.e. to get a
##'     gradient vector and a Hessian matrix.
##' 
##' @author Yves Deville
##'
##' @examples
##' ## Illustrate the effect of recycling rule.
##' pexp1(1.0, scale = 1:4, lower.tail = FALSE) - exp(-1.0 / (1:4))
##' pexp1(1:4, scale = 1:4, lower.tail = FALSE) - exp(-1.0)
##'
##' ## With gradient and Hessian.
##' pexp1(c(1.1, 1.7), scale = 1, deriv = TRUE, hessian = TRUE)
dexp1 <- function(x, scale = 1.0, log = FALSE,
                  deriv = FALSE, hessian = FALSE) {

    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    res <- .Call(Call_dexp1,
                 as.double(x),
                 as.double(scale),
                 as.integer(log),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(rownames(x), nm1))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(rownames(x), nm1, nm1))
        }
        
    }

    return(res)
    
}

pexp1 <- function (q, scale = 1.0, lower.tail = TRUE,
                   deriv = FALSE, hessian = FALSE) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }

    res <- .Call(Call_pexp1,
                 as.double(q),
                 as.double(scale),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(rownames(q), c("scale")))
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(rownames(q), nm1, nm1))
        }
    }
    return(res)
    
}

qexp1 <- function (p, scale = 1.0, lower.tail = TRUE,
                   deriv = FALSE, hessian = FALSE) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
    
    res <- .Call(Call_qexp1,
                 as.double(p),
                 as.double(scale),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(rownames(p), nm1))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(rownames(p), nm1, nm1))
        }
    }
    return(res)   

}

rexp1 <- function (n, scale = 1.0) {
    
    if (any(is.na(scale)) || any(scale <= 0) || !all(is.finite(scale))) {
        stop("'scale' must contain non-NA finite and positive numeric values")  
    }

    pexp1(runif(n), scale = scale)
    
}

