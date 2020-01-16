
## ***********************************************************************
##' Density, distribution function, quantile function and random
##' generation for the two-parameter Generalized Pareto Distribution
##' (GPD) distribution with \code{scale} and \code{shape}.
##'
##' @name GPD2
##' @rdname GPD2
##'
##' @aliases pGPD2 qGPD2 rGPD2
##' 
##' @title Density, Distribution Function, Quantile Function and
##' Random Generation for the Two-Parameter Generalized Pareto
##' Distribution (GPD)
##'
##' @usage
##' dGPD2(x, scale = 1.0, shape = 0.0, log = FALSE,
##'       deriv = FALSE, hessian = FALSE)
##' 
##' pGPD2(q, scale = 1.0, shape = 0, lower.tail = TRUE,
##'       deriv = FALSE, hessian = FALSE) 
##'
##' qGPD2(p, scale = 1.0, shape = 0.0, lower.tail = TRUE,
##'       deriv = FALSE, hessian = FALSE)
##' 
##' rGPD2(n, scale = 1.0, shape = 0.0)
##' 
##' @param scale Scale parameter. Numeric vector with suitable length,
##' see \bold{Details}.
##'
##' @param shape Shape parameter. Numeric vector with suitable length,
##' see \bold{Details}.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##' returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##' computed value w.r.t. the parameter vector is computed, and
##' returned as a \code{"gradient"} attribute of the result. This is a
##' numeric array with dimension \code{c(n, 2)} where \code{n} is the
##' length of the first argument, i.e. \code{x}, \code{p} or \code{q},
##' depending on the function.
##'
##' @param hessian Logical. If \code{TRUE}, the Hessian of each
##' computed value w.r.t. the parameter vector is computed, and
##' returned as a \code{"hessian"} attribute of the result. This is a
##' numeric array with dimension \code{c(n, 2, 2)} where \code{n} is
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
##' four lengths: that of the first argument and that of the two
##' parameters \code{scale} and \code{shape}. When \code{deriv} is
##' \code{TRUE}, the returned value has an attribute named
##' \code{"gradient"} which is a matrix with \eqn{n} lines and \eqn{2}
##' columns containing the derivatives. A row contains the partial
##' derivatives of the corresponding element w.r.t. the two parameters
##' \code{"scale"} and \code{"shape"} in that order.
##'
##' @details The probability functions \code{d}, \code{p} and \code{q}
##' all allow each of the two GP parameters to be a vector. Then the
##' recycling rule is used to get three vectors of the same length,
##' corresponding to the first argument and to the two GP
##' parameters. This behaviour is the standard one for the probability
##' functions of the \strong{stats} package but is unusual in R
##' packages devoted to Extreme Value in which the parameters must
##' generally have length one. Note that these functions can be used
##' e.g.  to evaluate the quantile with a given probability for a
##' large number of values of the parameter vector \code{c(shape,
##' scale)}. This is frequently required in he Bayesian framework with
##' MCMC inference.
##'
##' @note The attributes \code{"gradient"} and \code{"hessian"} have
##' dimension \code{c(n, 2)} and \code{c(n, 2, 2)}, even when \code{n}
##' equals \code{1}. Use the \code{drop} method on these objects to
##' drop the extra dimension if wanted i.e. to get a gradient vector
##' and a Hessian matrix.
##' 
##' @author Yves Deville
##'
##' @examples
##' ## Illustrate the effect of recycling rule.
##' pGPD2(1.0, scale = 1:4, shape = 0.0, lower.tail = FALSE) - exp(-1.0 / (1:4))
##' pGPD2(1:4, scale = 1:4, shape = 0.0, lower.tail = FALSE) - exp(-1.0)
##'
##' ## With gradient and Hessian.
##' pGPD2(c(1.1, 1.7), scale = 1, shape = 0, deriv = TRUE, hessian = TRUE)
dGPD2 <- function(x, scale = 1.0, shape = 0.0, log = FALSE,
                  deriv = FALSE, hessian = FALSE) {

    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    res <- .Call(Call_dGPD2,
                 as.double(x),
                 as.double(scale),
                 as.double(shape),
                 as.integer(log),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(rownames(x), nm2))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(rownames(x), nm2, nm2))
        }
        
    }

    return(res)
    
}

pGPD2 <- function (q, scale = 1.0, shape = 0, lower.tail = TRUE,
                   deriv = FALSE, hessian = FALSE) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }

    res <- .Call(Call_pGPD2,
                 as.double(q),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(rownames(q), c("scale", "shape")))
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(rownames(q), nm2, nm2))
        }
    }
    return(res)
    
}

qGPD2 <- function (p, scale = 1.0, shape = 0.0, lower.tail = TRUE,
                   deriv = FALSE, hessian = FALSE) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
    
    res <- .Call(Call_qGPD2,
                 as.double(p),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    n <- length(res)
    if (deriv) {
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(rownames(p), nm2))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(rownames(p), nm2, nm2))
        }
    }
    return(res)   

}

rGPD2 <- function (n, scale = 1.0, shape = 0.0) {
    
    if (any(is.na(scale)) || any(scale <= 0) || !all(is.finite(scale))) {
        stop("'scale' must contain non-NA finite and positive numeric values")  
    }
    if (any(is.na(shape)) || !all(is.finite(shape))) {
        stop("'shape' must contain non-NA finite numeric values")  
    }

    pGPD2(runif(n), scale = scale, shape = shape)
    
}

