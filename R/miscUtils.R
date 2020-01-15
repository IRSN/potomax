## *****************************************************************************
##' Distance between two lines.
##' 
##' This function is used to build diagnostics during the
##' construction of profile-likelihood confidence intervals. So the
##' function is not intended to be used by itself.
##' 
##' @title Distance Between Two Lines
##' 
##' @param x1 A numeric vector with positive norm.
##' 
##' @param x2 A numeric vector with positive norm and with the same
##' length as \code{x1}.
##' 
##' @return The distance.
##'
##' @details The vectors \code{x1} and \code{x2} are scaled to have an
##' unit Euclidean norm and to have the same direction. The distance
##' between the lines is simply the distance between the unit vectors
##' on the unit sphere.
##'
##' @note When a constrained optimisation is used for the
##' determination of the lower or upper bound of the interval, we have
##' to check that the gradient of the objective and that of the
##' constraint are colinear, which is performed with this function.
##' 
distLines <- function(x1, x2) {
    if (length(x1) != length(x2)) {
        stop("'x1' and 'x2' must have the same length")
    }
    if (any(is.na(x1)) || any(is.na(x2))) {
        return(NA)
    }
    n1 <- sqrt(sum(x1^2))
    n2 <- sqrt(sum(x2^2))
    if ((n1 <= 1e-6) || (n2 <= 1e-6)) {
        ## warning("'x1' and 'x2' should have Euclidean norm > 1e-6")
        return(NA)
    }
    x1 <- x1 / n1
    x2 <- x2 / n2
    i <- which.max(abs(x1))
    if (x1[i] * x2[i] < 0.0) x2 <- -x2
    s <- sum(x1 * x2)
    if (s > 1 - 1e-9) return(0.0)
    acos(s)
}


formatLevel <- function(level) {
    paste(gsub("\\.0$", "", sprintf("%4.1f", 100 * level)), "%", sep = "")
}

