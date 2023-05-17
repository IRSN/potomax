##*****************************************************************************
##' Confidence intervals for the parameters of a \code{poisGP} object.
##'
##' This method finds confidence intervals for each of the three
##' parameters of the chosen parameterisation: \code{"poisGP"}
##' (default) or \code{"PP"}. The recommended (and default) method
##' uses profile-likelihood as implemented in \code{proflik.default}.
##' The determination of the intervals relies on a new method based on
##' constrained optimisation: thus the profiled likelihood is not
##' computed as such, as opposed to what is usually done. The profiled
##' likelihood can be computed to check the results by using
##' \code{check = TRUE}, but no zero finding is used to find the
##' confidence limits following the classical method. The check is
##' thus entirely based on the graphics which must be carefully
##' inspected.
##' 
##' @title Confidence Intervals for a \code{poisGP} Object
##' 
##' @usage 
##' \method{confint}{poisGP}(object,
##'         parm = NULL,
##'         level = 0.95,
##'         type = c("poisGP", "PP"),
##'         method = c("proflik", "delta"),
##'         nSigma = 4,
##'         trace = 0,
##'         round = TRUE,
##'         out = c("array", "data.frame"),
##'         check = FALSE, nCheck = 50,
##'         ...) 
##' 
##' @param object An object with class \code{"poisGP"}.
##'
##' @param parm Not used yet. Confidence intervals are computed for
##' each of the three parameters.
##' 
##' @param level Confidence level(s).
##' 
##' @param type The type of parameterisation wanted: \emph{Poisson-GP}
##' or \emph{Point-Process}.
##'
##' @param method Character: \code{"delta"} leads to the simplistic
##' \emph{delta} method and \code{"proflik"} to the
##' \emph{profile-likelihood}.
##' 
##' @param nSigma Used only when \code{check} is \code{TRUE}. It
##' defines the interval around the ML estimate where the profile
##' log-likelihood will be evaluated in order to build a curve
##' allowing a check of the results. The value of \code{nSigma}
##' defines the number of standard deviation for the current parameter
##' that will be used. If needed, an an asymmetric interval can be
##' defined by using two numbers e.g. \code{c(3, 5)} if it is expected
##' that the confidence intervals spread more on the right side of the
##' estimates than they do on the left side.
##' 
##' @param trace Integer level of verbosity.
##' 
##' @param round Logical. If \code{TRUE} the confidence limits will be
##' rounded to a same small number of digits. This number is chosen
##' using the smallest of the standard deviations for the estimated
##' parameters.
##'
##' @param out Character giving the class of the output. By default
##' this is a three-dimensional array with dimensions:
##' \emph{parameter}, \emph{lower/upper} limit, and \emph{level}. If
##' \code{level} has length \code{1}, using \code{drop} on the output
##' will remove the third dimension and produce a matrix. The
##' \code{"data.frame"} gives the same results in "long
##' format".
##' 
##' @param check Logical. Used only when \code{method} is
##' \code{"proflik"}. If \code{TRUE} the function return results
##' intended to be used in a graphical check of the confidence limits
##' and taking the form of a list of two data frames. The first data
##' frame contains evaluations of the profile-negative log-likelihood
##' for each of the three parameters in order to draw
##' profile-likelihood curves. The second one contains the confidence
##' bounds in "long format".
##'
##' @param nCheck Number of evaluations of the profile log-likelihood if
##' \code{check} is \code{TRUE}.
##'
##' @param ... Further arguments passed to the \code{\link{profLik}}
##' method for the class of \code{object}. The arguments
##' \code{ftol_abs} and \code{ftol_rel} can be modified when some
##' problems are met.
##'
##' @return When \code{check} is \code{FALSE}, an array or a data
##' frame containing the lower and upper bounds \code{"L"} and
##' \code{"U"} of the confidence intervals. When \code{check} is
##' \code{TRUE} a list of two data frames which is given the class
##' \code{"confintCheck"} is order to use the \code{autoplot} method
##' that is implemented for this class.
##'
##' @seealso \code{\link{RL.poisGP}} for the computation of the return
##' levels with confidence intervals,
##' \code{\link{autoplot.confintCheck.poisGP}} for the graphical check
##' of the results. The \code{profLik.default} method is used by this
##' function.
##'
##' @note Remind that the \code{"PP"} parameterisation does not depend
##' on the threshold, as opposed to the \code{"poisGP"}
##' parameterisation. So \code{type = "PP"} should be used to
##' investigate \emph{threshold stability} for the full parameter
##' vector. However the confidence intervals on the two shape
##' parameter: Poisson-GP \eqn{\xi} and PP \eqn{\xi^\star}{\xi*} are
##' (or should be) identical.
##' 
##' @section Caution: The determination of the profile-likelihood
##' intervals can fail, so it is wise to set \code{check = TRUE} and
##' use the \code{autoplot} method on the returned object. Problems
##' seem to be more frequently met with \code{type = "PP"}, i.e. when
##' the Point-Process parameterisation is used.
##'
##' @importFrom stats qnorm qchisq
##' @import nloptr 
##' @method confint poisGP
##' @export
##' 
##' @examples
##' ## fit from the object Garonne from Renext (class "Rendata")
##' fit <- poisGP(Garonne, threshold = 2900)
##'
##' ci <- confint(fit, lev = c(0.70, 0.95), trace = 1)
##'
##' ## Check the results: this is quite time-consuming.
##' \dontrun{
##' cic <- confint(fit, lev = c(0.95, 0.70), check = TRUE)
##' autoplot(cic) + theme_gray() +
##'     ggtitle("Poisson-GP parameterisation")
##' 
##' cicPP <- confint(fit, type = "PP", lev = c(0.95, 0.70), check = TRUE)
##' autoplot(cicPP) + theme_gray() +
##'     ggtitle("Point-Process (PP) parameterisation")
##' }
confint.poisGP <- function(object,
                           parm = NULL,
                           level = 0.95,
                           type = c("poisGP", "PP"),
                           method = c("proflik", "delta"),
                           nSigma = 4,
                           trace = 0,
                           round = TRUE,
                           out = c("array", "data.frame"),
                           check = FALSE,
                           nCheck = 50,
                           ...) {
    
    out <- match.arg(out)
    method <- match.arg(method)
    type <- match.arg(type)
    p <- object$p
    parNames <- object$parNames

    if (method == "delta" && check) {
        warning("Since method is \"delta\", no check is needed and 'check' is ",
                "set to FALSE")
        check <- FALSE
    }
    
    if (check) {
        message("Use the 'autoplot' method on the result to check the results")
    }
    
    if (check && out == "array") {
        warning("Since 'check' is TRUE, 'out' set to \"data.frame\"")
        out <- "data.frame"
    } 

    ## ========================================================================
    ## Take into account the order of the levels and format them
    ## suitably.
    ## ========================================================================
    
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    ## ========================================================================
    ## Use lists to alternatively play with the two parameterisations
    ## "poisGP" and "PP".
    ## ========================================================================
    
    thetaHat <- sigHat <- pNames <- list()
    
    for (typi in c("poisGP", "PP")) {
        thetaHat[[typi]] <- coef(object, type = typi)
        sigHat[[typi]] <- sqrt(diag(vcov(object, type = typi)))
        pNames[[typi]] <- names(thetaHat[[typi]])
    }

    ## ========================================================================
    ## No great difficulty with the 'delta' method: select the
    ## suitable element in lists 'thetaHat', 'sigHat' and 'pNames'.
    ## ========================================================================

    if (method == "delta") {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        
        q <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)
        
        ci <- array(thetaHat[[type]], dim = c(p, 2L, nLevel),
                    dimnames = list(pNames[[type]], c("L", "U"), fLevel))
        
        cw <-  array(sigHat[[type]], dim = c(p, 2L, nLevel),
                     dimnames = list(pNames[[type]], c("L", "U"), fLevel)) 
        
        cw <- sweep(cw, MARGIN = c(p, 2), STATS = q, FUN = "*")
        ci <- ci + cw
        
    } else if (method == "proflik") {
                
        nSigma <- rep(nSigma, length.out = 2L)
        prob <- 1 - level
        
        cipoisGP <- array(NA, dim = c(p, 2L, nLevel),
                          dimnames = list(pNames[["poisGP"]], c("L", "U"),
                              fLevel)) 
        
        ci <- array(NA, dim = c(p, 2L, nLevel),
                    dimnames = list(pNames[[type]], c("L", "U"), fLevel)) 
        
        ## ====================================================================
        ## First use the standard 'profLik' method for each parameter
        ## with the functions dedicated to "poisGP". If type ==
        ## "poisGP" we are done! Else, we will use the results to the
        ## bounds on the 'poiGP' parameters, because it semms to help
        ## much for the profLik with the "PP" parameters.
        ## =====================================================================

        if (trace) {
            cat("\no Perform profile-likelihood for the \"poisGP\" ",
                "parameterisation.\n\n") }
        
        for (k in 1:p) {

            if (trace > 0) {
                cat("\no Parameter ", pNames[[type]][k],
                    "\n*******************\n")
            }

            myfun <- function(theta, object) {
                res <- theta[k]
                grad <- rep(0.0, p)
                grad[k] <- 1.0
                attr(res, "gradient") <- grad
                res
            }
            
            pl <- profLik(object, fun = myfun, level = level,
                          trace = trace, ...)
            attr(pl, "diagno") <- attr(pl, "theta") <- NULL
            cipoisGP[k, , ] <- pl[c("L", "U"), ]   
        }

        if (trace) {
            cat("\no Results for \"poisGP\"\n\n")
            print(cipoisGP)
        }

        
        if (type == "PP") {

            if (trace) {
                cat("\no Perform profile-likelihood for the \"PP\" ",
                    "parameterisation.\n")
            }
            
            for (k in 1:2) {
                
                if (trace > 0) {
                    cat("\no Parameter ", pNames[[type]][k],
                        "\n*******************\n")
                }
                
                ## =============================================================
                ## Define the function to be profiled. Note that this
                ## is a function of the "poisGP" parameter as is
                ## always true with our 'profLik' method.
                ## =============================================================
                
                myfun <- function(theta, object) {
                    
                    if ((theta[1] <= 0.0) || (theta[2] <= 0.0)) {
                        cat("PB\n")
                        print(theta) 
                    }

                    shape <- ifelse(p == 3, theta[3], 0.0)
                                  
                    ## XXXY
                    thetaStar <-
                        nieve::poisGP2PP(lambda = theta[1],
                                         scale = theta[2],
                                         shape = shape,
                                         loc = object$threshold, deriv = TRUE)
                    res <- thetaStar[k]
                    
                    ## In 3-rd dim we always have to remove the the GPD location,
                    ## an possibly the GPD shape
                    if (p == 3) {
                        grad <- attr(thetaStar, "gradient")[1, k, c(1, 3, 4)]
                    } else  {
                        grad <- attr(thetaStar, "gradient")[1, k, c(1, 3)] 
                    }
                    attr(res, "gradient") <- grad
                    res
                }

                ## =============================================================
                ## Perform profile-likelihood for the "PP" fun with
                ## the box bounds on the "poisGP" parameters set to
                ## their confidence limits found previously. This
                ## obviously must be done level by level. For
                ## instance, we know that by maxi/minimising the "PP"
                ## location 'muStar' in the 95% confidence region, the
                ## "poisGP" parameters will remain in their respective
                ## 95% 'marginal' confidence interval.
                ## =============================================================

                for (ilev in seq_along(level)) {
                    
                    object2 <- object
                    
                    object2$lb <- cipoisGP[ , "L", ilev]
                    object2$ub <- cipoisGP[ , "U", ilev]
                    plilev <- profLik(object = object2, fun = myfun,
                                      level = level[ilev], trace = trace)         
                    attr(plilev, "diagno") <- attr(plilev, "theta") <- NULL
                    ci[k, , ilev] <- plilev[c("L", "U"), 1]
                }
                
            }

            ## =================================================================
            ## For the shape parameter, no profiling is needed because
            ## the parameter is identical to the 'poisGP' shape!!!
            ## =================================================================

            if (p == 3) {
                ci[3L, , ] <- cipoisGP[3L, , ]
            }
            
        } else {
            ci <- cipoisGP
        }
        
        ## ====================================================================
        ## If a check is needed, then perform a series of
        ## optimisation(s)
        ## ====================================================================
        
        if (check) {
            
            if (type == "PP") {
                
                ## ============================================================
                ## Define the negLogLik to be minimised for the
                ## parameter number 'k'. This is a function of a
                ## vector of length 2: the PP parameter with its k-th
                ## element removed.
                ## ============================================================
                
                negLogLikNok <- function(thetaNok, k, thetak, deriv = TRUE) {
                    
                    theta <- rep(NA, p)
                    theta[-k] <- thetaNok
                    theta[k] <- thetak
                    if (p == 3) {
                        theta <-
                            try(nieve::PP2poisGP(locStar = theta[1L],
                                                 scaleStar = theta[2L],
                                                 shapeStar = theta[3L],
                                                 threshold = object$threshold,
                                                 deriv = deriv), silent = TRUE)
                    } else {
                        theta <-
                            try(nieve::PP2poisGP(locStar = theta[1L],
                                                 scaleStar = theta[2L],
                                                 shapeStar = 0.0,
                                                 threshold = object$threshold,
                                                 deriv = deriv), silent = TRUE)
                        
                        ## Caution here. The attributes are lost when indexing []
                        if (deriv) {
                            g <- attr(theta, "gradient")[ , -3, -3, drop = FALSE]
                            theta <- theta[-3]
                            attr(theta, "gradient") <- g
                        } else {
                            theta <- theta[-3]
                        }
                    }
   
                    ## ========================================================
                    ## Remind that an error occurs when the support of
                    ## the GEV distribution does not contain the
                    ## threshold as an interior point.
                    ## ========================================================
                    
                    if (inherits(theta, "try-error")) {
                        if (deriv) {
                            res <- list("objective" = Inf, 
                                        "gradient" = rep(NaN, 2))
                        } else res <- Inf
                    }
                    
                    if (all(is.finite(theta))) {
                        
                        nLL <-  object$negLogLikFun(theta, object = object,
                                                    deriv = deriv)
                        if (deriv) {
                            grad <-  attr(nLL, "gradient") %*%
                                attr(theta, "gradient")[1, , -k]
                            attr(nLL, "gradient") <- NULL
                            attr(nLL, "Cst") <- NULL
                            res <- list("objective" = nLL, "gradient" = grad)
                        } else {
                            attr(nLL, "gradient") <- NULL
                            attr(nLL, "Cst") <- NULL
                            res <- nLL
                        }
                    } else {
                        if (deriv) {
                            res <- list("objective" = NaN,
                                        "gradient" = rep(NaN, p - 1L))
                        } else res <- NaN
                    }
                    res
                    
                }
                
            } else if (type == "poisGP") {
                
                ## =============================================================
                ## Define the negLogLik to be minimised for the
                ## parameter numer 'k'. This is a function of a vector
                ## of length 2: the poisGP parameter with its k-th
                ## element removed.
                ## =============================================================
                
                negLogLikNok <- function(thetaNok, k, thetak, deriv = TRUE) {
                    
                    theta <- rep(NA, p)
                    theta[-k] <- thetaNok
                    theta[k] <- thetak
                    
                    if (all(is.finite(theta))) {
                        nLL <-  object$negLogLikFun(theta, object = object,
                                                    deriv = deriv)
                        if (deriv) {
                            grad <- attr(nLL, "gradient")
                            grad <- grad[1L, -k, drop = FALSE]
                            attr(nLL, "gradient") <- attr(nLL, "Cst") <- NULL
                            res <- list("objective" = nLL, "gradient" = grad)
                        } else {
                            attr(nLL, "gradient") <- NULL
                            attr(nLL, "Cst") <- NULL
                            res <- nLL
                        }
                    } else {
                        if (deriv) {
                            res <- list("objective" = NaN,
                                        "gradient" = rep(NaN, p - 1L))
                        } else res <- NaN
                    }
                    res
                    
                }
            }
            
            ## ================================================================
            ## Parameters to tune the optimisation. This settings have
            ## a considerable impact on the availability and even on
            ## the precision of the profile-likelihood value. 
            ##
            ## o 'ftol_rel' seems to improve but must be set to a very
            ## small value.
            ##
            ## Finding the best settings will need some work in the
            ## future...
            ## ================================================================
            
            optsNok <- list("algorithm" = "NLOPT_LN_COBYLA",
                            "xtol_rel" = 1.0e-8,
                            ## "xtol_abs" = 1.0e-8,
                            "ftol_abs" = 1e-8,
                            "ftol_rel" = 1e-8,
                            "maxeval" = 3000, "print_level" = 0,
                            "check_derivatives" = FALSE)
            
            optsNokDeriv <- list("algorithm" = "NLOPT_LD_LBFGS",
                                 "xtol_rel" = 1.0e-8,
                                 ## "xtol_abs" = 1.0e-8,
                                 "ftol_abs" = 1e-12,
                                 "ftol_rel" = 1e-12,
                                 "maxeval" = 8000, "print_level" = 0,
                                 "check_derivatives" = FALSE)
            
            for (k in 1:p) {
                
                thetaNok <- thetaHat[[type]][-k]
                
                thetaGridk <- seq(from = thetaHat[[type]][k] -
                                      nSigma[1L] * sigHat[[type]][k],
                                  to = thetaHat[[type]][k] +
                                      nSigma[2L] * sigHat[[type]][k],
                                  length.out = nCheck)
                negLogLikCk <- rep(NA, length(thetaGridk))
                
                status <- rep(NA,  length(thetaGridk))
                useGrad  <- rep(FALSE,  length(thetaGridk))
                
                ## ============================================================
                ## Manage bounds on the parameters. When 'type' is
                ## "poisGP" we can use the parameters set by the user
                ## at the creation, but this not possible in the "PP"
                ## case, except for the shape parameter. Note that we
                ## typically consider here parameter values that are
                ## outside of the confidence intervals.
                ## ============================================================
                
                if (type == "poisGP") {
                    lbk <- object$lb[-k]
                    ubk <- object$ub[-k]
                } else {
                    lbk <- c("loc" = -Inf, "scale" = 0.0)
                    ubk <- c("loc" = Inf, "scale" = Inf)
                    if (p == 3) {
                        lbk <- c(lbk, "shape" = unname(object$lb)[3])
                        ubk <- c(ubk, "shape" = unname(object$ub)[3])
                    } 
                    lbk <- lbk[-k]
                    ubk <- ubk[-k]
                }
                
                ## ============================================================
                ## Loop over all grid values. For each value we first
                ## try an optimisation with derivative. If this fail
                ## we try a new optimisation without derivative.
                ## ============================================================
                
                for (i in seq_along(thetaGridk)) {
                    
                    thetak <- thetaGridk[i]
                    thetaNok <- thetaHat[[type]][-k]
                    
                    resii <-  try(nloptr(x0 = thetaNok,
                                         eval_f = negLogLikNok,
                                         opts = optsNokDeriv,
                                         lb = lbk,
                                         ub = ubk,
                                         k = k,
                                         thetak = thetak,
                                         deriv = TRUE), silent = TRUE)
                    
                    if (inherits(resii, "try-error")) print(resii)
                    
                    if (!inherits(resii, "try-error") &&
                        (resii$status %in% 1:4)) {
                        negLogLikCk[i] <- resii$objective
                        useGrad[i] <- TRUE
                    } else {
                        resii <-  try(nloptr(x0 = thetaNok,
                                             eval_f = negLogLikNok,
                                             opts = optsNokDeriv,
                                             lb = lbk,
                                             ub = ubk,
                                             k = k,
                                             thetak = thetak,
                                             deriv = FALSE), silent = TRUE)
                        
                        if (!inherits(resii, "try-error") &&
                            (resii$status %in% 1:4)) {
                            negLogLikCk[i] <- resii$objective
                        }
                    }
                    
                    if (!inherits(resii, "try-error")) {
                        if (trace > 1) cat("status = ", resii$status, "\n")
                        status[i] <-  resii$status
                    } else {
                        if (trace > 1) cat("Optim error\n")
                        ## print(resii)
                    }
                }
                
                ## ============================================================
                ## define the data frame 'negLogLikC' containing the
                ## value of the profiled negative log-likelihood or
                ## add rows to this dataframe if k > 1
                ## ============================================================
                
                if (k == 1) {
                    negLogLikC <- data.frame(Name = rep(pNames[[type]][k],
                                                 nCheck),
                                             Par = thetaGridk,
                                             Value = negLogLikCk,
                                             stringsAsFactors = FALSE)
                    
                } else {
                    negLogLikC <-
                        dplyr::bind_rows(negLogLikC,
                                         data.frame(Name = rep(pNames[[type]][k],
                                                        nCheck),
                                                    Par = thetaGridk,
                                                    Value = negLogLikCk,
                                                    stringsAsFactors = FALSE)) 
                }      
            }
        }
    }


    if (round && (!any(is.na(sigHat[[type]])))) {
        d <- ceiling(-log(min(sigHat[[type]]), 10)) + 1
        ci <- round(ci, digits = d)
    }
    
    if (out == "data.frame") {
        
        ci <- reshape2::melt(ci, value.name = "Value",
                             varnames = c("Name", "LU", "Level"))
        nll <- object$negLogLik + qchisq(level, df = 1) / 2.0
        names(nll) <- fLevel
        ci <- cbind(ci, NegLogLik = nll[as.character(ci$Level)])
    }
    
    if (check) {
        
        ciPlus <- data.frame(Name = pNames[[type]],
                             LU = rep("est", p),
                             Level = rep("est", p),
                             Value = thetaHat[[type]],
                             NegLogLik = rep(-object$logLik, p))
        
        ci <- rbind(ci, ciPlus, deparse.level = 1)
            
        L <- list(ci = ci, negLogLikC = negLogLikC)
        class(L) <- "confintCheck.poisGP"
        L
    } else {
        ci
    }
    
}
