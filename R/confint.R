##*****************************************************************************
##' Confidence intervals for the parameters of a \code{poisGP} object.
##'
##' @title Confidence Intervals for a \code{poisGP} Object
##'
##' @method confint poisGP
##' 
##' @usage 
##' \method{confint}{poisGP}(object, parm, level = 0.95,
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
##' @param parm Not used yet.
##' 
##' @param level Confidence level(s).
##' 
##' @param method Character: \code{"delta"} leads to the delta method
##' and \code{"proflik"} to the profile-likelihood.
##' 
##' @param nSigma Used only when \code{check} is \code{TRUE}. It
##' defines the interval around the ML estimate where the profile
##' log-likelihood will be evaluated in order to build a curve
##' allowing a check of the results. The value of \code{nSigma}
##' defines the number of standard deviation for the current parameter
##' that will be used. If needed an an asymmetric interval can be
##' defined by using two numbers e.g. \code{c(3, 5)} if it is expected
##' that the confidence intervals spread more on the right side of the
##' estimates than they do on the left side.
##' 
##' @param trace Integer level of verbosity.
##' 
##' @param round Logical. If \code{TRUE} the confidence limits will be
##' rounded to a small number of digits. This number is chosen using the
##' smallest of the standard deviations for the estimated parameters.
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
##' @param ... Not used yet.
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
##' \code{\link{autoplot.confintCheck.poisGP}} for the graphical check of the
##' results.
##' 
##' @examples
##' ## Maybe could we use a 'potMax object here?
##' fit <- poisGP(data = Garonne$OTdata$Flow,
##'               effDuration = Garonne$OTinfo$effDuration,
##'               MAX.data = list("hist" = Garonne$MAXdata$Flow),
##'               MAX.effDuration = Garonne$MAXinfo$duration,
##'               threshold = 2800)
##' 
##' confint(fit, method = "prof", lev = c(0.70, 0.95), trace = 1)
##' cic <- confint(fit, method = "prof", lev = c(0.95, 0.70),
##'                nSigma = 3, check = TRUE, trace = 0)
##' autoplot(cic) + theme_gray()
##' 
confint.poisGP <- function(object,
                           parm, 
                           level = 0.95,
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

    ## take into account the order. 
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    if (method == "delta") {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        thetaHat <- object$estimate
        sigHat <- object$sd
        q <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)

        ci <- array(thetaHat, dim = c(object$p, 2L, nLevel),
                   dimnames = list(object$parNames, c("L", "U"), fLevel)) 
        cw <-  array(sigHat, dim = c(object$p, 2L, nLevel),
                     dimnames = list(object$parNames, c("L", "U"), fLevel)) 

        cw <- sweep(cw, MARGIN = c(3, 2), STATS = q, FUN = "*")
        ci <- ci + cw
        
    } else if (method == "proflik") {
        
        nSigma <- rep(nSigma, length.out = 2L)
        prob <- 1 - level
        thetaHat <- object$estimate
        sigHat <- object$sd
        ci <- array(NA, dim = c(object$p, 2L, nLevel),
                    dimnames = list(object$parNames, c("L", "U"), fLevel)) 

        ## ===================================================================
        ## For each parameter, we maximise/minimise it under the constraint
        ## that the logLik remains >= max logLik - delta where delta :=
        ## qchisq(1 - alpha) where alpha is given by the cofidence level.
        ##
        ## ===================================================================

        opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel" = 1.0e-12,
                      "ftol_abs" = 1.0e-12, "ftol_rel" = 1.0e-12,
                      "maxeval" = 8000,
                      "check_derivatives" = FALSE,
                      "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                          "xtol_rel" = 1.0e-12,
                          "maxeval" = 8000,
                          "ftol_abs" = 1.0e-12,
                          "ftol_rel" = 1.0e-12),
                      "print_level" = 0)

        if (trace >= 2) {
            opts1[["check_derivatives"]] <- TRUE
            opts1[["check_derivatives_print"]] <-  "all"
        }
        if (trace >= 3) {
              opts1[["print_level"]] <- 1
        }

        ## ==============================================================
        ## note that some arguments such as 'level' are unused but are
        ## required by the constraint
        ## ==============================================================
        
        f <- function(theta, k, chgSign = FALSE, level, object) {
            
            grad <- rep(0.0, object$p)
            
            if (chgSign) {
                grad[k] <- -1.0
                return(list("objective" = -theta[k], "gradient" = grad))
            } else {
                grad[k] <- 1.0
                return(list("objective" = theta[k], "gradient" = grad))
            }

        }

        g <- function(theta, k, chgSign = FALSE, level, object) {

            ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(theta = theta, object = object,
                                       deriv = TRUE)
            res2 <- list("constraints" = res - ellL,
                         "jacobian" = attr(res, "gradient"))
            res2 
            
        }
        
        ## =====================================================================
        ## note that although we recompute the gradient of the
        ## objective and the quantile of the chi-square distribution,
        ## this might be faster than re-defining the functions in the
        ## loop. Some experimentations would be needed to confirm
        ## this.
        ## =====================================================================

        
        for (k in 1L:object$p) {
            
            if (trace) cat(sprintf("\n\no Finding CI for \"%s\"\n", object$parNames[k]))

            ilevPrec <- 1L

            if (check) {

                thetaGridk <- seq(from = thetaHat[k] - nSigma[1L] * sigHat[k],
                                  to = thetaHat[k] + nSigma[2L] * sigHat[k],
                                  length.out = nCheck)
                
                negLogLikCk <- rep(NA, nCheck)
                
                negLogLikNok <- function(thetaNok, k, i) {
                    theta <- rep(NA, 3)
                    theta[-k] <- thetaNok
                    theta[k] <- thetaGridk[i]
                    ## cat("XXX", theta, "\n")
                    if (all(is.finite(theta))) {
                        negLogLikFun(theta, object = object, deriv = FALSE)
                    } else {
                        NaN
                    }
                }
                
                lb <- c(0, 0, -0.9)
                ub <- c(Inf, Inf, 2)
                
                optsNok <- list("algorithm" = "NLOPT_LN_COBYLA",
                                "xtol_rel" = 1.0e-8,
                                "xtol_abs" = 1.0e-8,
                                "ftol_abs" = 1e-5,
                                "maxeval" = 1000, "print_level" = 0,
                                "check_derivatives" = FALSE)
                
                for (ik in seq_along(thetaGridk)) {
                    resk <-  try(nloptr(x0 = thetaHat[-k],
                                        eval_f = negLogLikNok,
                                        lb = lb[-k],
                                        ub = ub[-k],
                                        opts = optsNok,
                                        k = k,
                                        i = ik))
                    
                    if (!inherits(resk, "try-error")) {
                        negLogLikCk[ik] <- resk$objective
                    }
                    
                }

                if (k == 1) {
                    negLogLikC <- data.frame(Name = rep(object$parNames[k], nCheck),
                                             Par = thetaGridk,
                                             Value = negLogLikCk,
                                             stringsAsFactors = FALSE)
                    
                } else {
                    negLogLikC <-
                        dplyr::bind_rows(negLogLikC,
                                         data.frame(Name = rep(object$parNames[k], nCheck),
                                                    Par = thetaGridk,
                                                    Value = negLogLikCk,
                                                    stringsAsFactors = FALSE)) 
                }
            }
            
            for (ilev in seq_along(level)) {
                
                lev <- level[ilev]
                
                if (trace) {
                    cat(sprintf("\n   o %s, lower bound: ", fLevel[ilev]))
                }

                ## =========================================================
                ## if we have successfully computed the result for a
                ## larger confidence level (and the same parameter),
                ## use it as initial guess
                ## ==========================================================
                
                if ((ilevPrec > 1L) && (!is.null(thetaLPrec))) {
                    theta0 <- thetaLPrec
                } else {
                    theta0 <- thetaHat
                }

                ## if (k == 3) opts1$print_level <- 3
                
                resL <- try(nloptr::nloptr(x0 = theta0,
                                           eval_f = f,
                                           eval_g_ineq = g,
                                           lb = c(0.0, 0, -0.90),
                                           ub = c(Inf, Inf, Inf),
                                           k = k, level = lev, chgSign = FALSE,
                                           opts = opts1,
                                           object = object),
                            silent = TRUE)

                ## if (k == 3) opts1$print_level <- 0
                
                if (!inherits(resL, "try-error")) {
                    if (trace == 1L) {
                        cat(sprintf("%7.2f\n", resL[["objective"]])) 
                    } else if (trace > 1L) {
                        cat("SOLUTION\n")
                        print(resL)
                    }
                } else {
                    print(resL)
                }
                
                ## the constraint must be active
                check1 <- g(resL$solution, object = object, k = k, level = lev)$constraints

                check2 <- object$negLogLikFun(theta = resL$solution,
                                              deriv = TRUE,
                                              object = object)

                if (trace) {
                    cat(sprintf("   Constraint check %10.7f, %10.4f\n", check1, check2))
                    grad0 <- drop(attr(check2, "gradient"))
                    cat("   Grad : ", sprintf("%8.6f", grad0 / norm(grad0, type = "2")),
                        "\n\n")
                }
                
                check1 <- (check1 > -1e-3)
                
                if (!inherits(resL, "try-error") && (resL$status %in% c(3, 4))) {
                    thetaLPrec <- resL[["solution"]]
                    ci[k, "L", ilev] <- thetaLPrec[k]
                } else {
                    thetaLPrec <- NULL
                }
                
                ## here we maximise will 'nloptr' only minimises things
                if (trace) {
                    cat(sprintf("   %s, upper bound: ", fLevel[ilev]))
                }

                ## =========================================================
                ## if we have successfully computed the result for a
                ## larger confidence level (and the same parameter),
                ## use it as initial guess
                ## ==========================================================
                
                if ((ilevPrec > 1L) && (!is.null(thetaUPrec))) {
                    theta0 <- thetaUPrec
                } else {
                    theta0 <- thetaHat
                }
            
                resU <- try(nloptr::nloptr(x0 = theta0,
                                           eval_f = f,
                                           eval_g_ineq = g,
                                           lb = c(0.0, 0, -0.90),
                                           ub = c(Inf, Inf, Inf),
                                           k = k, level = lev, chgSign = TRUE,
                                           opts = opts1,
                                           object = object))
                
                if (trace == 1L) {
                    cat(sprintf("%8.6f\n", -resU[["objective"]])) 
                } else  if (trace > 1L) {
                    cat("SOLUTION\n")
                    print(resU)
                }
                
                ## the constraint must be active
                check1 <- g(resU$solution, object = object, k = k, level = lev)$constraints

                check2 <- object$negLogLikFun(theta = resU$solution,
                                              deriv = TRUE, object = object)

                if (trace) {
                    cat(sprintf("   Constraint & value %10.7f, %10.4f\n", check1, check2))
                    grad0 <- drop(attr(check2, "gradient"))
                    cat("   Grad : ", sprintf("%7.5f", grad0 / norm(grad0, type = "2")),
                        "\n\n")
                }

                
                check1 <- (check1 > -1e-3)
                
                if (!inherits(resU, "try-error") && (resU$status %in% c(3, 4))) {
                    thetaUPrec <- resU[["solution"]]
                    ci[k, "U", ilev] <- thetaUPrec[k]
                } else {
                    thetaUPrec <- NULL
                }
                
                ilevPrec <- ilev
            }

        }
        
    }

    if (round && (!is.null(object$sd)) && (!any(is.na(object$sd)))) {
        d <- ceiling(-log(min(object$sd), 10)) + 1
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
        
        ciPlus <- data.frame(Name = object$parNames,
                             LU = rep("est", 3),
                             Level = rep("est", 3),
                             Value = object$estimate,
                             NegLogLik = rep(-object$logLik, 3))
        
        ci <- rbind(ci, ciPlus, deparse.level = 1)
            
        L <- list(ci = ci, negLogLikC = negLogLikC)
        class(L) <- "confintCheck.poisGP"
        L
    } else {
        ci
    }
}


