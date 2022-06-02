## *****************************************************************************
##' Profile-likelihood inference method. This method finds the bounds
##' of a confidence interval for the given function of the parameter
##' vector.
##'
##' Under suitable conditions such as the smoothness of the function
##' \eqn{f(\boldsymbol{\theta})}{f(\theta}) given in \code{fun},
##' \eqn{\eta := \boldsymbol{\theta}}{\eta := f(\theta)} can be
##' considered as a parameter of the model in a suitable
##' re-parameterisation of it. So it makes sense to use the
##' profile-likelihood method to derive confidence intervals on
##' it. Although different methods can be used for this the
##' \bold{potomax} package favours using an optimisation of
##' \eqn{f(\boldsymbol{\theta})}{f(\theta}) under a constraint of high
##' log-likelihood \eqn{\ell(\boldsymbol{\theta}) \geq
##' \ell_{\textrm{max}} - \delta}{l(\theta) >= lMax - \delta} where
##' \eqn{\delta} is a small positive value depending on the confidence
##' level.
##' 
##' @title Profile-Likelihood Inference Method
##'
##' @param object An object representing a fitted parametric model.
##'
##' @param fun A numeric function of the vector of parameters of the
##' model given in \code{object}.
##' 
##' @param ... Further arguments for methods.
##' 
##' @return The result, typically a numeric array containing
##' confidence bounds.
##'
##' @author Yves Deville
##' 
profLik <- function(object, fun, ...) {
    UseMethod("profLik")
}

## *****************************************************************************
##' Profile-likelihood inference for fitted model objects.
##'
##' @details Compute the lower and upper end-points of a
##' profile-likelihood based confidence interval. The (apparently new)
##' method used here relies on maximising and minimising the function
##' of interest, say \eqn{\eta(\boldsymbol(\theta)}{\eta(\theta)},
##' under the constraint that the log-likelihood is greater than the
##' maximal log-likelihood minus a positive quantity \eqn{\delta}
##' depending on the confidence level. This differs from the usual
##' method which relies on an univariate zero-finding for the
##' profile-likelihood function (minus a constant). Remind that each
##' evaluation of the profile requires a \eqn{p-1} dimensional
##' optimisation. As a major advantage, the new method does not
##' require a re-parameterisation of the model.
##'
##' The requirement for the method to work is that \code{object} has
##' suitable "slots", i.e. is a list with suitable elements. These
##' are: the number \code{p} of parameters, the vector \code{parNames}
##' of the parameter names, the vector \code{estimate} of ML
##' estimates, the closure \code{negLogLikFun} computing the negative
##' log-likelihood function and the value \code{negLogLik} of the
##' minimised negative log-likelihood. These requirements are
##' fulfilled for \code{\link{poisGP}} objects of the package. The
##' result returned by \code{object$negLogLikFun} should have a
##' \code{"gradient"} attribute to allow the use of derivatives.
##' 
##' @title Profile-Likelihood Inference for Fitted Parametric Model
##' Objects
##'
##' @usage 
##' \method{profLik}{default}(object,
##'         fun,
##'         level = 0.95,
##'         deriv = TRUE,
##'         trace = 0,
##'         diagno = TRUE,
##'         ftol_abs = 1e-12, ftol_rel = 1e-8,
##'         ...)
##' 
##' @param object A fitted model object. See \bold{Details}.
##'
##' @param fun A function of the parameter vector for which the
##' profile-likelihood will be carried over. This function must have
##' the arguments: \code{theta} for the vector of parameters and
##' \code{object} for the model object; so the function can use the
##' velue of some of slots of \code{object} see \bold{Details}. The
##' function must return a list with two elements with names
##' \code{"objective"} and \code{"gradient"} as required by
##' \code{\link[nloptr]{nloptr}} see \bold{Examples}. If needed, a
##' wrapper function can be used to use more arguments.
##' 
##' @param level Level of confidence. Can be of length \code{> 1}.
##'
##' @param deriv Logical. If \code{TRUE}, the function \code{fun} is
##' assumed to provide a gradient vector as an attribute named
##' \code{"gradient"} of the result. For now \code{deriv} can only be
##' \code{TRUE}, which implies that \code{fun} \emph{must} compute the
##' gradient.
##'
##' @param trace Level of verbosity; the value \code{0} prints
##' nothing.
##'
##' @param diagno Logical. When \code{TRUE} two diagnostics are
##' returned as two attributes \code{"diagno"} and \code{"theta"} of
##' the returned array. The array in attribute \code{"diagno"}
##' contains the return codes of \code{nloptr} (named
##' \code{"status"}), the values of the objective and of the
##' constraint functions as well as a value named \code{"gradDist"}
##' whcih measures the distance between the two directions of the two
##' gradients: objective and constraint.  The array in attribute
##' \code{"theta"} contains the value of the Poisson-GP parameter that
##' was found to maxi/minimise the function \code{fun} under the
##' constraint of a high log-likelihood. See section \bold{Note}.
##' 
##' @param ftol_abs,ftol_rel Absolute and relative tolerance to stop
##' the constrained optimisation \code{\link[nloptr]{nloptr}}. These
##' apply to the objective of the constrained optimisation that is to
##' the value of \code{fun}. Remind that \code{ftol_abs} is thus given
##' with the same unit as \code{fun}.
##' 
##' @param ... Not used yet. 
##'
##' @return An array with the value of the function and the
##' corresponding Lower and Upper end-points for the given confidence
##' levels. This array has two attributes with names \code{"diagno"}
##' and \code{"psi"} which both are arrays. The attributes provide
##' information about the numerical optimisation and the values of the
##' vector of parameter that maximised or minimised the function
##' \code{fun}.
##'
##' @author Yves Deville
##'
##' @references
##'
##' Deville Y. (2017) "Profile-Likelihood Using Constrained
##' Optimisation". Unpublished Tech. Report.
##' 
##' Johnson S.G. \emph{The NLopt Nonlinear-Optimization Package}.
##' \url{https://github.com/stevengj/nlopt}.
##'
##' Section \bold{Return values} in the manual
##' \href{https://nlopt.readthedocs.io/en/latest/NLopt_Reference}{NLOPT
##' Reference}.
##' 
##' @note For each confidence limit the numerical optimisation may
##' fail, in which case the limit will be \code{NA}. Using \code{trace
##' = 1} can be useful to further check the optimisation. We encourage
##' the user to inspect the set of important diagnostics returned when
##' \code{diagno = TRUE}, especially those in the attribute of the
##' result named \code{"diagno"}. The \emph{optimisation status} value
##' in \code{status} should be between \code{1} and \code{4} to
##' indicate that small changes on the parameter or on the objective
##' were eventually obtained. On the other hand, a value of \code{5}
##' for \code{status} indicates that the maximal number of iterations
##' was reached, which is considered here as a failure. Both the
##' \code{constraint} and \code{gradDist} values should be small
##' because the constraint must be active at the optimum and the two
##' gradients for the objective and the constraint must be colinear at
##' the optimum (Lagrange conditions). In other words the optimum
##' parameter vector must lie on the boundary of the high-likelihood
##' region corresponding to the chosen confidence level.
##'
##' @seealso \code{\link[nloptr]{nloptr}} for details on the
##' optimisation.
##' 
##' @examples
##' object <- poisGP(Garonne)
##'
##' ## =========================================================================
##' ## Define a function of the parameter vector: here the first component of
##' ## the "PP" parameter vector.
##' ## This is for code illustration only, since the the result can be obtained
##' ## using the 'confint' method with 'method = "proflik"', which
##' ## gives the profile-likelihood confidence intervals for each of the
##' ## three parameters.
##' ## =========================================================================
##'
##' numPP <- 2
##' myfun <- function(theta, object) {
##'     thetaStar <- poisGP2PP(lambda = theta[1], scale = theta[2], shape = theta[3],
##'                            loc = object$threshold, deriv = TRUE)
##'     res <- thetaStar[numPP]
##'     grad <- attr(thetaStar, "gradient")[1, numPP, c(1, 3, 4)]
##'     attr(res, "gradient") <- grad
##'     res
##' }
##'
##' pl <- profLik(object = object, fun = myfun, level = c(0.95, 0.70))
##' 
profLik.default <- function(object,
                            fun,
                            level = 0.95,
                            deriv = TRUE,
                            trace = 0,
                            diagno = TRUE,
                            ftol_abs = 1e-12,
                            ftol_rel = 1e-8,
                            ...) {

    .diagno <- diagno

    ## =========================================================================
    ## XXX TO IMPLEMENTED LATER. To allow the use of dots, we will
    ## need some 'do.call' and also the replacement method `formal<-`
    ## in order to add the formal arguments to both 'f' and 'g'. So
    ## some care will be required.
    ## =========================================================================
    
    if (FALSE) {
        dots <- match.call(expand.dots = FALSE)[["..."]]
        nmOk <- names(dots) %in% names(formals(fun))
        
        if (!all(nmOk)) {
            stop("all formals passed through the dots '...' must be ",
                 "formals of the function given in 'fun'")
        }
    }

    if (!("p" %in% names(object)) || !("estimate" %in% names(object)) ||
        !("parNames" %in% names(object)) || !("negLogLik" %in% names(object)) ||
        !("negLogLikFun" %in% names(object))) {
        stop("Invalid 'object'. Should be a list with elements named ",
             "'p', 'estimate', 'parNames' and 'negLogLikFun'")
    }
    
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    thetaHat <- object$estimate
    constrCheck <- -5e-3
    
    res <- array(NA, dim = c(Lim = 3L, Level = nLevel),
                 dimnames = list(Type = c("est", "L", "U"), Level = fLevel))

    ## =========================================================================
    ## We maximise/minimise 'fun' under the constraint that the logLik
    ## remains >= max logLik - delta where delta := qchisq(1 - alpha) / 2
    ## where alpha is given by the confidence level.
    ##
    ## The constrained optim is performed using an augmented
    ## Lagrangian method which requires a local companion algorithm
    ## with its own settings. So a sublist is used to tune the local
    ## optimisation.
    ## =========================================================================

    opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                  "xtol_rel" = 1.0e-12,
                  "ftol_abs" = ftol_abs,
                  "ftol_rel" = ftol_rel,
                  "maxeval" = 8000,
                  "check_derivatives" = FALSE,
                  "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-12,
                      "maxeval" = 8000,
                      "ftol_abs" = ftol_abs,
                      "ftol_rel" = ftol_rel),
                  "print_level" = 0)
    
    if (trace >= 2) {
        opts1[["check_derivatives"]] <- TRUE
        opts1[["check_derivatives_print"]] <- "all"
    }
    
    if (deriv) {
        
        f <- function(theta, object, level, chgSign) {
            
            res <- fun(theta, object)
            
            ## 'nloptr' fails on NA and NaN!
            if (is.na(res)) {
                if (chgSign) {
                    return(list("objective" = Inf,
                                "gradient" = rep(NaN, object$p)))
                 } else {
                     return(list("objective" = Inf,
                                 "gradient" = rep(NaN, object$p)))
                 }
            }

            gradtheta <-  attr(res, "gradient")
            
            if (chgSign) {
                 return(list("objective" = -res, "gradient" = -gradtheta))
             } else {
                 return(list("objective" = res, "gradient" = gradtheta))
             }
        }

        g <- function(theta, object, level, chgSign) {

            ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(theta = theta, object = object, deriv = TRUE)
            
            res2 <- list("constraints" = res - ellL,
                         "jacobian" = attr(res, "gradient"))
            res2 
        }
        
    } else {
        stop("for now, the method is only implemented for deriv = TRUE")
    }

    ## the confidence level is not used here
    val <- f(thetaHat, object, level = 0.95, chgSign = FALSE)$objective
    res["est",  ] <- val
    
    ## =========================================================================
    ## keep some information about optimisation: diagnostics and value
    ## of the parameter vector which lead to the min or max of the
    ## profiled function.
    ## =========================================================================
    
    diagno <-
        array(NA,
              dim = c(Lim = 2L,
                  Level = nLevel,
                  Diag = 4L),
              dimnames = list(Type = c("L", "U"),
                  Level = fLevel,
                  Diag = c("status", "objective", "constraint", "gradDist")))
    
    Theta <-
        array(NA,
              dim = c(Lim = 2L,
                  Level = nLevel,
                  theta = object$p),
              dimnames = list(Type = c("L", "U"),
                  Level = fLevel,
                  theta = object$parNames))
    
    labs <- c("L" = "Lower", "U" = "Upper")
    sign <- c("L" = 1.0, "U" = -1.0)
    chgSign <- c("L" = 0.0, "U" = 1.0)
    
    for (LU in c("U", "L")) {
        
        for (iLev in seq_along(level)) {
            
            lev <- level[iLev]
            
            if (trace) {
                cat(sprintf("\n\no %s bound for level %s\n",
                            labs[LU], fLevel[iLev]))
            }
            
            ## =================================================================
            ## if we have successfully computed the result for a
            ## larger confidence level (and the same function of the
            ## parameter vector), use the corresponding value of the
            ## parameter as initial guess.
            ## =================================================================
            
            if ((iLev > 1L) && !is.null(thetaPrec)) {
                theta0 <- thetaPrec
                if (trace > 1) {
                    cat("\nInitialising with the previous conf. level\n")
                }
            } else {
                theta0 <- thetaHat
                if (trace > 1) {
                    cat("\nInitialising with the MLE\n")
                }
            }

            ## print(f(theta0, object, level = lev, chgSign[LU]))
            
            resOpt <- try(nloptr::nloptr(x0 = theta0,
                                         eval_f = f,
                                         eval_g_ineq = g,
                                         lb = object$lb,
                                         ub = object$ub,
                                         level = lev,
                                         chgSign = chgSign[LU],
                                         opts = opts1,
                                         object = object))
            
            if (!inherits(resOpt, "try-error")) {
                diagno[LU, iLev, "status"] <- resOpt$status
                if (trace == 1L) {
                    cat(sprintf("    Optimisation status: %d\n", resOpt$status))
                    cat(sprintf("    Iterations:          %d\n", resOpt$iterations))
                }
                
                if (trace > 1L) {
                    cat("\nSOLUTION\n")
                    print(resOpt)
                }
            } else {
                if (trace) {
                   cat("    Optimisation error!\n")
                }
            }
            
            ## =================================================================
            ## compute value of the constaint as well as the distance
            ## between the two directions gradient of objective 'f'
            ## and gradient of constraint 'g' to see if the constaint
            ## is active at solution
            ## =================================================================

            if (!inherits(resOpt, "try-error") && (resOpt$status %in% 1:4)) {
                
                checkg <- g(theta = resOpt$solution,
                            level = lev,
                            chgSign = chgSign[LU],
                            object = object)
                
                checkf <- object$negLogLikFun(theta = resOpt$solution,
                                              deriv = TRUE,
                                              object = object)
                
                diagno[LU, iLev, "objective"] <- checkf[[1]]
                diagno[LU, iLev, "constraint"] <- checkg$constraints
                
                if (trace == 1L) {
                    cat(sprintf("    Objective value:  %10.7f\n", checkf[[1]]))
                    cat(sprintf("    Constraint check: %10.7f\n", checkg$constraints))
                }
                
                gradDist <- distLines(x1 = checkg$jacobian,
                                      x2 = attr(checkf, "gradient"))
                
                diagno[LU, iLev, "gradDist"] <- gradDist
                
                if (trace == 1L) {
                    cat(sprintf("    gradDist:        %10.7f\n", gradDist))
                    ## print(rbind("    f " = checkf$gradient,
                    ##             "    g " = checkg$jacobian))
                    
                }

                if ( (!is.na(gradDist)) && (gradDist < 0.05)) {
                
                    optDone <- TRUE
                    thetaPrec <- resOpt[["solution"]]
                    res[LU, iLev] <- sign[LU] * resOpt[["objective"]]
                    Theta[LU, iLev, ] <- resOpt[["solution"]]
                    
                } else {
                    thetaPrec <- NULL
                }
      
                
            } else {
                thetaPrec <- NULL
            }
      
        }
        
    }

    
    ## attach diagnostic and parameter values as attributes.
    if (.diagno) {
        attr(res, "diagno") <- diagno
        attr(res, "theta") <- Theta
    }
    
    invisible(res)

}
