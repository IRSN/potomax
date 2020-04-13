
## ****************************************************************************
##' Compute the ML estimate of the rate of the Poisson process given
##' the GPD parameters.
##'
##' @noRd
##' 
##' @title Compute the ML Estimate of the Rate Given the 
##'
##' @param thetaGP Parameter vector (GP part only).
##'
##' @param object The \code{poisGP} object corresponding to the
##' estimation.
##'
##' @param log Logical. If \code{TRUE} the function returns
##' \eqn{\log \hat{lambda}}{log(lambdaHat)} possibly with its gradient, else it
##' returns \eqn{\hat{\lambda}}{lambdaHat}.
##'
##' @param deriv Logical. Should the gradient be computed and be
##' returned as the \code{"gradient"} attribute of the result?
##'
##' @return The value \eqn{\hat{\lambda}}{lambdaHat} of the estimated
##' rate with an attribute named \code{"r"} corresponding to the
##' number of observations. When \code{deriv} is \code{TRUE} the
##' result has an attribute named \code{"gradient"} the value of which
##' is a row matrix with two columns.
##' 
lambdaHat <- function(thetaGP, object, log = TRUE, deriv = TRUE) {
    
    fd <-   object$fitData
    
    scale <- thetaGP[1]
    shape <- thetaGP[2]
    
    if (deriv) grad <- array(0.0, dim = c(1L, 2L))
    
    r <- 0
    w <- 0.0
    
    ## compute number of obs and discounted durations by blocks
    if (fd[["OT"]]$flag) {
        w <- fd[["OT"]]$effDuration
        r <- r + fd[["OT"]]$n
    }
    
    for (nm in c("MAX", "OTS")) {
        if (fd[[nm]]$flag) {
            S <- pGPD2(q = fd[[nm]]$threshold, scale = scale, shape = shape,
                       lower.tail = FALSE, deriv = deriv)
            w <- w + sum(fd[[nm]]$effDuration * S)
            r <- r + sum(fd[[nm]]$r)
            if (deriv) {
                grad <- grad + crossprod(fd[[nm]]$effDuration,
                                         attr(S, "gradient")) 
            }
        }
    }
    
    lambdaHat <- r / w
    if (log) lambdaHat <- log(lambdaHat)
    
    if (deriv) {
        if (log) {
            attr(lambdaHat, "gradient") <- - grad / w
        } else {
            attr(lambdaHat, "gradient") <- - (lambdaHat / w) * grad
        } 
    }
    
    attr(lambdaHat, "r") <- r
    
    lambdaHat
    
}

## *****************************************************************************
##' Concentrated (or profile) negative log-likelihood function for a
##' \code{poisGP} object.
##'
##' @noRd
##' 
##' @title Concentrated or Profile Negative Log-Likelihood Function
##' for a \code{poisGP} Object
##'
##' @param thetaGP Parameter vector (GP part only).
##'
##' @param object The \code{poisGP} object corresponding to the
##' estimation.
##'
##' @param deriv Logical. Should the gradient be computed and be
##' returned as the \code{"gradient"} attribute of the result?
##'
##' @return The value of the negative log-likelihood. When
##' \code{deriv} is \code{TRUE} the result has an attribute named
##' \code{"gradient"} the value of which is a row matrix with two
##' columns.
##' 
##'
negLogLikFunC <- function(thetaGP, object, deriv = TRUE) {

    scale <- thetaGP[1L]
    shape <- thetaGP[2L]

    if (scale < 0.0) return(NA)
    
    if (deriv) grad <- array(0.0, dim = c(1L, 2L))
    
    fd <-   object$fitData
    
    logLambdaHat <- lambdaHat(thetaGP = thetaGP, object = object,
                              log = TRUE, deriv = deriv)

    ## cat("lambdaHat = ", exp(logLambdaHat), "\n")
    
    ## Ordinary OT part caution all blocks OT are grouped as one
    ## here
    negLogL  <- -attr(logLambdaHat, "r") * logLambdaHat
    
    if (deriv) {
        grad <- grad - attr(logLambdaHat, "r") *
            attr(logLambdaHat, "gradient")
    }
    
    if (fd[["OT"]]$flag) {
        ldens <- dGPD2(x = fd[["OT"]]$data,
                       scale = scale, shape = shape,
                       log = TRUE, deriv = deriv)
        negLogL <- negLogL - sum(ldens)
        if (deriv) {
            grad <- grad - apply(attr(ldens, "gradient"), 2, sum)
        }
    }
    
    for (nm in c("MAX", "OTS")){
        if (fd[[nm]]$flag) {
            if (sum(fd[[nm]]$r) > 0) {
                ldens <- dGPD2(x = unlist(fd[[nm]]$data),
                               scale = scale, shape = shape,
                               log = TRUE, deriv = deriv)
                negLogL <- negLogL - sum(ldens)
                if (deriv) {
                    grad <- grad - apply(attr(ldens, "gradient"), 2, sum)
                }
            }
        }
    }

    if (deriv) attr(negLogL, "gradient") <- grad
    
    negLogL
    
}

negLogLikFunCD <- function(thetaGP, object) {

    nL <- negLogLikFunC(thetaGP = thetaGP, object = object, deriv = TRUE)
    
    list("objective" = nL, "gradient" = attr(nL, "gradient"))
    
}

## *****************************************************************************
##' Negative log-likelihood function for a \code{poisGP} object.
##'
##' @noRd
##' 
##' @title Negative Log-Likelihood Function for a \code{poisGP} Object
##'
##' @param theta Parameter vector Poisson rate and GP parameters.
##'
##' @param object The \code{poisGP} object corresponding to the
##' estimation.
##'
##' @param deriv Logical. Should the gradient be computed and be
##' returned as the \code{"gradient"} attribute of the result?
##'
##' @param hessian Logical. If \code{TRUE} the Hessian is computed and
##' returned as the \code{"hessian"} attribute of the result. Note that
##' \code{hessian} can be \code{TRUE} only when \code{deriv} is
##' \code{TRUE}.
##' 
##' @return The value of the negative log-likelihood. When
##' \code{deriv} is \code{TRUE} the result has an attribute named
##' \code{"gradient"} the value of which is a row matrix with three
##' columns corresponding to the parameters \code{lambda},
##' \code{scale} and \code{shape}.
##'
##' @section Caution: The negative log-likelihood is computed \emph{up
##' to a constant} which is unimportant in the optimisation and is
##' taken to be zero in the \emph{concentrated} version as computed by
##' \code{\link{negLogLikFunC}}, which actually is the workhorse of the
##' estimation. This maintains a compatibility with the log-likelihood
##' as computed by \code{Renext::Renouv} but not with that arising
##' from other packages. So please be careful when comparing
##' log-likelihoods across packages. This is also true for AIC and BIC
##' criteria.
##'
negLogLikFun <- function(theta, object, deriv = TRUE, hessian = FALSE) {
    
    lambda <- theta[1]
    scale <- theta[2]
    shape <- theta[3]

    if (is.na(lambda)) stop("'lambda' is NA")
    if (lambda < 0.0) stop("negative value of 'lambda'")
    if (scale < 0.0) stop("negative value of 'scale'")
    
    negLogL <- 0
    if (deriv) {
        grad <- array(0.0, dim = c(1L, 3L))
        if (hessian) {
            hess <- array(0.0, dim = c(1L, 3L, 3L))
        }
    }
    
    fd <-   object$fitData
    Cst <- 0
    
    if (fd[["OT"]]$flag) {

        n <- fd[["OT"]]$n
        w <- fd[["OT"]]$effDuration
        lw <- lambda * w
        negLogL  <- negLogL - n * log(lw) + lw
        Cst <- Cst - fd[["OT"]]$n * (1 - log(w))
        
        if (deriv) {
            grad[1L, 1L] <- grad[1L, 1L] - n / lambda + w
            if (hessian) {
                hess[1L, 1L, 1L] <- hess[1L, 1L, 1L] + 
                    n / lambda / lambda
            }
        }
        
        ## ====================================================================
        ## 2nd term: GP density. Note that it can be the case that
        ## OT$flag is TRUE but no observation exist over the
        ## threshold!
        ## ====================================================================
        
        if (fd[["OT"]]$n > 0) {

            ldens <- dGPD2(x = fd[["OT"]]$data, scale = scale, shape = shape,
                           log = TRUE, deriv = deriv, hessian = hessian)
            negLogL  <- negLogL - sum(ldens)
            
            if (deriv) {
                grad[1L, 2L:3L] <- grad[1L, 2L:3L] -
                    apply(attr(ldens, "gradient"), 2L, sum)
                if (hessian) {
                    hess[1L, 2L:3L, 2L:3L] <-  hess[1L, 2L:3L, 2L:3L] -
                        apply(attr(ldens, "hessian"), 2L:3L, sum)
                }
            }
        }
    }

    ## ========================================================================
    ## Note that this could be coded more efficiently because log(w)
    ## and log(lw) are both computed 
    ## ========================================================================
    
    for (nm in c("MAX", "OTS")) {
        
        if (fd[[nm]]$flag) {
            
            ## 1-st term
            r <- fd[[nm]]$r
            rSum <- sum(r)
            w <- fd[[nm]]$effDuration
            lw <- lambda * w
            negLogL <- negLogL - sum(r * log(lw))
            if (deriv) {
                grad[1L] <- grad[1L] - rSum / lambda 
                if (hessian) {
                    hess[1L, 1L, 1L] <- hess[1L, 1L, 1L] + 
                        rSum / lambda / lambda
                }
            }
            Cst <- Cst - sum(fd[[nm]]$r * (1 - log(w)))
            
            ## 2-nd term: GP survival
            S <- pGPD2(q = fd[[nm]]$threshold, scale = scale, shape = shape,
                       deriv = TRUE, hessian = hessian, lower.tail = FALSE)
            negLogL <- negLogL + sum(lw * S)

            if (deriv) {
                grad[1L] <- grad[1L] + sum(w * S) 
                grad[2L:3L] <- grad[2L:3L] + crossprod(lw, attr(S, "gradient"))
                if (hessian) {
                    ##  2nd order derivative w.r.t. 'lambda' and 'thetaGP'
                    hess[1L, 1L, 2L:3L] <- hess[1L, 1L, 2L:3L] +
                        crossprod(w, attr(S, "gradient"))
                    hess[1L, 2L:3L, 1L] <- hess[1L, 1L, 2L:3L]
                    ## 2nd order derivatives w.r.t. 'thetaGP' and 'thetaGP'
                    hess[1L, 2L:3L, 2L:3L] <- hess[1L, 2L:3L, 2L:3L] + lambda *
                        apply(attr(S, "hessian"), MARGIN = 2L:3L,
                              function(x) crossprod(w, x))
                }
            }
            
            ## 3-rd term: GP density
            if (rSum > 0) { 
                ldens <- dGPD2(x = unlist(fd[[nm]]$data),
                               scale = scale, shape = shape,
                               deriv = TRUE, hessian = hessian, log = TRUE)
                negLogL <- negLogL - sum(ldens)
                if (deriv) {
                    grad[1L, 2L:3L] <- grad[1L, 2L:3L] -
                        apply(attr(ldens, "gradient"), MARGIN = 2L, sum)
                    if (hessian) {
                        hess[1L, 2L:3L, 2L:3L] <- hess[1L, 2L:3L, 2L:3L] -
                            apply(attr(ldens, "hessian"), MARGIN = 2L:3L, sum)
                    }
                }
            }
        }
        
    }
    
    if (deriv) {
        attr(negLogL, "gradient") <- grad
        if (hessian) {
            attr(negLogL, "hessian") <- hess
        }
    }

    attr(negLogL, "Cst") <- Cst
     
    negLogL ## + Cst ???
     
}

## ****************************************************************************
##' Maximum-Likelihood Estimation of a Poisson-GP model using heterogeneous
##' data.
##'
##' The estimation proceeds by minimising a concentrated (or profile)
##' negative log-likelihood which depends on the two GPD parameters,
##' but not on the Poisson rate. So provided bound on this parameter
##' will have no effect on the estimation and the estimates can fail
##' to have its element \code{"lambda"} within the bounds when these
##' are not \code{0.0} and \code{Inf}. However the standard
##' (non-profile) negative log-likelihood function is built and
##' returned because it will be used to derive profile-likelihood
##' inference results.
##' 
##' @title Maximum-Likelihood Estimation of a Poisson-GP Model
##'
##' @method MLE poisGP
##' 
##' @usage
##'
##' \method{MLE}{poisGP}(object = NULL,
##'     parIni = NULL,
##'     estim = c("optim", "nloptr", "eval", "none"),
##'     coefLower, coefUpper,
##'     parTrack =  FALSE,
##'     scale = FALSE,
##'     trace = 0)  
##'
##' @param object A \code{poisGP} object that needs to be estimated.
##'
##' @param parIni Initial values for the parameter vector. This
##' is must be a named vector of length \eqn{2} with elements
##' names \code{"scale"} and \code{"shape"}.
##'
##' @param estim Type or method chosen for the estimation.
##'
##' @param coefLower,coefUpper Lower and Upper bounds for the
##' parameters. The should be numeric vectors with names in
##' \code{c("lambda", "scale", "shape")}. Only the bounds on the GP
##' parameters \code{"scale"} and \code{"shape"} can be used during
##' the estimation and they will only be used when \code{estim} is
##' \code{"nloptr"}. However the bounds are used in the inference
##' \code{\link{confint.poisGP}} and \code{\link{RL.poisGP}}.
##'
##' @param parTrack Not used yet.
##'
##' @param scale Logical. If \code{TRUE} the data used in the
##' optimisation are scaled , see \code{\link{threshData}}. NOT
##' IMPLEMENTED YET.
##'
##' @param trace Integer Level of verbosity.
##'
##' @return A list with the results of the likelihood
##' maximisation. The content of the list depends on the method as
##' given by \code{estim}, yet it should contain an element
##' \code{logLik} giving the maximised log-likelihood.
##' 
##' @author Yves Deville
##' 
MLE.poisGP <- function(object = NULL, 
                       parIni = NULL,
                       estim = c("optim", "nloptr", "eval", "none"),
                       coefLower,
                       coefUpper,
                       parTrack =  FALSE,
                       scale = FALSE,
                       trace = 0) {

    estim <-  match.arg(estim)
       
    res <- list(cvg = TRUE)

    ## ========================================================================
    ## Manage the bounds
    ## ========================================================================
    
    lb <- rep(c("lambda" = 0.0, "scale" = 0.0, "shape" = -0.99))
    
    if (!missing(coefLower)) {
        lmatch <- match(names(coefLower), names(lb))
        if ((length(lmatch) != length(coefLower)) || any(is.na(lmatch))) {
            stop("when given, 'coefLower' must be a named vector ",
                 "with suitable element names")
        }
        if (ind <- all(coefLower >= lb[lmatch])) {
            lb[lmatch] <- coefLower
        } else {
            stop("Bad value in 'coefLower', element ",
                 sprintf("\"%s\"", names(coefLower[!ind])),
                 ". Should be >= ", lb[lmatch][!ind])
        }           
    }

    if (scale) lb["scale"] <-  lb["scale"] / object$scaleData
    
    ub <- rep(c("lambda" = Inf, "scale" = Inf, "shape" = Inf))
    
    if (!missing(coefUpper)) {
        umatch <- match(names(coefUpper), names(ub))
        if ((length(umatch) != length(coefUpper)) || any(is.na(umatch))) {
            stop("when given, 'coefUpper' must be a named vector ",
                 "with suitable element names")
        }
        if (ind <- all(coefUpper <= ub[umatch])) {
            ub[umatch] <- coefUpper
        } else {
            stop("Bad value in 'coefUpper', element ",
                 sprintf("\"%s\"", names(coefUpper[!ind])),
                 ". Should be <= ", ub[umatch][!ind]) 
        }
    }

    if ((lb["lambda"] > 0.0) || (ub["lambda"] < Inf)) {
        stop("Bounds on 'lambda' must for now be '0.0' (lower) and 'Inf' (upper)")
    }
    
    
    if (scale) ub["scale"] <-  ub["scale"] / object$scaleData
    
    res$lb <- lb
    res$ub <- ub
    
    if (estim == "optim") {
        
        ## ====================================================================
        ## When `estim` is "optim" we use the standard BFGS algorithm,
        ## with no gradient, so `deriv` is FALSE.
        ## ====================================================================

        res$df <- 3
        
        ctrl <- list(maxit = 3000, trace = trace)
        if (!scale)  ctrl[["parscale"]] <- c("scale" = object$scaleData, "shape" = 1)
        
        res$fit <- try(optim(par = parIni[-1],
                             fn = negLogLikFunC,
                             deriv = FALSE,
                             method = "BFGS", control = ctrl,
                             object = object))

        if (!inherits(res$fit, "try-error")) {
            if (res$fit$convergence == 0) {
                estimate <- res$fit$par
                res$estimate <- estimate
                res$negLogLik <- res$fit$value
            } else {
                res$cvg <- FALSE
            }
        }

        eval <- res$cvg
        
    } else if (estim == "nloptr") {
        
        ## ====================================================================
        ## When `estim` is "nloptr" we use the BFGS algorithm with
        ## gradient. The result returned by the function must be a
        ## list. We can use box contraints on the parameters.
        ## ====================================================================

        if (trace) {
            cat("\nUsing the \"NLOPT_LD_BFGS\" algorithm with derivatives.\n")
        }
   
        opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                     "xtol_rel" = 1.0e-8,
                     "xtol_abs" = 1.0e-8,
                     "ftol_abs" = 1e-5,
                     "maxeval" = 1000, "print_level" = 0,
                     "check_derivatives" = FALSE)

        ## XXX caution! this works when the shape is constant only!!!
        p <- object$p - 1
   
        dfred <- sum(lb == ub)
        res$df <- 3 - dfred
        
        if (trace) {
            cat("\nInitial values 'x0' and bounds 'lb' and 'ub'\n")
            print(cbind(parIni = parIni, lb = lb, ub = ub))
            if (dfred) {
                cat("\nNumber of degree of freedom: ", object$df, "\n") 
            }
        }
        
        res$fit <- try(nloptr(x0 = parIni[-1L],
                              eval_f = negLogLikFunCD,
                              lb = lb[-1L],
                              ub = ub[-1L],
                              opts = opts,
                              object = object))
        
        if (!inherits(res$fit, "try-error")) {
            if (res$fit$status > 0) {
                estimate <- res$fit$solution
                
                if (any(estimate <= lb[-1L]) || any(estimate >= ub[-1L])) {
                    warning("some estimated parameters at ",
                            "the bounds, inference results are misleading")
                }

                names(estimate) <- c("scale", "shape") ## object$parNames
                res$estimate <- estimate
                res$negLogLik <- res$fit$objective
            } else {
                res$cvg <- FALSE
            }
        }
        eval <- res$cvg
        
    } else if (estim == "eval") {

        res$df <- 0
        
        ## ====================================================================
        ## When `estim` is "eval": if a valid vector of GP parameters
        ## is provided in 'parIni', we will find the corresponding
        ## Poisson rate, compute the log-likelihood and its
        ## derivatives. This is nearly what would be done by
        ## specifying two identical vectors in 'coefLower' and
        ## 'coefUpper'.
        ## ====================================================================

        if (is.null(parIni)) {
            eval <- FALSE
        } else if ((length(parIni) == 2) && (all.equal(names(parIni),
                             c("scale", "shape")))) {
            eval <- TRUE
            res$estimate <- parIni
        }

        warning("\nNo optimisation performed. Inference results can be misleading.\n")
        res$cvg <- FALSE
        
    }

    if (estim != "eval") {
        if (!res$cvg) {
            warning("\nConvergence not reached in optimisation.\n")
            estimate <- rep(NA, object$p)
            names(estimate) <- object$parNames
            res$negLogLik <- NA
            res$estimate <- estimate
            res$logLik <- NA
            
        } else {
            if (trace) {
                cat("Optimisation results\n")
                print(res)
            }       
        }
    }

    if (eval) {
        
        ## ====================================================================
        ## Find the estimated rate and compute the minimised negative
        ## log-likelihood with its gradient and hessian.
        ## ====================================================================
        
        .lambdaHat <- lambdaHat(thetaGP = res$estimate, object = object,
                                log = FALSE)
        res$estimate <- c("lambda" = .lambdaHat, res$estimate)
        
        if (FALSE) {
            res$hessianCheck <- optimHess(par = res$estimate,
                                          fn = negLogLikFun,
                                          deriv = FALSE,
                                          object = object)
            
            res$gradCheck <- numDeriv::grad(func = negLogLikFun,
                                       x = res$estimate,
                                       deriv = FALSE,
                                       object = object)

        }

        finalEval <- negLogLikFun(theta = res$estimate, object = object,
                                  deriv = TRUE, hessian = TRUE)
        
        res$hessian <- drop(attr(finalEval, "hessian"))
        res$gradNegLogLik <- drop(attr(finalEval, "gradient"))
        res$Cst <- attr(finalEval, "Cst")
        colnames(res$hessian) <- rownames(res$hessian) <-
            names(res$gradNegLogLik) <- object$parNames
        
        attr(finalEval, "gradient") <- attr(finalEval, "hessian") <- NULL
        attr(finalEval, "names") <- attr(finalEval, "Cst") <- NULL
        
        res$logLik <- -finalEval
        res$negLogLik <- finalEval
        
        res$cov <- try(solve(res$hessian))
        if (inherits(res$cov, "try-error") || (rev(eigen(res$cov)$value)[1] < 0.0)) { 
            res$cov <- array(NA, dim = c(3, 3))
        } 
        rownames(res$cov) <- colnames(res$cov) <- object$parNames
        res$sd <- sqrt(diag(res$cov))
        
        ## ====================================================================
        ## Find the 'PP' parameters using the POT2GEV transformation
        ## XXX. We could now use the 'poisGP2PP' function.
        ## ====================================================================

        estPP <- Renext::Ren2gev(res$estimate, threshold = object$threshold)
        covPP <- attr(estPP, "jacobian") %*% res$cov %*% t(attr(estPP, "jacobian"))
        attr(estPP, "jacobian") <- attr(estPP, "threshold") <- NULL
        res$PP <- list(estimate = estPP,
                       cov = covPP,
                       sd = sqrt(diag(covPP)))
        
        
    }
    
    ## if (parTrack) {
    ##     tpsi <-  matrix(trackEnv$psi, ncol = object$p,
    ##                             byrow = TRUE)
    ##     colnames(tpsi) <- object$parNames
    ##     res$tracked <-
    ##         list(psi = tpsi,
    ##              negLogLik = apply(tpsi, 1, negLogLikFun, deriv = FALSE,
    ##                  object = object))
    ## }
        
    res
    
    
    

}
