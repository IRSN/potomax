
## ****************************************************************************
##' Compute the ML estimate of the rate of the Poisson process given
##' the GPD parameters.
##'
##' @title Compute the ML Estimate of the Rate
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
##' Negative log-likelihood function for a
##' \code{poisGP} object.
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
##' @return The value of the negative log-likelihood. When
##' \code{deriv} is \code{TRUE} the result has an attribute named
##' \code{"gradient"} the value of which is a row matrix with three
##' columns corresponding to the parameters \code{lambda},
##' \code{scale} and \code{shape}.
##'
##' @section Caution: The negative log-likelihood is computed \emph{up
##' to a constant} which is unimportant in the optimisation and is
##' taken to be zero in the \emph{concentrated} version as computed by
##' \code{\link{negLogLikFunC}}, which actually is the worhorse of the
##' estimation. This maintains a compatibility with the log-likelihood
##' as computed by \code{Renext::Renouv} but not with that arising
##' from other packages. So please be careful when comparing
##' log-likelihoods across packages. This is also true for AIC and BIC
##' criteria.
##'
negLogLikFun <- function(theta, object, deriv = TRUE) {
    
    lambda <- theta[1]
    scale <- theta[2]
    shape <- theta[3]
    
    negLogL <- 0
    if (deriv) grad <- array(0.0, dim = c(1L, 3L))
    
    fd <-   object$fitData

    Cst <- 0
    
    if (fd[["OT"]]$flag) {

        lw <- lambda * fd[["OT"]]$effDuration
        negLogL  <- negLogL - fd[["OT"]]$n * log(lw) + lw
        Cst <- Cst - fd[["OT"]]$n * (1 - log(fd[["OT"]]$effDuration))
        
        if (deriv) {
            grad[1L] <- grad[1L] + fd[["OT"]]$effDuration - fd[["OT"]]$n / lambda
        }
        ## second term GP density. Note that it can be the case that
        ## no observation exist over the threshold.
        if (fd[["OT"]]$n > 0) {
            ldens <- dGPD2(x = fd[["OT"]]$data, scale = scale, shape = shape,
                           log = TRUE, deriv = deriv)
            negLogL  <- negLogL - sum(ldens)
            if (deriv) {
                grad[2:3] <- grad[2:3] - apply(attr(ldens, "gradient"), 2, sum)
            }
        }
    }
    
    for (nm in c("MAX", "OTS")) {
        
        if (fd[[nm]]$flag) {
            
            ## first term
            lw <- lambda * fd[[nm]]$effDuration
            negLogL <- negLogL - sum(fd[[nm]]$r * log(lw))
            if (deriv) grad[1] <- grad[1] - sum(fd[[nm]]$r) / lambda
            Cst <- Cst - fd[[nm]]$r * (1 - log(fd[[nm]]$effDuration))
            
            ## second term: GP survival
            S <- pGPD2(q = fd[[nm]]$threshold, scale = scale, shape = shape,
                       deriv = TRUE, lower.tail = FALSE)
            negLogL <- negLogL + sum(lw * S)

            if (deriv) {
                grad[1] <- grad[1] + sum(fd[[nm]]$effDuration * S)
                grad[2:3] <- grad[2:3] + crossprod(lw, attr(S, "gradient"))
            }
            
            ## third term GP density
            if (sum(fd[[nm]]$r) > 0) { 
                ldens <- dGPD2(x = unlist(fd[[nm]]$data),
                               scale = scale, shape = shape,
                               deriv = TRUE, log = TRUE)
                negLogL <- negLogL - sum(ldens)
                if (deriv) {
                    grad[2:3] <- grad[2:3] -
                        apply(attr(ldens, "gradient"), 2, sum) 
                }
            }
        }
        
    }
    
    if (deriv) attr(negLogL, "gradient") <- grad
    negLogL + Cst
     
}

## ****************************************************************************
##' Maximum-Likelihood Estimation of a Poisson-GP model using heterogeneous
##' data.
##'
##' The estimation proceeds by minimising a concentrated (or profile)
##' negative log-likelihood which depends on the two GPD parameters
##' but not on the Poisson rate. However the negative log-likelihood
##' function is built and returned because it will be used to derive
##' profile-likelihood inference results.
##' 
##' @title Maximum-Likelihood Estimation of a Poisson-GP Model
##'
##' @param object A \code{poisGP} object that needs to be estimated.
##'
##' @param parIni Initial values for the parameter vector.
##'
##' @param estim Type or method chosen for the estimation.
##'
##' @param coefLower,coefUpper Lower and Upper bounds for the
##' parameters.
##'
##' @param parTrack Not used yet.
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
MLE.poisGP <- function(object= NULL, 
                       parIni = NULL,
                       estim = c("optim", "nloptr", "none"),
                       coefLower = c("scale" = 0.0, "shape" = -0.99),
                       coefUpper = c("scale" = Inf, "shape" = 2.0),
                       parTrack =  FALSE,
                       trace = 0) {

    estim <-  match.arg(estim)
       
    cvg <- TRUE
    res <- list()

    if (estim == "optim") {
        
        res$fit <- try(optim(par = parIni,
                             fn = negLogLikFunC,
                             deriv = FALSE,
                             method = "BFGS",
                             control = list(maxit = 3000, trace = trace,
                                 parscale = c(1000, 1)),
                             ## hessian = TRUE,
                             object = object))

        if (!inherits(res$fit, "try-error")) {
            if (res$fit$convergence == 0) {
                estimate <- res$fit$par
                res$estimate <- estimate
                res$negLogLik <- res$fit$value
            } else {
                cvg <- FALSE
            }
        }
            
    } else if (estim == "nloptr") {
        
        opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                     "xtol_rel" = 1.0e-8,
                     "xtol_abs" = 1.0e-8,
                     "ftol_abs" = 1e-5,
                     "maxeval" = 1000, "print_level" = 0,
                     "check_derivatives" = FALSE)

        ## XXX caution! this works when the shape is constant only!!!
        p <- object$p - 1
        
        lb <- rep(c("scale" = -Inf, "shape" = -Inf))
           
        if (length(coefLower)) {
            lm <- match(names(coefLower), names(lb))
            if ((length(lm) != length(coefLower)) ||
                any(is.na(lm))) {
                stop("when given, 'coefLower' must be a named vector ",
                     "with suitable element names")
            }
            lb[lm] <- coefLower
        } 
        
        ub <- rep(c("scale" = Inf, "shape" = Inf))
        
        if (length(coefUpper)) {
            um <- match(names(coefUpper), names(ub))
            if ((length(um) != length(coefUpper)) ||
                any(is.na(lm))) {
                stop("when given, 'coefUpper' must be a named vector ",
                     "with suitable element names")
            }
            ub[lm] <- coefUpper
        }

        
        if (trace) {
            cat("Initial values 'parIni', bounds 'lb' and 'ub'\n")
            print(cbind(parIni = parIni, lb = lb, ub = ub))
        }
        
        res$fit <- try(nloptr(x0 = parIni,
                              eval_f = negLogLikFunCD,
                              lb = lb,
                              ub = ub,
                              opts = opts,
                              object = object))
        
        if (!inherits(res$fit, "try-error")) {
            if (res$fit$status > 0) {
                estimate <- res$fit$solution
                
                if (any(estimate <= lb) || any(estimate >= ub)) {
                    warning("some estimated parameters at ",
                            "the bounds, inference results are misleading")
                }

                names(estimate) <- c("scale", "shape") ## object$parNames
                res$estimate <- estimate
                res$negLogLik <- res$fit$objective
            } else {
                cvg <- FALSE
            }
        }
        
    }
    
    if (!cvg) {
        warning("convergence not reached in optimisation")
        estimate <- rep(NA, object$p)
        names(estimate) <- object$parNames
        res$negLogLik <- NA
        res$estimate <- estimate
        res$logLik <- NA
    } else {

        res$logLik <- -res$negLogLik

        if (trace) {
            cat("Optimisation results\n")
            print(res)
        }

        
        ## ## compute Hessian. 
        ## res$hessian <- numDeriv::hessian(func = negLogLikFun,
        ##                                  x = res$estimate, deriv = FALSE,
        ##                                  object = object)

        .lambdaHat <- lambdaHat(thetaGP = res$estimate, object = object, log = FALSE)
        res$estimate <- c("lambda" = .lambdaHat, res$estimate) 
        res$hessian <- optimHess(par = res$estimate,
                                 fn = negLogLikFun,
                                 deriv = FALSE,
                                 object = object)
        
        vcov <- try(solve(res$hessian), silent = TRUE)
        
        if (!inherits(vcov, "try-error")) {
            rownames(vcov) <- colnames(vcov) <- object$parNames
            res$vcov <- vcov
            res$sd <- sqrt(diag(vcov))
        }
    }
    
    if (parTrack) {
        tpsi <-  matrix(trackEnv$psi, ncol = object$p,
                                byrow = TRUE)
        colnames(tpsi) <- object$parNames
        res$tracked <-
            list(psi = tpsi,
                 negLogLik = apply(tpsi, 1, negLogLikFun, deriv = FALSE,
                     object = object))
    }
    
    res
    
    
    

}
