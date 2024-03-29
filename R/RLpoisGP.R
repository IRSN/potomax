## *****************************************************************************
##' Compute return levels along with confidence bounds for a Poisson-GP
##' model with ML inference results.
##'
##' @details
##' The return-level curve corresponding to the column or to the
##' dimension named \code{"Quant"} is obtained by plugging the ML
##' estimate of the Poisson-GP parameters in the quantile. The
##' confidence limits can be obtained by profile-likelihood or by the
##' standard 'delta' method.
##'
##' @usage
##' \method{RL}{poisGP}(object, period = NULL, level = 0.70,
##'    confintMethod = c("proflik", "delta", "none"),
##'    out = c("data.frame", "array"),
##'    trace = 0,
##'    check = FALSE, nCheck = 50, nSigma = 4,
##'    ftol_abs = 1e-12, ftol_rel = 1e-9,
##'    ...) 
##' 
##' @title Return Levels and Confidence Intervals for a Poisson-GP Model
##'
##' @param object An object with class \code{"poisGP"} representing
##' the inference results for a Poisson-GP model.
##'
##' @param period A vector of periods for which the return levels will
##' be computed. By default "round" periods covering the range from
##' \code{1} to \code{1000} are chosen.
##'
##' @param level Level of confidence. Can be a vector.
##'
##' @param confintMethod The method used to compute the confidence
##' intervals. The value \code{"proflik"} corresponds to the
##' profile-likelihood. The value \code{"delta"} corresponds to the
##' delta method and the value \code{"none"} can be used to obtain the
##' return levels without confidence limits on them.
##' 
##' @param out The type of outpout wanted, see the \strong{Value}
##' section.
##'
##' @param trace Integer level of verbosity. The default value is
##' \code{0}, but when \code{confintMethod} is \code{"proflik"}, it is
##' a good practice to use \code{trace = 1}, so this may change in the
##' future.
##'
##' @param check Logical If. \code{FALSE} the results are intended to
##' be used for a return-level plot, while the value \code{TRUE}
##' produce results for a graphical check of the computations.
##'
##' @param nCheck Number of points on a profile-likelihood curve for
##' each value \eqn{T} of the return period when \code{check} is
##' \code{TRUE}. These points are taken in an interval containing the
##' ML estimate \eqn{\hat{\rho}(T)}{rhoHat(T)} of the return level
##' \eqn{\rho(T)}{rho(T)}, with a range controled by \code{nSigma}.
##'
##' @param nSigma Range the of the return levels used to build the
##' curve when \code{check} is \code{TRUE}. Vector of length two with
##' its values \eqn{n_1}{n1} and \eqn{n_2}{n2} defining the range for the
##' values of \eqn{\rho(T)}{rho(T)} from \eqn{\hat{\rho}(T) - n_1
##' s(T)}{rhoHat(T) - n1 * s(T)} to \eqn{\hat{\rho}(T) + n_2
##' s(T)}{rhoHat(T) + n2 * s(T)} where \eqn{s(T)} is the estimated
##' standard deviation for \eqn{\hat{\rho}(T)}{rhoHat(T)}. The default
##' choice usually covers the profile-likelihood interval for all
##' periods. A vector of length one is recycled.
##'
##' @param ftol_abs,ftol_rel Absolute and relative tolerance to stop
##' the constrained optimisation \code{\link[nloptr]{nloptr}}. These
##' apply to the objective of the constrained optimisation, which is
##' here the return level as a function of the Poisson-GP parameter
##' vector. The smallest possible values reaching convergence should
##' be chosen. By increasing either of thes values the convergence
##' will be easier to get but the results may not be as precise as
##' wanted. This can/should be checked by using \code{check = TRUE}
##' and a subsequent call to \code{autoplot}.
##' 
##' @param ... Not used yet.
##'
##' @return When \code{check} is \code{FALSE} (default), an object of
##' class \code{"RL.poisGP"} inheriting either from
##' \code{data.frame}. It can also be a simple \code{array} containing
##' the results. In the first case, the object can be used with the
##' \code{autoplot} method to produce the \emph{return level} plot
##' without recomputing the results.
##'
##' When \code{check} is \code{TRUE} the result is an object with
##' class \code{RLCheck.poisGP} and can be used with \code{autoplot}
##' to build a \emph{graphical check of the results} of the
##' profile-likelihood method. This is a list containing two data
##' frames: \code{RL} contains the confidence limits and estimates,
##' while \code{negLogLik} contains a grid of values for the the
##' return level \eqn{\rho(T)}{rho(T)} along with the corresponding
##' values of the profile negative log-likelihood. This allows to plot
##' the curve.
##'
##' @note The check for \code{check = TRUE} is built by computing the
##' value of the profile-likelihood for each period \eqn{T} and each
##' candidate value of \eqn{\rho(T)}{rho(T)}. This is done by using a
##' two-parameter optimisation: maximise on the vector \eqn{[\lambda,
##' \, \xi]}{[\lambda, \xi]} of the Poisson rate \eqn{\lambda} and the
##' and GP shape \eqn{\xi}, the return level \eqn{\rho(T)}{rho(T)}
##' being fixed. This optimisation can fail to converge, in which case
##' the result is \code{NA}. For now a derivative-free optimisation is
##' used (COBYLA): the computations can be quite long.
##' 
##' @author Yves Deville
##'
##' @seealso \code{\link{poisGP}} for an example
##' \code{\link{autoplot.RLCheck.poisGP}}. The \code{confint} method
##' has a similar check possibility, see
##' \code{\link{confint.poisGP}}.
##'
##' @importFrom reshape2 melt
##' 
##' @method RL poisGP
##' @export
##' 
##' @examples
##' ## ================================================================
##' ## Use the 'Garonne' data from Renext, which embeds both OT data
##' ## and MAX data
##' ## ================================================================
##' fitp <- poisGP(data = Garonne, threshold = 2900, trace = 2)
##'
##' ## ================================================================
##' ## RL plot with profile-likelihood confidence levels
##' ## ================================================================
##' RL <- RL(fitp, out = "data", level = c(0.70, 0.95))
##' autoplot(RL)
##'
##' \dontrun{
##'    ## ================================================================
##'    ## CHECK the results. Quite slow!
##'    ## ================================================================
##'    RLc <- RL(fitp, out = "data", level = c(0.70, 0.95),
##'              check = TRUE)
##'    autoplot(RLc) + ylim(c(NA, 540)) +
##'        ggtitle("negative profile log-likelihood for rho(T)")
##'
##'    ## ================================================================
##'    ## Using values for ftol_abs or ftol_rel that are not small enough
##'    ## we get problems in the precision of the result.
##'    ## ================================================================
##'
##'    RLc <- RL(fitp, out = "data", level = c(0.70, 0.95), ftol_rel = 1e-5,
##'              check = TRUE)
##'    autoplot(RLc) + ylim(c(NA, 540)) +
##'        ggtitle(paste("negative profile log-likelihood for rho(T) ",
##'                      " with 'ftol_rel' too large"))
##' }
##' 
RL.poisGP <- function(object,
                      period = NULL,
                      level = 0.70,
                      confintMethod = c("proflik", "delta", "none"),
                      out = c("data.frame", "array"),
                      trace = 0,
                      check = FALSE,
                      nCheck = 50,
                      nSigma = 4,
                      ftol_abs = 1e-12,
                      ftol_rel = 1e-9,
                      ...) {

    Level <- NULL ## avoid warning in checks
    
    ## XXX CHANGE THE DEFAULT???
    ## if (missing(trace)) {
    ##     if (confintMethod = "proflik") trace <- 1
    ##     else trace <- 0
    ## }

    eps <-  1e-4
      
    out <- match.arg(out)
    parNames <- object$parNames
    distName <- object$distName
    p <- object$p
    thetaHat <- object$estimate
    
    if (is.null(period)) {
        period <- c(0.1, 0.2, 0.5, 1.0,  1.1, 1.5, 2, 5, 10, 20, 50, 75, 100,
                    125, 150, 175, 200, 250, 300, 500, 700, 1000) 
    } else {
        period <- sort(period)
    }

    ## =========================================================================
    ## The period 'T' must be such that the product p := 1 / lambdaHat
    ## / T is < 1.0. A value < 1.0 but too close to 1.0 leads to
    ## problems.
    ## =========================================================================
    
    ind <- (period * thetaHat[1]) > 1.5
    period <- period[ind]
    
    thetaHat <- object$estimate
    nPeriod <- length(period)
    prob <- 1.0 - 1.0 / period / thetaHat[1]
    fPeriod <- format(period)
    
    method <- match.arg(confintMethod)

    ## =========================================================================
    ## Cope with the check argument
    ## =========================================================================
   
    if (method == "delta" && check) {
        warning("Since method is \"delta\", no check is needed and 'check' is ",
                "set to FALSE")
        check <- FALSE
    }

    if (check) {
        message("Use the 'autoplot' method on the resulting object to ",
                "check the profile-likelihood results.")  }
    
    if (check && out == "array") {
        warning("Since 'check' is TRUE, 'out' set to \"data.frame\"")
        out <- "data.frame"
    }

    ## take into account the order. 
    indLevel <- order(level, decreasing = TRUE)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    if (method == "none") {
        
        RL <- array(NA, dim = c(Period = nPeriod, 1),
                    dimnames = list(Period = fPeriod, "Quant"))
        diagno <- NULL
        
        for (iPer in seq_along(period)) {
            RL[iPer, 1] <- object$threshold +
                Excd[[distName]]$qFun(p = prob[iPer],
                                      theta = thetaHat[-1])
        }
        
        
    } else if (method == "delta" || check) {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        covTheta <- vcov(object)
        p <- object$p
        
        if (is.null(covTheta) || any(is.na(covTheta))) {
            stop("vcov(object) must be a matrix with no NA. ",
                 "Consider changing 'confintMet'")
        }
        
        q <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)

        Quant <- rep(NA, nPeriod)
        
        RL <- array(NA,
                    dim = c(Period = nPeriod, Lim = 3L, Level = nLevel),
                    dimnames = list(Period = fPeriod,
                        Type = c("Quant", "L", "U"), Level = fLevel)) 

        diagno <- NULL
        
        for (iPer in seq_along(period)) {
            
            Quant[iPer] <- quant <- object$threshold +
                Excd[[distName]]$qFun(p = prob[iPer],
                                      theta = thetaHat[-1],
                                      deriv = TRUE)
            
            ## this is a 1 x p matrix
            gradTheta <- attr(quant, "gradient")

            ## =================================================================
            ## CAUTION: the density must be evaluated at a GPD quantile!
            ## =================================================================
            
            f <- Excd[[distName]]$dFun(x = Quant[iPer] - object$threshold,
                                       theta = thetaHat[-1],
                                       deriv = TRUE)
            
            gradTheta <- c("lambda" = 1 / f / thetaHat[1]^2 / period[iPer],
                           drop(gradTheta))
            sdRL <- sqrt(t(gradTheta) %*% covTheta %*% gradTheta) 
          
            RL[iPer, "Quant",  ] <- Quant[iPer] 
            RL[iPer, "L",  ] <- Quant[iPer] + outer(sdRL, q[ , 1L]) 
            RL[iPer, "U",  ] <- Quant[iPer] + outer(sdRL, q[ , 2L])     
            
        }
          
    }

    if (method == "proflik") {

        ## if check store the results of the delta method
        if (check) RLdelta <- RL
        
        constrCheck <- -5e-3
        constrCheckAbs <- 5e-3
        
        thetaHat <- object$estimate
        
        RL <- array(NA,
                    dim = c(Period = nPeriod, Lim = 3L, Level = nLevel),
                    dimnames = list(Period = fPeriod,
                        Type = c("Quant", "L", "U"), Level = fLevel))
        
        diagno <- array(NA,
                        dim = c(Period = nPeriod, Lim = 2L, Level = nLevel,
                            Diag = 3L),
                        dimnames = list(Period = fPeriod,
                            Type = c("L", "U"), Level = fLevel,
                            Diag = c("status", "constraint", "gradDist")))
        
        ## =====================================================================
        ## For each parameter, we maximise / minimise it under the
        ## constraint that the logLik remains >= max logLik - delta
        ## where delta := qchisq(1 - alpha) where 'alpha' is given by
        ## the confidence level.
        ## =====================================================================

        opts <- list()
        
        opts[[2]] <- list("algorithm" = "NLOPT_LD_AUGLAG",
                          "xtol_rel" = 1.0e-12,
                          "ftol_abs" = ftol_abs,
                          "ftol_rel" = ftol_rel,
                          "maxeval" = 8000,
                          "check_derivatives" = FALSE,
                          "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                              "xtol_rel" = 1.0e-12,
                              "maxeval" = 5000,
                              "ftol_abs" = ftol_abs,
                         "ftol_rel" = ftol_rel),
                          "print_level" = 0)
        
        opts[[1]] <- list("algorithm" = "NLOPT_LD_MMA",
                          "xtol_rel" = 1.0e-12,
                          "ftol_abs" = ftol_abs,
                          "ftol_rel" = ftol_rel,
                          "maxeval" = 5000,
                          "check_derivatives" = FALSE,
                          "print_level" = 0)
        
        if (trace >= 2) {
            opts[["check_derivatives"]] <- TRUE
            opts[["check_derivatives_print"]] <- "all"
        }
         
         ## ====================================================================
         ## note that some arguments such as 'level' are unused but are
         ## required by the constraint
         ## ====================================================================
         
        f <- function(theta, level, period, chgSign = FALSE, object) {

            prob <- 1.0 - 1.0 / theta[1] / period

            if (is.na(prob) || prob < 0.0 || prob > 1.0) {
                RL <- NA
            } else {
                RL <- object$threshold +
                    Excd[[distName]]$qFun(prob,
                                          theta = theta[-1],
                                          deriv = TRUE)
            }
            
            ## 'nloptr' fails on NA and NaN!
            if (is.na(RL)) {
                 if (chgSign) {
                     return(list("objective" = -Inf,
                                 "gradient" = rep(NaN, object$p)))
                 } else {
                     return(list("objective" = Inf,
                                 "gradient" = rep(NaN, object$p)))
                 }
             }
            
            gradTheta <- attr(RL, "gradient")
            
            ## add the derivative w.r.t. to the rate 'lambda', which
            ## comes by chain rule.
            
            f1 <- Excd[[distName]]$dFun(x = RL - object$threshold,
                                        theta = theta[-1],
                                        deriv = TRUE)
            gradTheta <- c("lambda" = 1.0 / f1 / theta[1]^2 / period,
                           drop(gradTheta))
            
            ## if (trace > 2) {
            ##     val <- RL
            ##     if (chgSign) val <- -val
            ##     ftol <- abs(val - valPrec) / valPrec
            ##     cat(sprintf("%d %7.4f %7.2f %7.4f, f = %16.14f ftol = %16.14f \n",
            ##                 count, theta[1], theta[2], theta[3], val, ftol))
            ##     count <<- count + 1
            ##     valPrec <<- val
            ## }
            
             if (chgSign) {
                 return(list("objective" = -RL, "gradient" = -gradTheta))
             } else {
                 return(list("objective" = RL, "gradient" = gradTheta))
             }
        }
         
        ## =====================================================================
        ## Up to (possibly unused) arguments, the constraint function
        ## is the same as that used in the 'confint' method.
        ## =====================================================================
         
        g <- function(theta, level, period, chgSign = FALSE, object) {
            
            ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(theta = theta, object = object,
                                       deriv = TRUE)
            grad <- attr(res, "gradient")
            list("constraints" = res[[1]] - ellL,
                 "jacobian" = grad)
            
        }
        
        ## =====================================================================
        ## note that although we recompute the gradient of the
        ## objective and the quantile of the chi-square distribution,
        ## this might be faster than re-defining the functions in the
        ## loop. Some experimentations would be needed to confirm
        ## this.
        ## =====================================================================
           
        for (iPer in seq_along(period)) {
            
            probi <- 1.0 - 1.0 / thetaHat[1] / period[iPer]
            
            if (!is.na(probi) && probi >= 0.0 && probi <= 1.0) {
                quant <- object$threshold +
                    Excd[[distName]]$qFun(p = probi,
                                          theta = thetaHat[-1],
                                          deriv = FALSE)
            } else {
                quant <- NA
            }
            
            for (iLev in rev(seq_along(level))) {
                lev <- level[iLev]    
                RL[iPer, "Quant", iLev] <- quant
            }
        }        
        
        ## =====================================================================
        ## Lower and upper bounds
        ## =====================================================================
        
        labs <- c("L" = "Lower", "U" = "Upper")
        sign <- c("L" = 1.0, "U" = -1.0)
        chgSign <- c("L" = 0.0, "U" = 1.0)
        
        if (trace) cat("\no Finding CI for Return Levels\n")        

        ## =====================================================================
        ## When check is TRUE, we evaluate the profil-likelihood on a
        ## grid of value of the Return Level rho(T) for each period T.
        ## We store the corresponding value in a vector and
        ## concatenate all vectors.
        ## =====================================================================

        if (check) {
            
            nSigma <-  rep(nSigma, length.out = 2)

            ## =================================================================
            ## Initialise vectors for the results and define a
            ## function fo be minimised.
            ## =================================================================
            
            rhoGrid <- numeric(0)
            periodGrid <- numeric(0)
            negLogLikCRho <- numeric(0)

            ## =================================================================
            ## Version with a gradient to be used with a "_LN_"
            ## algorithm such as 'COBYLA'
            ## 
            ## 'V' is the quantile of the distribution of a shape
            ## taken equal to 1.0
            ## =================================================================
            
            negLogLikNoRho <- function(thetaNoScale, period, iRho) {
                rho  <- rhoGridPer[iRho] - object$threshold
                theta <- rep(1.0, p)
                theta[-2] <- thetaNoScale
                q  <- 1.0 / theta[1] / period
                if (q > 1.0)  return(NaN)
                V <- Excd[[distName]]$qFun(q,
                                           theta[-1], lower.tail = FALSE)
                theta[2] <- rho / V
                if (all(is.finite(theta))) {
                    negLogLikFun(theta, object = object, deriv = FALSE)
                } else {
                    NaN
                }
            }

            ## =================================================================
            ## Version with a gradient to be used with a "_LD_"
            ## algorithm such as 'LBFGS'
            ##
            ## =================================================================
            
            negLogLikNoRhoGrad <- function(thetaNoScale, period, iRho) {
                
                rho  <- rhoGridPer[iRho] - object$threshold
                theta <- rep(1.0, p)
                theta[-2] <- thetaNoScale

                ## q := prob of exceedance
                q <- 1.0 / theta[1] / period
                ## print(c(theta[1], q))
                if (q > 1.0)  {
                    list("objective" = NaN,
                         "gradient" = rep(NaN, p - 1L))
                }

                ## =============================================================
                ## V is the quantile for a scale parameter equal to one, so
                ## the return level is rho = scale * V
                ## =============================================================
                
                V <- Excd[[distName]]$qFun(q, theta[-1],
                                           lower.tail = FALSE, deriv = TRUE)

                ## the first derivative is w.r.t. the scale parameter
                dVdthetaNoScale <- attr(V, "gradient")[-1]
                
                ## =============================================================
                ## use chain rule derivatition and implicit function derivation
                ##
                ##   1) dV / dlambda = (dV / dq)  * (dq / dlambda)
                ##   2) dV / dq = - 1 / f(q)
                ## =============================================================
                
                fV <- Excd[[distName]]$dFun(V, theta[-1])
                dVdLambda <- q / theta[1] / fV
                
                theta[2] <- rho / V
                
                negEll <- negLogLikFun(theta, object = object, deriv = TRUE)
                grad <- attr(negEll, "gradient")
                gradNoRho <- grad[1L, -2L]
                gradNoRho[1L] <- gradNoRho[1L] -
                    grad[1L, 2L] * rho * dVdLambda / V^2
                
                if (p > 2) {
                    ind <- 2:(p - 1L)
                    gradNoRho[ind] <- gradNoRho[ind] -
                        grad[1L, ind] * rho * dVdthetaNoScale / V^2
                }

                attributes(negEll) <- NULL
                
                if (all(is.finite(theta))) {
                    list("objective" = negEll,
                         "gradient" = gradNoRho)
                } else {
                    list("objective" = NaN,
                         "gradient" = rep(NaN, p -1L))
                }
            }

            ## =================================================================
            ## Tune the optimisation for the determination of the
            ## profile logLik. Note that we do not use gradients here.
            ## Recompile the package with check_derivatives = TRUE in
            ## case of doubt!!!
            ## =================================================================
            
            optsNoRho <- list("algorithm" = "NLOPT_LN_COBYLA",
                            "xtol_rel" = 1.0e-4,
                            "xtol_abs" = 1.0e-5,
                            "ftol_abs" = 1e-3,
                            "maxeval" = 3000, "print_level" = 0,
                              "check_derivatives" = FALSE)
            
            optsNoRhoGrad <- list("algorithm" = "NLOPT_LD_LBFGS",
                                  "xtol_rel" = 1.0e-8,
                                  "xtol_abs" = 1.0e-8,
                                  "ftol_abs" = 1e-5,
                                  "maxeval" = 3000, "print_level" = 0,
                                  "check_derivatives" = FALSE)
            
            for (iPer in seq_along(period)) {
                
                Lper <- RLdelta[iPer, "Quant", 1] - nSigma[1] *
                    (RLdelta[iPer, "U", 1] - RLdelta[iPer, "Quant", 1])
                Uper <- RLdelta[iPer, "Quant", 1] + nSigma[2] *
                    (RLdelta[iPer, "U", 1] - RLdelta[iPer, "Quant", 1])
                
                rhoGridPer <- seq(from = Lper, to = Uper, length.out = nCheck)
                periodGridPer <- rep(period[iPer], nCheck)
                negLogLikCRhoPer <- rep(NA, nCheck)
                
                for (iRho in seq_along(rhoGridPer)) {
                    resii <-  try(nloptr(x0 = thetaHat[-2],
                                         eval_f = negLogLikNoRhoGrad,
                                         lb = object$lb[-2],
                                         ub = object$ub[-2],
                                         opts = optsNoRhoGrad,
                                         period = period[iPer],
                                         iRho = iRho), silent = TRUE)
                    
                    if (!inherits(resii, "try-error") &&
                        (resii$status %in% 1:4)) {
                        if (trace > 1) {
                            cat("Optim #1 succesful!\n")
                        }
                        
                        negLogLikCRhoPer[iRho] <- resii$objective
                    } else {
                        if (trace > 1) {
                            ## print(resii)
                            cat(sprintf(paste0("Retrying optim -> no deriv. ",
                                               "T = %4.0f rho = %6.2f\n"),
                                        period[iPer],
                                        rhoGridPer[iRho]))
                        }
                        
                        resii <-  try(nloptr(x0 = thetaHat[-2],
                                             eval_f = negLogLikNoRho,
                                             lb = object$lb[-2],
                                             ub = object$ub[-2],
                                             opts = optsNoRho,
                                             period = period[iPer],
                                             iRho = iRho), silent = TRUE)
                        if (!inherits(resii, "try-error") &&
                            (resii$status %in% 1:4)) {
                            negLogLikCRhoPer[iRho] <- resii$objective
                        }
                    }
                }
                
                rhoGrid <- c(rhoGrid, rhoGridPer)
                periodGrid <- c(periodGrid, periodGridPer)
                negLogLikCRho <- c(negLogLikCRho, negLogLikCRhoPer)

            }

            negLogLikC <- data.frame(Period = periodGrid,
                                     rho = rhoGrid,
                                     Value = negLogLikCRho)
            
            
        }
        
        for (LU in c("L", "U")) {
            
            if (trace) {
                cat(sprintf("\n**************\n %s bounds \n**************\n",
                            labs[LU]))
            }
            
            for (iLev in rev(seq_along(level))) {
                
                lev <- level[iLev]

                if (trace) {
                    cat(sprintf("\nConfidence Level: %5.2f\n", lev))
                }
                
                ## =============================================================
                ## 2017-08-20 It was found that the convergence is much
                ## easier when the periods 'T' are taken in reverse order
                ## =============================================================
                for (iPer in rev(seq_along(period))) {
                    
                    if (trace) {
                        cat(sprintf("\n    - Period:  %5.1f\n", period[iPer]))
                    }

                    ## 2020-04-02 This turns out to be worst than using the
                    ## estimate!
                    ## 
                    ## if ((iPer < length(period)) && !is.null(thetaIniPrec)) {
                    ##     ## if ((iPer > 1L) && !is.null(psiIniPrec)) {
                    ##     theta0 <- thetaIniPrec
                    ## } else {
                    ##     theta0 <- thetaHat
                    ## }
                    
                    optDone <- FALSE
                    optNum <- 1

                    theta0 <- thetaHat
                    
                    while (!optDone && (optNum <= 2)) {
                        
                        if (trace && (optNum > 1)) {
                            cat("        <retrying optimisation!>\n")
                        }

                        ## if (trace > 2) {
                        ##     count <- 0
                        ##     valPrec <- NA
                        ## }
                        resOpt <- try(nloptr::nloptr(x0 = theta0,
                                                     eval_f = f,
                                                     eval_g_ineq = g,
                                                     lb = object$lb,
                                                     ub = object$ub,
                                                     level = lev,
                                                     period = period[iPer],
                                                     chgSign = chgSign[LU],
                                                     opts = opts[[optNum]],
                                                     object = object),
                                      silent = TRUE)
                        
                        if (!inherits(resOpt, "try-error")) {
                            
                            diagno[iPer, LU, iLev, "status"] <- resOpt$status
                        
                            if (trace == 1L) {
                                cat(sprintf("        Optimisation status: %d\n",
                                            resOpt[["status"]]))
                                cat(sprintf("        Iterations: %d\n",
                                            resOpt[["iterations"]]))
                                names(resOpt$solution) <- object$parNames
                                cat(sprintf("        Objective: %7.2f\n",
                                            resOpt[["objective"]]))
                            } else  if (trace > 2L) {
                                cat("\nSOLUTION\n")
                                print(resOpt)
                            }
                            
                            ## The constraint must be active. We have to check
                            ## that!
                            checkg <- g(theta = resOpt$solution,
                                        level = lev,
                                        period = period[iPer],
                                        chgSign = chgSign[LU],
                                        object = object)
                            
                            diagno[iPer, LU, iLev, "constraint"] <-
                                checkg$constraints
                            
                            ## The gradient of the objective must be colinear to the
                            ## jacobian of the constraint. We have to check that!
                            checkf <- f(theta = resOpt$solution,
                                        level = lev,
                                        period = period[iPer],
                                        chgSign = chgSign[LU],
                                        object = object)
                            
                            gradDist <- distLines(x1 = checkg$jacobian,
                                                  x2 = checkf$gradient)
                            
                            diagno[iPer, LU, iLev, "gradDist"] <- gradDist
                            
                            if (trace) {
                                cat(sprintf("        Optim status : %d\n",
                                            resOpt$status))
                                cat(sprintf("        Constraint check %10.7f\n",
                                            checkg$constraints))
                                cat(sprintf("        Gradient directions: %7.4f\n",
                                            gradDist))
                                if (trace == 2) {
                                    cat("        Gradients\n")
                                    print(rbind("        g" = checkg$jacobian,
                                                "        f" = checkf$gradient))
                                    cat("\n")
                                }
                            }
                            
                            ## It seems that the distance reached is smaller when
                            ## 'T' is large.

                        } else {
                            if (trace > 1) {
                                cat("'resOpt' is an error!!!\n")
                            }
                            ## test <- g(theta0, level = lev, period = period[iPer],
                            ##           chgSign = chgSign[LU],
                            ##           object = object)
                            ## print(test)
                        }
                        
                        gradLim <- 1.0 / period^0.6
                         
                        if (!inherits(resOpt, "try-error") &&
                            (resOpt$status %in% 1:4) &&
                            (all(abs(checkg$constraints) < constrCheckAbs))) {
                            ## && (!is.na(gradDist)) &&
                            ## (gradDist < gradLim)) {
                            
                            optDone <- TRUE
                            thetaIniPrec <- resOpt[["solution"]]
                            RL[iPer, LU, iLev] <- sign[LU] * resOpt[["objective"]]
                            
                        } else {

                            ## try another initial value. We could use the value
                            ## of theta obtained for a smaller confidence level?
                            if ((iPer < length(period) &&
                                 diagno[iPer + 1, LU, iLev, "status"] %in% 1:4)) {
                                theta0 <- thetaIniPrec
                                thetaIniPrec <- NULL
                            } 
                            optNum <- optNum + 1
                        }
                        
                    }
                    
                }
                
            }
        }
        

    }
    
    if (out == "data.frame") {
        
        if (requireNamespace("reshape2", quietly = TRUE)) {
            
            if (method == "none") {
                RL <-
                    reshape2::melt(RL[ , 1, drop = FALSE],
                                   value.name = "Quant",
                                   varnames = c("Period"))[ , c("Period", "Quant")]
            } else {
                ## =============================================================
                ## UGGLY CODE: there must be a simpler and more
                ## efficient way of doing that. The problem is that we
                ## want to drop the " Type" dimension but not the
                ## "Level" dimension even when it has no extension
                ## =============================================================
                
                df <- list()
                for (nm in c("Quant", "L", "U")) {
                    RL1 <- RL[ , nm, , drop = FALSE]
                    df[[nm]]  <-
                        reshape2::melt(RL1,
                                       value.name = nm,
                                       varnames = c("Period", "Type", "Level"))
                }
                RL <- data.frame(df[["Quant"]][, c("Period", "Level", "Quant")],
                                 L = df[["L"]]$L, U = df[["U"]]$U)
                
                ind <- with(RL, order(Level, Period))
                RL <- RL[ind, , drop = FALSE]
            }
            
        } else {
            stop("the package 'reshape2' could not be used")
        }

        for (nm in c("nOT", "yOT", "estimate", "threshold")) {
            attr(RL, nm) <- object[[nm]]
        }
   
        class(RL) <- c("RL.poisGP", "data.frame")
    }

    if (check) {
        
        RL <- reshape2::melt(RL,
                             value.name = "Value",
                             id.vars = c("Period", "Level"),
                             variable.name = c("LU"))

        RL1 <- subset(RL, LU != "Quant")
        nll <- object$negLogLik + qchisq(level, df = 1) / 2.0
        names(nll) <- fLevel
        RL1 <- cbind(RL1, NegLogLik = nll[as.character(RL1$Level)])
        
        RL2 <- subset(RL, (LU == "Quant") & (Level == fLevel[1])) 
        RL2$Level <- "est"
        RL2 <- cbind(RL2, NegLogLik = rep(-object$logLik, nPeriod))
                       
        RL <- rbind(RL1, RL2, deparse.level = 1)
        ## XXXY attributes neede here
        
        L <- list(RL = RL,
                  negLogLikC = negLogLikC,
                  ylim = -object$logLik + 3* (max(nll) + object$logLik))
        
        class(L) <- "RLCheck.poisGP"
        return(L)
        
    } else {
        
        attr(RL, "title") <- "Return Levels"
        attr(RL, "threshold") <- object$threshold
        
        return(RL)

    }
        
}


