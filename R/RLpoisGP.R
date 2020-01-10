## ***************************************************************************
##' Compute return levels along with confidence bounds for a Poisson-GP
##' model with ML inference results.
##'
##' The return level curve corresponding to the column or to the
##' dimension named \code{"Quant"} is obtained by plugging the ML
##' estimate of the poisson-GP parameters in the quantile. To a
##' certain extension, this can be compared to a Bayesian MAP when a
##' flat prior is used. However this is not exactly the case because
##' the ML Estimate is invariant under reparameterising, while a flat
##' prior refers to a specific parameterisation, hence so does the
##' corresponding MAP. The profile-likelihood is obtained by using a
##' specific algorithm based on constrained optimisation, and not on
##' the usual method of re-parameterisartion with quantiles.
##' 
##' @title Return Levels and Confidence Intervals for Poisson-GP Model
##'
##' @param object An object with class \code{poisGPML} representing
##' the inference results for a Poisson-GP model.
##'
##' @param period A vector of periods for which the return levels will
##' be computed.
##'
##' @param level Level of confidence.
##'
##' @param confintMethod The method used to compute the confidence
##' intervals.
##'
##' @param out The type of outpout wanted, see the \strong{Value}
##' section.
##'
##' @param biasCorrect Logical. Should bias correction be applied in
##' bootstrap?
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Not used yet.
##'
##' @return An object inheriting either from \code{data.frame} or from
##' \code{array} containing the results. In the first case, the object
##' has class \code{"RL.poisGP"}, which allows the use of
##' \code{autoplot} method without recomputing the results.
##'
##' @author Yves Deville
##' 
RL.poisGP <- function(object,
                      period = NULL,
                      level = 0.70,
                      confintMethod = c("delta", "none", "boot", "proflik"),
                      out = c("data.frame", "array"),
                      biasCorrect = FALSE,
                        trace = 1L,
                      ...) {
    
    eps <-  1e-4
      
    out <- match.arg(out)
    parNames <- c("lambda", "scale", "shape")
    thetaHat <- object$estimate
    
    if (is.null(period)) {
        period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 75, 100,
                    125, 150, 175, 200, 250, 300, 500, 700, 1000) 
    } else {
        period <- sort(period)
    }

    ind <- (period * thetaHat[1]) > 1.0
    period <- period[ind]
    
    thetaHat <- object$estimate
    nPeriod <- length(period)
    prob <- 1.0 - 1.0 / period / thetaHat[1]
    fPeriod <- format(period)
    
    method <- match.arg(confintMethod)
    if (!missing(biasCorrect) && (method != "boot")) {
        warning("the argument 'biasCorrect' provided will not ",
                "be used since 'confintMethod' != \"boot\"")
    }
        
    ## take into account the order. 
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    if (method == "none") {
        
        RL <- array(NA, dim = c(Period = nPeriod),
                    dimnames = list(Period = fPeriod))
        diagno <- NULL
        
        for (iPer in seq_along(period)) {
            RL[ , iPer] <- object$threshold +
                qGPD2(p = prob[iPer], scale = thetaHat[2], shape = thetaHat[3])
        }
        
        
    } else if (method == "delta") {
        
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
                qGPD2(p = prob[iPer], scale = thetaHat[2], shape = thetaHat[3],
                      deriv = TRUE)
            ## this is a 1 x p matrix
            gradTheta <- attr(quant, "gradient")

            ## CAUTION: the density must be evaluated at a GPD2 quantile!
            f <- dGPD2(x = Quant[iPer] - object$threshold,
                       scale = thetaHat[2], shape = thetaHat[3], deriv = TRUE)
            
            gradTheta <- c("lambda" = 1 / f / thetaHat[2]^2 / period[iPer],
                           drop(gradTheta))
            sdRL <- sqrt(t(gradTheta) %*% covTheta %*% gradTheta) 
          
            RL[iPer, "Quant",  ] <- Quant[iPer] 
            RL[iPer, "L",  ] <- Quant[iPer] + outer(sdRL, q[ , 1L]) 
            RL[iPer, "U",  ] <- Quant[iPer] + outer(sdRL, q[ , 2L])
            
        }
        
        ## change dim order: "Date", then "Level" then "Period" then "L or U"
        ## RL <- aperm(a = RL, perm = c(1, 4, 2, 3))
        
    } else if (method == "boot") {

        stop("method = \"boot\" is not implemented yet")
        
    } else if (method == "proflik") {
        
        constrCheck <- -5e-3
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
        
        ## ===================================================================
        ## For each parameter, we maximise/minimise it under the constraint
        ## that the logLik remains >= max logLik - delta where delta :=
        ## qchisq(1 - alpha) where alpha is given by the cofidence level.
        ##
        ## ===================================================================

        opts <- list()
        opts[[1]] <- list("algorithm" = "NLOPT_LD_MMA",
                          "xtol_rel" = 1.0e-5,
                          "ftol_abs" = 1.0e-7, "ftol_rel" = 1.0e-5,
                          "maxeval" = 200,
                          "check_derivatives" = FALSE,
                          "print_level" = 0)
        
        opts[[2]] <- list("algorithm" = "NLOPT_LD_AUGLAG",
                          "xtol_rel" = 1.0e-5,
                          "ftol_abs" = 1.0e-7, "ftol_rel" = 1.0e-5,
                          "maxeval" = 500,
                          "check_derivatives" = FALSE,
                          "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                              "xtol_rel" = 1.0e-5,
                              "maxeval" = 1000,
                              "ftol_abs" = 1.0e-7,
                              "ftol_rel" = 1.0e-5),
                          "print_level" = 0)
        
        
        if (trace >= 2) {
            opts[["check_derivatives"]] <- TRUE
            opts[["check_derivatives_print"]] <- "all"
        }
         
         ## ==============================================================
         ## note that some arguments such as 'level' are unused but are
         ## required by the constraint
         ## ===============================================================
         
        f <- function(theta, level, period, chgSign = FALSE, object) {

            
            prob <- 1.0 - 1.0 / theta[1] / period

            if (is.na(prob) || prob < 0.0 || prob > 1.0) {
                RL <- NA
            } else {
                RL <- object$threshold +
                    qGPD2(prob,
                          scale = theta[2], shape = theta[3],
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
            
            f1 <- dGPD2(x = RL - object$threshold,
                        scale = theta[2], shape = theta[3], deriv = TRUE)
            
            gradTheta <- c("lambda" = 1.0 / f1 / theta[2]^2 / period,
                           drop(gradTheta))
            
             if (chgSign) {
                 return(list("objective" = -RL, "gradient" = -gradTheta))
             } else {
                 return(list("objective" = RL, "gradient" = gradTheta))
             }
        }
         
        ## ==============================================================
        ## Up to (possibly unused) arguments, the constraint function
        ## is the same as that used in the 'confint' method.
        ## ===============================================================
         
        g <- function(theta, level, period, chgSign = FALSE, object) {
            
            ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(theta = theta, object = object, deriv = TRUE)
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
                    qGPD2(p = probi,
                      scale = thetaHat[2], shape = thetaHat[3],
                          deriv = FALSE)
            } else {
                quant <- NA
            }
            
            for (iLev in rev(seq_along(level))) {
                lev <- level[iLev]    
                RL[iPer, "Quant", iLev] <- quant
            }
        }        
        
        ## =========================================================
        ## Lower and upper bounds
        ## ==========================================================
        
        labs <- c("L" = "Lower", "U" = "Upper")
        sign <- c("L" = 1.0, "U" = -1.0)
        chgSign <- c("L" = 0.0, "U" = 1.0)
        
        if (trace) cat("\no Finding CI for Return Levels\n")        
        
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
                
                ## ========================================================
                ## 2017-08-20 It was found that the convergence is much
                ## easier when the periods 'T' are taken in reverse order
                ## ========================================================
                for (iPer in rev(seq_along(period))) {
                    
                    if (trace) {
                        cat(sprintf("\n    - Period:  %5.1f\n", period[iPer]))
                    }
                    
                    if ((iPer < length(period)) && !is.null(thetaIniPrec)) {
                        ## if ((iPer > 1L) && !is.null(psiIniPrec)) {
                        theta0 <- thetaIniPrec
                    } else {
                        theta0 <- thetaHat
                    }
                    
                    optDone <- FALSE
                    optNum <- 1
                    
                    while (!optDone && (optNum <= 2)) {
                        
                        if (trace && (optNum > 1)) {
                            cat("        <retrying optimisation!>\n")
                        }
                        
                        resOpt <- try(nloptr::nloptr(x0 = theta0,
                                                     eval_f = f,
                                                     eval_g_ineq = g,
                                                     level = lev,
                                                     period = period[iPer],
                                                     chgSign = chgSign[LU],
                                                     opts = opts[[optNum]],
                                                     object = object))
                        
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
                            
                            ## The constraint must be active. We have to check that!
                            checkg <- g(theta = resOpt$solution,
                                        level = lev,
                                        period = period[iPer],
                                        chgSign = chgSign[LU],
                                        object = object)
                            
                            diagno[iPer, LU, iLev, "constraint"] <- checkg$constraints
                            
                            ## The gradient of the objective must be colinear to the
                            ## jacobian of the constraint. We have to check that!
                            checkf <- f(theta = resOpt$solution,
                                        level = lev,
                                        period = period[iPer],
                                        chgSign = chgSign[LU],
                                        object = object)
                            
                            gradDist <- NSGEV::distLines(x1 = checkg$jacobian,
                                                         x2 = checkf$gradient)
                            
                            diagno[iPer, LU, iLev, "gradDist"] <- gradDist
                            
                            if (trace) {
                                cat(sprintf("        Constraint check %10.7f\n",
                                            checkg$constraints))
                                cat(sprintf("        Gradient directions: %7.4f\n", gradDist))
                                if (trace == 2) {
                                    cat("        Gradients\n")
                                    print(rbind("        g" = checkg$jacobian,
                                                "        f" = checkf$gradient))
                                    cat("\n")
                                }
                            }
                            
                            ## It seems that the distance reached is smaller when
                            ## 'T' is large.

                        }
                        
                        gradLim <- 1.0 / period^0.6
                        
                        if (!inherits(resOpt, "try-error") &&
                            (resOpt$status %in% c(3, 4)) &&
                            (all(checkg$constraints > constrCheck))) {
                                ## && (!is.na(gradDist)) &&
                                ## (gradDist < gradLim)) {
                            
                            optDone <- TRUE
                            thetaIniPrec <- resOpt[["solution"]]
                            RL[iPer, LU, iLev] <- sign[LU] * resOpt[["objective"]]
                            
                        } else {
                            optNum <- optNum + 1
                            thetaIniPrec <- NULL
                        }
                        
                    }
                    
                }
                
            }
        }
        

    }
    
    if (out == "data.frame") {
        
        if (requireNamespace("reshape2", quietly = TRUE)) {
            
            if (method == "none") {
                RL <- reshape2::melt(RL, value.name = "Quant",
                                     varnames = c("Date", "Period"))
            } else {
                ## ============================================================
                ## UGGLY CODE: there must be a simpler and more
                ## efficient way of doing that. The problem is that we
                ## want to drop the " Type" dimension but not the
                ## "Level" dimension even when it has no extension
                ## ============================================================
                
                df <- list()
                for (nm in c("Quant", "L", "U")) {
                    RL1 <- RL[ , nm, , drop = FALSE]
                    df[[nm]]  <- reshape2::melt(RL1,
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

    attr(RL, "title") <- "Return Levels"
    attr(RL, "threshold") <- object$threshold
    
    return(RL)
    
}


