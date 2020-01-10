
coefIni <- function(object, ...) {
    UseMethod("coefIni", object) 
}


## ****************************************************************************
##' Compute initial estimates for the parameters of a \code{poisGP}
##' model.
##'
##' The estimates are obtained as a weighted mean of up to three cheap
##' estimates corresponding to the different parts of the data: main
##' sample, MAX blocks and OTS blocks. In each cases, it is assumed
##' that the shape is zero. For the main sample the estimates of the
##' rate \code{lambda} and \code{sigma} are simply the exceedance rate
##' and the mean of the excesses. For MAX and OTS blocks, the
##' estimates are based on regression approximations as used in
##' \strong{Renext} package.
##' 
##' @title Initial Estimates for a \code{poisGP} Model
##'
##' @param object A \code{poisGP} object the parameters of which need
##' to be estimated.
##'
##' @param ... Not used yet.
##'
##' @return A vector of parameters for the \code{poisGP} object.
##'
##' 
coefIni.poisGP <- function(object, trace = 0, ...) {

    fd <- object$fitData
    
    if (fd$OT$flag && (fd$OT$n > 0)) {
        r <- fd$OT$n
        theta <- c("lambda" = r / fd$OT$effDuration,
                   "scale" = mean(fd$OT$data),
                   "shape" = 0.0)
        return(theta)
    } else {
        r <-  0L
        theta <- c("lambda" = 0.0, "scale" = 0.0, "shape" = 0.0)
    }
    
    if (fd$MAX$flag && (rMAX <- sum(fd$MAX$r) / 2)) {

        thetaMAX <- Renext::parIni.MAX(MAX = fd$MAX, threshold = 0.0,
                                       distname.y = "GPD")
        if (trace) {
            cat("\nInitial value using MAX data only\n")
            print(thetaMAX)
        }
        
        theta <- (r * theta + rMAX * thetaMAX) / (r + rMAX)
        r <- r + rMAX
    } 
    
    if (fd$OTS$flag && (rOTS <- sum(fd$OTS$r) / 2)) {
        
        thetaOTS <- Renext::parIni.OTS(OTS = fd$OTS, threshold = 0.0,
                                       distname.y = "GPD")
        if (trace) {
            cat("\nInitial value using OTS data only\n")
            print(thetaOTS)
        }
        
        theta <- (r * theta + thetaOTS) / (r + rOTS)
    } 

    theta
    
}

summary.poisGP <- function(object, ...) {
    out <- object
    class(out) <- "summary.poisGP"
    out
}

print.summary.poisGP <- function(x, ...) {

    cat("Poisson-GP model fitted Maximum-Likelihood\n")
    cat("\no Threshold: ", x$threshold, "\n")
    cat("\no Poisson-GP parameters and standard error\n")
    print(cbind(estimate = x$estimate, sd = x$sd))
    cat("\no Correlation between Poisson-GP estimates\n")
    print(cov2cor(x$cov))
    
}

logLik.poisGP <- function(object, ...) {
    res <- object$logLik
    attr(res, "df") <- object$df
    attr(res, "nobs") <- object$nobs
    res
}

##====================================================================
## AIC methods
##====================================================================

AIC.poisGP <- function(object, ..., k = 2) {
    return(NextMethod())
}

BIC.poisGP <- function(object, ...) {
    return(NextMethod())
}


## ****************************************************************************
##' Extract model coefficients for a Poisson-GP object.
##' 
##' @title Extract Model Coefficients for a Poisson-GP Object
##'
##' @param object A model with class \code{"poisGP"}.
##'
##' @param type The parameterisation to use.
##' 
##' @param ... Not used yet.
##'
##' @return A numeric vector of three coefficients.
##' 
coef.poisGP <- function(object, type = c("poisGP", "PP"), ...) {
    type <- match.arg(type)
    if (type == "poisGP") {
        return(object$estimate)
    } else {
        return(object$PP$estimate)
    }
}

## ****************************************************************************
##' Variance-covariance matrix for a Poisson-GP object.
##' 
##' @title Variance-Covariance Matrix  for a Poisson-GP Object
##'
##' @param object A model with class \code{"poisGP"}.
##'
##' @param type The parameterisation to use.
##' 
##' @param ... Not used yet.
##'
##' @return A covariance matrix of size three.
##' 
vcov.poisGP <- function(object, type = c("poisGP", "PP"), ...) {
    type <- match.arg(type)
    if (type == "poisGP") {
        return(object$cov)
    } else {
        return(object$PP$cov)
    }
}

## ****************************************************************************
##' Create a Poisson-GP model object, usually by ML estimation.
##'
##' The
##' \itemize{
##' \item{\code{estim = "optim"}}{
##'
##' The classical \code{stats::optim} function is used with
##' \code{method ="BFGS"}. The derivatives are not used and nor do the
##' bounds on the parameters given in \code{coefLower} and
##' \code{coefUpper}.
##'
##' }
##' \item{\code{estim = "nloptr"}}{
##'
##' The \code{nloptr::nloptr} function is used with the
##' \code{"NLOPT_LD_BFGS"} algorithm option. The derivatives are used
##' as well as the bounds on the parameters, leading to "box
##' constraints". The bounds can be used to fix the value of a GP
##' parameter by using the same value in \code{coefLower} and
##' \code{coefUpper}. For instance an exponential distribution can be
##' fitted by using a zero value for the shape both in
##' \code{coefLower} and \code{coefUpper}.
##' 
##' }
##' \item{\code{estim = "eval"}}{
##'
##' No optimisation is performed: the rate \code{lambda} corresponding
##' to the provided GP parameters is computed and the negative
##' log-likelihood and its first derivatives are evaluated, allowing
##' the determination of a (putative) covariance matrix for the
##' estimates. The named vector \code{parIni} should then contain
##' values for the Poisson-GP parameters, and valid values for the GP
##' parameters \code{"scale"} and \code{"shape"}. This possibility can
##' be used to check the results provided by other packages, e.g. to
##' recompute return levels. Note however that the provided parameters
##' may not be approximately maximising the likelihood and the
##' corresponding results will then be misleading.
##'
##' }
##' 
##' \item{\code{estim = "none"}}{
##'
##' No optimisation is performed. The log-likelihood and negative
##' log-likelihoods remain NA, and the initial values are ignored.
##' 
##' }}
##' 
##' @title Create a Poisson-GP Model Object
##' 
##' @param data A vector containing the observations for the "main" sample.
##'
##' @param threshold The \emph{main threshold} as commonly understood
##' in POT.
##' 
##' @param effDuration The effective duration of the observation
##' period corresponding to the \code{data}.
##'
##' @param MAX.data A list of numeric vectors corresponding to periods
##' or \emph{blocks}. Each vector contains the \eqn{r}-largest
##' observations for the block, where \eqn{r>0} can vary across
##' blocks. When a numeric vector is passed instead of a list, it is
##' understood that there is only one MAX block.
##'
##' @param MAX.effDuration A vector of positive numbers giving the
##' durations of the MAX blocks (in the same order). 
##' 
##' @param OTS.data A list of numeric vectors corresponding to periods
##' or \emph{blocks}. Each vector contains all observations for the
##' block which exceeded the corresponding threshold as given in
##' \code{OTS.threshold}. So the number \eqn{r \eg 0} vary across
##' blocks. When a numeric vector is passed instead of a list, it is
##' understood that there is only one OTS block.
##'
##' @param OTS.threshold A vector of positive numbers giving the
##' thresholds of the OTS blocks (in the same order). By construction
##' all the elements of \code{OTS.data[i]} are larger than
##' \code{OTS.threshol[i]} for each block index \code{i}.
##' 
##' @param OTS.effDuration A vector of positive numbers giving the
##' durations of the OTS blocks (in the same order). 
##' 
##' @param parIni A named parameter vector. This will be used to set
##' the values if \code{estim} is \code{"none"} or to provide initial
##' values in the other cases. When \code{parIni} is \code{NULL} the
##' parameter vector stored in the object as will contain \code{NA}
##' when \code{estim} is \code{"none"} or "good" initial values found
##' by devoted function.
##' 
##' @param estim Method or XXX
##'
##' @param coefLower,coefUpper Named vectors of bounds for the
##' parameters. The values can be infinite \code{Inf} or \code{-Inf}.
##' However, note that without bounds on the shape parameter \eqn{xi}
##' the maximum likelihood is infinite and obtained for \eqn{xi =
##' -1}. Note also that the bounds are ignored when \code{estim}
##' is set to \code{"optim"}.
##' 
##' @param scale Logical. If \code{TRUE} the observations in
##' \code{data}, \code{MAX.data} and \code{OTS.data} are all divided
##' by a common positive number in order to avoid numerical
##' problems. This number is returned as the \code{scaleData} element
##' of the returned list. Except from the numerical problems, the
##' value of the scale does not impact the results such as the
##' estimates of their covariance.
##' 
##' @param trace Integer level of verbosity.
##'
##' @return A list with among which the following objects are found
##'
##' \item{data}{
##'
##' A copy of the data provided in \code{data} and the other formals
##' \code{MAX.*} or \code{OTS.*}
##'
##' }
##' \item{fitData}{
##'
##' A modified version of \code{data} where the observations that do
##' not exceed the main threshold \eqn{u} are discarded and each
##' remaining observation \eqn{y_i} is replaced by the corresponding
##' excess \eqn{y_i - u}. So only positive observations are found in
##' the data vectors.
##'
##' }
##' 
##' @author Yves Deville
##'
##' @examples
##' fit <- poisGP(data = Garonne$OTdata$Flow, threshold = 2900,
##'               effDuration = 65,
##'               MAX.data = Garonne$MAXdata$Flow,
##'               MAX.effDuration = 143,
##'               estim = "none")
##' 
poisGP <- function(data = NULL,
                   threshold,
                   effDuration,
                   MAX.data = NULL,
                   MAX.effDuration = NULL,
                   OTS.data = NULL,
                   OTS.threshold = NULL,
                   OTS.effDuration = NULL,
                   ## distname.y = "GPD",
                   parIni = NULL,
                   estim = c("optim", "nloptr", "eval", "none"),
                   coefLower = c("scale" = 0.0, "shape" = -0.90),
                   coefUpper = c("scale" = Inf, "shape" = Inf),
                   scale = FALSE,
                   trace = 0) {

    eps <- 1e-10
    estim <-  match.arg(estim)

    if (is(data, "Rendata")) {
        stop("Using 'data' with class \"Rendata\" is not allled yet")
    }

    ## Note that we do not pass 'threshold' here.
    data <- potData(data = data,
                    effDuration = effDuration,
                    MAX.data = MAX.data,
                    MAX.effDuration = MAX.effDuration,
                    OTS.data = OTS.data,
                    OTS.threshold = OTS.threshold,
                    OTS.effDuration = OTS.effDuration)
    
    grandMin <- min(unlist(sapply(data, function(x) unlist(x$data))))

    scaleData <- 1.0
    
    if (missing(threshold)) {
        if (data$OT$flag || data$OTS$flag) {
            warning("'threshold' is missing and set just below the ",
                    "smallest obsevation")
        }
        threshold <- grandMin - eps
        fitData <- data
    } else if (threshold < grandMin) {
        warning("'threshold' is smaller than the smallest observation") 
    }

    ## ========================================================================
    ## Note that even when 'scale' is FALSE, 'scaleData' contains a
    ## value that can be used in 'control$sparscale' when estim is set
    ## to "optim". All in all the scaling is useful when estim is
    ## "nloptr".
    ## ========================================================================
    
    fitData <- threshData(threshold = threshold, data, scale = scale)
    scaleData <- attr(fitData, "scale")

    if (trace) {
        if (scale) {
            cat("\nThe data will be scaled by dividing then by ", scaleData, "\n")
        } else {
            cat("\nThe data will not be scaled\n")
        }
    }


    if ((!missing(coefLower) || !missing(coefUpper)) && (estim == "optim")) {
        warning("Provided values for 'coefLower' and 'coefUpper' will be ignored ",
                "since estim is set to \"optim\"")
    }
    
    
    ## ========================================================================
    ## Create a temporary object with class "poisGP" that can be used
    ## in some methods, not all be cause the object is not complete.
    ## ========================================================================
    
    thisPoisGP <- list(call = match.call(),
                       threshold = threshold,
                       data = data,
                       fitData = fitData,
                       p = 3L,
                       df = 3L,
                       parNames = c("lambda", "scale", "shape"),
                       scale = scale,
                       scaleData = scaleData)

    class(thisPoisGP) <- "poisGP"
    
    if (estim != "none") {

        if (is.null(parIni)) {
            parIni <- coefIni(thisPoisGP, trace = pmax(0, trace - 1))[-1]
        }
        
        if (trace) {
            cat("\nInitial values for the parameter\n")
            print(parIni)
        }
        
        res <- MLE.poisGP(object = thisPoisGP,
                          parIni = parIni,
                          estim = estim,
                          coefLower = coefLower,
                          coefUpper = coefUpper,
                          scale = scale,
                          trace = trace)
        
        for (nm1 in names(res)) {
            thisPoisGP[[nm1]] <- res[[nm1]]
        }
        
    } else {
        thisPoisGP$estimate <- c("lambda" = NA, "scale" = NA, "shape" = NA)
        thisPoisGP$negLogLik <- NA
        thisPoisGP$logLik <- NA
    }

    thisPoisGP$negLogLikFun <- negLogLikFun

    class(thisPoisGP) <- "poisGP"
    thisPoisGP
    
}
