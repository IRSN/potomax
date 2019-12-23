
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
coefIni.poisGP <- function(object, ...) {

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
        theta <- (r * theta +
            rMAX * Renext::parIni.MAX(MAX = fd$MAX, threshold = 0.0,
                                      distname.y = "GPD") ) / (r + rMAX)
        r <- r + rMAX
    } 
    
    if (fd$OTS$flag && (rOTS <- sum(fd$OTS$r) / 2)) {
        theta <- r * theta +
            rOTS * Renext::parIni.OTS(OTS = fd$OTS, threshold = 0.0,
                                      distname.y = "GPD") / (r + OTS)
    } 

    theta
    
}


## ****************************************************************************
##' Create a Poisson-GP model object, usually by ML estimation.
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
##' @param MAX.duration A vector of positive numbers giving the
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
##' @param OTS.duration A vector of positive numbers giving the
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
##' However, note that without bounds on the shape parameter \eqn{xi} the
##' maximum likelihood is infinite and obtained for \eqn{xi = -1}.
##'
##' @param trace Integer level of verbosity.
##'
##' @return A list with among which
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
##' 
##' fit <- poisGP(data = Garonne$OTdata$Flow, threshold = 2900,
##'               effDuration = 65,
##'               MAX.data = Garonne$MAXdata$Flow,
##'               MAX.effDuration = 143,
##'               estim = "none")
##' 
##'
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
                   estim = c("optim", "nloptr", "none"),
                   coefLower = c("scale" = 0.0, "shape" = -0.90),
                   coefUpper = c("scale" = Inf, "shape" = Inf),
                   trace = 0) {

    eps <- 1e-10
    estim <-  match.arg(estim)

    if (is(data, "Rendata")) {
        stop("Using 'data' with class \"Rendata\" is not allled yet")
    }

    ## Note that we do not pass 'threshold' here.
    data <- checkPoisGPData(data = data,
                            effDuration = effDuration,
                            MAX.data = MAX.data,
                            MAX.effDuration = MAX.effDuration,
                            OTS.data = OTS.data,
                            OTS.threshold = OTS.threshold,
                            OTS.effDuration = OTS.effDuration)

    grandMin <- min(unlist(sapply(data, function(x) unlist(x$data))))

    if (missing(threshold)) {
        if (data$OT$flag || data$OTS$flag) {
            warning("'threshold' is missing and set just below the ",
                    "smallest obsevation")
        }
        threshold <- gradMin - eps
        fitData <- data
    } else {
        if (threshold < grandMin) {
            warning("'threshold' is smaller than the smallest observation")
            fitData <- data
        } else {
            fitData <- threshData(threshold = threshold, data, scale = FALSE)
        }
    }
        
    thisPoisGP <- list(call = match.call(),
                       threshold = threshold,
                       data = data,
                       fitData = fitData,
                       p = 3L,
                       df = 3L,
                       parNames = c("lambda", "scale", "shape"))

    class(thisPoisGP) <- "poisGP"
    
    
    if (estim != "none") {

        if (is.null(parIni)) {
            parIni <- coefIni(thisPoisGP)[-1]
        }
        if (trace) {
            cat("Initial values for the parameter\n")
            print(parIni)
        }
        
        res <- MLE.poisGP(object = thisPoisGP,
                          parIni = parIni,
                          estim = estim,
                          coefLower = coefLower,
                          coefUpper = coefUpper,
                          trace = trace)
        
        for (nm1 in names(res)) {
            thisPoisGP[[nm1]] <- res[[nm1]]
        }
        
    } else {
        thisPoisGP$estimate <- c("lambda" = NA, "scale" = NA, "shape" = NA)
        thisPoisGP$negLogLik <- NA
        thisPoisGP$logLik <- NA
    }

    class(thisPoisGP) <- "poisGP"
    thisPoisGP
    
}
