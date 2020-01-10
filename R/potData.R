## ****************************************************************************
##' Check and structure data to define a Poisson-GP model.
##'
##' As opposed to what is done in \strong{Renext}, no 'main' threshold
##' is used here, because the data need to be checked \emph{before}
##' being censured using the threshold. So the MAX or OTS blocks can
##' later be used with any 'main' threshold, even if they contain
##' observations and thresholds that are smaller than the main
##' threshold.
##'
##' @title Check and Structure the data Provided to Define a
##' Poisson-GP Model.
##'
##' @param data A numeric vector containing the observations for the
##' main sample. If \code{NULL}, the main sample is assumed to be
##' absent.
##'
##' @param effDuration  Duration of the main sample. 
##'
##' @param MAX.data A numeric vector or list of numeric vectors
##' containing the observations for the \code{MAX} blocks.
##'
##' @param MAX.effDuration A numeric vector containing the durations
##' for the MAX blocks.
##'
##' @param OTS.data  A numeric vector or list of numeric vectors
##' containing the observations for the \code{OTS} blocks.
##'
##' @param OTS.threshold A numeric vector containing the thresholds
##' for the OTS blocks.
##'
##' @param OTS.effDuration A numeric vector containing the durations for
##' the OTS blocks.
##' 
##' @return A list with the element given on input checked and
##' suitably named.
##'
##' @examples
##' 
##' potData(data = Garonne$OTdata$Flow,
##'         effDuration = 65,
##'         MAX.data = Garonne$MAXdata$Flow,
##'         MAX.effDuration = 143)
##' 
potData <- function(data = NULL, effDuration = NULL,
                    MAX.data = NULL,
                    MAX.effDuration = NULL,
                    OTS.data = NULL,
                    OTS.threshold = NULL,
                    OTS.effDuration = NULL) {
    
    if (!is.null(data)) {
        data <- as.numeric(data)
        if ((length(effDuration) != 1L) || (effDuration <= 0.0)) {
            stop("When 'data' is provided, 'effDuration' must be a ",
                 "positive number")
        }
        OT <- list(flag = TRUE,
                   effDuration = effDuration,
                   n = length(data),
                   data = data)
    } else {
        OT <- list(flag = FALSE,
                   effDuration = NA,
                   n = NA,
                   data = NULL)
    }
    
    if (!is.null(MAX.data)) {
        if (is.numeric(MAX.data)) {
            if (length(MAX.data) == 0) {
                stop("A block MAX must contain at least an observation")
            }
            if ((length(MAX.effDuration) != 1L) || (MAX.effDuration <= 0.0)) {
                stop("When 'MAX.data' is numeric, 'MAX.effDuration' must be a ",
                     "positive number")
            }
            MAX <- list(flag = TRUE,
                        effDuration = MAX.effDuration,
                        r = length(MAX.data),
                        data = list("MAX block#1" = MAX.data))
        } else {
            if (!is.list(MAX.data) || !all(sapply(MAX.data, is.numeric))) {
                stop("'MAX.data' must be a numeric vector or a list of numeric ",
                     "vectors, one for each MAX block")
            }
            if (length(MAX.effDuration) != length(MAX.data) ||
                any(MAX.effDuration <= 0.0)) {
                stop("When 'MAX.data' is a list, 'MAX.effDuration' must be a ",
                     "numeric vector with length length(MAX.data)) and ",
                     "positive elements")
            }
            MAX.r <- sapply(MAX.data, length)
            if (any(MAX.r == 0L)) {
                stop("A MAX block can not be empty")
            }
            MAX <- list(flag = TRUE,
                        effDuration = MAX.effDuration,
                        r = MAX.r,
                        data = MAX.data)
        }            
    } else {
        MAX <- list(flag = FALSE,
                    effDuration = NA,
                    r = NA,
                    data = NULL)
        
    }
    
    if (!is.null(OTS.data)) {
        
        if (is.numeric(OTS.data)) {

            if (length(OTS.data) == 0) {
                stop("A block OTS must contain at least an observation")
            }
            if ((length(OTS.effDuration) != 1L) || (OTS.effDuration <= 0.0)) {
                stop("When 'OTS.data' is numeric, 'OTS.effDuration' must be a ",
                     "positive number")
            }
            if ((length(OTS.threshold) != 1L) ||
                OTS.threshold >= min(OTS.data)) {
                stop("When 'OTS.data' is numeric, 'OTS.threshold' must be a ",
                     "number < min(OTS.data)")
            }
           
            OTS <- list(flag = TRUE,
                        effDuration = OTS.effDuration,
                        r = length(MAX.data),
                        data = list("OTS block#1" = OTS.data))
        } else {
            if (!is.list(OTS.data) || !all(sapply(OTS.data, is.numeric))) {
                stop("'OTS.data' must be a numeric vector or a list of numeric ",
                     "vectors, one for each OTS block")
            }
            if (length(OTS.effDuration) != length(OTS.data) ||
                any(OTS.effDuration <= 0.0)) {
                stop("When 'OTS.data' is a list, 'OTS.effDuration' must be a ",
                     "numeric vector with length length(OTS.data)) and positive ",
                     "elements")
            }
            OTS.r <- sapply(OTS.data, length)
            OTS.min <- sapply(OTS.data, min)
            if (any(OTS.threshold >= OTS.min)) {
                stop("An OTS block must have its observations greater than the ",
                     "corresponding threshold")
            }
            
            OTS <- list(flag = TRUE,
                        effDuration = OTS.effDuration,
                        r = OTS.r,
                        data = OTS.data)
        }
        
    } else {
        OTS <- list(flag = FALSE,
                    effDuration = NA,
                    r = NA,
                    threshold = OTS.threshold,
                    data = NULL)
        
    }
    
    res <-  list(OT = OT, MAX = MAX, OTS = OTS)

    class(res) <- "potData"
    res
}

## ****************************************************************************
##' Return Periods for a \code{potData} object.
##' 
##' @title Return Periods for a \code{potData} object.
##'
##' @param object A \code{potData} object.
##' 
##' @param points Type of plotting positions to use: \code{ppoints} or
##' code Nelson's positions.
##'
##' @param a Passed to \code{link{ppoints}} when the value
##' \code{points} is \code{"p"}. Ignored else.
##'
##' @param defNames Not used.
##'
##' @param ... 
##'
##' @return A data frame.
##'
##' @examples
##' 
##' pdat <- potData(data = Garonne$OTdata$Flow,
##'                 effDuration = 65,
##'                 MAX.data = Garonne$MAXdata$Flow,
##'                 MAX.effDuration = 143)
##' L <- RP(pdat)
##' 
RP.potData <- function(object,
                       points = c("p", "H"),
                       a = 0.5,
                       defNames = FALSE,
                       ...) {

    points <- match.arg(points)
    
    if (!is.null(object$OT) && length(object$OT$data) > 0) {
        vals <- object$OT$data
        samp <- rep(1L, length(object$OT$data))
        nms <- names(object$OT$data)
        if (is.null(nms)) {
            nms <- rep("OT", length(object$OT$data))
        }
        names(samp) <- names(vals) <- nms
        mins <- min(object$OT$data)
        durs <- object$OT$effDuration
        groupNames <- "OT"
        sourceNames <- "OT"
        names(mins) <- names(durs) <- groupNames[1L]
        nBlock <- 1L
    }  else {
        vals <- numeric(0)
        samp <- numeric(0)
        mins <- numeric(0)
        durs <- numeric(0)
        groupNames <- character(0)
        sourceNames <- character(0)
        nBlock <- 0L
    }
    
    ## manage in the same way 'OTS' periods
    if (!is.null(object$OTS) && object$OTS$flag){
        r <- object$OTS$r
        nms <- names(object$OTS$data)
        if (is.null(nms)) {
            nms <- paste("OTS block", 1L:length(r), sep = "")
        } 
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("OTS", length(r)))
        newDurs <- object$OTS$effDuration
        newMins <- object$OTS$threshold
        names(newDurs) <- names(newMins) <- nms
        durs <- c(durs, newDurs)
        mins <- c(mins, newMins)
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(object$OTS$data)
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        nBlock <- nBlock + length(r)
    }
    
    ## manage 'MAX' periods
    if (!is.null(object$MAX) && object$MAX$flag){
        r <- object$MAX$r
        nms <- names(object$MAX$data)
        if (is.null(nms)) {
            nms <- paste("MAX block", 1L:length(r), sep = "")
        }
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("MAX", length(r)))
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(object$MAX$data)
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        nBlock <- nBlock + length(r)
        newDurs <- object$MAX$effDuration
        d <- diff(sort(vals))
        smallVal<- min(d[d > 0]) * 0.5
        newMins <- unlist(lapply(object$MAX$data, min)) - smallVal
        names(newDurs) <- names(newMins) <- nms
        durs <- c(durs, newDurs)
        mins <- c(mins, newMins)
    }
    
    ## order the thresholds
    ind <- order(mins)
    mins <- mins[ind]
    durs <- durs[ind]

    ## ========================================================================
    ## use more algorithmic notations as usesd in the 'Renext
    ## Computind Details' document, and remove ex aequo if necessary
    ## ========================================================================
    
    u <- mins
    W <- cumsum(durs)
    J <- length(u)

    ## for objects such as those created with 'RenouvNoEst'
    if (J == 0) {
        return(list(x = numeric(0),
                    group = numeric(0),
                    groupNames = character(0),
                    sourceNames = character(0),
                    lambda = NA,
                    S = numeric(0),
                    T = numeric(0),
                    thresh = c(object$threshold, Inf),
                    lambda.thresh = numeric(0),
                    S.thresh = numeric(0),
                    T.thresh = numeric(0)))
    }
    
    ind <- order(vals)
    vals <- vals[ind]
    samp <- samp[ind]
    
    ##=========================================================================
    ## First pass: compute 'A', 'D', 'rate' and the inverse return
    ## period AT THE THRESHOLDS. Store the indices of 'vals' falling
    ## between the thresholds.
    ## ========================================================================
    
    interv <- findInterval(x = vals, vec = u, rightmost.closed = TRUE)
    A <- D <- invT.thresh <- lambda.thresh <- rep(0, J + 1L)
    ind <- list()
    for (j in J:1L) {
        ind[[j]] <- (interv == j)
        A[j] <- sum(ind[[j]])
        lambda.thresh[j] <- A[j] / W[j] 
        D[j] <- D[j + 1L] + A[j] + lambda.thresh[j] * (W[J] - W[j])
        invT.thresh[j] <- D[j] / W[J] 
    }
    
    lambda <- sum(lambda.thresh[1:J])
    ## cat("lambda = ", lambda, "\n")
    
    ##===========================================================================
    ## Second pass: compute the inverse retun period at the data by
    ## interpolation, as well as the survival.  The 'Hpoints' options
    ## can be used only (at the time) for the upper slice of values.
    ## ===========================================================================
    
    T <- rep(NA, length(vals))
    S <- rep(NA, length(vals))
    Aprec <- FALSE
    
    for (j in J:1L) {
        if (A[j]) {
            ## reverse order
            prov <- invT.thresh[j + 1L] +
                (invT.thresh[j] - invT.thresh[j + 1L]) * ((A[j]:1L) - a) / (A[j] -2*a + 1)
            S[ind[[j]]] <- prov / lambda
            if (!Aprec && points == "H") {
                Sprov <- invT.thresh[j] / lambda
                Hprov <- -log(Sprov)
                Hprov <- Hprov + Hpoints(A[j])
                T[ind[[j]]]  <- exp(Hprov) / lambda
            } else {
                ## cat(sprintf("j = %d\n", j))
                T[ind[[j]]]  <- 1 / S[ind[[j]]] / lambda
            }
            Aprec <- TRUE
        }
    }

    if (TRUE) {

        lSourceNames <- factor(sourceNames, levels = c("OT", "MAX", "OTS"))
        ind <- order(lSourceNames, groupNames)
        
        fSourceNames <- factor(sourceNames[samp], levels = c("OT", "MAX", "OTS"))
        fGroupNames <- factor(groupNames[samp], levels = groupNames[ind])
        
        df <- data.frame(block = fGroupNames,          ## groupNames[samp],
                         source = fSourceNames,        ## sourceNames[samp],
                         S = S,
                         T = T,
                         x = vals)

        
    }
    
    list(data = df,
         x = vals,
         group = samp,
         groupNames = groupNames,
         sourceNames = sourceNames,
         lambda = lambda,
         S = S,
         T = T,
         thresh = c(u, Inf),
         lambda.thresh = lambda.thresh,
         S.thresh = invT.thresh / lambda,
         T.thresh = 1 / invT.thresh)
    
}


