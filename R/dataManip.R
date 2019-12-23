## ****************************************************************************
##' Check the provided data to define a Poisson-GP model.
##'
##' As opposed to what is done in \strong{Renext}, no 'main' threshold
##' is used here, because the data need to be checked \emph{before}
##' being censured using the threshold. So the MAX or OTS blocks can
##' later be used with any 'main' threshold, even if they contain
##' observations and thresholds that are smaller than the main
##' threshold.
##'
##' @title Check the data Provided to Define a Poisson-GP Model.
##'
##' @param data A numeric vector containing the observations for the
##' main sample. If \code{NULL}, the main sample is assumed to be
##' absent.
##'
##' @param effDuration  Duration of the main sample 
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
##' checkPoisGPData(data = Garonne$OTdata$Flow,
##'                 effDuration = 65,
##'                 MAX.data = Garonne$MAXdata$Flow,
##'                 MAX.effDuration = 143)
##' 
checkPoisGPData <- function(data = NULL, effDuration = NULL,
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
                        r = 1L,
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
                        r = 1L,
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
    
    list(OT = OT,
         MAX = MAX,
         OTS = OTS)
    
}

## ****************************************************************************
##' Select the observations over the given threshold within
##' heterogeneous data. The data possibly contain OT, MAX and OTS
##' blocks. The data can optionally be scaled using a scale that is
##' attached to the result as an attribute.
##'
##' When \code{scale} is \code{TRUE}, a suitable scale (a positive
##' number) is chosen as a power of \code{10} and is used to divide
##' the exceedances over \code{threshold}. This can in some cases
##' avoid numerical problems.
##' 
##' @title Select the Observations Over a Threshod within
##' Heterogeneous Data.
##'
##' @param threshold A "main" threshold used to select the
##' observations in each block.
##'
##' @param data A list with elements "OT", "MAX" and "OTS". Each
##' sublist contains a \code{flag} logical element a \code{data}
##' vector or list and a numeric \code{duration}.
##'
##' @return A list which is comparable to \code{data} but with the
##' observations below the threshold removed, and the related
##' information changed. For instance if \code{threshold} is greater
##' that some observations in a \code{"MAX"} block, these are
##' discarded and the number \code{r} is changed accordingly.
##'
##' @seealso \code{link{check}}
##' 
##' @examples
##' set.seed(123)
##' myData <-
##'     checkPoisGPData(data = rexp(50), duration = 50,
##'                     MAX.effDuration = c(25, 10),
##'                     MAX.data = list("MAX1" = tail(sort(rexp(20))),
##'                                     "MAX2" = tail(sort(rexp(12))),
##'                     OTS.effDuration = c(30, 65, 200),
##'                     OTS.threshold = c(1.8, 3.2, 5.0))
##'                     OTS.data = list("OTS1" = 1.8 + rexp(5),
##'                                     "OTS2" = 3.2 + rexp(3),
##'                                     "OTS3" = 5.0 + rexp(2)))
##' myData2 <- threshData(threshold = 4, data = myData, scale = FALSE)
##' 
threshData <- function(threshold = NULL, data, scale = FALSE) {

    eps <- 1e-4
    out <- data
    
    grandMin <- min(unlist(sapply(data, function(x) unlist(x$data))))
    grandMax <- max(unlist(sapply(data, function(x) unlist(x$data))))

    if (length(threshold) == 0) {
        warning("'threshold' is missing and set just below the ",
                "smallest obsevation")
        out$threshold <- gradMin - eps
    } else {
        if (threshold < grandMin) {
            warning("'threshold' is smaller than the smallest observation")
        }
        out$threshold <- threshold
    }

    sc <- 10^ceiling(log(grandMax - grandMin, 10))
    
    if (data$OT$flag) {
        out$OT$data <- (data$OT$data[data$OT$data > threshold] - threshold)
        if (scale) out$OT$data <- out$OT$data / sc
        out$OT$n <- length(out$OT$data)
    }

    ## =====================================================================
    ## Mind that for 'MAX' data there is no threshold, so this must be
    ## set on a by-block basis using the smallest observation.
    ## =====================================================================
    
    if (data$MAX$flag) {
        for (b in seq_along(data$MAX$data)) {
            out$MAX$threshold[b] <-
                pmax(threshold, min(data$MAX$data[[b]] - eps)) - out$threshold
            out$MAX$data[[b]] <-
                (data$MAX$data[[b]][data$MAX$data[[b]] > threshold] -
                     out$threshold)
            if (scale) out$MAX$data[[b]] <- out$MAX$data[[b]] /sc
            out$MAX$r[b] <- length(out$MAX$data[[b]])
        }
    }
    
    if (data$OTS$flag) {
        for (b in seq_along(data$OTS$data)) {
            out$OTS$threshold[b] <-
                pmax(threshold, data$OTS$threshold[[b]]) - out$threshold
            out$OTS$data[[b]] <-
                (data$OTS$data[[b]][data$OTS$data[[b]] > threshold] -
                     out$threshold) 
            if (scale) out$MAX$data[[b]] <- out$MAX$data[[b]] /sc
            out$OTS$r[b] <- length(out$OTS$data[[b]])
        }
    }
    
    if (scale) attr(out, "scale") <- sc
    out
    
}


                
