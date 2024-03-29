## ****************************************************************************
##' Select the observations over the given threshold within
##' heterogeneous data provided as a \code{potData} object. The data
##' possibly contain OT, MAX and OTS blocks. The data can optionally
##' be scaled using a scale that is attached to the result as an
##' attribute. This function is rather technical and should normally
##' to be needed by the user.
##'
##' When \code{scale} is \code{TRUE}, a suitable scale (a positive
##' number) is chosen as a power of \code{10} and is used to divide
##' the exceedances over \code{threshold}. This can in some cases
##' avoid numerical problems.
##' 
##' @title Select the Observations Over a Threshold within
##' Heterogeneous Data provided as a \code{potData} object
##'
##' @param threshold The "main" threshold used to select the
##' observations in each block.
##'
##' @param data A list with elements "OT", "MAX" and "OTS". Each
##' sublist contains a \code{flag} logical element a \code{data}
##' vector or list and a numeric \code{duration}.
##'
##' @param exceed Logical. If \code{TRUE} the data returned contains
##' the exceedances over the threshold, i.e. the threshold is
##' substracted from the data.
##' 
##' @param scale Logical. If \code{TRUE} the excesses over the
##' threshold will all be scaled by dividing them by a common "round"
##' positive number. This is intended to avoid numerical problems
##' during optimisation. Note that the round scaling number (a power of
##' 10) is always computed and returned as the \code{"scale"}
##' attribute of the result. Even though the data are not scaled this
##' number can be used in the \code{optim} function to set the
##' \code{parScale} element of the \code{control} list.
##'
##' @param warn Logical. If \code{TRUE} the function warns about a
##' threshold which is smaller than all the observation in 'data'.
##' 
##' @return A list which is comparable to \code{data} but with the
##' observations below the threshold removed, and the related
##' information changed. For instance if \code{threshold} is greater
##' that some observations in a \code{"MAX"} block, these are
##' discarded and the number \code{r} is changed accordingly.
##'
##' @note The 'main' threshold can exceed the threshold of some OTS
##' blocks and it can also exceed some observations in a MAX block of
##' \code{data}. In the later case the MAX block will be turned into
##' an OTS block with its threshold set to the main threshold; it can
##' then non longer have any observations if the main threshold
##' exceeds all the observations of the original MAX block.
##' 
##' @seealso \code{\link{potData}}.
##'
##' @export
##' 
##' @examples
##' set.seed(123)
##' myData <-
##'     potData(data = rexp(50), effDuration = 50,
##'             MAX.effDuration = c(25, 10),
##'             MAX.data = list("MAX1" = tail(sort(rexp(20))),
##'                             "MAX2" = tail(sort(rexp(12)))),
##'             OTS.effDuration = c(30, 65, 200),
##'             OTS.threshold = c(1.8, 3.2, 5.0),
##'             OTS.data = list("OTS1" = 1.8 + rexp(5),
##'                             "OTS2" = 3.2 + rexp(3),
##'                             "OTS3" = 5.0 + rexp(2)))
##' 
##' myData2 <- threshData(threshold = 4, data = myData)
##' autoplot(myData)
threshData <- function(threshold = NULL, data, exceed = FALSE, scale = FALSE,
                       warn = FALSE) {
    
    eps <- 1e-4
    out <- list()

    grandMin <- min(unlist(sapply(data, function(x) unlist(x$data))))
    grandMax <- max(unlist(sapply(data, function(x) unlist(x$data))))

    if (length(threshold) == 0) {
        if (warn) {
            warning("'threshold' is missing and set just below the ",
                    "smallest obsevation")
        }
        out$threshold <- threshold <- grandMin - eps
    } else {
        if (threshold < grandMin) {
            warning("'threshold' is smaller than the smallest observation")
        }
        out$threshold <- threshold
    }

    sc <- 10^floor(log(grandMax - grandMin, 10))

    out$OT <- data$OT

    if (data$OT$flag) {
        ind <- data$OT$data > threshold
        out$OT$data <- data$OT$data[ind]
        if (exceed)  out$OT$data <- out$OT$data - threshold
        
        if (scale) out$OT$data <- out$OT$data / sc
        out$OT$n <- length(out$OT$data)
    } 
  

    ## =====================================================================
    ## Mind that for 'MAX' data there is no threshold, so this must be
    ## set on a by-block basis using the smallest observation.
    ## =====================================================================
        
    out$OTS <- data$OTS
    
    if (data$OTS$flag) {
        for (b in seq_along(data$OTS$data)) {
            out$OTS$threshold[b] <- pmax(threshold, data$OTS$threshold[[b]])
            if (exceed) {
                out$OTS$threshold[b] <- out$OTS$threshold[b] - out$threshold
            }
            ind <- data$OTS$data[[b]] > threshold
            out$OTS$data[[b]] <- data$OTS$data[[b]][ind]
            if (exceed) {
                out$OTS$data[[b]] <- out$OTS$data[[b]] - out$threshold
            }
            if (scale) out$MAX$data[[b]] <- out$MAX$data[[b]] / sc
            out$OTS$r[b] <- length(out$OTS$data[[b]])
        }
        names(out$OTS$data) <- names(data$OTS$data)
    }
    
    out$MAX <- list(flag = FALSE)

    n.OTS <- n.OTS0 <- length(out$OTS$data)
    n.MAX <- 0L
    
    nms <- names(data$MAX$data)
    nmsFlag <- !is.null(nms)
    nms.OTS <- character(0)
    nms.MAX <- character(0)
    
    if (data$MAX$flag) {
        
        for (b in seq_along(data$MAX$data)) {

            ## ================================================================
            ## When the min of the MAX block is < threshold, the block
            ## must be turned into an OTS block
            ## ================================================================

            dur <- data$MAX$effDuration[b]
            dat <- data$MAX$data[[b]]
            ind <- (dat > threshold)
            dat <- dat[ind]

            if (exceed) dat <- dat - threshold
            if (scale) dat <- dat / sc

            ## cat(sprintf("b = %d, threshold = %5.3f\n", b, threshold))
            ## print(ind)
                
            if (!all(ind)) {

                n.OTS <- n.OTS + 1L

                if (!out$OTS$flag) {
                    out$OTS$flag <- TRUE
                    out$OTS$data <- list(dat)
                    out$OTS$r <- length(dat)
                    out$OTS$effDuration <- dur
                    out$OTS$threshold <- 1e-8
                } else {
                    out$OTS$data[[n.OTS]] <- dat
                    out$OTS$r <- c(out$OTS$r, length(dat))
                    out$OTS$effDuration <- c(out$OTS$effDuration, dur)
                    out$OTS$threshold <- c(out$OTS$threshold, 1e-8)
                 }
                
                if (nmsFlag) nms.OTS <- c(nms.OTS, nms[b]) 

            } else {
                
                n.MAX <- n.MAX + 1L
                
                if (!out$MAX$flag) {
                    out$MAX$flag <- TRUE
                    out$MAX$data <- list(dat)
                    out$MAX$r <- length(dat)
                    out$MAX$effDuration <- dur
                    out$MAX$threshold <- min(dat)
                } else {
                    out$MAX$data[[n.MAX]] <- dat
                    out$MAX$r <- c(out$MAX$r, length(dat))
                    out$MAX$effDuration <- c(out$MAX$effDuration, dur) 
                    out$MAX$threshold <- c(out$MAX$threshold, min(dat))
                }
                
                if (nmsFlag) nms.MAX <- c(nms.MAX, nms[b])
                
            }
        }

        
        if (nmsFlag) {
            if (out$MAX$flag) {
                names(out$MAX$data) <- nms.MAX
            }
            if (out$OTS$flag && (!data$OTS$flag || !is.null(names(data$OTS$data)))) {
                names(out$OTS$data)[(n.OTS0 + 1):n.OTS] <- paste("[MAX2OTS]", nms.OTS)
            }
        }
        
    }
    
    attr(out, "scale") <- sc

    if (!exceed) {
        out$threshold <- NULL
        class(out) <- "potData"
    }
    
    out
    
}


                
