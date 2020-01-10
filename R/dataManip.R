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
##' @title Select the Observations Over a Threshold within
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
    out <- list()
    
    grandMin <- min(unlist(sapply(data, function(x) unlist(x$data))))
    grandMax <- max(unlist(sapply(data, function(x) unlist(x$data))))

    if (length(threshold) == 0) {
        warning("'threshold' is missing and set just below the ",
                "smallest obsevation")
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
        out$OT$data <- (data$OT$data[data$OT$data > threshold] - threshold)
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
            out$OTS$threshold[b] <-
                pmax(threshold, data$OTS$threshold[[b]]) - out$threshold
            out$OTS$data[[b]] <-
                (data$OTS$data[[b]][data$OTS$data[[b]] > threshold] -
                     out$threshold) 
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
            dat <- dat[ind] - threshold
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
    out
    
}


                
