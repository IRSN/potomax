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
##' @noRd
##' 
##' @title Initial Estimates for a \code{poisGP} Model
##'
##' @param object A \code{poisGP} object the parameters of which need
##' to be estimated.
##'
##' @param trace Integer level of verbosity.
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

## ****************************************************************************
##' Log-likelihood of a Poisson-GP model object. 
##'
##' Although models for block maxima and \eqn{r}-largest models can be
##' considered as special cases of a Poisson-GP model, the
##' log-likelihood as computed in classical Extreme-Value packages can
##' differ by a constant related to the observation scheme. The
##' argument \code{type} aims at finding the same results in the two
##' cases where the results can be obtained with classical packages,
##' namely: \emph{"POT" fits with Poisson-GP representation} and
##' \emph{block maxima or r-largest}. In the first case, this is
##' hopefully attained by choosing \code{type} to be \code{"poisGP"}
##' and in the second case by choosing \code{type} to be \code{"PP"}.
##'
##' @method logLik poisGP
##'
##' @usage
##' \method{logLik}{poisGP}(object, type = c("poisGP", "PP"), ...) 
##' 
##' @title Log-Likelihood of a Poisson-GP Object
##'
##' @param object An object with class \code{"poisGP"}.
##'
##' @param type An experimental argument to make the log-likelihood
##' comparable to that of other packages/functions. See \strong{Details}.
##'
##' @param ... Not used yet.
##' 
##' @return Value of the log-likelihood.
##'
##' @section Caution: Comparing log-likelihoods or related indicators
##' such AIC, BIC across R packages can be irrelevant due to the use
##' of different constants.
##' 
logLik.poisGP <- function(object,
                          type = c("poisGP", "PP"),
                          ...) {
    type <- match.arg(type)
    
    if (type == "poisGP") {
        res <- object$logLik - object$Cst
    } else {
        res <-  object$logLik
    }

    attr(res, "df") <- object$df
    attr(res, "nobs") <- object$nobs
    res
    
}

## ****************************************************************************
##' Akaike's Information Criterion and Schwarz's Bayesian Information
##' Criterion for a Poisson-GP model object.
##' 
##' @method AIC poisGP
##' @aliases BIC.poisGP
##'  
##' @usage
##' \method{AIC}{poisGP}(object, ..., k)
##' 
##' \method{BIC}{poisGP}(object, ...)
##' 
##' @title Akaike's Information Criterion and Bayesian Information
##' Criterion for a Poisson-GP Object
##'
##' @param object An object with class \code{"poisGP"}.
##'
##' @param ... Not used yet.
##'
##' @param k See \code{\link[stats]{AIC}}.
##' 
##' @return Value of the criterion
##'
##' @note For technical reasons these methods do not have a
##' \code{type} argument as does \code{\link{logLik.poisGP}} and
##' consequently for a model fitted from block maxima or
##' \eqn{r}-largest the computed criteria will differ from those
##' computed by other packages because the log-likelihoods differ by a
##' constant.
##' 
##' @section Caution: Comparing log-likelihoods or related indicators
##' such AIC, BIC across R packages can be irrelevant due to the use
##' of different constants. Moreover the concept of \emph{number of
##' observations} is unclear when heterogeneous data are used. An
##' historical MAX or OTS information can have a very strong influence
##' on the estimation hence can not be compared to an ordinary OT
##' observation.
##'
##' @seealso \code{\link{logLik.poisGP}}
##' 
AIC.poisGP <- function(object, ..., k = 2) {
    return(NextMethod())
}

BIC.poisGP <- function(object, ...) {
    return(NextMethod())
}


## ****************************************************************************
##' Extract model coefficients for a Poisson-GP object.
##' 
##' When \code{type} is \code{"poisGP"} the three parameters are in
##' that order: \code{"lambda"} for the Poisson rate \eqn{\lambda},
##' \code{"scale"} for the GP scale \eqn{\sigma} and \code{"shape"}
##' for the GP shape \eqn{\xi}.
##'
##' When \code{type} is \code{"PP"} the \emph{Point Process}
##' parameterisation is used, in relation with a block duration. The
##' three parameters are in that order: \code{"loc"} for the location
##' \eqn{\mu}, \code{"scale"} for the scale \eqn{\sigma} and
##' \code{"shape"} for the shape \eqn{\xi}. They define a Generalised
##' Extreme Value (GEV) distribution the tail of which applies to the
##' maximum of the marks on a time block with the given duration. Note
##' that the \code{"PP"} scale differs from the \code{"poisGP"} scale
##' and that the location is not the threshold, and is smaller than
##' it.
##'
##' @title Extract Model Coefficients for a Poisson-GP Object
##'
##' @param object A model with class \code{"poisGP"}.
##'
##' @param type The parameterisation to use. With \code{"poisGP"} the
##' standard Poisson-GP parameters are used
##'
##' \code{"lambda"} for the
##' Poisson rate \eqn{\lambda}, \code{"scale"}
##' 
##' @param ... Not used yet.
##'
##' @return A numeric vector of three coefficients.
##'
##' @note The \code{"poisGP"} parameters can be transformed into
##' \code{"PP"} parameters by using \code{\link[Renext]{Ren2gev}}
##' function of the \strong{Renext} package. This requires giving both
##' the threshold and the block duration \code{w}.
##'
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
##' @note The correlation matrix can be obtained by using
##' \code{\link{cov2cor}} on the result.
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
##' The functions and methods used to maximise the likelihood are as follows.
##' \itemize{
##' \item{\code{estim = "optim"}}{
##' The classical \code{stats::optim} function is used with
##' \code{method ="BFGS"}. The derivatives are not used and nor do the
##' bounds on the parameters given in \code{coefLower} and
##' \code{coefUpper}.
##' }
##' \item{\code{estim = "nloptr"}}{
##' The \code{nloptr::nloptr} function is used with the
##' \code{"NLOPT_LD_BFGS"} algorithm option. The derivatives are used
##' as well as the bounds on the parameters, leading to "box
##' constraints". The bounds can be used to fix the value of a GP
##' parameter by using the same value in \code{coefLower} and
##' \code{coefUpper}. For instance an exponential distribution can be
##' fitted by using a zero value for the shape both in
##' \code{coefLower} and \code{coefUpper}.
##' }
##' \item{\code{estim = "eval"}}{
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
##' }
##' \item{\code{estim = "none"}}{
##' No optimisation is performed. The log-likelihood and negative
##' log-likelihoods remain NA, and the initial values are ignored.
##' 
##' }}
##' 
##' @title Create a Poisson-GP Model Object
##'
##' @usage
##' poisGP(data = NULL, threshold, effDuration,
##'        MAX.data = NULL, MAX.effDuration = NULL,
##'        OTS.data = NULL, OTS.threshold = NULL, OTS.effDuration = NULL,
##'        parIni = NULL,
##'        estim = c("optim", "nloptr", "eval", "none"),
##'        coefLower = c("scale" = 0.0, "shape" = -0.90),
##'        coefUpper = c("scale" = Inf, "shape" = Inf),
##'        scale = FALSE,
##'        trace = 0)
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
##' \code{OTS.threshold}. So the number \eqn{r \geq 0} vary across
##' blocks. When a numeric vector is passed instead of a list, it is
##' understood that there is only one OTS block.
##'
##' @param OTS.threshold A vector of positive numbers giving the
##' thresholds of the OTS blocks (in the same order). By construction
##' all the elements of \code{OTS.data[i]} are larger than
##' \code{OTS.threshold[i]} for each block index \code{i}.
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
##' @param estim \code{Character} defining the function and the method
##' that will be used to maximise the likelhood. See \strong{Details}.
##'
##' @param coefLower,coefUpper Named vectors of bounds for the
##' parameters. The values can be infinite \code{Inf} or \code{-Inf}.
##' However, note that without bounds on the shape parameter \eqn{\xi}
##' the maximum likelihood is infinite and obtained for \eqn{\xi =
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
##' ## =====================================================================
##' ## Use Garonne data from Renext
##' ## =====================================================================
##'
##' fit1p <- poisGP(data = Garonne$OTdata$Flow, threshold = 2900,
##'               effDuration = 65,
##'               MAX.data = Garonne$MAXdata$Flow,
##'               MAX.effDuration = 143)
##' fit1R <- Renouv(Garonne, threshold = 2900, distname.y = "GPD",
##'                 plot = FALSE)
##'
##' cbind(Renext = coef(fit1R), potomax = coef(fit1p))
##' 
##' ## CAUTION when comparing log-likelihoods, see ?logLik.poisGP
##' cbind(Renext = logLik(fit1R), potomax = logLik(fit1p))
##' 
##' ## ==============================================================
##' ## Use the 'venice' data from the 'ismev' package. Contains
##' ## r-largest observations as a matrix with one row by year and NA.
##' ## So some transformations are needed. Note that first
##' ## of 'venice' must be removed, and that the 'venice' data from
##' ## the evd package may misleadingly be used instead.
##' ## ==============================================================
##'
##' rm(venice)
##' data(venice, package = "ismev")
##' MAX.data <- as.list(as.data.frame(t(venice[ , -1])))
##' MAX.data <- lapply(MAX.data, function(x) x[!is.na(x)])
##' MAX.effDuration <- rep(1, length(MAX.data))
##'
##' fit2i <- ismev::rlarg.fit(venice[ , -1])
##' fit2R <- Renext::fGEV.MAX(MAX.data = MAX.data,
##'                           MAX.effDuration = MAX.effDuration)
##' fit2p <- poisGP(MAX.data = MAX.data,
##'                 MAX.effDuration = MAX.effDuration)
##'
##' ## To compare the coefficients we must use the "PP" coefficients
##' ## of the poisGP object rather than the standard "poisGP"
##' ## coefficients.
##' 
##' cbind("ismev" = fit2i$mle,
##'       "Renext" = fit2R$estimate,
##'       "potomax" = coef(fit2p, type = "PP"))
##'
##' ## CAUTION when comparing log-likelihoods, see ?logLik.poisGP
##' ## We choose here the "PP" type which usually makes the result
##' ## comparable to those based on block maxima or on r-largest.
##' 
##' cbind("ismev" = -fit2i$nllh,
##'        "Renext" = fit2R$loglik,
##'        "potomax" = logLik(fit2p, type = "PP"))
##'
##' ## profile-likelihood confidence intervals on parameters
##' confint(fit2p)
##' 
##' autoplot(fit2p) + ggtitle("Venice r-largest")
##'
##' ## Now censor the MAX data. This can not be done with the
##' ## other packages 
##' fit3p <- poisGP(MAX.data = MAX.data,
##'                 MAX.effDuration = MAX.effDuration,
##'                 threshold = 100)
##' coef(fit3p)
##' autoplot(fit3p) + ggtitle("Venice r-largest with threshold 100 cm")
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

    ## XXXX
    if (scale) {
        stop("'scale = TRUE' is not fully implemented for now. ",
             "This needs work to rescale the estimated parameters ",
             "but also to change the log-likelihood values.")
    }


    
    if (is(data, "Rendata")) {
        stop("Using 'data' with class \"Rendata\" is not allowed yet")
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
            cat("\nThe data will be scaled by dividing then by ",
                scaleData, "\n")
        } else {
            cat("\nThe data will not be scaled\n")
        }
    }

    ## ========================================================================
    ## Compute a possible number of observations. Since in the general
    ## case the observations can have a very different impact on the
    ## estimation, the definition of 'nobs' is unclear here. Maybe we
    ## could consider that an empty OTS block is worth an observation,
    ## but OTS blocks with different durations should not be
    ## equivalent.
    ## ========================================================================
    
    nobs <- 0
    if (fitData$OT$flag) nobs <- nobs + fitData$OT$n
    if (fitData$MAX$flag) nobs <- nobs + sum(fitData$MAX$r)
    if (fitData$OTS$flag) nobs <- nobs + sum(fitData$OTS$r)
    

    if ((!missing(coefLower) || !missing(coefUpper)) && (estim == "optim")) {
        warning("Provided values for 'coefLower' and 'coefUpper' will be ",
                "ignored since 'estim' is set to \"optim\"")
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
                       nobs = nobs,
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
