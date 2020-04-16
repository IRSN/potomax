

## ****************************************************************************
##' Draw a ggplot layer for a \code{potData} object.
##'
##' @title Draw a ggplot Layer for a \code{potData} Object
##'
##' @method autolayer potData
##'
##' @usage
##' \method{autolayer}{potData}(object,
##'           type = c("RLplot", "timeplot"),
##'           aes = FALSE,
##'           xVar = c("T", "p"),
##'           group,
##'           points = c("p", "H"), a = 0.5,
##'           ...)
##' 
##' @param object A \code{potData} object.
##'
##' @param type See \code{\link{autoplot.potData}}.
##'
##' @param aes Logical. If \code{TRUE} the colour, fill and shape of
##' the points are used within the aesthetic function \code{aes}
##' wo they are registred to appear in legends.
##' 
##' @param xVar See \code{\link{autoplot.potData}}.
##'
##' @param group See \code{\link{autoplot.potData}}.
##' 
##' @param points,a  See \code{\link{autoplot.potData}}.
##' ' \code{"p"}. See \code{\link{RP.potData}}.
##'
##' @param ... Other arguments passed to \code{geom_point}.
##'
##' @return A ggplot object.
##'
##' @seealso \code{\link{potData}}, \code{\link{autoplot.potData}}.
##'
##' @note With \code{aes = TRUE} it can be needed to reset the scales
##' to avoid the creation of confusing legends. The package
##' \strong{ggnewscale} can be of some help then.
##'
##' @examples
##' fit <- poisGP(data = Garonne, threshold = 2800)
##' RL <- RL(fit)
##' g <- autoplot(RL)
##' g
##' pdat <- as.potData(Garonne)
##' g1 <- g + autolayer(pdat, aes = TRUE)
##' g1
##' 
autolayer.potData <- function(object,
                              type = c("RLplot", "timeplot"),
                              aes = FALSE,
                              xVar = c("T", "p"),
                              group, ## = c("block", "source"),
                              points = c("p", "H"),
                              a = 0.5,
                              ...) {
    block <- NULL ## to avoid warning
    
    pColour <-  c("black", "orangered", "SpringGreen3",  "SteelBlue3")
    pShape <-  c(21, 21, 23, 24, 25)
    pFill <-  c("black", "gold", "SpringGreen1", "SteelBlue1")
   
    type <- match.arg(type)
    xVar <- match.arg(xVar)
    ## group <- match.arg(group)
    
    if (type == "RLplot") {
        
        df <- RP(object, points = points, a = a)$data
        if (nlevels(df$block) < 6 ) {
            group <- "block"
            df <- within(df, Type <- as.factor(block))
        } else {
            group <- "source"
            df <- within(df, Type <- as.factor(source))
        }
        
    } else {
        stop("Not implemented yet!")
    }
    
    out <- list()

    if (aes) {
        out[[1]]  <- geom_point(data = df,
                                mapping = aes_string(x = "T", y = "x",
                                    group = "Type",
                                    shape = "Type", colour = "Type", fill = "Type"),
                                stroke = 0.9, alpha = 0.8)
    } else {
        i <- 1
        for (lev in levels(df$Type)) {
            out[[i]]  <- geom_point(data = subset(df, Type == lev),
                                    mapping = aes_string(x = "T", y = "x",
                                        group = "Type"),
                                    shape = pShape[i], colour = pColour[i],
                                    fill = pFill[i],
                                    stroke = 0.9, alpha = 0.8)
            i <- i + 1
        }

    }
 
    out

}

## ****************************************************************************
##' Build ggplot graphics for a \code{potData} object.
##'
##' @title Build ggplot Graphics for a \code{potData} Object
##'
##' @method autoplot potData
##'
##' @usage
##' \method{autoplot}{potData}(object,
##'           type = c("RLplot", "timeplot"),
##'           aes = FALSE,
##'           xVar = c("T", "p"),
##'           group,
##'           points = c("p", "H"), a = 0.5,
##'           blockDuration, 
##'           ...)
##' 
##' @param object A \code{potData} object.
##'
##' @param type Type of plot wanted.
##'
##' @param aes Logical. If \code{TRUE} the colour, fill and shape of
##' the points are used within the aesthetic function \code{aes} wo
##' they are registred to appear in legends.
##' 
##' @param xVar Used when \code{type} is \code{"RLplot"}. The variable
##' to map the absissa \eqn{x} in the return level plot. When
##' \code{xVar} is \code{"T"} the plotting positions are computed as
##' explained in \emph{Renext Computing Details}. For \code{xVar =
##' "p"} a further step is required to probability of exceedance
##' \eqn{p} in reference with a given block duration. 
##'
##' @param group Character with value in \code{"block"} or
##' \code{"source"}. In the first case the color and the shape of the
##' points depend on the block. In the second case the color and the
##' shape depend on the type or source of the block: \code{"OT"},
##' \code{"MAX"} and \code{"OTS"}. By default the value of
##' \code{group} is chosen to be \code{"block"} when the number of
##' block is small enough and \code{"source"} otherwise.
##' 
##' @param points Type of plotting positions to be used for a RL
##' plot. See \code{\link{RP.potData}}.
##'
##' @param a Parameter to used in \code{ppoints} when \code{points} is
##' \code{"p"}. See \code{\link{RP.potData}}.
##'
##' @param blockDuration A block duration used to find the plotting
##' positions when \code{xVar} is \code{"p"}.
##' 
##' @param ... Other arguments passed to \code{geom_point}.
##'
##' @return A ggplot layer.
##'
##' @seealso \code{\link{potData}}, \code{\link[ggplot2]{scale_manual}}.
##'
##' @note The user might have change the colour and the shape of the
##' points by using \code{scale_colour_manual} and
##' \code{scale_shape_manual}. 
##'
##' @references
##' 
##' Yves Deville (2020). \emph{Renext Computing Details}. Technical
##' Report.
##' 
##' 
##' @examples
##' pdat <- as.potData(Garonne)
##' autoplot(pdat) + ggtitle("Garonne potData")
##' g <- autoplot(pdat, aes = TRUE) + ggtitle("Garonne potData")
##' g
##' g <- g +  scale_colour_manual(values = c("firebrick", "SpringGreen")) +
##'        scale_shape_manual(values = c(21, 22))
##' 
##' ## use autolayer
##' ggplot() + autolayer(pdat) + scale_x_log10() + theme_gray()
##' 
##' ## Change the label of the historical block
##' pdat <- potData(data = Garonne$OTdata$Flow, effDuration = 65,
##'                 MAX.data = list("hist" = Garonne$MAXdata$Flow),
##'                 MAX.effDuration = 143)
##' autoplot(pdat, aes = TRUE)
##' 
autoplot.potData <-  function(object,
                              type = c("RLplot", "timeplot"),
                              aes = FALSE,
                              xVar = c("T", "p"),
                              group, ## = c("block", "source"),
                              points = c("p", "H"),
                              a = 0.5,
                              blockDuration, 
                              ...) {

    type <- match.arg(type)
    if (type == "p") {
        stop("'type = \"p\"' not allowed yet for potData") 
    }
    
    gg <- ggplot()
    gg <- gg + autolayer(object = object,
                         type = type,
                         aes = aes,
                         xVar = xVar,
                         points = points,
                         a = a,
                         ...) +
        scale_x_log10()  + xlab("Period (years)") + ylab("Quantile")
    gg
}
        
## ****************************************************************************
##' Autoplot method for objects of class \code{"RL.poisGP"}
##' representing a table of return levels.
##'
##' @title Autoplot Method for Return Levels
##'
##' @method autoplot poisGP
##'
##' @usage
##' \method{autoplot}{RL.poisGP}(object, ...)
##' 
##' @param object An object of class \code{"RL.poisGP"}. It inherits
##' from \code{data.frame} and has columns containing the quantiles
##' and confidence limits.
##'
##' @param ... Not used yet.
##' 
##' @return A ggplot graphics object.
##'
##' 
autoplot.RL.poisGP <- function(object, ...) {

    gg <- ggplot()

    if (all(c("L", "U") %in% names(object))) {
        gg <- gg +
            geom_ribbon(
                data = object,
                mapping = aes_string(x = "Period", ymin = "L", ymax = "U",
                    group = "Level"),
                fill = "SteelBlue2", alpha = 0.3)
        
        gg <- gg +
            geom_line(
                data = object,
                mapping = aes_string(x = "Period", y = "L", group = "Level",
                    linetype = "Level"),
                colour = "SteelBlue4")
        
        gg <- gg +
            geom_line(
                data = object,
                mapping = aes_string(x = "Period", y = "U", group = "Level",
                    linetype = "Level"),
                colour = "SteelBlue4")
    }
    
    gg <- gg +
        geom_line(
            data = object,
            mapping = aes_string(x = "Period", y = "Quant"),
            colour = "orangered", size = 0.8)

    u <- attr(object, "threshold")
    if (!is.null(u) && !is.na(u)) {
        gg <- gg + geom_hline(yintercept = u, 
                              color = "SpringGreen3", alpha = 0.6)
    }
    
    gg <- gg + scale_x_log10() + ylab("Quantile")
 
    gg

}

## ****************************************************************************
##' Autoplot method for objects representing fitted Poisson-GP models. 
##'
##' @title Autoplot Method for \code{poisGP} Objects
##'
##' @method autoplot poisGP
##' 
##' @usage
##' \method{autoplot}{poisGP}(object, 
##'          which = c("RL", "pp"),
##'          level = 0.70,
##'          points = c("H", "p", "none"), a = 0.5,
##'          allPoints = FALSE,
##'          trace = 0,
##'          ...)
##' 
##' @param object An object with class \code{"poisGP"}.
##' 
##' @param which The type of plot to be built: \code{"RL"} is for a
##' Return Level plot and \code{"pp"} is for a probability plot.
##'
##' @param level Used when \code{which} is \code{"RL"} and then passed
##' to \code{\link{RL.poisGP}}.
##' 
##' @param points The plotting position system to be used to show the
##' data attached to the object. When \code{points} is \code{"none"},
##' the data are not shown. With the value \code{"H"} the data are
##' shown using the "H-points" of Nelson's plotting positions. With
##' the value \code{"p"}, the data are shown using the "p-points" are
##' used and the argument \code{a} can be passed to the \code{ppoints}
##' function.
##'
##' @param a The parameter used to compute the plotting positions when
##' \code{points} is \code{"p"}. See \code{\link{RP.potData}}.
##'
##' @param allPoints Logical. With the default \code{FALSE}, the
##' limits for the \code{y} axis are set to show only the points over
##' the threshold. The value \code{TRUE} leads to showing all points
##' in the data including those which are below the threshold.
##'
##' @param trace Integer level of verbosity. Passed to
##' \code{\link{RL.poisGP}}.
##' 
##' @param ... Not used.
##'
##' @return A ggplot graphics object.
##'
##' @seealso \code{\link[Renext]{Hpoints}} in \strong{Renext} and
##' \code{\link[stats]{ppoints}} in \strong{stats}.
##' 
autoplot.poisGP <- function(object,
                            which = c("RL", "pp"),
                            level = 0.70,
                            points = c("H", "p", "none"),
                            a = 0.5,
                            allPoints = FALSE,
                            trace = 0,
                            ...) {

    which <- match.arg(which)
    points <- match.arg(points)

    if (which == "RL") {

        RLs <- RL(object, level = level, trace = trace)
        gg <- autoplot(RLs, ...)
        if (!allPoints)  gg <- gg + ylim(object$threshold, NA)
        
        if (points != "none") {
            gg <- gg + autolayer(object$data, points = points, a = a,
                                 stroke = 0.9, alpha = 0.7)
            gg <- gg +
                scale_colour_manual(name = "Type",
                                    values = c("black", "orangered", "SpringGreen3",
                                               "SteelBlue3"))
            gg <- gg + scale_shape_manual(name = "Type",
                                          values = c(21, 21, 22, 24))
            gg <- gg +
                scale_fill_manual(name = "Type",
                                    values = c("black", "gold", "SpringGreen1", "SteemBlue1"))
        }

    } else {
        cat("Not implemented yet!")
    }
   
    
    gg

}
## ****************************************************************************
##' Autoplot a list of RL tables, generally using the same data or
##' using different data but sharing some observations.
##'
##' @title Autoplot Method for a List of Return Levels Tables
##'
##' @method autoplot RL.poisGPList
##' 
##' @param object A list of Return Level tables, with class
##' \code{"RL.poisGPList"}.
##' 
##' @param ... Not used yet.
##' 
##' @return A ggplot graphics object.
##'
##' @examples
##' u <- seq(from = 2610, to = 3000, length.out = 4)
##' fits <- RLs <-  list()
##' for (iu in seq_along(u)) {
##'     fits[[iu]] <- poisGP(data = Garonne$OTdata$Flow,
##'                          effDuration = Garonne$OTinfo$effDuration,
##'                          MAX.data = list("hist" = Garonne$MAXdata$Flow),
##'                          MAX.effDuration = Garonne$MAXinfo$duration,
##'                          threshold = u[iu])
##'     RLs[[iu]] <- RL(fits[[iu]], conf = "prof", out = "data")
##' }
##' names(fits) <- names(RLs) <- paste("thresh =",  u)
##'
##' ## force the S3 class of the result 
##' class(RLs) <- "RL.poisGPList"
##' autoplot(RLs) + theme_gray() +
##'    ggtitle("Return level plots for Garonne and several thresholds")
##' 
autoplot.RL.poisGPList <- function(object, ...) {

    chck <-  sapply(object, function(x) inherits(x, "RL.poisGP"))
    if (!all(chck)) {
        stop("All elements of `object` must inherit from \"RL.poisGP\"")
    }
    
    nms <- names(object)
    
    for (i in seq_along(object)) {

        if (i == 1) {
            df <- data.frame(".name" = nms[i], object[[i]],
                             stringsAsFactors = FALSE)
        } else {
            df <- dplyr::bind_rows(df,
                                   data.frame(".name" = nms[i], object[[i]],
                                              stringsAsFactors = FALSE))
        }
        
    }
    df <- within(df, .name <- factor(.name))

    ## uggly hack, but avoids writing code :)
    class(df) <- c("RL.poisGP", "data.frame")
    autoplot(df) + facet_wrap(.name ~ .)

}
## ****************************************************************************
##' Autoplot method for objects devoted to the check of the
##' profile-likelihood confidence intervals on the parameters for
##' Poisson-GP models.
##'
##' @details
##' The object is a list containing two data frames named \code{ci}
##' and \code{negLogLikC}. The \code{ci} data frame contains the
##' confidence intervals in a long format suitable for plotting; it
##' will be used to draw horizontal and vertical lines showing the
##' confidence limits. The \code{negLogLikC} data frame contains
##' evaluations of the profile log-likelihood for a grid of value of
##' each parameter; it will be used to display a curve for each
##' parameter.
##'
##' For each parameter, the two verical lines corresponding to a given
##' confidence level represent the lower and upper bounds of the
##' confidence interval for the parameter. The curve represent the
##' negative profiled log-likelihood for the considered parameter. For
##' each confidence level the two vertical lines should intersect the
##' curve at points where the value of the profiled negative
##' log-likelhood corrresponds to the horizontal line for the the
##' level. The horizontal and vertical lines for the estimate should
##' intersect at the point where the value of the negative
##' log-likelhood is minimal. If the later condition is not true then
##' the convergence of the ML estimation was not reached for some
##' reason.
##' 
##' @title Autoplot Method for \code{confintCheck.poisGP} Objects
##'
##' @method autoplot confintCheck.poisGP
##' 
##' @param object An object with class \code{"confIntCheck.poisGP"}
##' generated by the \code{confint} method see
##' \code{\link{confint.poisGP}}.
##'
##' @param ... Not used yet.
##'
##' @return A ggplot graphics object.
##'
##' example(confint.poisGP)
##' 
autoplot.confintCheck.poisGP <- function(object, ...) {
    
    gg <- ggplot()
    gg <- gg + geom_line(data = object$negLogLikC,
                         mapping = aes_string(x = "Par", y = "Value", group = "Name"))
    
    gg <- gg + geom_point(data = object$negLogLikC,
                          mapping = aes_string(x = "Par", y = "Value", group = "Name"),
                          size = 0.8)
    
    gg <- gg + geom_vline(data = object$ci,
                          mapping = aes_string(xintercept = "Value", group = "Name",
                              linetype = "Level", colour = "Level"))
    
    gg <- gg + geom_hline(data = object$ci,
                          mapping = aes_string(yintercept = "NegLogLik", group = "Name",
                              linetype = "Level", colour = "Level"))
    
    gg <- gg + facet_wrap(Name ~ ., scales = "free_x") + ylab("negative log lik.") +
        xlab("param. value")
    
    gg

}


## ****************************************************************************
##' Autoplot method for objects devoted to the check of the
##' profile-likelihood confidence intervals on the return levels for
##' Poisson-GP models.
##'
##' @details The object is a list containing two data frames named
##' \code{RL} and \code{negLogLikC}. The \code{RL} data frame contains
##' the confidence intervals in a long format suitable for plotting;
##' it will be used to draw horizontal and vertical lines showing the
##' confidence limits. The \code{negLogLikC} data frame contains
##' evaluations of the profile log-likelihood for a grid of value of
##' each return period; it will be used to display a curve for each
##' parameter.
##'
##' The check is similar to that used for the confidence intervals as
##' implemented in \code{\link{autoplot.confintCheck.poisGP}}. Each
##' return level is here used as a parameter.
##' 
##' @title Autoplot Method for \code{RLCheck.poisGP} Objects
##'
##' @method autoplot RLCheck.poisGP
##' 
##' @param object An object with class \code{"RLCheck.poisGP"}
##' generated by the \code{confint} method see
##' \code{\link{confint.poisGP}}.
##'
##' @param ... Not used yet.
##'
##' @return A ggplot graphics object.
##'
##' @seealso \code{\link{autoplot.confintCheck.poisGP}}
##' 
autoplot.RLCheck.poisGP <- function(object, ...) {
    
    gg <- ggplot()
    gg <- gg + geom_line(data = object$negLogLikC,
                         mapping = aes_string(x = "rho", y = "Value", group = "Period"))
    
    gg <- gg + geom_point(data = object$negLogLikC,
                          mapping = aes_string(x = "rho", y = "Value", group = "Period"),
                          size = 0.8)
    
    gg <- gg + geom_vline(data = object$RL,
                          mapping = aes_string(xintercept = "Value", group = "Period",
                              linetype = "Level", colour = "Level"))
    
    gg <- gg + geom_hline(data = object$RL,
                          mapping = aes_string(yintercept = "NegLogLik", group = "Period",
                              linetype = "Level", colour = "Level"))
    
    gg <- gg + facet_wrap(Period ~ ., scales = "free_x", labeller = label_both)
    
    gg <- gg + ylim(c(NA, object$ylim)) + ylab("negative log lik.") +
        xlab("Return Level")

    gg
}


