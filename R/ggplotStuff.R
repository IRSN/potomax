## ****************************************************************************
##' Draw a ggplot layer for a \code{potData} object.
##'
##' @title Draw a ggplot Layer for a \code{potData} object
##'
##' @method autolayer potData
##'
##' @usage
##' \method{autolayer}{potData}(object,
##'           type = c("RLplot", "timeplot"),
##'           group,
##'           points = c("p", "H"), a = 0.5,
##'           ...)
##' 
##' @param object A \code{potData} object.
##'
##' @param type Type of plot wanted.
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
##' @param ... Other arguments passed to \code{geom_point}.
##'
##' @return A ggplot layer.
##'
##' @seealso \code{\link{potData}}, \code{\link[ggplot2]{scale_colour_manual}} and
##' \code{\link[ggplot2]{scale_shape_manual}}.
##'
##' @note The user might have change the colour and the shape of the
##' points by using \code{scale_colour_manual} and
##' \code{scale_shape_manual}.
##' 
##' @examples
##' pdat <- potData(data = Garonne$OTdata$Flow, effDuration = 65,
##'                 MAX.data = Garonne$MAXdata$Flow, MAX.effDuration = 143)
##' ggplot() + autolayer(pdat) + scale_x_log10() + theme_gray()
##' ## Change the label of the historical block
##' pdat <- potData(data = Garonne$OTdata$Flow, effDuration = 65,
##'                 MAX.data = list("hist" = Garonne$MAXdata$Flow),
##'                 MAX.effDuration = 143)
##' ggplot() + autolayer(pdat) + scale_x_log10() + theme_gray()
autolayer.potData <- function(object,
                              type = c("RLplot", "timeplot"),
                              group, ## = c("block", "source"),
                              points = c("p", "H"),
                              a = 0.5,
                              ...) {
    type <- match.arg(type)
    ## group <- match.arg(group)
    
    if (type == "RLplot") {

        df <- RP.potData(object, points = points, a = a)$data
        if (nlevels(df$block) < 6 ) {
            group <- "block"
        } else {
            group <- "source"
        }

        
    } else {
        stop("Not implemented yet!")
    }

    if (group == "block") {
        geom_point(data = df,
                   mapping = aes_string(x = "T", y = "x", group = "block",
                       shape = "block", colour = "block", fill = "block"),
                   ...) 
    } else {
        geom_point(data = df,
                   mapping = aes_string(x = "T", y = "x",
                       shape = "source", group = "source", colour = "source",
                                        fill = "source"),
                   ...)
    }

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
    
    gg <- gg + scale_x_log10() 
 
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
##'          points = c("H", "p", "none"), a = 0.5,
##'          allPoints = FALSE,
##'          ...)
##' 
##' @param object An object with class \code{"poisGP"}.
##'
##' @param which The type of plot to be built: \code{"RL"} is for a
##' Return Level plot and \code{"pp"} is for a probability plot.
##'
##' @param points The plotting position system to be used to show the
##' data attached to the object. When \code{points} is \code{"none"}, 
##'
##' @param a The parameter used to compute the plotting positions when
##' \code{points} is \code{"p"}. See \code{\link{RP.potData}}.
##'
##' @param allPoints Logical. With the default \code{FALSE}, the
##' limits for the \code{y} axis are set to show only the points over
##' the threshold. The value \code{TRUE} leads to showing all points
##' in the data including those which are below the threshold.
##' 
##' @param ... Not used.
##'
##' @return A ggplot graphics object.
##' 
autoplot.poisGP <- function(object,
                            which = c("RL", "pp"),
                            points = c("H", "p", "none"),
                            a = 0.5,
                            allPoints = FALSE,
                            ...) {

    which <- match.arg(which)
    points <- match.arg(points)

    if (which == "RL") {

        RLs <- RL(object)
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
##' @title Autoplot method for a list of Return Levels Tables
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
##' ##' Autoplot method for objects representing fitted Poisson-GP models. 
##'
##' @title Autoplot Method for \code{poisGP} Objects
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
    
    gg <- gg + facet_wrap(Name ~ ., scales = "free_x")
    
    gg

}


