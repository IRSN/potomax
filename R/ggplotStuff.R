## ****************************************************************************
##' Draw a ggplo layer for a \code{potData} object.
##'
##' @title Draw a ggplot Layer for a \code{potData} object.
##'
##' @param object A \code{potData} object.
##'
##' @param type Type of plot wanted.
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
##' @param object An object of class \code{"RL.poisGP"}. It inherits
##' from \code{data.frame} and has columns c
##'
##' @param ... Not used yet.
##' 
##' @return A ggplot graphic object.
##'
##' 
autoplot.RL.poisGP <- function(object, ...) {

    gg <- ggplot()
    
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
##' Autoplot method for objects representing fitted Poisson-Gp models. 
##'
##' @title Autoplot Method for \code{poiGP} Objects
##'
##' @param object An object with class \code{"poisGP"}.
##'
##' @param points The plotting position system to be used to show the
##' data attached to the object. When \code{points} is \code{"none"}, 
##'
##' @param a See \code{\link{RP.potData}}.
##'
##' @param ... 
##'
##' @return
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


    }
   
    
    gg

}

autoplot.RL.poisGPList <- function(object, facets = TRUE,
                                   ...) {

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
