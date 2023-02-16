context("PP2poisGP")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the transformation PP -> Poisson-GP and its derivative
## ***************************************************************************
library(numDeriv)
library(testthat)

set.seed(1234)

n <- 4
locStar <- rnorm(n, sd = 20)
scaleStar <- rgamma(n, shape = 2)
shapeStar <- c(rnorm(n -1, sd = 0.1), 0.0)

p <- runif(n)
threshold <- rep(NA, n)
for (i in 1:n) {
    threshold[i] <- qgev(p = p[i], loc = locStar[i], scale = scaleStar[i],
                         shape = shapeStar[i])
}

cases <- c(rep("non-zero", n - 1), "non-zero")

## ===========================================================================
## Check the n - 1 cases with non-zero shape
## ===========================================================================

for (i in n - 1) {
    
    res <- PP2poisGP(locStar = locStar[i],
                     scaleStar = scaleStar[i],
                     shapeStar = shapeStar[i],
                     threshold = threshold[i],
                     deriv = TRUE)
    
    thetaStar <- c("loc" = locStar[i], "scale" = scaleStar[i],
                   "shape" = shapeStar[i])
    
    ## Compute with Renext
    res0 <- Renext::gev2Ren(thetaStar, threshold = threshold[i],
                            distname.y = "GPD")
    
    g <- drop(attr(res, "gradient"))
    g0 <- attr(res0, "jacobian")
    
    attr(res, "gradient") <- NULL
    attr(res0, "jacobian") <- attr(res0, "threshold") <- NULL
    
    e <- drop(res) - res0
    test_that(desc = sprintf("Case shape %s, consistency of value with Renext",
                  cases[i]),
              expect_lt(max(abs(e)), 1e-10))

    e <- g - g0[-2, ]

    test_that(desc = sprintf("Case shape %s, consistency of jacobian with Renext",
                  cases[i]),
              expect_lt(max(abs(e)), 1e-10))

}

## ===========================================================================
## Check case 'n' zero shape. Use the Jacobian
## ===========================================================================

PP2poisGPFun <- function(thetaStar, threshold) {
    PP2poisGP(locStar = thetaStar[1],
              scaleStar = thetaStar[2],
              shapeStar = thetaStar[3],
              threshold = threshold,
              deriv = FALSE)
    
}

thetaStar <- c("loc" = locStar[n], "scale" = scaleStar[n],
               "shape" = shapeStar[n])

jacNum <- jacobian(func = PP2poisGPFun, x = thetaStar, threshold = threshold[n])

res <- PP2poisGP(locStar = locStar[n],
                 scaleStar = scaleStar[n],
                 shapeStar = shapeStar[n],
                 threshold = threshold[n],
                 deriv = TRUE)
g <- drop(attr(res, "gradient"))

e <- g - jacNum

test_that(desc = sprintf("Case shape %s, Jacobian and numeric Jacobian",
              cases[n]),
          expect_lt(max(abs(e)), 1e-7))
