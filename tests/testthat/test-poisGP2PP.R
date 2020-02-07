context("poisGP2PP")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the transformation poisGP -> PP and its derivative
## *************G**************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

n <- 1
lambda <- rexp(n)
loc <- rnorm(n, sd = 20)
scale <- rgamma(n, shape = 2)
shape <- c(rnorm(n, sd = 0.1) , 0.0)
cases <- c("non-zero", "zero")

for (i in seq_along(shape)) {
    
    res <- poisGP2PP(lambda = lambda,
                     loc = loc, scale = scale, shape = shape[i],
                     deriv = TRUE)
    
    theta <- c("lambda" = lambda, "scale" = scale, "shape" = shape[i])

    ## Compute with Renext
    res0 <- Renext::Ren2gev(theta, threshold = loc)
    
    g <- drop(attr(res, "gradient"))
    g0 <- attr(res0, "jacobian")
    
    attr(res, "gradient") <- NULL
    attr(res0, "jacobian") <- attr(res0, "threshold") <- NULL
    
    e <- drop(res) - res0
    test_that(desc = sprintf("Case shape %s, consistency of value with Renext",
                  cases[i]),
              expect_lt(max(abs(e)), 1e-10))

    e <- g[ , -2] - g0

    cat("e in test\n")
    print(e)
    test_that(desc = sprintf("Case shape %s, consistency of jacobian with Renext",
                  cases[i]),
              expect_lt(max(abs(e)), 1e-10))

}
