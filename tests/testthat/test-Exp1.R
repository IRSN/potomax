context("GPD2")

## *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the implementation of the Exp1 distribution (C code used via
## .Call)
## *****************************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

## =============================================================================
## check that the Exp1 quantile and distribution functions are consistent
## =============================================================================

n <- 100
sigma <- rexp(1)

x <- as.vector(rexp1(n = n, scale = sigma))
F <- pexp1(x, scale = sigma)
e <- x - qexp1(p = F, scale = sigma)

test_that(desc = "Consistency of the Exp1 cdf and quantile funs #1",
          expect_lt(max(abs(e)), 1e-10))

p <- runif(n)
q <- qexp1(p, scale = sigma)
e <- p - pexp1(q, scale = sigma)
test_that(desc = "Consistency of the Exp1 cdf and quantile funs #2",
          expect_lt(max(abs(e)), 1e-10))

## =============================================================================
## check that the Exp1 density and distribution functions are
## consistent. Note that we do not check that the computation is 
## precise but rather that there is no error in the code.
## =============================================================================

n <- 100
sigma <- rexp(1)

x <- rexp1(n, scale = sigma)
    
F <- function(x) {
    pexp1(x, scale = sigma)
}

fval <- dexp1(x, scale = sigma)
eps <- 1e-6
e <- fval - (F(x + eps) - F(x - eps)) / 2 / eps

test_that(desc = "Consistency of the exp1 cdf density funs",
          expect_lt(max(abs(e)), 1e-6))

## =============================================================================
## check the gradient and the Hessian of the log-density 
## =============================================================================

n <- 1

sigma <- rexp(1)
x <- as.vector(rexp1(n = n, scale = sigma))

f <- function(theta) {
    dexp1(x, scale = theta[1], log = TRUE, deriv = TRUE,
          hessian = TRUE)
}

fval <- dexp1(x = x, scale = sigma, log = TRUE,
              deriv = TRUE, hessian = TRUE)

theta0 <- c(scale = sigma)
Jnum <- jacobian(func = f, x = theta0)
Hnum <- hessian(func = f, x = theta0)

Jcomp <- attr(fval, "gradient")
errGrad <- abs(Jcomp - Jnum)
errGradRel <- abs(errGrad / Jnum) 
test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)

test_that(desc = "Gradient of the exp1 log-density",
          expect_true(all(test)))

Hcomp <- drop(attr(fval, "gradient"))
errHess <- abs(Hcomp - Hnum)
errHessRel <- abs(errHess / Hnum) 
test <- (errHess <- 1e-4) | (errHessRel < 1e-3)

test_that(desc = "Hessian of the exp1 log-density",
          expect_true(all(test)))

    
## =============================================================================
## check that the Hessian of the log-likelihood is OK
## =============================================================================

n <- 300
sigma <- rexp(1)
y <- as.vector(rexp1(n = n, scale = sigma))
theta <- c(scale = sigma)

## this is for numerical differentiation
f2 <- function(theta, y) {
    sum(dexp1(x = y, scale = theta[1], log = TRUE))
}

## this is to get the Hessian as an attribute
fval <- dexp1(x = y, scale = sigma, log = TRUE,
             deriv = TRUE, hessian = TRUE)

H2 <- apply(attr(fval, "hessian"), MARGIN = 2, FUN = sum)
H2num <- hessian(func = f2, x = theta, y = y)

test_that(desc = "Hessian of the exp1 log-lik",
            expect_lt(max(abs(H2 - H2num)), 1e-3))


## ==========================================================================
## Check the gradient and the Hessian of the quantile function
## ==========================================================================

n <- 1

f <- function(theta) {
    qexp1(p, scale = theta[1])
}

sigma <- rexp(1)
p <- runif(1)

theta0 <- c(scale = sigma)
fval <- qexp1(p = p, scale = sigma)

Jnum <- jacobian(func = f, x = theta0)
Hnum <- hessian(func = f, x = theta0) 

Jcomp <- attr(fval, "gradient")
errGrad <- abs(Jcomp - Jnum)
errGradRel <- abs(errGrad / Jnum) 
test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)

test_that(desc = "Gradient of the Exp1 quantile",
          expect_true(all(test)))

Hcomp <- drop(attr(fval, "gradient"))
errHess <- abs(Hcomp - Hnum)
errHessRel <- abs(errHess / Hnum) 
test <- (errHess <- 1e-4) | (errHessRel < 1e-3)

test_that(desc = "Hessian of the Exp1 quantile",
          expect_true(all(test)))

    
## ==========================================================================
## Check the gradient and the Hessian of the distribution function
## ==========================================================================

f <- function(theta) {
    pexp1(x, scale = theta[1])
}

n <- 1

sigma <- rexp(1)
x <- rexp1(n, scale = sigma) 

theta0 <- c(scale = sigma)
fval <- pexp1(x, scale = sigma)

Jnum <- jacobian(func = f, x = theta0)
Hnum <- hessian(func = f, x = theta0)

errGrad <- abs(Jcomp - Jnum)
errGradRel <- abs(errGrad / Jnum) 
test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)

test_that(desc = "Gradient of the Exp1 dist function",
          expect_true(all(test)))

Hcomp <- drop(attr(fval, "gradient"))
errHess <- abs(Hcomp - Hnum)
errHessRel <- abs(errHess / Hnum) 
test <- (errHess <- 1e-4) | (errHessRel < 1e-3)

test_that(desc = "Hessian of the Exp1 dist function",
          expect_true(all(test)))
