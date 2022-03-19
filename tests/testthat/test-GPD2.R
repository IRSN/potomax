context("GPD2")

## *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the implementation of the GPD2 distribution (C code used via
## .Call)
## *****************************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

## =============================================================================
## check that the GPD2 quantile and distribution functions are consistent
## =============================================================================

n <- 100
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)

x <- as.vector(rGPD2(n = n, scale = sigma, shape = xi))
F <- pGPD2(x, scale = sigma, shape = xi)
e <- x - qGPD2(p = F, scale = sigma, shape = xi)

test_that(desc = "Consistency of the GPD2 cdf and quantile funs #1",
          expect_lt(max(abs(e)), 1e-10))

p <- runif(n)
q <- qGPD2(p, scale = sigma, shape = xi)
e <- p - pGPD2(q, scale = sigma, shape = xi)
test_that(desc = "Consistency of the GPD2 cdf and quantile funs #2",
          expect_lt(max(abs(e)), 1e-10))

## =============================================================================
## check that the GPD2 density and distribution functions are
## consistent. Note that we do not check that the computation is 
## precise but rather that there is no error in the code.
## =============================================================================

n <- 100
sigma <- rexp(1)

for (xi in c(-0.1, -1e-4, 0.0, 1e-4, 0.1)) {
    
    x <- rGPD2(n, scale = sigma, shape = xi)
    
    F <- function(x) {
        pGPD2(x, scale = sigma, shape = xi)
    }
    
    fval <- dGPD2(x, scale = sigma, shape = xi)
    eps <- 1e-6
    e <- fval - (F(x + eps) - F(x - eps)) / 2 / eps
    
    test_that(desc = "Consistency of the GPD2 cdf density funs",
              expect_lt(max(abs(e)), 1e-6))
}


## =============================================================================
## check the gradient and the Hessian of the log-density Note that the
## numerical error on the Hessian is larger than that on the gradient
## and that the computation is more difficult when xi is close to
## zero.
## =============================================================================

n <- 1

for (xi in c(-0.2, -0.1, -1e-3, -1e-4, 0.0, 1e-4, 0.1, 0.2)) {

    sigma <- rexp(1)
    x <- as.vector(rGPD2(n = n, scale = sigma, shape = xi))
    
    f <- function(theta) {
        dGPD2(x, scale = theta[1], shape = theta[2], log = TRUE, deriv = TRUE,
              hessian = TRUE)
    }
    
    fval <- dGPD2(x = x, scale = sigma, shape = xi, log = TRUE,
                  deriv = TRUE, hessian = TRUE)
    
    theta0 <- c(scale = sigma, shape = xi)
    Jnum <- jacobian(func = f, x = theta0)
    Hnum <- hessian(func = f, x = theta0)

    Jcomp <- attr(fval, "gradient")
    errGrad <- abs(Jcomp - Jnum)
    errGradRel <- abs(errGrad / Jnum) 
    test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)
    
    test_that(desc = sprintf("Gradient of the GPD2 log-density for xi = %7.5f",
                             xi),
              expect_true(all(test)))

    Hcomp <- drop(attr(fval, "gradient"))
    errHess <- abs(Hcomp - Hnum)
    errHessRel <- abs(errHess / Hnum) 
    test <- (errHess <- 1e-4) | (errHessRel < 1e-3)
    
    test_that(desc = sprintf("Hessian of the GPD2 log-density for xi = %7.5f",
                             xi),
              expect_true(all(test)))
}
    
## =============================================================================
## check that the Hessian of the log-likelihood is OK
## =============================================================================

n <- 300
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)
y <- as.vector(rGPD2(n = n, scale = sigma, shape = xi))
theta <- c(scale = sigma, shape = xi)

## this is for numerical differentiation
f2 <- function(theta, y) {
    sum(dGPD2(x = y, scale = theta[1],
             shape = theta[2], log = TRUE))
}

## this is to get the Hessian as an attribute
fval <- dGPD2(x = y, scale = sigma, shape = xi, log = TRUE,
             deriv = TRUE, hessian = TRUE)

H2 <- apply(attr(fval, "hessian"), MARGIN = c(2, 3), FUN = sum)
H2num <- hessian(func = f2, x = theta, y = y)

test_that(desc = "Hessian of the GPD2 log-lik",
            expect_lt(max(abs(H2 - H2num)), 1e-3))


## ==========================================================================
## Check the gradient and the Hessian of the quantile function
## ==========================================================================

n <- 1

f <- function(theta) {
    qGPD2(p, scale = theta[1], shape = theta[2])
}

for (xi in c(-0.2, -0.1, -1e-3, -1e-4, 0.0, 1e-4, 0.1, 0.2)) {
    sigma <- rexp(1)
    p <- runif(1)
    
    
    theta0 <- c(scale = sigma, shape = xi)
    fval <- qGPD2(p = p, scale = sigma, shape = xi,
                  deriv = TRUE, hessian = TRUE)
    
    Jnum <- jacobian(func = f, x = theta0)
    Hnum <- hessian(func = f, x = theta0) 
    
    Jcomp <- attr(fval, "gradient")
    errGrad <- abs(Jcomp - Jnum)
    errGradRel <- abs(errGrad / Jnum) 
    test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)
    
    test_that(desc = sprintf("Gradient of the GPD2 quantile for xi = %7.5f",
                             xi),
              expect_true(all(test)))

    Hcomp <- drop(attr(fval, "gradient"))
    errHess <- abs(Hcomp - Hnum)
    errHessRel <- abs(errHess / Hnum) 
    test <- (errHess <- 1e-4) | (errHessRel < 1e-3)
    
    test_that(desc = sprintf("Hessian of the GPD2 quantile for xi = %7.5f",
                             xi),
              expect_true(all(test)))
}
    
## ==========================================================================
## Check the gradient and the Hessian of the distribution function
## ==========================================================================

f <- function(theta) {
    pGPD2(x, scale = theta[1], shape = theta[2])
}

n <- 1

for (xi in c(-0.2, -0.1, -1e-3, -1e-4, 0.0, 1e-4, 0.1, 0.2)) {

    sigma <- rexp(1)
    x <- rGPD2(n, scale = sigma, shape = xi) 
    
    theta0 <- c(scale = sigma, shape = xi)
    fval <- pGPD2(x, scale = sigma, shape = xi, deriv = TRUE, hessian = TRUE)
    
    Jnum <- jacobian(func = f, x = theta0)
    Hnum <- hessian(func = f, x = theta0)

    errGrad <- abs(Jcomp - Jnum)
    errGradRel <- abs(errGrad / Jnum) 
    test <- (errGrad <- 1e-5) | (errGradRel < 1e-4)
    
    test_that(desc = sprintf("Gradient of the GPD2 dist function for xi = %7.5f",
                             xi),
              expect_true(all(test)))

    Hcomp <- drop(attr(fval, "gradient"))
    errHess <- abs(Hcomp - Hnum)
    errHessRel <- abs(errHess / Hnum) 
    test <- (errHess <- 1e-4) | (errHessRel < 1e-3)
    
    test_that(desc = sprintf("Hessian of the GPD2 dist function for xi = %7.5f",
                             xi),
              expect_true(all(test)))
    
}

## Note that the numerical Hessian is very poor when 'sigma'
## is > 1
sigma <- runif(1)
x <- rGPD2(n, scale = sigma, shape = 0.0) 

theta0 <- c(scale = sigma, shape = xi)
fval <- pGPD2(x, scale = sigma, shape = xi, deriv = TRUE, hessian = TRUE)

Jnum <- jacobian(func = f, x = theta0)
Hnum <- hessian(func = f, x = theta0)

## print(theta0)
## print(Hnum)
## print(attr(fval, "hessian")[1, , ])

test_that(desc = "Gradient of the GPD2 distribution function shape == 0",
            expect_lt(max(abs(attr(fval, "gradient") - Jnum)), 1e-6))

test_that(desc = "Hessian of the GPD2 distribution function shape == 0",
            expect_lt(max(abs(drop(attr(fval, "hessian")) - Hnum)),
                      1e-3))
