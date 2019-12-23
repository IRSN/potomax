library(numDeriv)
library(testthat)

context("logLikFuns")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com> GOAL: Test the
## implementation of the log-likelihood functions and companions:
## 'lambdaHat', 'logLikFunC' and 'logLikFun'.
## ***************************************************************************

## ============================================================================
## Check that the estimated rate is computed by 'lambdaHat' is close
## to that found by Renext::Renouv
## ============================================================================

set.seed(1234)
u <- seq(from = 2550, to = 3000, by = 50)
testLambdaHat <- testLogLikFun <- testLogLikFunC <- rep(NA, length(u))

for (iu in seq_along(u)) {
    fitp <- poisGP(data = Garonne$OTdata$Flow,
                   threshold = u[iu],
                   effDuration = Garonne$OTinfo$effDuration,
                   MAX.data = Garonne$MAXdata$Flow,
                   MAX.effDuration = Garonne$MAXinfo$duration,
                   estim = "none")
    
    fitR <- Renouv(Garonne, dist = "GPD", threshold = u[iu], plot = FALSE)

    llamb <- lambdaHat(thetaGP = coef(fitR)[-1], object = fitp)
    testLambdaHat[iu] <- abs(exp(llamb) - coef(fitR)[1])
    
    coGPD <- coef(fitR)[-1]
    testLogLikFunC[iu] <-
        max(abs(negLogLikFunC(thetaGP = coGPD, object = fitp) +
                    logLik(fitR)))

    coAll <- coef(fitR)
    testLogLikFun[iu] <-
        max(abs(negLogLikFun(theta = coAll, object = fitp) +
                    logLik(fitR)))
    
}

## ============================================================================
## Consistency of the logLik (at the optimum only) with Renext
## ============================================================================

test_that(desc = "Consistency of 'LambdaHat' and Renext",
          expect_lt(max(testLambdaHat), 1e-5))

test_that(desc = "Consistency of 'logLikFunC' and Renext",
          expect_lt(max(testLogLikFunC), 1e-4))

test_that(desc = "Consistency of 'logLikFun' and Renext",
          expect_lt(max(testLogLikFunC), 1e-4))

## ============================================================================
## Now let us check some gradients
## ============================================================================

coGPD <- coef(fitR)[-1]  + sqrt(diag(vcov(fitR))[-1]) * runif(2)
logLamb <- lambdaHat(thetaGP = coGPD, object = fitp)
gradNum <- grad(func = lambdaHat, x = coGPD, object = fitp, deriv = FALSE)
test_that(desc = "Check the gradient of 'lambdaHat' with log = TRUE",
          expect_lt(max(abs(gradNum - attr(logLamb, "gradient"))), 1e-4))

lamb <- lambdaHat(thetaGP = coGPD, object = fitp, log = FALSE)
gradNum <- grad(func = lambdaHat, x = coGPD, object = fitp,
                log = FALSE, deriv = FALSE)
test_that(desc = "Check the gradient of 'lambdaHat' with log = FALSE",
          expect_lt(max(abs(gradNum - attr(lamb, "gradient"))), 1e-4))

nllC <- negLogLikFunC(thetaGP = coGPD, object = fitp)
gradNum <- grad(func = negLogLikFunC, x = coGPD, object = fitp,
                deriv = FALSE)
test_that(desc = "Check the gradient of 'logLikFunC",
          expect_lt(max(abs(gradNum - attr(nllC, "gradient"))), 1e-4))

coAll <- coef(fitR)  + sqrt(diag(vcov(fitR))) * runif(3)
nll <- negLogLikFun(theta = coAll, object = fitp)
gradNum <- grad(func = negLogLikFun, x = coAll, object = fitp,
                deriv = FALSE)
test_that(desc = "Check the gradient of 'logLikFun",
          expect_lt(max(abs(gradNum - attr(nll, "gradient"))), 1e-4))

