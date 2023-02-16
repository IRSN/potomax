**potomax** Package News
===========================

# News in version 0.2.3

## Changes

- The GPD probability functions `dGPD2`, `pGPD2`, `qGPD2` and `rGPD2`
  and the utility functions `poisGP2PP` and `PP2poisGP` have been
  moved (with some changes) to the **nieve** package, which is
  available both [on   CRAN](https://CRAN.R-project.org/package=nieve)
  and
  [on GitHub](https://github.com/yvesdeville/nieve/). So **potomax**
  now imports these functions from ** nieve**. As a consequence the
  compiled code used by this functions has been discarded, making the
  package easier to install at least for Windows users. 
  
# News in version 0.2.3

## Changes 

- Changed the dependencies.
  
# New in version 0.2.1

- New vignette: *R Package potomax: Overview*, including the Venice
  example.

## Version 0.2.0

- The one-parameter exponential distribution can now be used as a
  distribution for the excesses. This leads to Gumbel-distributed
  block maxima.
  
- Bug fix. The quantile function `qGPD2` did not work as expected with
 `lower.tail = FALSE`.

- The approximation of the GPD2 distribution functions for $\xi
  \approx 0$ has been improved. For each of the three functions
  `dGPD2`, `pGPD2` and `qGPD2`, a quadratic function of $\xi$ is used
  for small $\xi$ providing a Taylor series approximation. The
  derivatives w.r.t. the parameters are consistently taken as
  functions of $\xi$ that are linear (gradient) and constant
  (hessian). The distribution functions are now smooth functions of
  $\xi$ even for very small $\xi$, say $\xi \approx 10^{-6}$ or
  less. Some tests devoted to this "small $\xi$" approximation have
  been added to the package.
