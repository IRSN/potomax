# TODO list for the `potmax` package


## Improvements in the existing functions and doc

- Implement the `scale = TRUE` possibility in `poisGP`. This needs
  some care because the log-likelhood will depend on the scale and is
  used in several functions or method, not just by the creator
  `poisGP`.

- Write or improve the `summary` and `print` methods for the `"poisGP"` and
  `"potData"` classes.

- **[X]** Warn about AIC and BIC. In the general context of
  heterogeneous data, the concept of *number of observations* can be
  misleading, because a single "MAX" observation can have more weight
  than dozens of other observations. So the AIC and BIC criterion
  should not be used when there is some historical information. A
  possible solution is setting the number of observation to `NA` when
  there is some historical information.

- Improve the exchanges between the class `"potData"` and the class
  `"Rendata"` of **Renext**. We should be able to simulate data, which
  is very useful for testing.

## New functions or methods

- Write functions for the classical diagnostics in *threshold
  choice*. These should include the trace of the estimated shape for a
  varying threshold, along with (profile-likelihood) confidence
  intervals.

- Define *generalised residuals* for the context. Note that these
  should not necessarily be attached to the observations because an
  OTS block with no observation deserves a residual and can have a
  huge weight in the estimation. So some theoretical work is required
  here.

## Changes in the classes and creator functions

- Enhance the `potData` S3 class by allowing it to contain *optional*
  information about time and events. For each block there could be a
  start and an end indication leading to two vectors `start` and `end`
  for each of the three parts `"OT"`, `"MAX"` and `"OTS"`. Similarly
  there could be vectors of class `"Date"` or `"POSIXct"` giving the
  time at which the events occurred. This would allow a time plot for
  the data as implemented by `plot.Rendata` in **Renext**. A
  `simulate` method could also be implemented for the class.

- Consistently with **Renext**, allow the `data` argument of `poisGP`
  to be of class `"potData"`, implying that the other arguments
  `MAX*` and `OTS*` are left missing as well as `effDuration`.

## Data and examples

- Attach to the package several datasets in `potData` format:
  `Garonne`, `venice`, and more. It would be fun to have data for the
  Potomac river:)

## Performance

- Write the log-likelihood functions in C in order to make faster the
  functions/methods using profile-likelihood such as `confint` and `RL`.

- Re-factor the `Ren2gev` function of **Renext** into a C function of
**potomax** accepting vectors for the Poisson-GP parameters `lambda`,
`scale` and `shape` as well as for the block duration `w` and returning
a matrix of GEV or `"PP"` parameters with a suitable number of
rows. This function would be useful in a Bayesian framework where the
transformation sometimes has to be made on a full sequence of MCMC
iterates, the number of rows being then typically of several
thousands.