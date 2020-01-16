
## Improvements in the existing functions and doc

- Warn about AIC and BIC. In the general context of heterogeneous
  data, the concept of *number of observations* can be misleading,
  because a single "MAX" observation can have more weight tand dozens
  of other observations. So the AIC and BIC criterion should not be
  used when there is some historical information. A possible solution
  is setting the number of observation to `NA` when there is some
  historical information.

## New functions or methods

- Write functions for the classical diagnostics in *threshold
  choice*. These should include the trace of the estimated shape for a
  varying threshold, along with (profile-likelihood) conficence
  intervals.

- Write the `confint` method for the class `poisGPList` and a
  subsequent `autoplot` method for the results. These should be
  similar to what exists for linear models in **nlme** package.

- Define generalised residuals for the context. Note that these should
  not necessarily be attached to the observations because an OTS block
  with no observation deserves a residual and can have a huge weight
  in the estimation. So some theoretical work is required here.


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
  `Garonne`, `venice`, and more.