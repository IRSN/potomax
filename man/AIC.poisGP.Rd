% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poisGP.R
\name{AIC.poisGP}
\alias{AIC.poisGP}
\alias{BIC.poisGP}
\title{Akaike's Information Criterion and Bayesian Information
Criterion for a Poisson-GP Object}
\usage{
\method{AIC}{poisGP}(object, ..., k)

\method{BIC}{poisGP}(object, ...)

\method{BIC}{poisGP}(object, ...)
}
\arguments{
\item{object}{An object with class \code{"poisGP"}.}

\item{...}{Not used yet.}

\item{k}{See \code{\link[stats]{AIC}}.}
}
\value{
Value of the criterion
}
\description{
Akaike's Information Criterion and Schwarz's Bayesian Information
Criterion for a Poisson-GP model object.
}
\note{
For technical reasons these methods do not have a
\code{type} argument as does \code{\link{logLik.poisGP}} and
consequently for a model fitted from block maxima or
\eqn{r}-largest the computed criteria will differ from those
computed by other packages because the log-likelihoods differ by a
constant.
}
\section{Caution}{
 Comparing log-likelihoods or related indicators
such AIC, BIC across R packages can be irrelevant due to the use
of different constants. Moreover the concept of \emph{number of
observations} is unclear when heterogeneous data are used. An
historical MAX or OTS information can have a very strong influence
on the estimation hence can not be compared to an ordinary OT
observation.
}

\seealso{
\code{\link{logLik.poisGP}}
}
