% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/potData.R
\name{RP.potData}
\alias{RP.potData}
\title{Return Periods for a \code{potData} Object}
\usage{
\method{RP}{potData}(object, points = c("p", "H"), a = 0.5, ...)
}
\arguments{
\item{object}{A \code{potData} object, usually created by using
\code{\link{potData}} or by coercing an object of class
\code{"Rendata"}.}

\item{points}{Type of plotting positions to use: \code{ppoints} or
code Nelson's positions.}

\item{a}{Passed to \code{\link{ppoints}} when the value
\code{points} is \code{"p"}. Ignored else.}

\item{...}{Not used yet.}
}
\value{
A list with several elements. The element \code{data} is a
    data frame which contains so-called \emph{plotting positions}
    as required to display the observations as points on en
    empirical return level plot. The columns of \code{data}
    include \code{T} giving the return period in years, \code{S}
    giving the survival or probability of exceedance and \code{x}
    return level. The columns \code{OT} and \code{source} keep
    trace of the blocks from which the observation was extracted.
}
\description{
Return Periods for a \code{potData} object.
}
\details{
This is a refactoring of the \code{\link[Renext]{SandT}} function
of the \strong{Renext} package.
}
\examples{
pdat <- potData(data = Garonne$OTdata$Flow,
                effDuration = 65,
                MAX.data = Garonne$MAXdata$Flow,
                MAX.effDuration = 143)
L <- RP(pdat)

}
\references{
Chap. 4 in Yves Deville(2020) \emph{Renext Computing Details}. Tech. Report.
}
\seealso{
\code{\link{potData}}, the \code{\link[Renext]{SandT}}
function of \strong{Renext}.
}
