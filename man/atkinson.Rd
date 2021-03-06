% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atkinson.R
\name{atkinson}
\alias{atkinson}
\title{Atkinson's inequality index}
\usage{
atkinson(
  x,
  w = rep(1, length(x)),
  eta = 1/2,
  data = NULL,
  na.rm = TRUE,
  no.negatives = FALSE,
  ...
)
}
\arguments{
\item{x}{a numerical vector whose index is to be estimated}

\item{w}{an optional vector of non-negative integer values weights.}

\item{eta}{the inequality aversion parameter. Larger value indicate higher
aversion, zero inequality neutrality and negative values inequality
preference.}

\item{data}{the data frame where the variables are found.}

\item{na.rm}{a logical indicating whether NA's should removed.}

\item{no.negatives}{a logical indicating wether to exclude negative values
of x}
}
\value{
The inequality index.
}
\description{
This function estimates Atkinson's inequality index.
}
\details{
None.
}
\examples{

atkinson(runif(100), eta = 1/2)

}
\references{
\insertRef{lambert1993}{incdist}
}
\seealso{
\code{\link{gini}} \code{\link{ge}} \code{\link{prospect}}
}
\author{
Markus Jantti \email{markus.jantti@iki.fi}
}
