% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mobility.R
\name{mob_mat}
\alias{mob_mat}
\alias{mtr}
\alias{mdet}
\alias{ml}
\alias{mf}
\alias{mb}
\alias{mmtr}
\alias{mmdet}
\alias{mml}
\alias{mmf}
\alias{nnb}
\title{A mobility matrix}
\usage{
mob_mat(x1, x2, w, Q1 = NULL, Q2 = NULL, q = 5, labels = FALSE)
}
\arguments{
\item{x1}{A vector of income data}

\item{x2}{A vector of income data}

\item{w}{A vector of weights}

\item{Q1}{Vector of cut offs to be applied to x1}

\item{Q2}{Vector of cut offs to be applied to x2}

\item{q}{Number of classes}

\item{data}{A data frame in which to look for the data}

\item{na.rm}{A logical, indicating whether or not to remove NAs}

\item{index}{An inequality index.}
}
\value{
A list including the conditional and unconditional matrices (EDIT).
}
\description{
This function a mobility matrix for two-dimensional data.
}
\references{
Cite something suitable. (EDIT)
}
\author{
Markus Jantti \email{markus.jantti@iki.fi}
}
