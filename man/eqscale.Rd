% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incdist.R
\name{eqscale}
\alias{eqscale}
\title{Transform an incdist object with an equivalence scale}
\usage{
eqscale(object)
}
\arguments{
\item{object}{an incdist object.}
}
\value{
Returns an object where the income variables in frm have been
replaced the equivalised version and the old unequavalised form has been
copied old.frm.
}
\description{
A utility function used by summary.incdist to transform "raw" income
variables by an equivalence scale
}
\details{
Presuposes the present in object of a list called eqscale, with (at least)
components formula, coef. Variables in formula must be found in the
object\$frm.
}
\examples{

x1 <- runif(50)
x2 <- rexp(50)
ad <- ifelse(ad <- rpois(50, 1), ad, 1)
ch <- rpois(50, 2)
## sum to total income
y1 <- x1 + x2
## generate a grouping
g <- factor(c(rep(1, 10), rep(2, 10), rep(1, 15), rep(2, 15)),
            labels = letters[1:2]) # 2 groups
## generate a partitioning variable
p <- factor(c(rep(1, 20), rep(2, 30)), labels = LETTERS[1:2])
## generate some weights
w <- rpois(50, 5)
## put it all into a data frame
test.d <- data.frame(x1, x2, y1, g, p, ad, ch, w)

id.0 <- incdist(y1 ~ x1 +  x2, part = ~ p, group = ~ g, weights = w,
data = test.d)
id.0$idvarnames <- c("a", "c")
id.0$eqscale <- list(formula = ~ (ad + k*ch)^d, coef=list(k=.8, d=.75),
             type="citro-michael")

eqscale(id.0)

}
\author{
Markus Jantti \email{markus.jantti@iki.fi}
}
