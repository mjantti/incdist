# poverty measurement



#' The Foster Greer Thorbecke poverty index
#'
#' This function estimates the Foster Greer Thorbecke poverty index for a
#' resource variable x and optional weights w.
#'
#' This function estimates the Foster-Greer-Thorbecke poverty index for a
#' poverty line z. The poverty line z is taken by default to be 1/2 of the
#' weighted median, but can be any fraction of the median or mean or a value
#' given as an argument. For \eqn{\alpha=0}, FGT is the proportion of units
#' with income below z \eqn{\Pr(x < z)}. For \eqn{\alpha=1}, FGT equals
#' \eqn{\Pr(x < z)(1 - \bar x/z )} and for \eqn{\alpha = 2} FGT equals
#' \eqn{\Pr(x < z)(1 - \bar x/z )CV^2}, where \eqn{CV^2} is the squared
#' coefficient of variation.
#'
#' The poverty line can be set explicitly or can be taken to be a fraction of
#' either the median (default) or the mean.
#'
#' @param x the resource variable
#' @param w the sampling weights
#' @param alpha the FGT parameter
#' @param mean a logical indicating whether the default poverty line relate to
#' the mean rather than the median (which is the default)
#' @param fraction the fraction of the mean or median to use as the poverty
#' line (default is 0.5)
#' @param z the poverty line (defaults to the product of median and fraction))
#' @param data a data frame holding the data that contain x and w
#' @param na.rm a logical indicating whether NA's should removed
#' @return A list with elements \item{fgt}{the index value} \item{n}{sample
#' size} \item{z}{poverty line} \item{sumq}{sum of weights}
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references
#' @examples
#'
#' fgt(runif(100))
#' n <- 100
#' income <- exp(10 + rnorm(n)*sqrt(2))
#' weight <- rpois(n,7)
#' fgt(income, weight)
#' sapply(seq(.3, .6, by = .1), function(x) fgt(income, weight, fraction = x)$fgt)
#'
#' @export fgt
fgt <-
  function (x, w = rep(1, length(x)), alpha = 0,
            mean = FALSE, fraction = 0.5,
            z = ifelse(mean, fraction * weighted.mean(x, w),
              fraction * weighted.median(x, w)),
            data = NULL,
            ranked = NULL,
            na.rm = TRUE)
{
   ## attach a possible data frame. Remember to detach it!
  if(!is.null(data) & is.data.frame(data))
      {
          attach(data)
          on.exit(detach(data))
      }
  if (fraction < 0 || fraction > 1)
    stop("Fraction outside unit interval!")
  # moved treatment of NA's, missing values and others to utility function
  # "clean.income"
  incmat <- clean.income(x, w)
  x <- incmat[, 1]
  w <- incmat[, 2]
  retval <- list()
#  if (missing(z))
#  {
#    if (mean)
#      z <- fraction * weighted.mean(x, w)
#    else
#      z <- fraction * weighted.median(x, w)
#  }
  fgt <- weighted.mean((x < z) * (1 - x/z)^alpha, w)
  retval$fgt <- fgt
  retval$n <- length(x)
#  if(length(z) != 1) warning("Poverty line is a vector.")
  retval$z <- z
  retval$sumw <- sum(w)
  #structure(retval, class = "poverty")
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}

fgt0 <- function(...)
    {
        fgt(...)$fgt
    }

fgt1 <- function(...)
    {
        fgt(alpha=1, ...)$fgt
    }

fgt2 <- function(...)
    {
        fgt(alpha=2, ...)$fgt
    }
