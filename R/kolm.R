#' Kolm's inequality index
#'
#' This function estimates the inequality indices (relative or absolute) of
#' Kolm.
#'
#' The inequality index.
#'
#' @param x a numerical vector whose index is to be estimated
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the inequality aversion parameter (>0). Larger value indicate
#' higher aversion, zero inequality neutrality.
#' @param data the data frame where the variables are found.
#' @param na.rm a logical indicating whether NA's should removed.
#' @param absolute a logical, indicating if the absolute rather than relative
#' measure be given.
#' @return The index value
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{gini}} \code{\link{ge}}, \code{\link{atkinson}}
#' @references
#' \insertRef{lambert1993}{incdist}
#'
#' @examples
#'
#' kolm(runif(100), eta = 1)
#'
#' @export
kolm <-
  function(x, w = rep(1,length(x)), eta = 1, data = NULL,
           absolute = FALSE, na.rm = TRUE)
{
  ## attach a possible data frame. Remember to detach it!
  if(!is.null(data) & is.data.frame(data)) {
          attach(data)
          on.exit(detach(data))
      }
  # moved treatment of NA's, missing values and others to utility function
  # "clean_income"
  incmat <- clean_income(x, w)
  x <- incmat[,1]
  w <- incmat[,2]
  #retval <- list()
  m <- weighted_mean(x,w)
  km <- eta*(m - x)
  km <- weighted_mean(exp(km), w)
  retval <- 1/eta*log(km)
  if(absolute) retval <- m - retval*m ## this can not be right
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}
