#' Atkinson's inequality index
#'
#' This function estimates Atkinson's inequality index.
#'
#' None.
#'
#' @aliases atkinson ede
#' @param x a numerical vector whose index is to be estimated
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the inequality aversion parameter. Larger value indicate higher
#' aversion, zero inequality neutrality and negative values inequality
#' preference.
#' @param data the data frame where the variables are found.
#' @param na.rm a logical indicating whether NA's should removed.
#' @param no.negatives a logical indicating wether to exclude negative values
#' of x
#' @return The inequality index.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{gini}} \code{\link{ge}}
#' @references
#' \insertRef{lambert1993}{incdist}
#'
#' @examples
#'
#' atkinson(runif(100), eta = 1/2)
#'
#' @export atkinson
atkinson <-
  function(x, w = rep(1,length(x)), eta = 1/2, data = NULL, na.rm = TRUE,
            no.negatives = FALSE, ...){
  ## attach a possible data frame. Remember to detach it!
  if(!is.null(data) & is.data.frame(data))
      {
          attach(data)
          on.exit(detach(data))
      }
  # moved treatment of NA's, missing values and others to utility function
  # "clean_income"
  incmat <- clean_income(x, w, no.negatives = no.negatives, na.rm = na.rm)
  x <- incmat[,1]
  w <- incmat[,2]
  #retval <- list()
  m <- weighted_mean(x,w)
  if (eta == 1)
    {
      ## geom.mean <- function(x, w){
      ##   sumw <- sum(w)
      ##   gm <- prod(w * x)^(1/sumw)
      ##   gm}
      ## edei <- geom.mean(x, w)
      edei <- exp(weighted_mean(log(x), w))
      retval <- 1 - edei/m
    }
  else
    {
      edei <- weighted_mean(x^(1-eta),w)
      retval <- 1 - edei^(1/(1-eta))/m
    }
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}
