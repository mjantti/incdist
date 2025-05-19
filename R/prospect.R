#' Inequality index based on prospect theory
#'
#' This function prospect-theory-based estimates of Atkinson's inequality index for some
#' reference value.
#'
#' None.
#'
#' @aliases prospect
#' @param x a numerical vector whose index is to be estimated, consisting of
#' income less some reference value
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the inequality aversion parameter. Larger value indicate higher
#' aversion, zero inequality neutrality and negative values inequality
#' preference.
#' @param data the data frame where the variables are found.
#' @param na.rm a logical indicating whether NA's should removed.
#' @return The inequality index.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{gini}} \code{\link{ge}}  \code{\link{atkinson}}
#' @references
#' \insertRef{janttikanburnyyssolaandpirttila2014}{incdist}
#'
#' @examples
#'
#' prospect(runif(100), eta = 1/2)
#'
#' @export prospect
prospect <-
  function(x, w = rep(1,length(x)), eta = 1, data = NULL, na.rm = TRUE){
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
  if (eta == 1){
      geom.mean <- function(x, w){
        sumw <- sum(w)
        gm <- prod(w * x)^(1/sumw)
        gm}
      edei <- geom.mean(x, w)
    }
  else
    edei <- weighted_mean(x^(1-eta),w)
  m <- weighted_mean(x ,w)
  retval <- 1 - edei^(1/(1-eta))/m
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}
