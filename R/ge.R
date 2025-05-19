## the generalized entropy class of index
## ranked does not really do anything here...
# check the special cases


#' Generalized entropy inequality index
#'
#' This function estimates members of the Generalized Entropy family of
#' inequality indices.
#'
#' If \eqn{\eta = 0}, ge() estimates Theils 0 index, for \eqn{\eta = 1} ge()
#' estimate Thiels log mean difference. If \eqn{\eta = 2}, ge is 1/2 of the
#' square coefficient of variation, \eqn{\frac{\sigma^2_x}{2\mu_x}}.
#'
#' @aliases ge ge0 ge1 theil cv2
#' @param x a numerical vector whose index is to be estimated
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the inequality aversion parameter, \eqn{\eta \ge 0}.  Larger
#' values place more weight on the inequality of high incomes.
#' @param data the data frame where the variables are found.
#' @param na.rm a logical indicating whether NA's should removed.
#' @return The inequality index.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{gini}} \code{\link{atkinson}}
#' @references
#' \insertRef{lambert1993}{incdist}
#'
#' @examples
#'
#' ge(runif(100), eta = 2)
#'
#' @export ge
ge <- function(x, w = rep(1,length(x)), eta = 2,
               data = NULL,  ranked=x, na.rm = TRUE,...){
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
  if(eta == 1) {
    # check how this business works
    Th <- x/weighted_mean(x,w)
    Th <- w %*% (x * log(Th))
    retval <- Th/(w %*% x)
  }
  else if (eta == 0){
    retval <- log(weighted_mean(x,w)) - weighted_mean(log(x), w)
  }
  else {
   edei <- (x/weighted_mean(x,w))^eta
   retval <- weighted_mean(edei - 1,w)/(eta * (eta - 1))
  }
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}
#' @export cv2
cv2 <-
    function(x, w = rep(1,length(x)), data = list(),
             ranked=x, na.rm = TRUE, ...)
  2*ge(x, w, eta = 2, data,  ranked, na.rm = na.rm, ...)
#' @export ge0
ge0 <-
    function(...)
        ge(eta=0, ...)
#' @export ge1
ge1 <-
    function(...)
    ge(eta=1, ...)
