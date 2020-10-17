#' Generalized gini or concentration coefficient
#'
#' This function estimates the generalized Gini coefficient or, if an
#' other variable is used for ordering, the concentration coefficient
#' of a variable with respect to the ordering variable. For the
#' generalizing parameter eta = 2, this yields the ordinary Gini
#' coefficient.
#'
#' The function gini.default works on vectors, whereas the
#' gini.incdist function works on incdist objects (and currently
#' produces lots of output, much of which is unlikely to be useful.)
#'
#' @aliases gini gini.default gini.incdist concentration_coefficient
#' @param x a numerical vector whose index is to be estimated.
#' @param object an incdist object.
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the generalization parameter. For the default value of
#'     eta = 2, gini() yields the ordinary gini coefficient.
#' @param na.rm a logical indicating whether NA's should removed.
#' @param data the data frame where the variables are found.
#' @param ranked an optional variable whch defines the ranking. If not
#'     the default, the resulting index is a concentration
#'     coefficient.
#' @return The generalized gini or concentration coefficient.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{inequality}}, \code{\link{incdist}}, \code{\link{lorenz}}
#' @references
#' \insertRef{lambert1993}{incdist}
#'
#' @examples
#'
#' gini(runif(20), eta = 2)
#' gini(runif(20), rpois(20, 5), eta = 2)
#'
#' @export
gini <- function(x, ...)
  {
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("gini", x)
}
#' @export
gini.incdist <- function(x, ...)
  {
    inequality.incdist(x, type = "gini", ...)
  }

                                        # this one now seems OK.
#' @export
gini.default <-
  function(x, w = rep(1,length(x)), eta = 2, na.rm = TRUE,
           data = NULL, ranked = x, no.negatives = FALSE,
           type="ord",
           absolute=ifelse(type=="absolute", TRUE, FALSE), ...)
{
  ## attach a possible data frame. Remember to detach it!
  if(!is.null(data) & is.data.frame(data))
      {
          attach(data)
          on.exit(detach(data))
      }
  ## moved treatment of NA's, missing values and others to utility function
  ## "clean_income"
  incmat <- clean_income(cbind(x, ranked), w,
                         no.negatives = no.negatives, na.rm = na.rm)
  x <- incmat[,1]
  ranked <- incmat[,2]
  w <- incmat[,3]
  n <- length(x)
  ## check if weights are probs or unit correspondences
  ## not implemented and not really needed now.
  ind <- order(ranked)
  ranked <- ranked[ind]
  x <- x[ind]
  w <- w[ind]
  w.tmp <- sum(w) - cumsum(w) + 1
  mu <- weighted_mean(ranked, w)
  ## or maybe I need to know here if weights are probs or freq weights?
  r <- (c(w.tmp[1:(length(w.tmp)-1)], 0) + 1/2*w)/sum(w)
  if(n > 100)
    {
      gini.m <- cov.wt(cbind(x, (1 - r)^(eta-1)),w)
      ## does a concentration coefficient divide by the mean of x or
      ## ranked? by ranked!
      ret <- 2*gini.m$cov[1,2]/mu
      names(ret) <- NULL
    }
  else
    {
      ret <- 1 + 1/sum(w) - 2*sum(w * w.tmp * x)/(sum(w)^2*mu)
    }
  ## there should probably be a finite sample correction here?
  if(absolute) ret <- ret*mu
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##  detach(data)
  ret
}
#' @export
concentration_coef <- function(x,  w = rep(1,length(x)), eta = 2, na.rm = TRUE,
           data = NULL, ranked = x, ...)
  gini(x,  w, eta, na.rm, data, ranked)
