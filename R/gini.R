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
#'     eta = 2, gini() yields the ordinary gini coefficient.(Makes little sense to type = "absolute".)
#' @param na.rm a logical indicating whether NA's should removed.
#' @param data the data frame where the variables are found.
#' @param ranked an optional variable whch defines the ranking. If not
#'     the default, the resulting index is a concentration
#'     coefficient.
#' @param no.negatives indicates whether or not negative incomes should be
#' allowed.
#' @param type whether to estimate an absolute or relative gini.
#' @param absolute alternative for type = "absolute".
#' @param small if set to TRUE estimates the gini coefficient by integrating a  micro-data-based Lorenz curve. (This ensures the resulting gini respects the Lorenz order)
#' @return The generalized gini or concentration coefficient.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{inequality}}, \code{\link{incdist}}, \code{\link{lorenz}}
#' @references
#' \insertRef{lambert1993}{incdist} \insertRef{davidson2009}{incdist}
#'
#' @examples
#'
#' gini(runif(20), eta = 2)
#' gini(runif(20), rpois(20, 5), eta = 2)
#'
#' @export
gini <- function(x, ...) {
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("gini", x)
}
#' @export
gini.incdist <- function(x, ...) {
    inequality.incdist(x, type = "gini", ...)
  }

## this one now seems OK.
#' @export
gini.default <-
    function(x, w = rep(1,length(x)), eta = 2, na.rm = TRUE, data = NULL, ranked = x, no.negatives = FALSE, type="ord", absolute=ifelse(type=="absolute", TRUE, FALSE), small = FALSE,  ...) {
        ## attach a possible data frame. Remember to detach it!
        if (!is.null(data) & is.data.frame(data)) {
            attach(data)
            on.exit(detach(data))
        }
        incmat <- clean_income(cbind(x, ranked), w, no.negatives = no.negatives, na.rm = na.rm)
        x <- incmat[,1]
        ranked <- incmat[,2]
        w <- incmat[,3]
        n <- length(x)
        ind <- order(ranked)
        ranked <- ranked[ind]
        x <- x[ind]
        w <- w[ind]
        mu <- weighted_mean(ranked, w)
        ## NB: here are some details to work out, esp wrt to Davidson 2009 eq 5
        w.tmp <- sum(w) - cumsum(w) + 1
        ## or maybe I need to know here if weights are probs or freq weights?
        ## NB: since I've reversed the weights, the correction is the same as in Davidson?
        r <- (c(w.tmp[1:(length(w.tmp)-1)], 0) + 1/2*w)/sum(w)
        if (small & eta == 2) { ## used to be n <= 100
            ## used to be the below. Now use a Lorenz curve estimated on microdata and integrate it
            ## ret <- 1 + 1/sum(w) - 2*sum(w * w.tmp * x)/(sum(w)^2*mu)
            lorcurve <- incdist::lorenz.default(x, w, q=FALSE)
            p <-c(0, lorcurve$p)
            l <- c(0, lorcurve$ordinates)
            n <- length(p)
            dp <- diff(p)
            lhat <- .5*(l[2:n] + l[1:n-1])
            ret <- 1 - 2*sum(lhat * dp)
        }
        if (small & eta != 2) { ## I do not really know if the below inclusion of eta is correct
            ret <- 1 + 1/sum(w) - 2*sum(w * w.tmp^(eta - 1) * x)/(sum(w)^2*mu)
        }
        if (!small) { ## NB: this needs to be thought through
            gini.m <- cov.wt(cbind(x, (1 - r)^(eta-1)), w)
            ret <- 2*gini.m$cov[1,2]/mu
            names(ret) <- NULL
        }
        ## there should probably be a finite sample correction here?
        if(absolute) ret <- ret*mu
        ret
}
#' @export
concentration_coef <- function(x,  w = rep(1,length(x)), eta = 2, na.rm = TRUE, data = NULL, ranked = x, ...)
  gini(x,  w, eta, na.rm, data, ranked)
