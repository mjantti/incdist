#' Decomposition of inequality index by population subgroups
#'
#' This function decomposes generalized entroy (GE) inequality
#' indices by mutually exclusive population subgroup into a
#' within- and between-group term.
#'
#' None.
#'
#' @param x a numerical vector whose index is to be estimated, consisting of
#' income less some reference value
#' @param w an optional vector of non-negative integer values weights.
#' @param eta the inequality aversion parameter. Larger value indicate higher
#' aversion, zero inequality neutrality and negative values inequality
#' preference.
#' @param data the data frame where the variables are found.
#' @param ranked the variable by which observations should be ranked (has no meaning for GE based within/between rankings)
#' @param na.rm a logical indicating whether NA's should removed.
#' @param group the (one-way) grouping varaiable
#' @param func the inequality index
#' @return a vector with overall inequality, the within and the between term
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{ge}}
#'  @references Lambert, P. (1993).  \emph{The distribution and redistribution
#' of income. A mathematical analysis.} Manchester University Press,
#' Manchester.
#' @examples
#'
#' td <- data.frame(y=exp(rnorm(500)), g=sample(letters[1:4], 500, replace=TRUE))
#' decomp_groups(y, group=g, data=td)
#' @export
decomp_groups <- function(x, w = rep(1,length(x)), eta = 2,
               data = NULL,  ranked=x, na.rm = TRUE, group=NULL, func="ge",
                          ...) {
      if(!is.null(data) & is.data.frame(data)) {
          attach(data)
          on.exit(detach(data))
      }
      ind <-
        do.call(func, list(x=x, w=w, eta=eta, na.rm=na.rm))
      ## now collapse
      xb <- tapply(x, list(g=group), mean, na.rm=TRUE)
      wb <- tapply(w, list(g=group), sum, na.rm=TRUE)
      ind_b <- do.call(func, list(x=xb, w=wb, eta=eta, na.rm=TRUE))
      ind_w <- ind - ind_b
      retval <- c(ind, ind_w, ind_b)
      retval
  }

## by income components (for Gini)
## CV^2
