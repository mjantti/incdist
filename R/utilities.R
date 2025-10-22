#
integrate.locfit <- function(x, y, ...) {
  n <- length(x)
  x <- sort(x)
  dx <- x[2:n] - x[1:n-1]
  xhat <- .5*(x[2:n] + x[1:n-1])
  fhat <- predict(y,xhat)
  t.i <- sum(fhat*dx)
  t.m <- sum(xhat*fhat*dx)
  t.v <- sum(xhat**2*fhat*dx) - t.m**2
  res <- c(n,t.i,t.m,t.v)
  names(res) <- c("Length","Integration", "Mean (lf)", "Variance")
  res
}

integrate_lorenz <- function(x, ...) {
    if(!is.lorenz(x)) stop("Not a lorenz curve object!")
    p <-c(0, x$p)
    l <- c(0, x$ordinates)
    n <- length(p)
    dp <- diff(p)
    lhat <- .5*(l[2:n] + l[1:n-1])
    g <- 1 - 2*sum(lhat * dp)
    g
}

## why on earth is this here?
"trace.plot" <-
    function(x, y, g, data = sys.parent(),
             add = FALSE, col = "black", ...) {
        g <- factor(g)
        dfr <- na.omit(data.frame(x=x,y=y,g=g))
        dfr <- dfr[order(dfr$g,dfr$x),]
        if(add == FALSE) plot(x, y, type="n", ...)
        for ( i in levels(g) )
          # lty = i too many line types make postscript pictures bad
          # col = i colors don't print well
          # try lty = i %/% 2
            evalq(lines(x, y, lty = as.numeric(i), col = col), subset(dfr,g==i))
    }

# for weighted data.

# weighted moments (mean, variance, std)
# am using  standard R functions for mean, var and std

#' Weighted statistic
#'
#'
#' The weighted_moment is a wrapper around R's weighted_mean which raises the
#' variable x to the power a before taking the weighted mean.  weighted_var and
#' weighted_std are based on weighted_moment.
#'
#' The weighted_quantile function yields a version of quantiles which does not
#' quite reach the sophistication of that built into R. In particular, for a
#' probability p this function yields the highest x for which cumsum(w)/sum(w)
#' < p, rather than some interpolated value.  weighted_median is
#' weighted_quantile(x,w, p = .5)
#'
#' @aliases weighted_moment weighted_var weighted_std weighted_quantile
#' weighted_median weighted_sum weighted_moment
#' @param x the variable whose moment is to be estimated
#' @param w the weights
#' @param probs the probabilities of the quantiles
#' @param a the degree of the moment
#' @param names no idea what this does
#' @param na.rm a logical indicating whether NA's should removed
#' @return The desired statistic. Weighted quantile produces a vector with
#' length(probs).
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{weighted_mean}}, \code{\link{cov.wt}}
#' @examples
#'
#' x <- rexp(100)
#' w <- rpois(100,5)
#' weighted_median(x, w)
#' weighted_moment(x,w, 1)
#' weighted_mean(x,w, 1)
#' weighted_var(x,w)
#'
#' @importFrom Hmisc wtd.quantile
#'
#' @export
weighted_moment <- function(x, w, a=1, ranked=x, na.rm = FALSE) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    i <- !is.na(x) & !is.na(w)
    w <- w[i]
    x <- x[i]
  }
  weighted.mean(x^a,w, na.rm = na.rm)
}
#' @export
weighted_mean <- function(x, w, na.rm = FALSE) {
    weighted.mean(x, w, na.rm)
}
#' @export
weighted_var <- function(x, w, ranked=x, na.rm = FALSE) {
    if(missing(w)) {
        if(is.vector(x)) as.vector(cov.wt(as.matrix(x))$cov)
        else diag(cov.wt(as.matrix(x))$cov)
      }
    else {
        if (na.rm) {
          i <- !is.na(x) & !is.na(w)
          w <- w[i]
          x <- x[i]
        }
        if(is.vector(x)) as.vector(cov.wt(as.matrix(x),w)$cov)
        else diag(cov.wt(as.matrix(x),w)$cov)
      }
  }
#' @export
weighted_std <- function(x, w, ranked=x, na.rm = FALSE)
  sqrt(weighted_var(x,w, na.rm = na.rm))


#
# weighted quantiles. various versions, none good.
#' @export
weighted_quantile <-
  function(x, w = rep(1, length(x)),
           probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE)  {
    ## call wtd.quantile
    wtd.quantile(x, w, probs, na.rm=na.rm)
    ##weighted_quantile.1(x, w, probs, na.rm=na.rm)
  }

# define this last.
#' @export
weighted_median <- function(x, w, ranked=x, na.rm = FALSE, ...) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
    }
  weighted_quantile(x, w, .5, ...)
}


#' @export
weighted_crosstable <- function(x1, x2, w1 = rep(1, length(x1))) {
  ## error checks: are factors, do weights exist, at least
  y <- tapply(seq(along=w1), list(x1, x2),
              function(i, w = w1) sum(w[i]))
  class(y) <- "table"
  y
}


#' Utility functions
#'
#' This utility us used to tidy up or provide warnings in the functions that
#' use weights
#'
#'
#' @param x an income vector or matrix
#' @param w the vector of weights
#' @param no.negatives a logical indicating whether negative values in the data
#' should be removed. The default value is false
#' @param no.nans a logical indicating whether NaN:s in the data should be
#' removed. The default value is TRUE
#' @param no.infinites a logical indicating whether infinite values in the data
#' should be removedd. The default value is TRUE
#' @param na.rm A logical indicating whether NA:s should be removed from the
#' data. The default value is TRUE
#' @return A matrix with as many columns equal to the number of columns in x
#' plus one for the weights, and rows with the inappropriate oobservations
#' removed.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{complete.cases}}, \code{\link{na.omit}}
#' @examples
#'
#' x <- c(NA, 1:5, Inf, -1, NaN, -Inf)
#' w <- rpois(length(x), 2)
#' clean_income(x, w)
#' x <- cbind(x, y = rpois(length(x), 10))
#' clean_income(x, w)
#'
#' @export clean_income
clean_income <- function(x, w = rep(1,length(x)),
                         no.negatives = FALSE,
                         no.nans = TRUE,
                         no.infinites = TRUE,
                         na.rm = TRUE) {
    if(is.matrix(x)) k <- ncol(x)
    else {
        k <- 1
        x <- as.matrix(x)
        colnames(x) <- "x"
      }
    ## convert inf:s, nan:s and negatives to NA but record numbers
    ##df <- data.frame(w=w, x=x)
    df <- cbind(w, x)
    nas <-  apply(is.na(df), 2, sum)
    ## replacements
    ## this is a little tricky, since the restriction should
    ## *really* be that either all(x) > 0 or all(x)<0
    ## should this assume x is a matrix or scalar, or be better tuned?
    if(no.negatives) {
        if(any(colMeans(df[,-1, drop=FALSE], na.rm=TRUE)>0) && any(df<0)) {
             negs <- apply(df < 0, 2, sum)
             df <- replace(df, df < 0, NA)
           }
        if(any(colMeans(df[,-1, drop=FALSE], na.rm=TRUE)<0) && any(df>0)) {
             negs <- apply(df > 0, 2, sum)
             df <- replace(df, df > 0, NA)
           }
      }
    if(no.infinites) {
        infs <- apply(is.infinite(df), 2, sum)
        df <- replace(df, is.infinite(df), NA)
      }
    if(no.nans) {
        nans <- apply(is.nan(df), 2, sum)
        df <- replace(df, is.nan(df), NA)
      }
    if(na.rm) {
        ind <- complete.cases(df)
        df <- df[ind,, drop = FALSE]
      }
    x <- df[,-1, drop=FALSE]
    w <- df[,1]
    cbind(x, w)
    ##df
  }



## the influence function a la Cowell & Flachaire

#' Variance of inequality index
#' based on influence
#'
#' @param x an income vector or matrix
#' @param w the vector of weights
#' @param ranked variable by which data are ranked
#' @param na.rm A logical indicating whether NA:s should be removed from the
#' data. The default value is TRUE
#' @param index the index in question, one of gini, ge2, ge0, ge1
#' @references
#' \insertRef{cowellandflachaire2015}{incdist}
#' @export
ineqvar <- function(x, w, ranked=x, na.rm = FALSE, index="gini",...) {
    if (missing(w))
        w <- rep(1, length(x))
    n <- length(x)
    if(index=="gini") {
        i.sorted <- order(x)
        x <- x[i.sorted]
        w <- w[i.sorted]
        i <- cumsum(w)
        }
    ## map all cv2 etc to ge
    ## get value of inequality index
    ineq <- do.call(index, list(x=x, w=w, ...))
    mu <- weighted_mean(x, w)
    ## calculate the Z
    ## Gini (p 406, eq 6.77
    eta <- switch(index, gini=2, ge=2, ge0=0, ge1=1)
    zgini <- -1*(ineq + 1)*x + (2*i - 1)/n*x - 2/n*cumsum(x)
    ## eta != 0,1
    zge <- (eta^2 - eta)^(-1)*(x/mu)^eta - eta*(x/mu)*(ineq + (eta^2 - eta)^(-1))
    ## eta=0
    zge0 <- (x/mu)*log(x)
    ## eta=1
    zge1 <- (x/mu)*(log(x/mu) - ineq - 1)
    zi <- switch(index,
                 gini = zgini,
                 ge = zge,
                 ge0 = zge0,
                 ge1 = zge1)
    switch(index,
           gini=var(zi)/(n*mu^2),
           ge = var(zi),
           ge0 = var(zi),
           ge1 = var(zi))
    }
