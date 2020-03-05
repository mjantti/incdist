#
integrate.locfit <- function(x, y, ...){
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

## why on earth is this here?
"trace.plot" <-
    function(x, y, g, data = sys.parent(),
             add = FALSE, col = "black", ...){
        g <- factor(g)
        dfr <- na.omit(data.frame(x=x,y=y,g=g))
        dfr <- dfr[order(dfr$g,dfr$x),]
        if(add == FALSE) plot(x, y, type="n", ...)
        for ( i in levels(g) )
          # lty = i too many line types make postscript pictures bad
          # col = i colors don't print well
          # try lty = i %/% 2
            evalq(lines(x,y, lty = as.numeric(i), col = col), subset(dfr,g==i))
    }

# for weighted data.

# weighted moments (mean, variance, std)
# am using  standard R functions for mean, var and std

#' Weighted statistic
#'
#'
#' The weighted.moment is a wrapper around R's weighted.mean which raises the
#' variable x to the power a before taking the weighted mean.  weighted.var and
#' weighted.std are based on weighted.moment.
#'
#' The weighted.quantile function yields a version of quantiles which does not
#' quite reach the sophistication of that built into R. In particular, for a
#' probability p this function yields the highest x for which cumsum(w)/sum(w)
#' < p, rather than some interpolated value.  weighted.median is
#' weighted.quantile(x,w, p = .5)
#'
#' @aliases weighted.moment weighted.var weighted.std weighted.quantile
#' weighted.median weighted.sum
#' @param x the variable whose moment is to be estimated
#' @param w the weights
#' @param probs the probabilities of the quantiles
#' @param a the degree of the moment
#' @param names no idea what this does
#' @param na.rm a logical indicating whether NA's should removed
#' @return The desired statistic. Weighted quantile produces a vector with
#' length(probs).
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{weighted.mean}}, \code{\link{cov.wt}}
#' @references
#' @examples
#'
#' x <- rexp(100)
#' w <- rpois(100,5)
#' weighted.median(x, w)
#' weighted.moment(x,w, 1)
#' weighted.mean(x,w, 1)
#' weighted.var(x,w)
#'
#'

weighted.moment <- function(x,w,a=1, ranked=x, na.rm = FALSE) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    i <- !is.na(x) & !is.na(w)
    w <- w[i]
    x <- x[i]
  }
  weighted.mean(x^a,w, na.rm = na.rm)
}

weighted.var <- function(x, w, ranked=x, na.rm = FALSE)
  {
    if(missing(w))
      {
        if(is.vector(x)) as.vector(cov.wt(as.matrix(x))$cov)
        else diag(cov.wt(as.matrix(x))$cov)
      }
    else
      {
        if (na.rm) {
          i <- !is.na(x) & !is.na(w)
          w <- w[i]
          x <- x[i]
        }
        if(is.vector(x)) as.vector(cov.wt(as.matrix(x),w)$cov)
        else diag(cov.wt(as.matrix(x),w)$cov)
      }
  }

weighted.std <- function(x, w, ranked=x, na.rm = FALSE)
  sqrt(weighted.var(x,w, na.rm = na.rm))


#
# weighted quantiles. various versions, none good.

weighted.quantile.2 <-
  function(x, w = rep(1, length(x)),
           probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE) {
    # n <- length(x)
    # check:
    # - if w exists (length(w)>0) => if not, then w <- rep(1,length(x))
    # - if w is the weight or the inclusion probability
    #   + a test sum(w) > length(x)??????
    if (length(x) != length(w))
      stop("Weights and variable vectors are of unequal length!")
    if (na.rm)
      { valid.obs <- !is.na(x) & !is.na(w)
        x <- x[valid.obs]
        w <- w[valid.obs]}
    else if (any(is.na(x)) | any(is.na(w)))
      stop("Missing values and NaN's not allowed if `na.rm' is FALSE")
    if (any((p.ok <- !is.na(probs)) & (probs < 0 | probs > 1)))
      stop("probs outside [0,1]")
    if (na.p <- any(!p.ok)) {
      o.pr <- probs
      probs <- probs[p.ok]
    }
    np <- length(probs)
    if (missing(w))
      w <- rep(1, length(x))
    if (na.rm) {
      w <- w[i <- !is.na(x)]
      x <- x[i]
    }
    # check if weights are probs or unit correspondences
    ind <- order(x)
    cumw <- cumsum(w[ind])/sum(w)
    x <- x[ind]
    med <- numeric(0)
    for(i in 1:length(probs)) {
      if(probs[i] > 0 & probs[i] < 1) med[i] <- max(x[cumw<=probs[i]])
      else if (probs[i] == 0) med[i] <- min(x)
      else if (probs[i] == 1) med[i] <- max(x)
    }
    # the treatment of names
    if (names && np > 0) {
      dig <- max(2, getOption("digits"))
      names(med) <- paste(if (np < 100)
                         formatC(100 * probs, format = "fg", wid = 1, dig = dig)
      else format(100 * probs, trim = TRUE, dig = dig), "%",
                         sep = "")
    }
    if (na.p) {
      o.pr[p.ok] <- med
      names(o.pr)[p.ok] <- names(med)
      o.pr
    }
    else med
  }

# this version uses R's quantile and repeats each obs by w.
# very inefficient!
#
# this is a version that departs from the unweighted quantile funtion in R
# incomplete

weighted.quantile.defunct <-
  function (x, w = rep(1,length(x)), ranked=x,
                 probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE)
{
  # if no weights are given, use the built-in function
  if (missing(w)) quantile.default(x, probs, na.rm, names)
  else if (length(x) != length(w))
    stop("Weights and variable vectors are of unequal length!")
  if (na.rm)
    { valid.obs <- !is.na(x) & !is.na(w)
      x <- x[valid.obs]
      w <- w[valid.obs]}
  else if (any(is.na(x)) | any(is.na(w)))
    stop("Missing values and NaN's not allowed if `na.rm' is FALSE")
  if (any((p.ok <- !is.na(probs)) & (probs < 0 | probs > 1)))
    stop("probs outside [0,1]")
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
  }
  np <- length(probs)
  n <- length(x)
  if (n > 0 && np > 0) {
    wsum <- sum(w)
    index <- 1 + (wsum - 1) * probs
    lo <- floor(index)
    hi <- ceiling(index)
    i.sorted <- order(x)
    x <- x[i.sorted]
    w <- w[i.sorted]
    cw <- cumsum(w)
    i <- cw[lo]
    qs <- x[i]
    i <- seq(along = i)[i & !is.na(i)][qs[i] > -Inf]
    .minus <- function(x, y) ifelse(x == y, 0, x - y)
    qs[i] <- qs[i] + .minus(x[hi[i]], x[lo[i]]) * (index[i] - lo[i])
  }
  else {
    qs <- rep(as.numeric(NA), np)
  }
 qs
}

weighted.quantile.1 <-
  function(x, w, probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE)
{
  if(missing(w))
    quantile.default(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE)
  else
    {
      # this is *incredibly* inefficient!
      # make w an integer, if not already
      # should test if any of the weights are not integers first
      # is.integer(w) does not work
      w <- round(w)
      # generate the population level x
      x <- rep(x,w)
      # call quantile.default
      quantile.default(x, probs, na.rm = FALSE, names = TRUE) }
}

weighted.quantile <-
  function(x, w = rep(1, length(x)),
           probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE)
  {
    ## call wtd.quantile
    wtd.quantile(x, w, probs, na.rm=na.rm)
    ##weighted.quantile.1(x, w, probs, na.rm=na.rm)
  }

# define this last.

weighted.median <- function(x, w, ranked=x, na.rm = FALSE, ...) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
    }
  weighted.quantile(x, w, .5, ...)
}

weighted.sum <- function(x, w, ranked=x, na.rm = FALSE)
  {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    i <- !is.na(x) & !is.na(w)
    w <- w[i]
    x <- x[i]
  }
  w %*% x
}


# some additional weight funtions
# stolen from definition of ecdf in library(stepfun)
# do I need to load stepfun first?
weighted.ecdf <-
  function (x, w = rep(1,length(x)), ranked=x)
{
  order.x <- order(x)
  x <- x[order.x]
  w <- w[order.x]
  n <- length(x)
  w.sum <- sum(w)
  rval <- approxfun(x, cumsum(w)/w.sum, method = "constant", yleft = 0,
                    yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  attr(rval, "call") <- sys.call()
  rval
}

weighted.rank <-
  function (x, w = rep(1,length(x)), ranked=x)
{
  if((n <- length(x)) != (n.w <- length(w))) stop("x and w are not of equal length")
  # what happens to ties when usign order?
  order.x <- order(x)
  x <- x[order.x]
  w <- w[order.x]
  n <- length(x)
  # need to offset this by one to get r_j = \sum_{i < j} w_i + 1/2*w_j
  tmp.w <- c(0,cumsum(w)[1:(n-1)])
  rank <- 1/sum(w)*(1/2*w + tmp.w)
  rank
}

weighted.table.mean <-
  function (x1, w1, l, ...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted.mean(x[i], w[i],...))
}

weighted.table.old <-
  function (x1, w1, l, ...) weighted.table.mean(x1=x1, w1=w1, l=l, ...)

weighted.table.sum <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted.sum(x[i], w[i],...))
}

weighted.table.median <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted.median(x[i], w[i],...))
}

weighted.table.var <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted.var(x[i], w[i],...))
}
weighted.table <-
  function (w1, l, na.rm=TRUE, relative=FALSE, ...)
{
  tm <- tapply(seq(along = w1), l, function(i, w = w1) sum(w[i],...))
  if(relative=="rows")
    tm <- tm/matrix(rep(rowSums(tm, na.rm=na.rm), ncol(tm)), nrow=nrow(tm))
  else if(relative=="cols")
    tm <- tm/matrix(rep(colSums(tm, na.rm=na.rm), nrow(tm)), ncol=ncol(tm))
  tm
}

##

weighted.crosstable <- function(x1, x2, w1 = rep(1, length(x1)))
{

  ## error checks: are factors, do weights exist, at least
  y <- tapply(seq(along=w1), list(x1, x2),
              function(i, w = w1) sum(w[i]))
  class(y) <- "table"
  y
}

weighted.crosstable.old <- function(l, w1 = rep(1, length(l[[1]])))
{

  ## error checks: are factors, do weights exist, at least
  y <- tapply(seq(along=w1), l,
              function(i, w = w1) sum(w[i]))
  class(y) <- "table"
  y
}

# missing and NA treatment

## write a new  function to handle x matrices


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
#' @references
#' @examples
#'
#' x <- c(NA, 1:5, Inf, -1, NaN, -Inf)
#' w <- rpois(length(x), 2)
#' clean.income(x, w)
#' x <- cbind(x, y = rpois(length(x), 10))
#' clean.income(x, w)
#'
#' @export clean.income
clean.income <- function(x, w = rep(1,length(x)),
                         no.negatives = FALSE,
                         no.nans = TRUE,
                         no.infinites = TRUE,
                         na.rm = TRUE)
  {
    if(is.matrix(x)) k <- ncol(x)
    else
      {
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
    if(no.negatives)
      {
        if(any(colMeans(df[,-1, drop=FALSE], na.rm=TRUE)>0) && any(df<0))
           {
             negs <- apply(df < 0, 2, sum)
             df <- replace(df, df < 0, NA)
           }
        if(any(colMeans(df[,-1, drop=FALSE], na.rm=TRUE)<0) && any(df>0))
           {
             negs <- apply(df > 0, 2, sum)
             df <- replace(df, df > 0, NA)
           }
      }
    if(no.infinites)
      {
        infs <- apply(is.infinite(df), 2, sum)
        df <- replace(df, is.infinite(df), NA)
      }
    if(no.nans)
      {
        nans <- apply(is.nan(df), 2, sum)
        df <- replace(df, is.nan(df), NA)
      }
    if(na.rm)
      {
        ind <- complete.cases(df)
        df <- df[ind,, drop = FALSE]
      }
    x <- df[,-1, drop=FALSE]
    w <- df[,1]
    cbind(x, w)
    ##df
  }

## two sample test

distdiff <- function(df1, df2, income, weight=NULL, stat, B=1000)
    {
        n1 <- dim(df1)[1]
        n2 <- dim(df2)[1]
        n <- n1 + n2
        df <- rbind(df1, df2)
        res <- vector(mode="numeric", length=B)
        for(b in 1:B)
            {
                ndf1 <- df[sample(1:n, n1, replace=TRUE),]
                ndf2 <- df[sample(1:n, n2, replace=TRUE),]
                if(!is.null(weight))
                    {
                        wgt1 <- ndf1[[weight]]
                        wgt2 <- ndf2[[weight]]
                    }
                else
                    {
                        wgt1 <- rep(1, n1)
                        wgt2 <- rep(1, n2)
                    }
        res[b] <-
            do.call(stat, list(x=ndf1[[income]], w=wgt1))-
                do.call(stat, list(x=ndf2[[income]], w=wgt2))
            }
        res
    }

## I have no idea why these are here!
## d1 <- data.frame(y=exp(rnorm(100)), wgt=1+rpois(100, 4))
## d2 <- data.frame(y=exp(1+2*rnorm(150)), wgt=1+rpois(150, 5))

## gini.default(d1$y, d1$wgt)
## gini.default(d2$y, d2$wgt)
## gini.default(d1$y, d1$wgt)-gini.default(d2$y, d2$wgt)
## dg <- distdiff(d1, d2, "y", weight="wgt", "gini.default", B=10000)
## plot(density(dg))
## quantile(dg, c(.025, .975))

## the influence function a la Cowell & Flachaire

ineqvar <- function(x, w, ranked=x, na.rm = FALSE, index="gini",...){
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
    mu <- weighted.mean(x, w)
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
