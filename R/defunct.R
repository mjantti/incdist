weighted_quantile_2 <-
  function(x, w = rep(1, length(x)),
           probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE) {
    # n <- length(x)
    # check:
    # - if w exists (length(w)>0) => if not, then w <- rep(1,length(x))
    # - if w is the weight or the inclusion probability
    #   + a test sum(w) > length(x)??????
    if (length(x) != length(w))
      stop("Weights and variable vectors are of unequal length!")
    if (na.rm) { valid.obs <- !is.na(x) & !is.na(w)
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

weighted_quantile_defunct <-
  function (x, w = rep(1,length(x)), ranked=x,
                 probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE)
{
  # if no weights are given, use the built-in function
  if (missing(w)) quantile(x, probs, na.rm, names)
  else if (length(x) != length(w))
    stop("Weights and variable vectors are of unequal length!")
  if (na.rm) { valid.obs <- !is.na(x) & !is.na(w)
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
#' @export
weighted_quantile_1 <-
  function(x, w, probs = seq(0, 1, 0.25), ranked=x, na.rm = FALSE, names = TRUE)
{
  if(missing(w))
    quantile(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE)
  else {
      # this is *incredibly* inefficient!
      # make w an integer, if not already
      # should test if any of the weights are not integers first
      # is.integer(w) does not work
      w <- round(w)
      # generate the population level x
      x <- rep(x,w)
      # call quantile.default
      quantile(x, probs, na.rm = FALSE, names = TRUE) }
}

weighted_sum <- function(x, w, ranked=x, na.rm = FALSE) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    i <- !is.na(x) & !is.na(w)
    w <- w[i]
    x <- x[i]
  }
  w %*% x
}


weighted_ecdf <-
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
# some additional weight funtions
# stolen from definition of ecdf in library(stepfun)
# do I need to load stepfun first?





weighted_rank <-
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

weighted_table_mean <-
  function (x1, w1, l, ...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted.mean(x[i], w[i],...))
}

weighted_table_old <-
  function (x1, w1, l, ...) weighted_table.mean(x1=x1, w1=w1, l=l, ...)
#' @export
weighted_table_sum <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted_sum(x[i], w[i],...))
}
#' @export
weighted_table_median <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted_median(x[i], w[i],...))
}
#' @export
weighted_table_var <-
  function (x1, w1, l,...)
{
  tapply(seq(along = x1), l, function(i, x = x1, w = w1)
         weighted_var(x[i], w[i],...))
}
#' @export
weighted_table <-
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

weighted_crosstable_old <- function(l, w1 = rep(1, length(l[[1]])))
{

  ## error checks: are factors, do weights exist, at least
  y <- tapply(seq(along=w1), l,
              function(i, w = w1) sum(w[i]))
  class(y) <- "table"
  y
}

# missing and NA treatment

## write a new  function to handle x matrices

## two sample test

distdiff <- function(df1, df2, income, weight=NULL, stat, B=1000) {
        n1 <- dim(df1)[1]
        n2 <- dim(df2)[1]
        n <- n1 + n2
        df <- rbind(df1, df2)
        res <- vector(mode="numeric", length=B)
        for(b in 1:B) {
                ndf1 <- df[sample(1:n, n1, replace=TRUE),]
                ndf2 <- df[sample(1:n, n2, replace=TRUE),]
                if(!is.null(weight)) {
                        wgt1 <- ndf1[[weight]]
                        wgt2 <- ndf2[[weight]]
                    }
                else {
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
