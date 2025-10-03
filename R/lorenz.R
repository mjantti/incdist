#'
#' A generalized lorenz curve
#'
#' \code{lorenz}, \code{lorenz.default}, \code{lorenz.locfit} and
#' \code{lorenz.incdist} compute a (possibly generalized) Lorenz curve. The
#' default method assumes the data are vectors, the locfit method assumes that
#' the function is applied to a locfit density object, and the incdist method
#' assumes the object has been created by the incdist methods. The curve is
#' estimated at a given number of points or at every point in the data.
#'
#' \code{var.cov} is a utility function, invoked when \code{cov = T} for
#' \code{lorenz.default}, and estimates the variance-covariance matrix of the
#' lorenz ordinates and income shares.
#'
#' print and summary methods are provided, as well as coercion to data fame
#' (can be useful when combined e.g. with lattice to produce good graphics.
#'
#' \code{dominates.lorenz} compares two lorenz curves and returns (1=x L
#' dominates y, 0 = curves cross or are identical, -1 = y L dominates x).
#'
#' \code{dominates.lorenz_list} takes a list of lorenz curves and computes an
#' upper triangular matrix with elements dominates.lorenz(row, column)
#'
#' The function computes an ordinary lorenz curve. The summary and plot methods
#' allow also a generalized (L*mean) be plotted and summarized.  The weights
#' are for the time being assumed to be positive integers, although this will
#' hopefully change soon.
#'
#' Later on this function should know about income distribution objects. I am
#' also considering extending it to be able to handle \code{\link{locfit}},
#' \code{\link{density}} and \code{\link{ecdf}} objects, at least.
#'
#' Simple summary, print and plot methods are provided.
#'
#' There is a substantial literature on calculating standard errors for Lorenz
#' curve ordinates. If \code{cov = TRUE} is chosen and \code{q = FALSE}, the
#' returned object includes the variance matrix of the ordinates and of the
#' income shares, estimated as shown in Beach and Davidson (1983) and Beach and
#' Kaliski (1986).  The R package \code{\link{boot}} is probably quite useful
#' for doing so here, not least because all interesting work involving Lorenz
#' curves concerns the comparison of two (or more) Lorenz curves, making it
#' most often a two-sample problem.
#'
#' @aliases lorenz.default lorenz lorenz.locfit lorenz.incdist is.lorenz
#' summary.lorenz print.summary.lorenz as.data.frame.lorenz dominates.lorenz
#' dominates.lorenzlist variance.lorenz
#' @param x a numerical vector whose Lorenz curve is to be estimated
#' @param w an optional vector of non-negative integer values weights
#' @param obj a locfit density object
#' @param ranked an optional vector of length(x) which defines the ordering. If
#' not x, lorenz returns the concentration curve of x for the variable given in
#' ranked.
#' @param q an integer or logical. If an integer, it is the number of equally
#' spaced points at which the Lorenz curve ordinates are estimated. If FALSE,
#' then the Lorenz curve is estimated at every point in x.
#' @param p the population points at which the Lorenz curve is estimated. Defaults to
#' equally spaced point on the unit interval. If given, then q is set to the length of
#' this (minus one).
#' @param data the data frame where the variables are found.
#' @param subset use this subset of the data
#' @param na.rm a logical value governing the treatment of NA:s. If true, NAs
#' are removed.
#' @param cov a logical indicating whether the variance covariance matrix of
#' the lorenz curve is to be estimated
#' @param no.negatives indicates whether or not negative incomes should be
#' allowed
#' @param cutoffs the values of the q-1 cutoffs for the quantile groups (may be
#' used in within-group estimation)
#' @param group.cutoffs a logical, indicating whether each group in the incdist
#' object should use their own quantile cutoffs or if the overall cutoffs
#' should be used.
#' @return A object of class lorenz with elements \item{p}{the population
#' shares.} \item{ordinates}{the Lorenz curve ordinates.} \item{mean}{the
#' (weighted) mean of \eqn{x}{x}.} \item{n}{the sample size.}
#' \item{sum.weights}{the sum of weights.} \item{quintile.means}{(if \eqn{q}{q}
#' is not false) the within-group means.} \item{quintile.shares}{(if \eqn{q}{q}
#' is not false) the within-group incomes shares.}
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{weighted_mean}}, \code{\link{tip}},
#' \code{\link{locfit}}
#' @references
#' \insertRef{lambert1993}{incdist}
#' \insertRef{beachanddavidson1983}{incdist}
#' \insertRef{beachandkaliski1986}{incdist}
#'
#'
#' @importFrom Rdpack reprompt
#'
#' @examples
#'
#' lorenz(runif(10), q = FALSE)
#' n <- 100
#' income <- exp(10 + rnorm(n)*sqrt(2))
#' weight <- rpois(n,7)
#' lor <- lorenz(income, weight, p = c(0, .5, .7, .9, .98, 1))
#' lor <- lorenz(income, weight, q = 5)
#' plot(lor, main = "Lorenz curve")
#' plot(lor, lor.type = "gen", main = "Generalized Lorenz curve")
#' lor.1 <- lorenz(runif(100))
#' lor.2 <- lorenz(rexp(100))
#' lor.3 <- lorenz(rexp(100), cov = TRUE)
#'
#' plot.lorenz_list(list("1"=lor.1, "2"=lor.2))
#' dominates.lorenz(lor.1, lor.2)
#' dominates.lorenz_list(list("1"=lor.1, "2"=lor.2))
#' library(locfit)
#' x <- runif(500)
#' y <- rexp(500)
#' lf.1 <- lorenz(locfit(~ x))
#' lf.2 <- lorenz(locfit(~ y))
#' dominates.lorenz(lor.1, lor.2)
#' dominates.lorenz_list(list("1"=lor.1, "2"=lor.2))
#' plot.lorenz_list(list("1"=lor.1, "2"=lor.2))
#'
#'
#' @export
lorenz <- function(x, ...) {
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("lorenz", x)
}

## I would like to add lorenz methods for
## - data frames
## - density objects
## - locfit objects
## - incdist objects
## - ecdf objects

#' @export lorenz.default
#' @export
lorenz.default <- function(x, w = rep(1,length(x)),
                           ranked = x,
                           q = 5,
                           p = NULL,
                           data = NULL,
                           subset,
                           na.rm = TRUE,
                           cov = FALSE,
                           no.negatives = FALSE,
                           cutoffs = NULL,
                           ...) {
    ## attach a possible data frame. Remember to detach it!
    if(!is.null(data) & is.data.frame(data)) {
          attach(data)
          on.exit(detach(data))
      }
    ## this does not work. How should the argument to detach be given?
    ## on.exit(detach(as.name("data"), character.only=TRUE))
    ## moved treatment of NA's, missing values and others to utility function
    ## "clean_income"
    if (length(x) != length(ranked))
      warning("The ordering variable length is unequal to x!")
    ## some changes needed since I introduced the possibility to give p as an argument
    if(is.null(p) & q)
        p <- seq(0, 1, 1/q)
    ## if p is given and we are not doing microdata (q is not FALSE), make q equal to its length - 1
    if(!is.null(p) & q)
        q <- length(p) - 1
    ## The microdata case (note that lorenz.default needs to deal with this case)
    ## something is not right here (2021-12-30)
    if(is.null(p) & !q)
        p <- q
    if(!is.null(p) & any(p < 0 | p > 1))
        stop("Population proportions must be within the unit interval!")
    incmat <- clean_income(cbind(x, ranked), w, no.negatives, na.rm)
    x <- incmat[,1]
    ranked <- incmat[,2]
    w <- incmat[,3]
    n <- length(x)
    retval <- list()
    mu <- weighted_mean(x, w)
    sx <- c(x %*% w)
    sigma2 <- weighted_var(x, w)
    if (!q) { # the microdataversion
        ind <- order(ranked)
        x <- x[ind]
        w <- w[ind]
        wsum <- sum(w)
        p <- cumsum(w)/wsum
        lorenz <- cumsum(w * x)/sx
        quantx <- NULL # needed later on
    }
    ## using this weighted quantile function is very inefficient!
    else {
        if(2*q > length(x))
            warning("Few obs per class. You should probably reduce q!")
        ##h.quantx <- wtd.quantile(x, w, probs=seq(0,1,1/q), na.rm = na.rm)
        if(is.null(cutoffs))
            quantx <- weighted_quantile(ranked, w, probs=p, na.rm = na.rm, names = FALSE)
        else  {
            if(length(cutoffs)!= q + 1)
                stop("cutoffs must have length equal to q + 1")
            quantx <-  cutoffs
        }
        ## need to add a check here that the cutoffs are unique!
        duplsx <- duplicated(quantx)
        ## modified to take the unique cut points. Need to modify the
        ## lorenz curve if  there are duplicates.
        cutsx <- cut(ranked, unique(quantx), labels = FALSE, right = FALSE,
                     include.lowest = TRUE)#, right=F)q
        ## sometimes, these have too few (no support for all quantx)
        ## and sometimes too many (weird!) elements
        ## try not having the "as.vector"
        qmeanx <-tapply(seq(along=x), cutsx, function(i, x1=x, w1=w)  weighted_mean(x1[i], w1[i]))
        nqmeanx <- dimnames(qmeanx)[[1]]
        qmeanx <- as.vector(qmeanx)
        qsumx <- tapply(seq(along=x), cutsx, function(i, x1=x, w1=w) sum(x1[i] * w1[i]))
        nqsumx <- dimnames(qsumx)[[1]]
        qsumx <- as.vector(qsumx)
        ## Need to generate qmeanx and qsumx vectors of length p
        ## Scenarios where too few elements occur
        ## 1. Lots of zeros (i.e., effective left censoring)
        ## 2. "flat" regions on the middle
        ## 3. Lots of max or near max values (effective right censoring)
        if(any(duplsx)) {
            warning("Some cut-off points duplicated!")
            ## start with scenario 1 (the most likely)
            ## dupl.index <- which(duplsx)
            ## how to guard against 2 and 3?
            ## note also that the coding is quite convoluted. This should be cleaned up!
            ## work backwards, i.e., sort
            tmp.qmeanx <- sort(qmeanx, decreasing = TRUE)
            tmp.qsumx <- sort(qsumx, decreasing = TRUE)
            tqmeanx <- tqsumx <- rep(0, q)
            ## use tmp.qmeansx to fill in qmeanx from back to front
            g <- 1
            for(h in seq(along = qmeanx)){
                tqmeanx[h] <- ifelse(h <= g, tmp.qmeanx[h], tmp.qmeanx[length(tmp.qmeanx)])
                tqsumx[h] <- ifelse(h <= g, tmp.qsumx[h], tmp.qsumx[length(tmp.qsumx)])
                g <- g + 1
            }
            ## now reverse again
            qmeanx <- sort(tqmeanx, decreasing = FALSE)
            qsumx <- sort(tqsumx, decreasing = FALSE)
            ## reduce clutter
            rm(tmp.qmeanx, tmp.qsumx, tqmeanx, tqsumx)
        }
        ## create the lorenz curve
        lorenz <- cumsum(qsumx)/sx
        ## the below is stupid
        ## p <- seq(0,1, len = q + 1)[2:(q+1)]
        ## is it OK to define these here?
        retval$quantile.means <- qmeanx
        retval$quantile.shares <- qsumx/sx
    }
    retval$q <- q
    retval$ordinates <- lorenz
    ## drop the lowest unless these are microdata
    if (!q) {
        retval$p <- p
    }
    else {
        retval$p <- p[-1]
    }
    retval$mean <- mu
    retval$sigma2 <- sigma2
    ##unnecessary. methods can get this.
    ##retval$g.ordinates <- lorenz*retval$mean
    retval$n <- n
    retval$sum.weights <- sum(w)
    ## only define cut.offs if not based on micro-data
    retval$cut.offs <- quantx
    ##retval$cut.offs.h <- h.quantx
    ## detach the data
    ##if(!is.null(data) & is.data.frame(data))
    ##  detach(data)
    ## start calculating the covariance matrices, if requested
    retval$variances <- NULL
    if(cov) {
        if(!q)
            warning("Cannot calculate the variance for microdata!")
        else {
            object <- structure(retval, class = "lorenz")
            retval$variances <- var_lor(object)
        }
    }
    ## and now return the object as a structure.
    structure(retval, class = "lorenz")
}
## lorenz for a locfit density object
#' @export lorenz.locfit
#' @export
lorenz.locfit <- function(x, ...) {
  if(!inherits(x, "locfit")) stop("Not a locfit object!")
  ## I want to generate the empirical CDF from the object
  ## should I use the data points or the
  ## predicted density at (evalution points, data, ?)
  ##x <- locfit.matrix(object)$x # <- microdata based
  ## reusing x is not good practice
  object <- x
  x <- lfknots(object, what = "x")
  x <- sort(x)
  n <- length(x)
  dx <- diff(x)
  pd <- predict(object, x)
  area.1 <-  pmin(pd[1:(n-1)], pd[2:n]) * dx #<- the rectangle
  area.2 <- 1/2 * abs(diff(pd)) * dx #<- the triangle
  p <- cumsum(area.1 + area.2)#/sum(area.1+area.2) # for to be one?
  p <- c(p, 1)
  lorenz <- cumsum(x)/sum(x)
  retval <- list()
  retval$ordinates <- lorenz
  retval$p <- p
  retval$mean <- c(x %*% pd)
  ##unnecessary. methods can get this.
  ##retval$g.ordinates <- lorenz*retval$mean
  retval$n <- n
  retval$sum.weights <- n
  structure(retval, class = "lorenz")
}
## lorenz for a density object



## the predicate function
#' @export
is.lorenz <- function(x, ...) inherits(x, "lorenz")


## summary and print methods
#' @export
summary.lorenz <- function(x, ...)
  structure(x, class = c("summary.lorenz", class(object)))

#' @export
print.summary.lorenz <- function(x, ...) {
  q <- length(x$p)
  cat("Lorenz curve at ", q, " points.\n")
  cat("Overall mean: ", format(x$mean, ...), "\n")
  cat("Sample size: ", format(x$n, ...), "\n")
  cat("Sum of weights: ", format(x$sum.weights, ...), "\n")
  mat <- cbind(x$p, x$ordinates, x$quantile.means)
  colnames(mat) <-
    c("Population prop", "Cum income share", "Group means")
  rownames(mat) <- as.character(1:length(x$ordinates))
  print(mat)
  invisible(x)
}

#' @export
print.lorenz <- function(x, ...) {
    object <- x
    if(!is.lorenz(object)) stop("Not a Lorenz curve!")
    cat("Lorenz curve at ", length(object$p), " points.\n")
    cat("Overall mean: ", object$mean, "\n")
    cat("Sample size: ", object$n, "\n")
    cat("Sum of weights: ", object$sum.weights, "\n")
    invisible(object)
  }
## an as.data.frame method
#' @export as.data.frame.lorenz
#' @export
as.data.frame.lorenz <- function(x, row.names, optional, ...) {
    object <- x
    if(!is.lorenz(object)) stop("Not a Lorenz curve!")
    ## a problem with very highly concentrated data, work around by select only [1:k]
    k <- length(object$p)
    ordinates <- c(0, object$ordinates[1:k])
    p <- c(0, object$p)
    ## change first element to a zero
    ## why? populate with the mean instead!
    ## populate the full vector
    mean <- rep(object$mean, k+1)
    n <- rep(object$n, k+1)
    sum.weights <- rep(object$sum.weights, k+1)
    ## add the GL and the AL ordinates
    gl.ordinates <- mean*ordinates
    abs.ordinates <- (ordinates - p)*mean
    if(!object$q) {
          quantile.means <-  quantile.shares <- cut.offs <- quantile.groups <- rep(NA, k+1)
      }
    else { ## if based on grouped data
        quantile.means <- c(NA, object$quantile.means[1:k])
        quantile.shares <- c(NA, object$quantile.shares[1:k])
        cut.offs <- object$cut.offs
        quantile.groups <- paste(c("0", 1:k))
      }
    ret <-
      data.frame(p, ordinates, quantile.means, quantile.shares, cut.offs,
                 quantile.groups,
                 mean, gl.ordinates, abs.ordinates,
                 n, sum.weights)
    ret
  }

## a function to convert a list of lorenz curves into a dataset
#' @export as.data.frame.lorenz_list
#' @export
as.data.frame.lorenz_list <- function(x, row.names, optional, ...) {
    obj <- x
    if(!is.list(obj))
      stop("Object is not a list!")
    if(!all(sapply(obj, is.lorenz)))
      stop("Some element  is not a Lorenz curve!")
    ret <- lapply(obj, as.data.frame)
    if(is.null(list.names <- names(obj)))
        names(ret) <- paste(seq(along=ret))
    else
      names(ret) <- list.names
    for(i in names(ret)) ret[[i]]$elname <- rep(i, dim(ret[[i]])[1])
    ret <- do.call("rbind", ret)
    ret
  }

## compare two lorenz curves

#' @export dominates.lorenz
#' @export
dominates.lorenz <- function(x, y, lor.type="ord", rep.num=TRUE, above.p=FALSE, ...) {
    if(!is.lorenz(x)||!is.lorenz(y)) stop("Not a Lorenz curve!")
    ## based on micro.data?
    microdata <- FALSE
    if(!x$q || !y$q) microdata <- TRUE
    if(x$q != y$q) stop("x and y of unequal lengths!")
    if(microdata) {
        if(lor.type=="ord") {
            lx <- x$ordinates
            ly <- y$ordinates
           }
        if(lor.type=="gen") {
            lx <- x$ordinates*x$mean
            ly <- y$ordinates*y$mean
           }
        if(lor.type=="abs") {
            lx <- (x$ordinates-x$p)*x$mean
            ly <- (y$ordinates-y$p)*y$mean
           }
        px <- x$p
        py <- y$p
        names(lx) <- names(px) <- paste("x", seq(along=lx), sep="")
        names(ly) <- names(py) <- paste("y", seq(along=ly), sep="")
        pxy <- c(px, py)
        pxy <- pxy[order(pxy)]
        llx <- lx[names(pxy)]
        ## weird thing with the names
        tc <- names(llx)
        tc[is.na(tc)] <- names(py)
        names(llx) <- tc
        lly <- ly[names(pxy)]
        ## weird thing with the names
        tc <- names(lly)
        tc[is.na(tc)] <- names(px)
        names(lly) <- tc
        ## interpolate
        for(i in seq(along=llx)) {
            if(i==1) {
                llx[i] <- ifelse(is.na(llx[i]), 0, llx[i])
                lly[i] <- ifelse(is.na(lly[i]), 0, lly[i])
              }
            if(i==length(llx)) {
                llx[i] <- ifelse(is.na(llx[i]), llx[i-1], 1)
                lly[i] <- ifelse(is.na(lly[i]), lly[i-1], 1)
              }
            if(i>1 && i<length(llx)) {
                llx[i] <- ifelse(is.na(llx[i]), llx[i-1], llx[i])
                lly[i] <- ifelse(is.na(lly[i]), lly[i-1], lly[i])
              }
          }
        diff <- llx-lly
        if(above.p) diff <- diff[pxy>above.p]
      }
    else {
          q <- length(x$ordinates)
          omit.these <- 1
          if(lor.type=="ord") {
                  diff <- x$ordinates - y$ordinates
                  omit.these <- c(1,q)
              }
          if(lor.type=="gen")
              diff <-   x$ordinates*x$mean - y$ordinates*y$mean
          if(lor.type=="abs")
              diff <- (x$ordinates - x$p)*x$mean - (y$ordinates - y$p)*y$mean
          diff <- diff[-omit.these]
      }
    if(rep.num==TRUE) {
        if(all(diff >= 0)) ret <- TRUE
        ##if((any(diff > 0) && any(diff < 0)) || all(diff == 0)) ret <- NA
        if(!all(diff >= 0)) ret <- FALSE
      }
    if(rep.num=="tex") {
        if(any(diff > 0) && !any(diff < 0)) ret <- "$<$"
        if((any(diff > 0) && any(diff < 0)) || all(diff == 0)) ret <- "$\\sim$"
        if(!any(diff > 0) && any(diff < 0)) ret <- "$>$"
      }
    ret
  }

#' @export dominates.lorenz_list
#' @export
dominates.lorenz_list <-
  function(x, symmetric=FALSE, lor.type="ord", rep.num=TRUE, above.p=FALSE, ...) {
      object <- x
      if(!is.list(object))
          stop("Object is not a list!")
    if(!all(sapply(object, is.lorenz)))
      stop("Some element  is not a Lorenz curve!")
    k <- length(object)
    ret <- matrix(FALSE, ncol = k, nrow = k)
    ## give rows and labels names
    if(!is.null(names(object))) {
        rownames(ret) <- names(object)
        colnames(ret) <- names(object)
      }
    else {
        rownames(ret) <- paste(1:k)
        colnames(ret) <- paste(1:k)
      }
    ## populate the matrix
    ## this does the lower half. redo the code
##     for(i in 1:(k-1))
##       for(j in (i+1):k)
##           ret[i,j] <- dominates.lorenz(object[[i]], object[[j]],
##                                        lor.type=lor.type, rep.num=rep.num)
##     if(symmetric)
##       {
##         for(i in 1:k)
##           for(j in 1:i)
##           {
##             if(j<i) ret[i,j] <- -1*ret[j,i]
##             else ret[i,j] <- 0
##           }
##       }
    for(i in 1:k) {
        for(j in 1:k) {
            ret[i,j] <-
              dominates.lorenz(object[[i]], object[[j]],
                               lor.type=lor.type, rep.num=rep.num,
                               above.p=above.p, ...)
          }
      }
    ret[is.na(ret)] <- FALSE
    diag(ret) <- TRUE
    ret
  }

## the variance using Beach and Davidson and Beach and Kaliski

#' @export
var_lor <- function(x, what = "lorenz", ...) {
    object <- x
    if(!is.lorenz(object)) stop("Object is not a Lorenz curve!")
    mu <- object$mean
    sigma2 <- object$sigma2
    q <- object$q
    muj <- object$quantile.means
    sigma2j <- object$quantile.shares
    pr <- object$p
    ksi <- object$cut.offs
    n <- object$n
    ## gamma and lambda as in 2a-sc in B&K
    gamma <- lambda2 <- numeric(q)
    gamma[1] <- muj[1]
    lambda2[1] <- sigma2j[1]
    for(i in 2:q) {
        gamma[i] <- (i - 1/i)*gamma[i] + 1/i*muj[i]
        lambda2[i] <- (i - 1/i)*lambda2[i-1] + 1/i*sigma2j[i] +
          (i - 1)*(gamma[i] - gamma[i-1])^2
      }
    ## be careful that the dimensions of cutsx etx are right!
    ## work first with theorem 1 in beach and davidson
    ## and now for eq. 8
    Omega <- matrix(0, ncol = q, nrow = q)
    for(j in 1:q)
      for(i in 1:j) {
          Omega[i, j] <-
            pr[i]*(lambda2[i] +
                   (1 - pr[j])*(ksi[i]-gamma[i])*
                   (ksi[j]-gamma[j]) +
                   (ksi[i]-gamma[i])*(gamma[j]-gamma[i]))
          ## the upper triangle for symmetry...
          if(i < j) Omega[j, i] <- Omega[i, j]
        }
    ## now for V with elements in eq 20 for theorem 2
    V.l <- matrix(0, ncol = q-1, nrow = q-1)
    ## check!! indexing
    for(j in 1:(q-1))
      for(i in 1:j) {
          V.l[i, j] <-
            1/mu^2*Omega[i, j] +
              pr[i]*gamma[i]/mu^2 *
                pr[j]*gamma[j]/mu^2*sigma2 -
                  pr[i]*gamma[i]/mu^3*Omega[j, q] -
                    pr[j]*gamma[j]/mu^3*Omega[i, q]
          if(i < j) V.l[j, i] <- V.l[i, j]
        }
    ## check!! this will not work, something is wrong with
    ## dimensions. Check V.l[q, ]!
    ## my q ==  K + 1 in Beach and Davidson
    t.V <- rbind(
                 cbind(
                       rbind(rep(0, q),
                             cbind(rep(0, q-1), V.l)),
                       rep(0, q)),
                 rep(0, q+1))
    ## now for the variance matrix of income shares
    V.s <- matrix(0, ncol = q, nrow = q)
    for(j in 1:q)
      for(i in 1:j) {
          V.s[i, j] <- t.V[i+1, j+1] +  t.V[i, j] +
            - t.V[i+1, j] - t.V[i, j+1]
          if(i < j) V.s[j, i] <- V.s[i, j]
        }
    retval <- list()
    dimnames(V.l) <- list(paste(1:(q-1)), paste(1:(q-1)))
    dimnames(V.s) <- list(paste(1:q), paste(1:q))
    ## V.l and V.s are the variances for sqrt(n)*(t - theta)!
    retval$ordinates <- V.l/n
    retval$quantile.shares <- V.s/n
    retval
  }


## lorenz curves for a incdist object
#' @export lorenz.incdist
#' @export
lorenz.incdist <- function(x, q=5, p = NULL, equivalise = FALSE, group.cutoffs=TRUE, concentration=FALSE, ...) {
    object <- x
    ## test if this is an incdist  object
    if (!inherits(object, "incdist")) {
      stop("First argument must be the incdist object!")
    }
    ## strategy:
    ## to proper Lorenz curves for "y" and concentration curves for x:s
    ## code stolen from summary.incdist
        ## 2. this is where the income inequality code should start
    ## y contains the income variable
    ## check if y == x (to some tolerance)
    ## i. for each part
    ## a. inequality & poverty indices, lorenz & tip objects
    ## b. concentration curves & indices
    ## c. figure out the group decompositions (i.e, do a-b for these)
    ## ii. for the full set of partitions
    ## a-c
    attach(object)
    on.exit(detach(object))
    ## this does not work. How should the argument to detach be given?
    ##on.exit(detach("object", character.only=TRUE))
    ## some changes needed since I introduced the possibility to give p as an argument
    if(is.null(p) & q)
        p <- seq(0, 1, 1/q)
    ## if p is given and we are not doing microdata (q is not FALSE), make q equal to its length - 1
    if(!is.null(p) & q)
        q <- length(p) - 1
    ## the microdata case (note that lorenz.default needs to deal with this case)
    if(!q)
        p <- q
    income <- terms(object$formula, data = frm)
    if(length(panames))
        pal <- levels(frm[[panames]])
    if(length(grnames))
        grl <- levels(frm[[grnames]])
    ## need to make some stuff explicitly null
    if(!length(grnames)) grl <- NULL
    if(!length(panames)) pal <- NULL
    eqscl <- object$eqscale
    ## save the old frame, will be needed if we want another e.s.
    old.frm <- frm
    if(!is.null(eqscl) & equivalise)
        frm <- eqscale(object)
    # 0. the top level statistics
    if(length(panames)) frm.list <- split(frm, frm[[panames]])
    else {
        frm.list <- list()
        frm.list[[1]] <- frm
      }
    ## create the lists (structures)
    ## ret.x and ret.y to hold the statistics
    ## n to hold sample sizes
    ## sumw to hold sum of weights
    ## must be done here
    ret.y <- n <- sumw <- list() ##ret.y{i:group}{l:partition}
    for(i in 1:(1+length(grl))) {
          ret.y[[i]] <- n[[i]] <- sumw[[i]] <- list()
          for(l in 1:length(frm.list)) {
                ret.y[[i]][[l]] <- n[[i]][[l]] <- sumw[[i]][[l]] <- NA
              }
          names(ret.y[[i]]) <- names(n[[i]]) <- names(sumw[[i]]) <- pal
        }
    names(ret.y) <- names(n) <- names(sumw) <- c("all", grl)
    ret.x <- list() ##ret.y{j:incomecomps}{i:group}{l:partition}
    ## skip the first incnames, which is in "ret.y"
    if(length(incnames)) {
        for(j in 1:length(incnames)) {
            ret.x[[j]] <- list()
            for(i in 1:(1+length(grl))) {
                ret.x[[j]][[i]] <- list()
                for(l in 1:length(frm.list)) {
                    ret.x[[j]][[i]][[l]] <- NA
                  }
                names(ret.x[[j]][[i]]) <- pal
              }
            names(ret.x[[j]]) <- c("all", grl)
          }
        names(ret.x) <- incnames
      }
    ## start doing the work
    for(l in 1:length(frm.list)) { ## l indexes partitions {
        if(length(grnames)) {
            count.grnames <- c(dim(frm.list[[l]])[1],
                               table(frm.list[[l]][[grnames]]))
            frm.list.l <- split(frm.list[[l]],
                                ## should I give a levels argument
                                ## to ensure that all partitions have
                                ## the same structure (no. of levels)?
                                as.factor(frm.list[[l]][[grnames]]))
            frm.list.l <- c(list(frm.list[[l]]), frm.list.l)
          }
        else {
            frm.list.l <- list(frm.list[[l]])
            count.grnames <- dim(frm.list[[l]])[1]
          }
        ## this used to be (??)
        ## i indexes groups present. The first is all.
        ## add this argument for the case in case you want to use common (overall)
        ## quantile cutoffs in each group
        cutoffs <- NULL
        for(i in 1:(1+length(grl))) {
            doing.name <- grnames[i]
            if(group.cutoffs && i>1)
                cutoffs <- ret.y[[1]][[l]][["cut.offs"]]
            ## must add some code here to
            ## a. check if every i has a dataframe in frm.list.l
            ## b. if not, create an empty data frame for those
            if(count.grnames[i] == 0) {
                frm.list.l[[i]] <- subset(frm.list[[l]], TRUE)
                if(attr(income, "response"))
                  y <- 0
                n[[i]][[l]] <- 0
              }
            if (count.grnames[i] != 0 && attr(income, "response")) {
              y <- model.extract(frm.list.l[[i]], response)
              n[[i]][[l]] <- length(y)
              ##yname <- deparse(formula[[2]])
            }
            else {
              y <- NULL
            }
            if(length(incnames)) {
                x <- as.matrix(frm.list.l[[i]][, incnames])
                if (length(incnames) == dim(x)[2])
                  dimnames(x) <- list(NULL, incnames)
              }
            weights <- model.extract(frm.list.l[[i]], weights)
            sumw[[i]][[l]] <- sum(weights)
            ## "func" holds the functions to be calculated
            ## these must accept two (an exactly two) arguments
            ## the data and the weights
            if (!count.grnames[i]) next
            if(!is.null(weights))
              ret.y[[i]][[l]] <-
                lorenz.default(as.vector(y),  w = weights, p = p , cutoffs = cutoffs)
            else
              ret.y[[i]][[l]] <-
                lorenz.default(as.vector(y), p = p, cutoffs = cutoffs)
            if(length(incnames)) {
                for(j in 1:length(incnames)) {
                    ## "func" holds the functions to be calculated
                    ## these must accept two (an exactly two) arguments
                    ## the data and the weights
                    if (!count.grnames[i]) next
                    if(!is.null(weights)) {
                        if(concentration)
                          ret.x[[j]][[i]][[l]] <-
                              lorenz.default(x[,j], w = weights, ranked = as.vector(y), p = p,
                                             cutoffs = cutoffs)
                        else
                          ret.x[[j]][[i]][[l]] <-
                            lorenz.default(x[,j], w = weights, p = p, cutoffs = cutoffs)
                      }
                    else {
                        if(concentration)
                          ret.x[[j]][[i]][[l]] <-
                            lorenz.default(x[,j], ranked=as.vector(y), p = p, cutoffs = cutoffs)
                        else
                          ret.x[[j]][[i]][[l]] <-
                            lorenz.default(x[,j], p = p, cutoffs = cutoffs)
                    }
                  } ## j
              } ## if j
        }
      } ## partitions (years, mostly, could be countries
    ## detach(object)
    ret <- list(ret.y, ret.x)
    ## this may not work! Moreover, maybe lorenz_incdist is sufficient?
    structure(ret, class = c("lorenz_incdist", "lorenz"))
  }


## and an as.data.frame method
#' @export as.data.frame.lorenz_incdist
#' @export
as.data.frame.lorenz_incdist <- function(x, ...) {
        obj <- x
        if(!is.list(obj))
            stop("Object is not a list!")
        ## this needs to change: does not function well!
        ## first the "all"        part
        gs <- names(obj[[1]])
        tdl.top <-
            lapply(gs,
                   function(x) {
                       ids <- names(obj[[1]][[x]])
                       ret <- lapply(obj[[1]][[x]], as.data.frame.lorenz)
                       for(i in seq(along=ids)) {
                           ret[[i]][["id"]] <- ids[i]
                       }
                       ret <- do.call("rbind", ret)
                       ret[["comp"]] <- "overall"
                       ret[["group"]] <- x
                       ret
                   })
        tdl.top <- do.call("rbind", tdl.top)
        ## now the components
        if(length(obj)<2)
            td <- tdl.top
        else {
                comps <- names(obj[[2]])
                tdl.comp <- vector(mode="list", length=length(comps))
                names(tdl.comp) <- comps
                for(j in seq(along=tdl.comp)) {
                        tdl.comp[[j]] <-
                            lapply(gs,
                                   function(x) {
                                       ids <- names(obj[[2]][[j]][[x]])
                                       ret <- lapply(obj[[2]][[j]][[x]], as.data.frame.lorenz)
                                       for(i in seq(along=ids)) {
                                           ret[[i]][["id"]] <- ids[i]
                                       }
                                       ret <- do.call("rbind", ret)
                                       ret[["comp"]] <- "overall"
                                       ret[["group"]] <- x
                                       ret
                                   })
                        tdl.comp[[j]] <- do.call("rbind", tdl.comp[[j]])
                        tdl.comp[[j]][["comp"]] <- comps[j]
                    }
                tdl.comp <- do.call("rbind", tdl.comp)
                td <- rbind(tdl.top, tdl.comp)
            }
        td
    }

## convert from a data.frame to a (single) lorenz curve
#' @export as.lorenz
#' @export
as.lorenz <- function(df, namel=list(p="p", ordinates="ordinates"), ...) {
    q <- dim(df)[1]-1
    obj <- list(q=q,
                quantile.means=rep(NA, q),
                quantile.shares=rep(NA, q),
                 cut.offs=rep(NA, q+1),
                 quantile.groups=rep(NA, q),
                 mean=NA,
                 n=NA,
                sum.weights=NA)
    ## stop()
    for(i in names(namel)) obj[[i]] <- df[[namel[[i]]]]
    structure(obj, class="lorenz")
}
