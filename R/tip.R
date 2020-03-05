# a tip curve



#' The "three I's of poverty" (TIP) or "Cumulated Poverty Gap" (CGP) curve
#'
#' This function returns the cumulated poverty gap.
#'
#' This function computes the `TIP`curve, that is, the cumulated poverty gap as
#' a function of the population proportion as \deqn{TIP(p) := \sum_{F^{-1}(x_i)
#' < p} I(x_i < z) (z - x_i),} or, if the normalized by the poverty line,
#' \deqn{TIP(p) := \sum_{F^{-1}(x_i) < p} I(x_i < z) (1 - x_i/z).} The TIP
#' curve is increasing for values of x below the poverty line z, at which point
#' it becomes horizontal. That point gives the proportion below the poverty
#' line (p) and the cumulative average poverty gap (the value of the TIP
#' curve).  If weights are given, the function returns the weighted TIP curve.
#' If \eqn{q = FALSE}, tip() returns an object whose elements \eqn{p} and
#' \eqn{ordinates} are of length(x).
#'
#' Simple summary, print and plot methods are provided.
#'
#' @aliases tip.default tip cpg is.tip summary.tip print.summary.tip
#' dominates.tip dominates.tiplist
#' @param x the resource variable
#' @param w the sampling weights
#' @param q an integer or logical. If an integer, it is the number of equally
#' spaced points at which the Lorenz curve ordinates are estimated. If FALSE,
#' then the Lorenz curve is estimated at every point in x.
#' @param z the poverty line
#' @param data the data frame where the variables are found.
#' @param alpha the FGT parameter
#' @param normalize a logical. If true, tip returns the normalized gap (i.e.,
#' the gap divided by the poverty line).
#' @param mean a logical indicating whether the default poverty line relate to
#' the mean rather than the median
#' @param fraction the fraction of the mean or median to use as the poverty
#' line
#' @param na.rm a logical indicating whether NA's should removed
#' @return A object of class tip with elements \item{p}{the population
#' proportion} \item{ordinates}{the TIP curve value at p} \item{hc}{the
#' proportion of units with x less than z} \item{relgap}{the average poverty
#' gap (across the poor population)} \item{n}{sample size} \item{z}{the poverty
#' line} \item{sum.weights}{the sum of weights}
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{lorenz}}, \code{\link{fgt}}
#' @references Jenkins, Stephen P & Peter J Lambert ( 1997), `Three `i's of
#' poverty curves, with an analysis of U.K. poverty trends', Oxford Economic
#' Papers 49(3), 317-327.
#' @examples
#'
#' income <- rexp(100)
#' weight <- rpois(100,5)
#' tip.ex <- tip(income, weight, q = FALSE)
#' plot(tip.ex, lwd = 2)
#'
#' @export tip
tip <- function(x, ...)
  {
    UseMethod("tip", x)
  }
#' @export cpg
cpg <- function(x, ...)
  {
    tip(x, ...)
  }
#' @export tip.default
tip.default <- function(x, w = rep(1,length(x)), q = 20, data = NULL,
                        alpha = 0,
                        fraction = 0.5,
                        normalize = TRUE,
                        mean = FALSE,
                        z = ifelse(mean,
                          fraction*weighted.mean(x,w),
                          fraction*weighted.median(x,w)),
                        na.rm = TRUE, ...)
{
  ## attach a possible data frame. Remember to detach it!
  if(!is.null(data) & is.data.frame(data))
      {
          attach(data)
          on.exit(detach(data))
        }
  # moved treatment of NA's, missing values and others to utility function
  # "clean.income"
  incmat <- clean.income(x, w, na.rm)
  x <- incmat[,1]
  w <- incmat[,2]
  retval <- list()
  if (!q) { # the microdataversion
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]
    wsum <- sum(w)
    p <- cumsum(w)/wsum
    if(normalize)
      gap <- (z - x)/z*(x < z)
    else
      gap <- (z - x)*(x < z)
    sg <- w %*% gap
    cg <- cumsum(w * gap)/sum(w)#/sg # Do not divide by total gap!
    p <- c(0,p)
    tip <- c(0,cg)
  }
  # using this weighted quantile function is very inefficient!
  else {
    if(2*q > length(x))
      warning("Few obs per class. You should probably reduce q!")
    ## this must be modified so it handles the case of non-unique quantiles
    quantx <- weighted.quantile(x, w, probs=seq(0,1,1/q), names = FALSE)
     ## need to add a check here that the cutoffs are unique!
    duplsx <- duplicated(quantx)
    cutsx <- cut(x, unique(quantx), labels = FALSE,, right = FALSE,
                 include.lowest = TRUE)#, right=F)

    if(normalize)
      gap <- (z - x)/z*(x < z)
    else
      gap <- (z - x)*(x < z)
    qmeanx <-
      as.vector(tapply(seq(along=gap),cutsx,
                       function(i, x1=gap, w1=w)  weighted.mean(x1[i], w1[i])))
    qsumx <-
      as.vector(tapply(seq(along=gap),cutsx,
                       function(i, x1=gap, w1=w) sum(x1[i] * w1[i])))
    tip <- c(0,cumsum(qsumx)/sum(w))
    ## see above: do not divide by /(gap %*% w))
    ## handle duplicates
     if(any(duplsx))
      {
        warning("Some cut-off points duplicated!")
        dupl.index <- which(duplsx)
        tmp.tip <- tip
        tmp.qmeanx <- qmeanx
        tmp.qsumx <- qsumx
        for(i in dupl.index)
          {
            j <- i+1 # for lorenz curve, i for the others
            tmp.tip <- c(tmp.tip[1:(i-1)],
                         tmp.tip[i-1],
                         tmp.tip[i:length(tmp.tip)])
            tmp.qmeanx <- c(tmp.qmeanx[1:(i-1)],
                            tmp.qmeanx[i-1],
                            tmp.qmeanx[i:length(tmp.qmeanx)])
            tmp.qsumx <- c(tmp.qsumx[1:(i-1)],
                           tmp.qsumx[i-1],
                           tmp.qsumx[i:length(tmp.qsumx)])
          }
        tip <- tmp.tip
        qmeanx <- tmp.qmeanx
        qsumx <- tmp.qsumx
      }
    p <- seq(0,1,1/q)
    # is it OK to define these here?
  }
  retval$ordinates <- tip
  retval$p <- p
  retval$q <- q
  retval$hc <- weighted.mean((x < z), w)
  retval$relgap <- weighted.mean(gap,w)
  retval$z <- z
  retval$n <- length(x)
  retval$sum.weights <- sum(w)
  class(retval) <- "tip"
  ## detach the data
  ## if(!is.null(data) & is.data.frame(data))
  ##   detach(data)
  retval
}

is.tip <- function(object) inherits(object, "tip")

# a print method


summary.tip <- function(object, ...) {
    micro <- object$n == length(object$p) - 1
    if(micro) {
      ind <- object$p < object$hc
      ind[sum(ind)+1] <- TRUE
      p <- object$p[ind]
      ordinates <- object$ordidates[ind]
    }
    structure(object, class = c("sum.tip", class(object)))
  }

print.sum.tip <- function(x, ...){
  # how to handle a micro curve? in summary?
  q <- length(x$p) - 1
  cat("TIP curve at ", q, " points.\n")
  cat("Poverty head count ratio: ", format(x$hc, ...), "\n")
  cat("Poverty relative poverty gap: ", format(x$relgap, ...), "\n")
  cat("Sample size: ", format(x$n, ...), "\n")
  cat("Sum of weights: ", format(x$sum.weights, ...), "\n")
  mat <- cbind(x$p, x$ordinates)
  ind <- x$p < x$hc
  ind[sum(ind) + 1] <- TRUE
  mat <- mat[ ind, ]
  colnames(mat) <-
    c("Population prop", "Cum poverty gap")
  print(mat)
  invisible(x)
}

## as.data.frame

as.data.frame.tip <- function(x, row.names, optional, ...)
  {
    object <- x
    if(!is.tip(object)) stop("Not a TIP curve!")
    k <- length(object$p)
    p <- object$p
    ordinates <- object$ordinates
    hc <- rep(object$hc, k)
    relgap <- rep(object$relgap, k)
    z <- rep(object$z, k)
    n <- rep(object$n, k)
    sum.weights <- rep(object$sum.weights, k)
    ret <-
      data.frame(p, ordinates, hc, relgap, z, n, sum.weights)
    ret
  }

## a function to convert a list of lorenz curves into a dataset
as.data.frame.tip.list <- function(x, row.names, optional, ...)
  {
    obj <- x
    if(!is.list(obj))
      stop("Object is not a list!")
    if(!all(sapply(obj, is.tip)))
      stop("Some element  is not a TIP curve!")
    ret <- lapply(obj, as.data.frame.tip)
    if(is.null(list.names <- names(obj)))
        names(ret) <- paste(seq(along=ret))
    else
      names(ret) <- list.names
    for(i in names(ret)) ret[[i]]$elname <- rep(i, dim(ret[[i]])[1])
    ret <- do.call("rbind", ret)
    ret
  }


## tip dominance

dominates.tip <-  function(x, y, rep.num=TRUE,
                           above.p=FALSE)
  {
    if(!is.tip(x)||!is.tip(y)) stop("Not a Tip curve!")
    ## based on micro.data?
    if(!x$q || !y$q) stop("Not implemented yet for micro data!")
    if(x$q != y$q) stop("x and y of unequal lengths!")
    diff <- x$ordinates - y$ordinates
    ## drop the first one
    diff <- diff[-1]
    ## copied from dominates.lorenz
    if(rep.num==TRUE)
      {
        if(all(diff >= 0)) ret <- TRUE
        ##if((any(diff > 0) && any(diff < 0)) || all(diff == 0)) ret <- NA
        if(!all(diff >= 0)) ret <- FALSE
      }
    if(rep.num=="tex")
      {
        if(any(diff > 0) && !any(diff < 0)) ret <- "$<$"
        if((any(diff > 0) && any(diff < 0)) || all(diff == 0)) ret <- "$\\sim$"
        if(!any(diff > 0) && any(diff < 0)) ret <- "$>$"
      }
    ret
  }


dominates.tip.list <- function(object, rep.num=TRUE, above.p=FALSE, ...)
  {
    if(!is.list(object))
      stop("Object is not a list!")
    if(!all(sapply(object, is.tip)))
      stop("Some element  is not a Tip curve!")
    k <- length(object)
    ret <- matrix(FALSE, ncol = k, nrow = k)
    if(!is.null(names(object)))
      {
        rownames(ret) <- names(object)
        colnames(ret) <- names(object)
      }
    else
      {
        rownames(ret) <- paste(1:k)
        colnames(ret) <- paste(1:k)
      }
    for(i in 1:k)
      for(j in 1:k)
          ret[i,j] <- dominates.tip(object[[i]], object[[j]],
                                    rep.num=rep.num,
                                    above.p=above.p, ...)
    ret[is.na(ret)] <- FALSE
    diag(ret) <- TRUE
    ret
  }

## tip curves for a incdist object

tip.incdist <- function(object, q=5, equivalise = FALSE, only.aggregate=TRUE,
                           concentration=TRUE, ...)
  {
    ## test if this is an incdist  object
    if (!inherits(object, "incdist")){
      stop("First argument must be the incdist object!")
    }
    ## strategy:
    ## to proper Tip curves for "y" and concentration curves for x:s
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
    income <- terms(formula, data = frm)
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
    else
      {
        frm.list <- list()
        frm.list[[1]] <- frm
      }
    ## create the lists (structures)
    ## ret.x and ret.y to hold the statistics
    ## n to hold sample sizes
    ## sumw to hold sum of weights
    ## must be done here
    ret.y <- n <- sumw <- list() ##ret.y{i:group}{l:partition}
    for(i in 1:(1+length(grl)))
        {
          ret.y[[i]] <- n[[i]] <- sumw[[i]] <- list()
          for(l in 1:length(frm.list))
              {
                ret.y[[i]][[l]] <- n[[i]][[l]] <- sumw[[i]][[l]] <- NA
              }
          names(ret.y[[i]]) <- names(n[[i]]) <- names(sumw[[i]]) <- pal
        }
    names(ret.y) <- names(n) <- names(sumw) <- c("all", grl)
    ## ret.x <- list() ##ret.y{j:incomecomps}{i:group}{l:partition}
    ## skip the first incnames, which is in "ret.y"
    ## if(length(incnames))
    ##   {
    ##     for(j in 1:length(incnames))
    ##       {
    ##         ret.x[[j]] <- list()
    ##         for(i in 1:(1+length(grl)))
    ##           {
    ##             ret.x[[j]][[i]] <- list()
    ##             for(l in 1:length(frm.list))
    ##               {
    ##                 ret.x[[j]][[i]][[l]] <- NA
    ##               }
    ##             names(ret.x[[j]][[i]]) <- pal
    ##           }
    ##         names(ret.x[[j]]) <- c("all", grl)
    ##       }
    ##     names(ret.x) <- incnames
    ##   }
    ## start doing the work
    for(l in 1:length(frm.list)) ## l indexes partitions
      {
        if(length(grnames))
          {
            count.grnames <- c(dim(frm.list[[l]])[1],
                               table(frm.list[[l]][[grnames]]))
            frm.list.l <- split(frm.list[[l]],
                                ## should I give a levels argument
                                ## to ensure that all partitions have
                                ## the same structure (no. of levels)?
                                as.factor(frm.list[[l]][[grnames]]))
            frm.list.l <- c(list(frm.list[[l]]), frm.list.l)
          }
        else
          {
            frm.list.l <- list(frm.list[[l]])
            count.grnames <- dim(frm.list[[l]])[1]
          }
        ## this used to be (??)
        ## i indexes groups present. The first is all.
        for(i in 1:(1+length(grl)))
          {
            doing.name <- grnames[i]
            ## must add some code here to
            ## a. check if every i has a dataframe in frm.list.l
            ## b. if not, create an empty data frame for those
            if(count.grnames[i] == 0)
              {
                frm.list.l[[i]] <- subset(frm.list[[l]], T)
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
            if(length(incnames))
              {
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
                tip.default(as.vector(y), weights, q=q)
            else
              ret.y[[i]][[l]] <-
                tip.default(as.vector(y), q=q)
            ## if(length(incnames))
            ##   {
            ##     for(j in 1:length(incnames))
            ##       {
            ##         ## "func" holds the functions to be calculated
            ##         ## these must accept two (an exactly two) arguments
            ##         ## the data and the weights
            ##         if (!count.grnames[i]) next
            ##         if(!is.null(weights))
            ##           {
            ##             if(concentration)
            ##               ret.x[[j]][[i]][[l]] <-
            ##                 tip.default(x[,j], w=weights,
            ##                                ranked=as.vector(y),
            ##                                q=q)
            ##             else
            ##               ret.x[[j]][[i]][[l]] <-
            ##                 tip.default(x[,j], w=weights,
            ##                                q=q)
            ##           }
            ##         else
            ##           {
            ##             if(concentration)
            ##               ret.x[[j]][[i]][[l]] <-
            ##                 tip.default(x[,j],
            ##                                ranked=as.vector(y),
            ##                                q=q)
            ##             else
            ##               ret.x[[j]][[i]][[l]] <-
            ##                 tip.default(x[,j],
            ##                                q=q)
            ##           }
            ##       } ## j
            ##   } ## if j
          }
      } ## partitions (years, mostly, could be countries
    ## ret <- list(ret.y, ret.x)
    ret <- ret.y
    structure(ret, class = c("tip.incdist", "tip"))
  }

## and an as.data.frame method

as.data.frame.tip.incdist <-
    function(x, ...)
    {
        obj <- x
        if(!is.list(obj))
          stop("Object is not a list!")
        ## this needs to change: does not function well!
        ## first the "all"        part
        gs <- names(obj)
        tdl.top <-
            lapply(gs,
                   function(x)
                     {
                       ids <- names(obj[[x]])
                       ret <- lapply(obj[[x]], as.data.frame.tip)
                       k <- sapply(ret, function(x) dim(x)[1])
                       ret <- do.call("rbind", ret)
                       ret[["id"]] <- rep(ids, each=k)
                       ret[["comp"]] <- "overall"
                       ret[["group"]] <- x
                       ret
                   }
                   )
        td <- do.call("rbind", tdl.top)
        td
    }
