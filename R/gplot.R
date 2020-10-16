## a file for ggplot2-based plotting of incdist objects

## lorenz

## plot the lorenz-curve for a data vector
## note that add will no longer do much?
## (should it just add a curve to an existing plot?)


#' Plot a Lorenz curve
#'
#' plot.lorenz plots a (generalized) Lorenz curve.
#'
#' plot.lorenz_list plots a list of lorenz curve on the same graph and optionally
#' places a legend using the names of the lorenz curves as labels.
#'
#' An ordinary Lorenz curve is bounded in the unit box, is positive and
#' increasing and lies everywhere below the 45 degree line. By default the
#' function draws the 45 degree line and sets ylim and xlim to c(0,1).
#'
#' If several generalized Lorenz curves are drawn, or ordinary Lorenx curves
#' are drawn with diff = FALSE, care need to be taken that the ylim gets
#' correectly set.
#'
#' @aliases plot.lorenz lines.lorenz plot.lorenz_list
#' @param x A Lorenz curve object.
#' @param add A logical. If TRUE, the Lorenz curve is added to a plot.
#' @param lor.type The type of Lorenz curve to be drawn. If "ord", the graph is
#' an ordinary Lorenz curve \eqn{L(p)}, if "gen", the graph is a generalized
#' Lorenz curve \eqn{\mu_x L(p)}.
#' @param diff A logical. If TRUE, graph is of \eqn{p - L(p)} rather than
#' \eqn{L(p)}.
#' @param ylab The y axis label. Defaults to "Income share".
#' @param xlab The x axis label. Defaults to "Population share".
#' @param ... Additional graphics arguments to pass to the plot function.
#' @return A plot.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{lorenz}}.
#' @references
#' @examples
#'
#' plot(lorenz(runif(100), q = 5))
#'
#' @importFrom graphics plot lines
#' @importFrom ggplot2 ggplot
#'
#' @export
plot.lorenz <- function(x, add = FALSE, lor.type = "ord",
                        diff = FALSE,
                        xlab = FALSE, ylab = FALSE,
                        ylim = FALSE, xlim = FALSE,
                        type = "l", lty= FALSE,
                        ...)
{
    object <- x
  if(!is.lorenz(object)) stop("Not a Lorenz curve!")
  ## use the as.data.frame
  td <- as.data.frame(object)
  ## graph parameters. Need to be checked
  tyend <- 1
  ## lor.type <- c("ord", "gen", "abs")
  if(missing(xlab))
    xlab <- "Population share p"
  ## the standard case
  if(lor.type=="ord" & !diff)
    {
      td$y <- td$ordinates
      td1 <- data.frame(p=c(0, 1), y=c(0, 1))
      if(missing(ylab))
        ylab <- "Cumulative income share L(p)"
    }
  if(lor.type=="ord" & diff)
    {
      td$y <- td$ordinates - td$p
      td1 <- data.frame(p=c(0, 1), y=c(0, 0))
      tyend <- 1
  if(missing(ylab) & diff)
    ylab <- "Cumulative income share - diagonal p - L(p)"
    }
  if(lor.type=="gen")
    {
      td$y <- td$gl.ordinates
      td1 <- data.frame(p=c(0, 1), y=c(0, td$mean[2]))
      tyend <- td$mean[2]
      if(missing(ylab))
        ylab <- "Cumulative income share time mean L(p)*mean"
    }
  if(lor.type=="abs") {
    td$y <- td$abs.ordinates
    td1 <- data.frame(p=c(0, 1), y=c(0, 0))
    if(missing(ylab))
      ylab <- "Absolute Lorenz curve ordinate (p - L(p))*mean"
  }
  if(!add)
    d <- ggplot(td, aes(x=p, y=y)) + geom_line() +
       geom_segment(x=0, y=0, xend=1, yend=tyend) + ## geom_line(data=td1) +
      ylab(ylab) + xlab(xlab)
  else
    d <- geom_line(data=td, aes(x=p, y=y))
  d
}

## check the above against the below

# a plot method
# this is missing lor.type at least. Why?
plot.lorenz.defunct <-
  function(x, add = FALSE,
           xlab = FALSE, ylab = FALSE,
           xlim = FALSE, ylim = FALSE, type = FALSE, ...)
{
        object <- x
  if(!is.lorenz(object)) stop("Not a Lorenz curve!")
  ## use the as.data.frame
  td <- as.data.frame(object)
  if(missing(xlab))
    xlab <- "Population share"
    ## the standard case
  if(lor.type=="ord" & !diff)
    {
      td$y <- td$ordinates
      td1 <- data.frame(p=c(0, 1), y=c(0, 1))
      if(missing(ylab))
        ylab <- "Cumulative income share L(p)"
    }
  if(lor.type=="ord" & diff)
    {
      td$y <- td$ordinates - td$p
      td1 <- data.frame(p=c(0, 1), y=c(0, 0))
      if(missing(ylab) & diff)
        ylab <- "Cumulative income share - diagonal p - L(p)"
    }
  if(lor.type=="gen")
    {
      td$y <- td$gl.ordinates
      td1 <- data.frame(p=c(0, 1), y=c(0, td$mean[2]))
      if(missing(ylab))
        ylab <- "Cumulative income share time mean L(p)*mean"
    }
  if(lor.type=="abs") {
    td$y <- td$abs.ordinates
    td1 <- data.frame(p=c(0, 1), y=c(0, 0))
    if(missing(ylab))
      ylab <- "Absolute Lorenz curve ordinate (p - L(p))*mean"
  }
  if(!add)
    d <- ggplot(td, aes(x=p, y=y)) + geom_line() +
      geom_line(data=td1) +
      ylab(ylab) + xlab(xlab)
  else
    d <- geom_line(data=td, aes(x=p, y=y))
  d
}

#' @export
lines.lorenz <-
  function(x, lor.type = "ord", diff = FALSE, ...)
{
  ## do the "..." get correctly passed?
  plot.lorenz(x=x, add=TRUE, lor.type=lor.type, diff=diff, ...)
}

## make a function to plot a list of Lorenz curves
#' @export
plot.lorenz_list <-  function(x, add = FALSE, lor.type = "ord",
                          diff = FALSE,
                          xlab = FALSE, ylab = FALSE,
                          ylim = FALSE, xlim = FALSE,
                          type = "l",
                          col = FALSE,
                          legend = TRUE,
                          ...)
{
            object <- x
    if(!is.list(object))
      stop("Object is is not a list!")
    if(!all(sapply(object, is.lorenz)))
        stop("Some element  is not a Lorenz curve!")
    ## invocation of " as.data.frame.lorenz_list" should be superfluous
    td <- as.data.frame.lorenz_list(object)
    k <- length(unique(td$elname))
    ## lor.type <- c("ord", "gen", "abs")
    if(missing(xlab))
      xlab <- "Population share p"
    ## the standard case
    if(lor.type=="ord" & !diff)
      {
        td$y <- td$ordinates
        td1 <- data.frame(p=c(0, 1), y=c(0, 1))
        if(missing(ylab))
          ylab <- "Cumulative income share L(p)"
      }
    if(lor.type=="ord" & diff)
      {
        td$y <- td$ordinates - td$p
        td1 <- data.frame(p=c(0, 1), y=c(0, 0))
        if(missing(ylab) & diff)
          ylab <- "Cumulative income share - diagonal p - L(p)"
      }
    if(lor.type=="gen")
      {
        td$y <- td$gl.ordinates
        td1 <- data.frame(p=c(0, 1), y=c(0, td$mean[2]))
        if(missing(ylab))
        ylab <- "Cumulative income share time mean L(p)*mean"
      }
    if(lor.type=="abs") {
      td$y <- td$abs.ordinates
      td1 <- data.frame(p=c(0, 1), y=c(0, 0))
      if(missing(ylab))
        ylab <- "Absolute Lorenz curve ordinate (p - L(p))*mean"
    }
    ## expand td1 and add elname
    td1$elname <- rep(unique(td$elname)[1], 2)
    if(!add){
        d <- ggplot(td, aes(x=p, y=y, group=elname, colour=elname)) + geom_line() +
            ylab(ylab) + xlab(xlab)
        if(lor.type!="gen")
            d <-  d + geom_line(data=td1, col="black")
        }
    else
      d <- geom_line(data=td, aes(x=p, y=y, colour=elname))
    d
  }

## the incdist method
#' @export
plot.lorenz_incdist <-  function(x, add = FALSE, lor.type = "ord",
                          diff = FALSE,
                          xlab = FALSE, ylab = FALSE,
                          ylim = FALSE, xlim = FALSE,
                          type = "l",
                          col = FALSE,
                          legend = TRUE,
                          ...)
{
            object <- x
    if(!inherits(object, "lorenz_incdist")) stop("Not an incdist summary object!")
    td <- as.data.frame(object)
    ## preserve the order of group, comp
    for(i in c("id", "comp", "group"))
      td[[i]] <- factor(td[[i]], labels=unique(td[[i]]))
    ## for later reference
    td2 <- unique(td[c("id", "comp", "group", "mean")])
    td2$p <- 1
    td2$y <-
      ifelse(lor.type=="gen", td2$mean,
             ifelse(lor.type=="ord", 1, 0))
    td2$mean <- NULL
    td1 <- td2[c("id", "comp", "group")]
    td1$p <- 0
    td1$y <- 0
    td1 <- rbind(td1, td2)
    rm(td2)
    ## lor.type <- c("ord", "gen", "abs")
    if(missing(xlab))
      xlab <- "Population share p"
    ## the standard case
    if(lor.type=="ord" & !diff)
      {
        td$y <- td$ordinates
        if(missing(ylab))
          ylab <- "Cumulative income share L(p)"
      }
    if(lor.type=="ord" & diff)
      {
        td$y <- td$ordinates - td$p
        if(missing(ylab) & diff)
          ylab <- "Cumulative income share - diagonal p - L(p)"
      }
    if(lor.type=="gen")
      {
        td$y <- td$gl.ordinates
        if(missing(ylab))
        ylab <- "Cumulative income share time mean L(p)*mean"
      }
    if(lor.type=="abs") {
      td$y <- td$abs.ordinates
      if(missing(ylab))
        ylab <- "Absolute Lorenz curve ordinate (p - L(p))*mean"
    }
    if(!add)
      d <- ggplot(td, aes(x=p, y=y, group=id, colour=id)) + geom_line() +
        ## this stuff is not working properly. Think more carefully!
        geom_line(data=td1, col="black") +
        ylab(ylab) + xlab(xlab) +
        facet_grid(comp ~ group)
    else
      d <- geom_line(data=td, aes(x=p, y=y, colour=id))
    d
  }

## tip

# a plot method



#' Plot a TIP or CPG curve
#'
#' plot.tip plots a so-called "Three I's of Poverty" (TIP) or Cumulative
#' Poverty Gap (CPG) object.
#'
#' lines.tip add lines to an existing tip plot.
#'
#' plot.tipl plots a list of tip curves.
#'
#' Plots the TIP curve.
#'
#' @aliases plot.tip lines.tip plot.tipl
#' @param object An object of class tip.
#' @param add A logical. If TRUE, the curve is added to an existing plot.
#' @return A plot.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{tip}}
#' @references Jenkins, Stephen P & Peter J Lambert ( 1997), `Three `i's of
#' poverty curves, with an analysis of U.K. poverty trends', Oxford Economic
#' Papers 49(3), 317-327.
#' @examples
#'
#' plot(tip(runif(100)))
#'
#' @export
plot.tip <-
  function(x, add = FALSE,
           xlab = FALSE, ylab = FALSE,
           xlim = FALSE, ylim = FALSE, type = FALSE, ...)
{
            object <- x
            if(!is.tip(object)) stop("Not a tip curve!")
    ## use the as.data.frame
  td <- as.data.frame(object)
  if(missing(xlab))
    xlab <- "Population share"
  if(missing(ylab))
    ylab <- "Cumulative poverty gap"
  td$y <- td$ordinates
  if(!add)
    d <- ggplot(td, aes(x=p, y=y)) + geom_line() +
      ylab(ylab) + xlab(xlab)
  else
    d <- geom_line(data=td, aes(x=p, y=y))
  d
}

#' @export
lines.tip <- function(x, ...)
{
    plot.tip(x, add = TRUE, ...)
  }


## make a method to plot a list of Tip curves
#' @export
plot.tip_list <-  function(x, add = FALSE,
                       col = FALSE,
                       xlab = FALSE, ylab = FALSE,
                       xlim = FALSE, ylim = FALSE,
                       type = "l", legend = TRUE,
                        ...)
{
    object <- x
    if(!is.list(object))
      stop("Object is not a list!")
    if(!all(sapply(object, is.tip)))
      stop("Some element  is not a TIP curve!")
    td <- as.data.frame.tip_list(object)
    k <- length(unique(td$elname))
    if(missing(xlab))
      xlab <- "Population share"
    if(missing(ylab))
      ylab <- "Cumulative poverty gap"
    td$y <- td$ordinates
    if(!add)
      d <- ggplot(td, aes(x=p, y=y, group=elname, colour=elname)) + geom_line() +
        ylab(ylab) + xlab(xlab)
    else
      d <- geom_line(data=td, aes(x=p, y=y, colour=elname))
    d
  }
## the incdist method
#' @export
plot.tip_incdist <-  function(x, add = FALSE, lor.type = "ord",
                          diff = FALSE,
                          xlab = FALSE, ylab = FALSE,
                          ylim = FALSE, xlim = FALSE,
                          type = "l",
                          col = FALSE,
                          legend = TRUE,
                          ...)
{
    object <- x
    if(!inherits(object, "tip_incdist")) stop("Not an incdist summary object!")
    td <- as.data.frame(object)
    ## preserve the order of group, comp
    for(i in c("id", "comp", "group"))
      td[[i]] <- factor(td[[i]], labels=unique(td[[i]]))
    if(missing(xlab))
      xlab <- "Population share p"
    if(missing(ylab))
      ylab <- "Cumulative poverty gap"
    td$y <- td$ordinates
    if(!add)
      d <- ggplot(td, aes(x=p, y=y, group=id, colour=id)) + geom_line() +
        ylab(ylab) + xlab(xlab) +
        facet_grid( ~ group)
    else
      d <- geom_line(data=td, aes(x=p, y=y, colour=id))
    d
  }
