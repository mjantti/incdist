# the fields and ok mobility measure



#' Fields mobility index
#'
#' This function estimates the Fields income mobility index.
#'
#'
#' @param x A matrix with income data
#' @param w A vector of weights
#' @param data A data frame in which to look for the data.
#' @param na.rm A logical, indicating whether or not to remove NAs
#' @return A scalar.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references Cite some suitable paper by Fields.
#' @examples
#'
#'
#' @export fields
fields <- function(x, w = rep(1,dim(x)[1]), data = list(), na.rm=TRUE) {
  # missing etc treatment
  n <- dim(x)[1]
  k <- dim(x)[2]
  incmat <- clean.income(x, w)
  x <- incmat[,-(k+1)]
  w <- incmat[,(k+1)]
  ## this must be wrong
  diff.xy <- abs(log(x[,1]) - log(x[,k]))
  retval <- weighted.mean(diff.xy, w)
  retval
}

# the shorrocks type measure (for various inequality indices)

#' Shorrocks mobility index
#'
#' This function estimates the Shorrocks income mobility index for a specific
#' inequality indext
#'
#' @param x A matrix with income data
#' @param w A vector of weights
#' @param data A data frame in which to look for the data
#' @param na.rm A logical, indicating whether or not to remove NAs
#' @param index An inequality index.
#' @return A scalar.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references Cite Shorrocks's Econometrica paper
#' @examples
#'
#'
#' @export shorrocks

shorrocks <- function(x, w = rep(1,dim(x)[1]),
                      data = NULL, na.rm = TRUE,
                      index = "gini", ...)
{
    ## attach a possible data frame. Remember to detach it!
    if(!is.null(data) & is.data.frame(data))
    {
        attach(data)
        on.exit(detach(data))
    }
    ## missing etc treatment
    n <- dim(x)[1]
    k <- dim(x)[2]
    if(!k>1) stop("x has only one column!")
    incmat <- clean.income(x, w, na.rm = na.rm)
    #XXX:  <- incmat[,-(k+1)]
    w <- incmat[,(k+1)]
    x.bar <- apply(x, 1, mean)
    ## check if weights are probs or unit correspondences
    ## not implemented and not really needed now.
    index.x <- apply(x, 2, function(y) do.call(index, list(y, w, ...)))
    mean.#XXX:  <- apply(x, 2, function(y) weighted.mean(y, w))
    ## x.bar
    index.xbar <- do.call(index, list(x.bar, w, ...))
    mean.xbar <- weighted.mean(x.bar, w)
    ## and for m
    m <- 1 - index.xbar/sum(mean.x/mean.xbar*index.x)
    ## if(!is.null(data) & is.data.frame(data))
    ##   detach(data)
    ret.val <- list()
    ret.val$index <- index.x
    ret.val$mean <- mean.x
    ret.val$index.xbar <- index.xbar
    ret.val$mean.xbar <- mean.xbar
    ret.val$m <- m
    ret.val$m
    ret.val
}


## matrix based stuff

#' A mobility matrix
#'
#' This function a mobility matrix for two-dimensional data.
#'
#' @aliases mtr mdet ml mf mb mmtr mmdet mml mmf nnb
#'
#' @param x1 A vector of income data
#' @param x2 A vector of income data
#' @param w A vector of weights
#' @param Q1 Vector of cut offs to be applied to x1
#' @param Q2 Vector of cut offs to be applied to x2
#' @param q Number of classes
#' @param data A data frame in which to look for the data
#' @param na.rm A logical, indicating whether or not to remove NAs
#' @param index An inequality index.
#' @return A list including the conditional and unconditional matrices (EDIT).
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references Cite something suitable. (EDIT)
#' @examples
#'
#'
#' @export mob.mat

mob.mat <-
  function(x1, x2, w, Q1 = NULL, Q2 = NULL, q = 5, labels = FALSE)
{
  n01 <- length(x1)
  n02 <- length(x2)
  if (n01 != n02)
    stop("x1 and x2 not of equal length!")
  if(missing(w)) w <- rep(1, length(x1))
  df <- data.frame(x1, x2, w)
  n.orig <- dim(df)[1]
  df <- subset(df,
               !is.nan(x1) & !is.nan(x2) & !is.nan(w) &
               is.finite(x1) & is.finite(x2) & is.finite(w) &
               !is.na(x2) & !is.na(x2) & !is.na(w))
  if(missing(Q1)) Q1 <- weighted.quantile(df$x1, df$w, 0:q/q)
  X1 <- cut(df$x1, Q1, include.lowest=TRUE, right = FALSE, labels = labels)
  if(missing(Q2)) Q2 <- weighted.quantile(df$x2, df$w, 0:q/q)
  X2 <- cut(df$x2, Q2, include.lowest=TRUE, right = FALSE, labels = labels)
  n1 <- wtd.table(X1, df$w, type = 'table')
  n2 <- wtd.table(X2, df$w, type = 'table')
  n <- dim(df)[1]
  ## how should this be done in a weighted way?
  m <- weighted.crosstable(X1,X2,df$w)
  rm <- m/sum(n1)
  p <- m/rep(n1,nrow(m))
  p1 <- n1/sum(n1)
  p2 <- n2/sum(n1)
  list(p=p, m=m, rm=rm, n01=n01, n02=n02, p1=p1, p2=p2,
       n1=n1, n2=n2, n=n, n.orig = n.orig, Q1 = Q1, Q2 = Q2)
}


## matrix mobility measures. move these soon (ha!) to package incdist
## these are for mobility (i.e., conditional) matrices
mtr <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    ret <- (k[1] - sum(diag(mat)))/(k[1]-1)
    ret
  }

mdet <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    ret <- 1 - abs(det(mat))^(1/(k[1]-1))
    ret
  }

ml <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    ret <- 1 - abs(eigen(mat, only.values = TRUE)$values[2])
    ret
  }
mf <-  function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    lim.dist <- rep(1/k[1], k[1])
    ret <- 1 - 1/k[1]^2*sum(abs(mat/(1/k[1]) - 1))
    ret
  }

mb <-  function(mat, vec, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    if(k[1] != length(vec))
      stop("Vector length not equal to matrix dim!")
    i <- 1:k[1]
    tmp.mat1 <- mat * vec
    tmp.mat2 <- abs(outer(i, i, "-"))
    ret <- sum(mat * vec * abs(outer(i, i, "-")))
    ret
  }
## rewrite in terms of the (unconditional) frequencies for the real data
mmtr <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    mat <- mat/rep(apply(mat, 1, sum), k[1])
    ret <- (k[1] - sum(diag(mat)))/(k[1]-1)
    ret
  }

mmdet <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    mat <- mat/rep(apply(mat, 1, sum), k[1])
    ret <- 1 - abs(det(mat))^(1/(k[1]-1))
    ret
  }

mml <- function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    mat <- mat/rep(apply(mat, 1, sum), k[1])
    ret <- 1 - abs(eigen(mat, only.values = TRUE)$values[2])
    ret
  }
mmf <-  function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    lim.dist <- rep(1/k[1], k[1])
    mat <- mat/rep(apply(mat, 1, sum), k[1])
    ret <- 1 - 1/k[1]^2*sum(abs(mat/(1/k[1]) - 1))
    ret
  }

mmb <-  function(mat, ...)
  {
    if(!is.matrix(mat)) stop("argument is not a matrix!")
    k <- dim(mat)
    if(k[1] != k[2]) stop("Not a square matrix!")
    i <- 1:k[1]
    mat <- mat/sum(mat)
    ret <- sum(mat * abs(outer(i, i, "-")))
    ret
  }
