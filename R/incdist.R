## test to see if both make it to NAMESPACE
# income -- R functions for income distribution analysis
# v 0.1 1998/01/07 mjantti@abo.fi
# Goals:
# This library will supply standard methods for income distribution
# research, including
# - (relative, absolute and generalized) inequality and poverty indices
# - (relative, absolute and generalized) inequality and poverty curves
# - statistical inference for the above
# - for microdata, grouped data and using kernel density objects
# - fitting distribution functions to grouped and microdata
# - robust methods
# v 0.2 2003/07/15 mjantti@abo.fi
## move to use new style methods
## .. abandoned
## v. 0.3 2012-12-03 markus.jantti@sofi.su.se

### an income distribution object


#' An income distribution object
#'
#' This function creates an object of class incdist that can be used to
#' facilitate income distribution analysis.
#'
#' The incdist object can be summarized, and the summary can be printed and
#' plotted (plotting is currently broken).
#'
#' @aliases incdist incdist.default incdist.formula summary.incdist
#' print.summary.incdist plot.summary.incdist
#' @param formula a one-or two-sided formula that defines the resource
#' variables of interest.
#' @param weights sampling weights.
#' @param data a data frame that holds the data.
#' @param subset defines a subset of the data to be used.
#' @param eqsale defines a transformation of each resource variable. Defaults
#' to NULL, that is, no transformation.
#' @param group a one-sided formula that defines the grouping structure.
#' @param part a one-sided formula that defines breaks in the data, for which
#' analyses are done separately.
#' @param idvar a one-sided formula that auxiliary variables to be kept in the
#' resulting object, for instance, variables needed for equivalence scale
#' transformation or deflating the analysis variables.
#' @param na.rm a logical indicating if only complete cases should be used.
#' @param func a character vector giving which income distribution functionals
#' to calculate.
#' @return An object of class incdist.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references
#' \insertRef{lambert1993}{incdist}
#'
#' @examples
#'
#' ## generate some data
#' ## income components
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- rexp(n)
#' a <- ifelse(a <- rpois(n, 1), a, 1)
#' c <- rpois(n, 2)
#' ## sum to total income
#' y1 <- x1 + x2
#' ## generate a grouping
#' g <- factor(rbinom(n, 1, .5), labels = letters[1:2])
#' ## generate a partitioning variable
#' p <- factor(sample(c(1, 2, 3), n, replace = TRUE),
#'             labels = paste("t", 1:3, sep = ""))
#' ## generate some weights
#' w <- rpois(n, 5)
#' ## put it all into a data frame
#' test.d <- data.frame(x1, x2, y1, g, p, a, c, w)
#'
#' summary(incdist(y1 ~ x1 +  x2, data=test.d, group = ~ g))
#' summary(incdist(y1 ~ x1 +  x2, data=test.d, part = ~ p))
#' id.0 <- incdist(y1 ~ x1 +  x2, data=test.d, part = ~ p, group = ~ g,
#'                 idvar = ~ a + c,  weights = w)
#' id.0$eqscale <- list(formula = ~ (a + k*c)^d, coef=list(k=.8, d=.75),
#'              type="citro-michael")
#'
#' print(summary(id.0, equivalise = FALSE), all = TRUE)
#' print(summary(id.0, equivalise = TRUE), all = TRUE) # the default
#' print(id.0)
#' ## there is not such function: plot(id.0)
#' @importFrom stats approxfun as.formula complete.cases cov.wt formula is.empty.model
#' model.extract model.weights na.omit predict terms var weighted.mean
#'


## the formula method, currently only this exists
## will this need to be modified
#' @export
incdist <- function(x, ...)
{
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("incdist", x)
}
#' @export
incdist.formula <-
  function(formula, weights, data = sys.frame(sys.parent()), subset,
           eqscale = NULL,
           group = ~ 1, part = ~ 1, idvar = ~ 0,
           na.rm = TRUE, ...)
  {
    if (!inherits(formula, "formula")){
      stop("First argument must be the income formula!")
    }
    # 1. Set up the data.
    # figure out the income part, the partitioning and the group structure
    # code copied  from locfit{locfit} &  sem{tsls.formula}
    income <- terms(formula, data = data)
    ## if there are no components, let the intercept be
    if(length(attr(income, "term.labels")))
      attr(income, "intercept") <- 0
    ## allow for no response (in the future...)
    grterms <- terms(group, data = data)
    if(length(attr(grterms, "term.labels")))
      attr(grterms, "intercept") <- 0
    paterms <- terms(part, data = data)
    ## here also, if paterms is ~1 (which should be the default),
    ## keep the intercept
    if(length(attr(paterms, "term.labels")))
      attr(paterms, "intercept") <- 0
    idterms <- terms(idvar, data = data)
    attr(idterms, "intercept") <- 0
    c1 <- as.character(formula)
    c2 <- as.character(group)
    c3 <- as.character(part)
    c4 <- as.character(idvar)
    ## how to get weights named correctly
    formula <- as.formula(paste(c1[2], c1[1], c1[3], "+", c2[2], "+", c3[2],
                                "+", c4[2]))
    m <- match.call()
    m$formula <- formula
    m$eqscale <- m$group <- m$part <- m$method <- m$idvar <- NULL
    m[[1]] <- as.name("model.frame")
    frm <- eval(m, sys.frame(sys.parent()))
    w <- model.weights(m)
    incnames <- as.character(attributes(income)$variables)[-1]
    ## initialize at NULL
    ## should panames and grnames be coerced to be factors here?
    yname <- panames <- grnames <- idvarnames <- NULL
    if (attr(income, "response"))
      {
        incnames <- incnames[-1]
        yname <- deparse(formula[[2]])
      }
    if(!is.empty.model(part))
      {
        panames <- as.character(attributes(paterms)$variables)[-1]
        for(i in panames) frm[[i]] <- as.factor(frm[[i]])
      }
    if(!is.empty.model(group))
      {
        grnames <- as.character(attributes(grterms)$variables)[-1]
        for(i in grnames) frm[[i]] <- as.factor(frm[[i]])
      }
    if(!is.empty.model(idvar))
      {
        idvarnames <- as.character(attributes(idterms)$variables)[-1]
      }
    ## start constructing the return object.
    ## return a structure with class incdist
    ## which contains the data etc
    ## data, formula, yname, xnames, part names, groupnames
    ret <- list(frm, formula, eqscale,
                     yname, incnames, panames, grnames, idvarnames,
                     group, part, idvar, m)
    names(ret) <- c("frm", "formula", "eqscale", "yname",
                    "incnames", "panames", "grnames",
                    "idvarnames", "group",
                    "part", "idvar", "m")
    class(ret) <- "incdist"
    ret
  }

print.incdist <- function(object, ...)
  {
    if(!is.incdist(object)) stop("Not an incdist object!")
    cat("An incdist object\n")
    cat("Overall variable: \t", object$yname, "\n")
    cat("Components: \t", object$incnames, "\n")
    cat("Grouping (within part): \t", object$grnames, "\n")
    cat("Partitioned by: \t", object$panames, "\n")
    invisible(object)
  }

## summary.incdist here
## part is the partitioning factor
## eqscale should be a formula? that is applied to every element in income
## method should specify an index
## is the first argument a formula?
## group if the grouping factor(s)
## NB: problem with Gini coefficient use
summary.incdist <- function(object,
                            equivalise = FALSE,
                            func = c("weighted_mean", "weighted_std",
                              "gini.default", "concentration.coef"),
                            poverty=FALSE, concentration = !poverty,
                            povertyline=NULL,
                            povertyline.function="weighted_median",
                            povertyline.fraction=0.5,
                            frame=TRUE, ...)
  {
    if (!inherits(object, "incdist")){
      stop("First argument must be the incdist object!")
    }
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
    ## see lorenz
    ##on.exit(## detach(object))
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
    if(frame)
        old.frm <- frm
    else
        old.frm <- ""
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
    ret.x <- list() ##ret.y{j:incomecomps}{i:group}{l:partition}
    ## skip the first incnames, which is in "ret.y"
    if(length(incnames))
      {
        for(j in 1:length(incnames))
          {
            ret.x[[j]] <- list()
            for(i in 1:(1+length(grl)))
              {
                ret.x[[j]][[i]] <- list()
                for(l in 1:length(frm.list))
                  {
                    ret.x[[j]][[i]][[l]] <- NA
                  }
            names(ret.x[[j]][[i]]) <- pal
              }
            names(ret.x[[j]]) <- c("all", grl)
          }
        names(ret.x) <- incnames
      }
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
        ## this used to be
        for(i in 1:(1+length(grl))) ## i indexes groups present. The first is all.
          {
            doing.name <- grnames[i]
            ## must add som code here to
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
            ## figure out the poverty line if any of the functions is poverty
            ## only for group 1 ("all")
            ## this and the next place it is used needs to be revisited!
            if(poverty & i==1)
              {
                if(is.null(povertyline))
                  {
                    if(!is.null(weights))
                      poverty.line <-
                        povertyline.fraction*
                          do.call(povertyline.function,
                                  list(x=as.numeric(y), w=weights))
                    else
                      poverty.line <-
                        povertyline.fraction*
                          do.call(povertyline.function, list(x=as.numeric(y)))
                  }
                else
                  poverty.line <- frm.list.l[[i]][[povertyline]]
              }
            ## "func" holds the functions to be calculated
            ## these must accept two (an exactly two) arguments
            ## the data and the weights
            tmp.res <- list()
            for(k in 1:length(func))
              {
                ## initialize to 0 for the case the group is empty
                tmp.res[[k]] <- 0
                if (!count.grnames[i]) next
                if(!is.null(weights))
                  {
                    if(length(grep("poverty", func[k]))>0 & poverty)
                      {
                        if(length(grep("absolute", func[k]))>0)
                          tmp.res[[k]] <-
                          do.call(func[k], list(as.vector(y), w = weights,
                                                ...))
                        else
                          tmp.res[[k]] <-
                            do.call(func[k], list(as.vector(y), w = weights,
                                                  z = as.numeric(poverty.line),
                                                  ...))
                      }

                    else
                      tmp.res[[k]] <-
                        do.call(func[k], list(as.vector(y), w = weights, ...))
                  }
                else
                  tmp.res[[k]] <- do.call(func[k], list(as.vector(y)), ...)
              }
            names(tmp.res) <- func
            ret.y[[i]][[l]] <- tmp.res
            if(length(incnames))
              {
                for(j in 1:length(incnames))
                  {
                    ## "func" holds the functions to be calculated
                    ## these must accept two (an exactly two) arguments
                    ## the data and the weights
                    tmp.res <- list()
                    for(k in 1:length(func))
                      {
                        ## initialize to zero and check if groups is empty
                        tmp.res[[k]] <- 0
                        if (!count.grnames[i]) next
                        if(!is.null(weights))
                          {
                            if(func[k] == "concentration.coef" || concentration)
                              {
                                tmp.res[[k]] <-
                                  do.call(func[k],
                                          list(x[,j], w = weights,
                                               ranked = as.vector(y), ...))
                              }
                            else
                              {
                                if(func[k] == "relative.poverty" || poverty)
                                  tmp.res[[k]] <-
                                    do.call(func[k],
                                            list(x[,j],
                                                 w=weights,
                                                z = as.numeric(poverty.line),
                                                 ...))
                                else
                                  tmp.res[[k]] <-
                                    do.call(func[k],
                                            list(x[,j], w=weights, ...))
                              }
                          }
                        else
                          {
                            if(func[k] == "concentration.coef" || concentration)
                              tmp.res[[k]] <-
                                do.call(func[k],
                                        list(x[,j], ranked = as.vector(y)),
                                        ...)
                            else
                              {
                                if(func[k] == "relative.poverty" || poverty)
                                  tmp.res[[k]] <-
                                    do.call(func[k], list(x[,j], ...))
                                else
                                  tmp.res[[k]] <-
                                    do.call(func[k], list(x[,j], ...))
                              }
                          }
                      }
                    names(tmp.res) <- func
                    ret.x[[j]][[i]][[l]] <- tmp.res
                  }
              }
          }
      }
    ## start constructing the return object.
    ## detach(object)
    if(frame)
      {
        ret <- list(ret.y, ret.x,
                    n, sumw,
                    object$yname, object$incnames,
                    object$grnames, object$panames,
                    func,
                    object$frm, old.frm)
        names(ret) <- c(object$yname, "components",
                        "sample.size", "sum.weights",
                        "overall", "comp",
                        "group", "partition",
                        "func", "frm", "old.frm")
      }
    else
      {
        ret <- list(ret.y, ret.x,
                    n, sumw,
                    object$yname, object$incnames,
                    object$grnames, object$panames,
                    func)
        names(ret) <- c(object$yname, "components",
                        "sample.size", "sum.weights",
                        "overall", "comp",
                        "group", "partition",
                        "func")
      }
    structure(ret, class = c("summary.incdist", class(object)))
  }



## do a print.summary function for incdist

print.summary.incdist <- function(object, all = FALSE, what = TRUE,
                                  relative = FALSE, ...)
  {
    if(!inherits(object, "summary.incdist")) stop("Not an incdist summary object!")
    ##
    comps <- object$comp
    if(what == TRUE)
      funcs <- object$func
    else
      funcs <- what
    if(length(object$partition))
      parts <- levels(object$frm[[object$partition]])
    else parts <- 1
    if(length(object$group))
      groups <- c("all", levels(object$frm[[object$group]]))
    else groups <- 1
    stats <- list()
    for(k in 1:length(funcs))
      {
        ## this assumes that the funcs have returned scalars for each part, comps
        stats[[k]] <- list()
        ##
        for(j in 1:(1+length(comps)))
          {

            stats[[k]][[j]] <- matrix(0, nrow = length(parts),
                                      ncol = length(groups))
            rownames(stats[[k]][[j]]) <- parts
            colnames(stats[[k]][[j]]) <- groups
            for(l in 1:length(parts))
              {
                for(i in 1:length(groups))
                  {
                    if(j==1)
                      stats[[k]][[j]][l, i] <- object[[1]][[i]][[l]][[k]]
                    else
                      stats[[k]][[j]][l, i] <- object[[2]][[j-1]][[i]][[l]][[k]]
                    if(relative == TRUE)
                      {
                        if(j==1)
                          stats[[k]][[j]][l, i] <-object[[1]][[i]][[l]][[k]]/
                            object[[1]][[1]][[l]][[k]]
                       else
                        stats[[k]][[j]][l, i] <- object[[2]][[j-1]][[i]][[l]][[k]]/
                          object[[2]][[j-1]][[1]][[l]][[k]]
                      }
                  }
              }
          }
        names(stats[[k]]) <- c(object$overall, comps)
      }
    ## this is not very elegant, but not doing it leads to problems
    stats[[k+1]] <- list()
    stats[[k+1]][[1]] <- matrix(0, nrow = length(parts),
                                ncol = length(groups))
    rownames(stats[[k+1]][[1]]) <- parts
    colnames(stats[[k+1]][[1]]) <- groups
    stats[[k+2]] <- list()
    stats[[k+2]][[1]] <- matrix(0, nrow = length(parts),
                                ncol = length(groups))
    rownames(stats[[k+2]][[1]]) <- parts
    colnames(stats[[k+2]][[1]]) <- groups
    for(l in 1:length(parts))
      {
        for(i in 1:length(groups))
          {
            stats[[k+1]][[1]][l, i] <- object[["sample.size"]][[i]][[l]]
            stats[[k+2]][[1]][l, i] <- object[["sum.weights"]][[i]][[l]]
            if(relative == TRUE)
              stats[[k+2]][[1]][l, i] <- object[["sum.weights"]][[i]][[l]]/
                object[["sum.weights"]][[1]][[l]]
          }
      }
    names(stats) <- c(funcs, "sample.size", "sum.weights")
    retval <- list(stats, funcs, parts, c(object$overall, comps), groups)
    names(retval) <-
      c("result.matrices", "functions", "partitions", "variables", "groups")
    #structure(retval, class = c("summary.incdist", class(x)))
    cat("An incdist.summary object\n")
    cat("Variables: ", retval$variables, "\n")
    cat("Partitions: ", retval$partitions, "\n")
    cat("Statistics: ", retval$functions, "\n")
    if(all)
      {
        cat("The matrices: \n")
        for(k in 1:(length(retval$result.matrices)-2)) ## indexes functions
          {
            for(j in 1:(1+length(comps)))
              {
                cat("Statistic: ", retval$functions[k],
                    "Income component: ", retval$variables[j], "\n")
                print(retval$result.matrices[[k]][[j]])
              }
          }
        cat("Observations: \n")
        print(retval$result.matrices[[k+1]][[1]])
        cat("Sum of weights: \n")
        print(retval$result.matrices[[k+2]][[1]])
        cat("\n")
      }
    invisible(retval)
  }

## I need better methods to display my results
## the summary object has results for the top-level variables in [[1]]
## and all components in [[2]]

incdist.as.array <- function(object, ...)
  {
    ##if(!inherits(object, "summary.incdist"))
    ##  stop("Not an incdist summary objectect!")
    fun <- object$functions
    var <- object$variables
    par <- object$partitions
    gro <- object$groups
    ## make a large array that holds all
    ret <- array(unlist(object[["result.matrices"]]),
                 dim=c(length(par),
                   length(gro),
                   length(var),
                   length(fun)),
                 dimnames=
                 list(partitions=par,
                      groups=gro,
                      variables=var,
                      functions=fun))
    ret
  }

## the predicate function

is.incdist <- function(object) inherits(object, "incdist")

#' Transform an incdist object with an equivalence scale
#'
#' A utility function used by summary.incdist to transform "raw" income
#' variables by an equivalence scale
#'
#' Presuposes the present in object of a list called eqscale, with (at least)
#' components formula, coef. Variables in formula must be found in the
#' object\$frm.
#'
#' @param object an incdist object.
#' @return Returns an object where the income variables in frm have been
#' replaced the equivalised version and the old unequavalised form has been
#' copied old.frm.
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso
#' @references
#' @examples
#'
#' x1 <- runif(50)
#' x2 <- rexp(50)
#' ad <- ifelse(ad <- rpois(50, 1), ad, 1)
#' ch <- rpois(50, 2)
#' ## sum to total income
#' y1 <- x1 + x2
#' ## generate a grouping
#' g <- factor(c(rep(1, 10), rep(2, 10), rep(1, 15), rep(2, 15)),
#'             labels = letters[1:2]) # 2 groups
#' ## generate a partitioning variable
#' p <- factor(c(rep(1, 20), rep(2, 30)), labels = LETTERS[1:2])
#' ## generate some weights
#' w <- rpois(50, 5)
#' ## put it all into a data frame
#' test.d <- data.frame(x1, x2, y1, g, p, ad, ch, w)
#'
#' id.0 <- incdist(y1 ~ x1 +  x2, part = ~ p, group = ~ g, weights = w,
#' data = test.d)
#' id.0$idvarnames <- c("a", "c")
#' id.0$eqscale <- list(formula = ~ (ad + k*ch)^d, coef=list(k=.8, d=.75),
#'              type="citro-michael")
#'
#' eqscale(id.0)
#'
#' @export

eqscale <- function(object)
  {
    .tmp.eqscale <-   function (object)
      {
        if(!is.incdist(object)) stop("Not an incdist objectect!")
        if(is.null(object$eqscale)) stop("No equivalence scale present in objectect!")
        form <- object$eqscale$form
        coef <- object$eqscale$coef
        data <- object$frm
        thisEnv <- environment()
        env <- new.env()
        for (i in names(data)) {
          assign(i, data[[i]], envir = env)
        }
        if(!is.null(coef))
          ind <- as.list(coef)
        else ind <- NULL
        parLength <- 0
        for (i in names(ind)) {
          temp <- coef[[i]]
          assign(i, temp, envir = env)
        }
        eval(form[[2]], envir = env)
      }

    if(!is.incdist(object)) stop("Not an incdist objectect!")
    ## this is probably not needed?
    ##attach(object)
    data <- object$frm
    tmp.es <- .tmp.eqscale(object)
    ##env <- new.env() # not needed, I think
    ## this gives a waring on R CMD check. Need to define the explicitly
    for(i in c(object$yname, object$incnames)) data[[i]] <- data[[i]]/tmp.es
    data
  }


## functions that "bind" to income distribution functionals
## an inequality index
#'
#' Estimate inequality indices
#'
#' This function estimates inequality indices.
#'
#'
#' @aliases inequality inequality.incdist inequality.default
#' @param x a numerical vector or incdist object
#' @param w an optional vector of non-negative integer values weights.
#' @param indices to estimate
#' @return the estimated indices
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{gini}}, \code{\link{ge}}, \code{\link{atkinson}}
#' @references Lambert, P. (1993).  \emph{The distribution and
#'     redistribution of income. A mathematical analysis.} Manchester
#'     University Press, Manchester.

#' @export

inequality <- function(x, ...)
{
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("inequality", x)
}
#' @export
inequality.default <- function(x, type = c("gini", "ge", "atkinson", "cv2"), ...)
{
  if(!is.numeric(x)) stop("argument x is non-numeric!")
  args <- list(...)
  ret <- numeric(length(type))
  names(ret) <- type
  for(i in type) ret[i] <- do.call(i, c(list(x), args))
  ret
}
#' @export
inequality.incdist <-
  function(x, what = 1, type = c("gini", "ge", "atkinson", "cv2"),
           frame=FALSE,
           equivalise=FALSE, ...)
{
    object <- x
    if(!is.incdist(object)) stop("Not an incdist object!")
    ## for what
    ## what to calculate  (change this so multiple measures can be used.)
    ## type <- match.arg(type)
    obj <- summary(object, func = type, frame=frame, equivalise=equivalise,...)
    ##print(obj, all = TRUE)
    obj[[what]][[1]]
}


## a poverty index
#'
#' Estimate poverty indices
#'
#' This function estimates poverty indices.
#'
#'
#' @aliases povery poverty.incdist poverty.default
#' @param x a numerical vector or incdist object
#' @param w an optional vector of non-negative integer values weights.
#' @param indices to estimate
#' @return the estimated indices
#' @author Markus Jantti \email{markus.jantti@@iki.fi}
#' @seealso \code{\link{inequality}}, \code{\link{fgt}}
#' @references Lambert, P. (1993).  \emph{The distribution and
#'     redistribution of income. A mathematical analysis.} Manchester
#'     University Press, Manchester.
#' @export
poverty <- function(x, ...)
{
  if(is.null(class(x))) class(x)  <- data.class(x)
  UseMethod("poverty", x)
}

#' @export
poverty.default <- function(x, type = "fgt", ...)
{
  if(!is.numeric(x)) stop("argument x is non-numeric!")
  args <- list(...)
  ret <- numeric(length(type))
  names(ret) <- type
  for(i in type) ret[i] <- do.call(i, c(list(x), args))
  ret
}

#' @export
poverty.incdist <-
  function(x, what = 1, type = "fgt", frame=FALSE, ...)
{
    object <- x
  if(!is.incdist(object)) stop("Not an incdist object!")
  ## for what
  ## what to calculate
  ## type <- match.arg(type)
  obj <- summary(object, func = type, poverty=TRUE, frame=frame, ...)
  ##print(obj, all = TRUE)
  obj[[what]][[1]]
}


## make a "as.data.frame" method
as.data.frame.incdist <- function(object, all=TRUE, what=TRUE,
                                  relative=FALSE, ...)
  {
    if(!inherits(object, "summary.incdist")) stop("Not an incdist summary object!")
    vars <- c(object$overall, object$comp)
    groupdim <- object$group
    grouplevs <- names(table(object$frm[[groupdim]]))
    ngroup <- length(grouplevs)
    partdim <- object$partition
    partlevs <- names(table(object$frm[[partdim]]))
    npart <- length(partlevs)
    func <- object$func
    nfunc <- length(func)
    ## which approach?
    ## constructing "by hand" may be unreliable

    #ret1 <- as.numeric(unlist(object[1:2]))
    #nstat <- length(ret1)
    #Variable <- rep(vars, each=ngroup*npart*nfunc)
    #Group <- rep()
    #Partition <- rep()
    #Function <- rep()

    #ret <-
    #  data.frame(val=ret1,
    #             Variable=rep(vars, each=ngroup*npart*nfunc),
    #             Group=rep(),
    #             Partition=rep(),
    #             Function=rep())

    ## or create a dataset by unlisting the 1st and 2nd compont
    ## Care needs to be taken to deal with different kind of objects
    ## 1. z ~ y + x (works)
    ## 2. z ~ 1 (works)
    ## 3. ~ x + y
    if(length(object[[1]])>0)
      {
        tx1 <- unlist(object[1])
        ne <- length(strsplit(names(tx1[1]), "\\.")[[1]])
        td1 <- data.frame(val=as.numeric(tx1))
        ## a temporary fix
        tn <- unlist(strsplit(names(tx1), "\\.", perl=TRUE))
        ##tn <- tn[-grep("%", tn)]
        td2 <-
          as.data.frame(matrix(tn, ncol=ne, byrow=T), stringsAsFactors = FALSE)
      }
    if(length(object[[2]])>0)
      {
        tx3 <- unlist(object[2])
        td3 <- data.frame(val=as.numeric(tx3))
        ## a temporary fix
        tn <- unlist(strsplit(names(tx3), "\\.", perl=TRUE))
        ##tn <- tn[-grep("%", tn)]
        td4 <-
          as.data.frame(matrix(tn, ncol=ne+1, byrow=T)[,-1], stringsAsFactors = FALSE)
      }
    if(length(object[[1]])>0 && length(object[[2]])>0)
      ret <- rbind(cbind(td1, td2), cbind(td3, td4))
    if(length(object[[1]])>0 && !length(object[[2]])>0)
      ret <- cbind(td1, td2)
    if(!length(object[[1]])>0 && length(object[[2]])>0)
      ret <- cbind(td3, td4)
    ## this *must* be generalised at some point
    ret$Variable <- ret$V1
    ret$Group <- ret$V2
    ret$Partition <- ret$V3
    ret$Function <- paste(ret$V4, ret$V5, sep=".")
    if(ne>5) {ret$Value <- ret$V6; ret$V6 <- NULL}
    ret$V1 <- ret$V2 <- ret$V3 <- ret$V4 <- ret$V5 <- NULL
    ret
  }

## design and write a "c()" method for two incdist objects
