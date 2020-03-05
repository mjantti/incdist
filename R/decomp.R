## functions for decomposing inequality

## the generalized entropy class for within/between

## rudimentary

decomp.groups <- function(x, w = rep(1,length(x)), eta = 2,
               data = NULL,  ranked=x, na.rm = TRUE, group=NULL, func="ge",
                          ...)
  {
      if(!is.null(data) & is.data.frame(data))
      {
          attach(data)
          on.exit(detach(data))
      }
      ind <-
        do.call(func, list(x=x, w=w, eta=eta, na.rm=na.rm))
      ## now collapse
      xb <- tapply(x, list(g=group), mean, na.rm=TRUE)
      wb <- tapply(w, list(g=group), sum, na.rm=TRUE)
      ind.b <- do.call(func, list(x=xb, w=wb, eta=eta, na.rm=TRUE))
      ind.w <- ind - ind.b
      retval <- c(ind, ind.w, ind.b)
      retval
  }

## by income components (for Gini)
## CV^2
