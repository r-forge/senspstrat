## Function to calculate w given by equation
##
##           1
## ----------------------
##   - beta y - alpha
## %e                 + 1
.calc.w <- function(alpha, beta.y)
  1L/(1L + exp(-alpha - beta.y))

## Need to solve equaition 1 for alpha
##
##  n
## ====            dF
## \                 i
##  >    ----------------------- - C                                         (1)
## /       - beta y  - alpha
## ====            i
## i = 1 %e                  + 1
##
## to do so find the minimum of the integral of the equation 1 shown in
## equation 2
##
##  n
## ====            beta y  + alpha
## \                     i
##  >    dF  log(%e                + 1) - alpha C                            (2)
## /       i
## ====
## i = 1
##
.alpha.est <- function(alpha, beta.y, dF, C)
  sum(log(1L + exp(alpha + beta.y)) * dF) - alpha*C

.calc.alphahat <- function(beta.y, dF, C) {
  alphahat <- optimize(f=.alpha.est, interval=c(-100L, 100L), beta.y=beta.y,
                       dF=dF, C=C)$minimum
  
  if(alphahat > 90 || alphahat < -90) {
    warning("optimize overflow alphahat value invalid")
  }
  
  alphahat
}

.calc.ecdf <- function(x) {
  n <- length(x)
  
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")

  vals <- sort(unique(x))

  index <- match(x, vals)
  val.count = tabulate(index)
  return(list(F=cumsum(val.count)/n, vals=vals, index=index))
}

.foldUpperTri <- function(x) {
  xrow <- row(x)
  xcol <- col(x)
  xnrow <- nrow(x)
  upper <- xrow < xcol

  x[(xrow[upper] - 1)*xnrow + xcol[upper]] <- x[upper]
  return(x)
}

".sumCrossUpperTri<-" <- function(x, na.rm=TRUE, value) {
  xrow <- row(x)
  xcol <- col(x)
  indx <- xrow <= xcol
  
  x[indx] <- rowSums(value[xcol[indx],]*value[xrow[indx],], na.rm=na.rm)
  x
}
