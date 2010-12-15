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

print.sensitivity.0d <- function(x, ...) {
  labs <- attr(x, "parameters")

  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("ACE:\t", paste("E(Y(", labs$z1, ") - Y(", labs$z0,") | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  ACE <- x$ACE
  print(ACE)

  ci.dim <- which(names(dimnames(x$ACE.ci)) == "ci.method")
  ci.method <- dimnames(x$ACE.ci)[[ci.dim]]
  ci.slice <- slice.index(x$ACE.ci, ci.dim)
  newdim <- dim(x$ACE.ci)[-ci.dim]
  newdimnames <- dimnames(x$ACE.ci)[-ci.dim]

  cat("\nACE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(array(x$ACE.ci[ci.slice == which(ci.method == "analytic")],
                dim=newdim, dimnames=newdimnames))
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
    cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(array(x$ACE.ci[ci.slice == which(ci.method == "bootstrap")],
                dim=newdim, dimnames=newdimnames))
  }

  invisible(NULL)
}

print.sensitivity.1d <- function(x, ...) {
  labs <- attr(x, "parameters")

  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("SCE:\t", paste("E(Y(", labs$z1, ") - Y(", labs$z0,") | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  SCE <- x$SCE
  print(SCE)

  ci.dim <- which(names(dimnames(x$SCE.ci)) == "ci.method")
  ci.method <- dimnames(x$SCE.ci)[[ci.dim]]
  ci.slice <- slice.index(x$SCE.ci, ci.dim)
  newdim <- dim(x$SCE.ci)[-ci.dim]
  newdimnames <- dimnames(x$SCE.ci)[-ci.dim]

  cat("\nSCE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(array(x$SCE.ci[ci.slice == which(ci.method == "analytic")],
                dim=newdim, dimnames=newdimnames))
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
    cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(array(x$SCE.ci[ci.slice == which(ci.method == "bootstrap")],
                dim=newdim, dimnames=newdimnames))
  }

  invisible(NULL)
}

plot.sensitivity.1.0d <- function(x, xlim, ylim,
                                  xlab=expression(beta), ylab='ACE',
                                  display = c("analytic", "bootstrap"),
                                  col='black', line.col=col, point.col=col,
                                  analytic.col="red",
                                  analytic.line.col=analytic.col,
                                  analytic.point.col=analytic.col,
                                  bootstrap.col="green",
                                  bootstrap.line.col=bootstrap.col,
                                  bootstrap.point.col=bootstrap.col,
                                  panel.last=NULL,
                                  type='l', ...) {

  display <- match.arg(display, several.ok=TRUE)

  sortIndx <- sort.list(x$beta)
  beta <- x$beta[sortIndx]
  ACE <- x$ACE[sortIndx]
  ACE.ci <- x$ACE.ci[sortIndx,,, drop=FALSE]
  
  finIndx <- is.finite(beta)
  infIndx <- is.infinite(beta)
  
  beta.fin <- beta[finIndx]
  ACE.fin <- ACE[finIndx]
  ACE.ci.fin <- ACE.ci[finIndx,,, drop=FALSE]

  beta.inf <- beta[infIndx]
  ACE.inf <- ACE[infIndx]
  ACE.ci.inf <- ACE.ci[infIndx,,, drop=FALSE]

  if(missing(ylim)) {
    ylim <- range(c(ACE, ACE.ci))
  }

  if(missing(xlim)) {
    xlim <- range(beta, finite=TRUE)
  }
  
  indx <- match(beta.inf, c(-Inf, Inf))

  plot.default(x=beta.fin, y=ACE.fin, xlim=xlim, ylim=ylim, type=type,
               xlab=xlab, ylab=ylab, ...)
  
  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
  points(x=inf.x, y=ACE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$ACE.var)) {
    for(i in seq_len(dim(ACE.ci.fin)[[2]])) {
      lines(x=beta.fin, y=ACE.ci.fin[, i, 'analytic'], lty=2,
            col=analytic.line.col)
    }
    points(x=rep.int(inf.x, times=dim(ACE.ci.fin)[[2]]),
           y=ACE.ci.inf[,indx, "analytic"], pch=3, col=analytic.point.col)
  }

  if('bootstrap' %in% colnames(x$ACE.var) && 'bootstrap' %in% display) {
    for(i in seq_len(dim(ACE.ci.fin)[[2]])) {
      lines(x=beta.fin, y=ACE.ci.fin[,i,"bootstrap"], lty=2,
            col=bootstrap.line.col)
    }
    points(x=rep(inf.x, times=dim(ACE.ci.fin)[[2]]),
           y=ACE.ci.inf[,indx,"bootstrap"], pch=3, col=bootstrap.point.col)
  }
}

plot.sensitivity.2.0d <- function(x, xlim, ylim,
                                  xlab=expression(beta[0]),
                                  ylab=expression(beta[1]),
                                  display = c("analytic", "bootstrap"),
                                  col=c(gray(.9),gray(1),gray(.8)),
                                  panel.last=NULL, ...) {

  display <- match.arg(display, several.ok=FALSE)

  sortIndx0 <- sort.list(x$beta0)
  sortIndx1 <- sort.list(x$beta1)
  
  beta0 <- x$beta0[sortIndx0]
  beta1 <- x$beta1[sortIndx1]
  
  ACE <- x$ACE[sortIndx0,sortIndx1,]
  
  ACE.ci <- x$ACE.ci[sortIndx0,sortIndx1,,,, drop=FALSE]

  reject <- ifelse(ACE < ACE.ci[,,,1,display],
                  -1,
                  ifelse(ACE > ACE.ci[,,,2,display], 
                         1,
                         0))

  finIndx0 <- is.finite(beta0)
  infIndx0 <- is.infinite(beta0)

  finIndx1 <- is.finite(beta1)
  infIndx1 <- is.infinite(beta1)

  beta0.fin <- beta0[finIndx0]
  beta1.fin <- beta1[finIndx1]

  reject.fin <- reject[finIndx0,finIndx1,, drop=FALSE]

  beta0.inf <- beta0[infIndx0]
  beta1.inf <- beta1[infIndx1]

  reject.inf <- reject[infIndx0,infIndx1,,drop=FALSE]

  if(missing(ylim)) {
    ylim <- range(beta1, finite=TRUE)
  }

  if(missing(xlim)) {
    xlim <- range(beta0, finite=TRUE)
  }
  

  for(i in seq_len(dim(reject)[3])) {
    image(x=beta0.fin, y=beta1.fin, z=reject.fin[,,i], xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, breaks=c(-1.5,-0.5,0.5,1.5), axes=FALSE,
          col=col,
          sub=bquote(.(as.symbol(names(dimnames(ACE))[3])) == .(format(as.numeric(dimnames(ACE)[[c(3,i)]]), digits=3))),
          ...)
  }
  
#  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
#  points(x=inf.x, y=ACE.inf, pch=1, col=point.col)
}

plot.sensitivity.1.1d <- function(x, xlim, ylim, xlab=expression(beta), ylab="SCE",
                                 t.point,
                                 display = c("analytic", "bootstrap"),
                                 col='black', line.col=col, point.col=col,
                                 analytic.col="red",
                                 analytic.line.col=analytic.col,
                                 analytic.point.col=analytic.col,
                                 bootstrap.col="green",
                                 bootstrap.line.col=bootstrap.col,
                                 bootstrap.point.col=bootstrap.col,
                                 panel.last=NULL,
                                 type='l', ...) {
  
  display <- match.arg(display, several.ok=TRUE)

  sortIndx <- sort.list(x$beta)
  beta <- x$beta[sortIndx]
  SCE <- x$SCE[sortIndx]
  SCE.ci <- x$SCE.ci[sortIndx, t.point,,, drop=FALSE]

  finIndx <- is.finite(beta)
  infIndx <- is.infinite(beta)
  
  beta.fin <- beta[finIndx]
  SCE.fin <- SCE[finIndx]
  SCE.ci.fin <- SCE.ci[finIndx, t.point,,, drop=FALSE]

  beta.inf <- beta[infIndx]
  SCE.inf <- SCE[infIndx]
  SCE.ci.inf <- SCE.ci[infIndx, t.point,,, drop=FALSE]
  
  if(missing(ylim)) {
    ylim <- range(c(SCE, SCE.ci))
  }

  if(missing(xlim)) {
    xlim <- range(beta, finite=TRUE)
  }
  
  indx <- match(beta.inf, c(-Inf, Inf))

  plot.default(x=beta.fin, y=SCE.fin, xlim=xlim, ylim=ylim, type=type,
               xlab=xlab, ylab=ylab, ...)
  inf.x <- par('usr')[indx] + strwidth("m")/c(1,-1)[indx]
  points(x=inf.x, y=SCE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$SCE.var)) {
    for(i in seq_len(dim(SCE.ci.fin)[[3]])) {
      lines(x=beta.fin, y=SCE.ci.fin[, t.point, i,'analytic'], lty=2,
            col=analytic.line.col)
      points(x=inf.x, y=SCE.ci.inf[,t.point, i, "analytic"], pch=3, col=analytic.point.col)
    }
  }

  if('bootstrap' %in% colnames(x$SCE.var) && 'bootstrap' %in% display) {
    for(i in seq_len(dim(SCE.ci.fin)[[3]])) {
      lines(x=beta.fin, y=SCE.ci.fin[,t.point,i,'bootstrap'], lty=2,
            col=analytic.line.col)
    }
    points(x=rep(inf.x, times=2), y=SCE.ci.inf[indx, t.point,,"bootstrap"], pch=3, col=bootstrap.point.col)
  }
}
