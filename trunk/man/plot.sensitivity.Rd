\name{plot.sensitivity}
\alias{plot.sensitivity}
\alias{plot.sensitivity.1.0d}
\alias{plot.sensitivity.2.0d}

\title{
  plots the results of calls to the sensitivity functions.
}
\description{
  Functions used to plot the objects created by the sensitivityPStrat
  family of functions.
}
\usage{
\method{plot}{sensitivity.1.0d}(x, xlim, ylim,
     xlab = expression(beta), ylab ="ACE", display = c("analytic", "bootstrap"),
     col = "black", line.col = col, point.col = col,
     analytic.col = "red", analytic.line.col = analytic.col,
     analytic.point.col = analytic.col,
     bootstrap.col = "green", bootstrap.line.col = bootstrap.col,
     bootstrap.point.col = bootstrap.col,
     panel.last = NULL, type = "l", ...)

\method{plot}{sensitivity.2.0d}(x, xlim, ylim, xlab = expression(beta[0]), ylab = expression(beta[1]),
     display = c("analytic", "bootstrap"), col = c(gray(.9), gray(1), gray(.8)),
     panel.last = NULL, ...)

}
\arguments{
  \item{x}{
    sensitivity object
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{display}{
%%     ~~Describe \code{display} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{line.col}{
%%     ~~Describe \code{line.col} here~~
}
  \item{point.col}{
%%     ~~Describe \code{point.col} here~~
}
  \item{analytic.col}{
%%     ~~Describe \code{analytic.col} here~~
}
  \item{analytic.line.col}{
%%     ~~Describe \code{analytic.line.col} here~~
}
  \item{analytic.point.col}{
%%     ~~Describe \code{analytic.point.col} here~~
}
  \item{bootstrap.col}{
%%     ~~Describe \code{bootstrap.col} here~~
}
  \item{bootstrap.line.col}{
%%     ~~Describe \code{bootstrap.line.col} here~~
}
  \item{bootstrap.point.col}{
%%     ~~Describe \code{bootstrap.point.col} here~~
}
  \item{panel.last}{
%%     ~~Describe \code{panel.last} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\author{
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}

\seealso{
  \code{\link{plot.default}}  
}
\examples{
data(vaccine.trial)

plot(with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    beta1=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    phi=c(0.95,0.90,0.80), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=100)
         ))

}
\keyword{ hplot }

