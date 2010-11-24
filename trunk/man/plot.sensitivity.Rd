\name{plot.sensitivity}
\alias{plot.sensitivity.1.0d}
\alias{plot.sensitivity.2.0d}
\alias{plot.sensitivity.1d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plots the results of calls to the sensitivity functions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.sensitivity.1.0d(x, xlim, ylim, xlab = expression(beta), ylab =
"ACE", display = c("analytic", "bootstrap"), col = "black", line.col =
col, point.col = col, analytic.col = "red", analytic.line.col =
analytic.col, analytic.point.col = analytic.col, bootstrap.col =
"green", bootstrap.line.col = bootstrap.col, bootstrap.point.col =
bootstrap.col, panel.last = NULL, type = "l", ...)
plot.sensitivity.2.0d(x, xlim, ylim,
                      xlab=expression(beta[0]),
                      ylab=expression(beta[1]),
                      display = c("analytic", "bootstrap"),
                      col=c(gray(.9),gray(1),gray(.8)),
                      panel.last=NULL, ...) {

}
%- maybe also 'usage' for other objects documented here.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
