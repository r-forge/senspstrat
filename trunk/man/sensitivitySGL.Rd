\name{sensitivitySGL}
\alias{sensitivitySGL}
\title{principal stratification sensitivity analysis with time to event data}
\description{
  principal stratification sensitivity analysis with time to event data
  using the method described by Shepherd, Gilbert, and Lumley (2007).
}
\usage{
sensitivitySGL(z, s, d, y, beta, tau, time.points, selection, trigger,
               groupings, empty.principle.stratum, ci=0.95,
               ci.method = c("analytic", "bootstrap"), na.rm = FALSE,
               N.boot = 100L, oneSidedTest = FALSE, twoSidedTest = TRUE,
               verbose = getOption("verbose"))
}
\arguments{
  \item{z}{
    vector; contains the grouping values (e.g., treatment assignment)
    for each record.
  }
  \item{s}{
    vector; indicates whether a record is selected.
  }
  \item{d}{
    vector; indicates whether a post-selection event has occurred. Can
    be \code{NA} for unselected records.
  }
  \item{y}{
    vector; the length of time from selection until event (d) or
    censoring.
    Can be \code{NA} for unselected records.
  }
  \item{beta}{
    vector; values of the sensitivity parameter.  \code{Inf} and
    \code{-Inf} are acceptable.
  }
  \item{tau}{
    maximum observed follow-up time after selection.  Selection weights
    are constant for t>tau.
  }
  \item{time.points}{
    vector; time points, t, at which SCE(t) will be estimated.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{trigger}{
    logical; the value of \code{d} that denotes the post-selection event.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, the first element \code{g0} being the
    value of \code{z} the delineates the first group, the last element
    \code{g1} being the value of \code{z} which delineates the second group.
  }
  \item{empty.principal.stratum}{
    vector of two elements \code{c(s0,s1)}; describes the \code{s}
    values that select the empty principal stratum.  If
    \code{empty.principal.stratum=c(s0,s1)}, then stratum defined by
    S(g0)==s0 and S(g1)==s1 is the empty stratum.  In this example s0
    and s1 refer to the two possible values of s. (Note: method only
    works if s0 != s1).
  }
  \item{ci}{
    numeric vector; confidence interval level, defaults to 0.95.
  }
  \item{ci.method}{
    character;  method by which the convidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. 
  }
  \item{na.rm}{
    logical; indicates whether records that are invalid due to \code{NA}
    values should be removed from the data set.
  }
  \item{N.boot}{
    integer; number of bootstrap repitions that will be run when running with
    \code{ci.method} equal to \dQuote{bootstrap}.
  }
  \item{oneSidedTest}{
    logical.  Return a one sided confidence interval for
    SCE(t). Defaults to \code{FALSE}.
  }
  \item{twoSidedTest}{
    logical.  Return a two sided confidence interval for SCE(t).
    Defaults to \code{TRUE}.
  }
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is happening.
  }
}
\details{
  
}
\value{
  object of class \code{sensitivity3d}

  \item{SCE}{
    array; Calculated values of SCE for all combinations of the values from
    \code{phi}/\code{Pi}/\code{psi}, and
    \code{time.points}.  Array dimentions are
    \code{length(time.points)}, \code{length(psi)}.
  }
  \item{beta}{
    beta values used.
  }
  \item{psi}{
    phi values used.
  }
  \item{SCE.ci}{
    array; Confindence interval of the SCE value.  Confidence interval
    determined by \code{quantile} if using  \code{ci.method}
    \dQuote{bootstrap}.  Other wise calculated using analytic
    varance. Array dimentions the same as element \code{SCE}.
  }
  \item{SCE.var}{
    array; Variance of the SCE value.  Array dimentions the same as
    element \code{SCE}.
  }
}
\references{
  Shepherd BE, Gilbert PB, Lumley T (2007), "Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization," Journal of the American Statistical Association
  102, 573-582.
}
\author{
  Bryan E. Shepherd \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr

  
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}

\seealso{
  \code{\link{sensitivitySGD}}, \code{\link{sensitivityJR}}
  \code{\link{Surv}}
}

\examples{

data(vaccine.trial)
sens.time<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(.25, 0,-.25,-.5), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          N.boot=1000)
               )
sens.time

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design }
