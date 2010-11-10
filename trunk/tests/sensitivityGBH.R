library(sensitivityPStrat)

data(vaccine.trial)

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(0,.25,.5,.75,1,1.25,1.5),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=1000)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf, 0,.25,.5,.75,1,Inf),                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=1000)
         )
ans
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=-Inf,                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=1000)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))
