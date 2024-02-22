\name{density }
\alias{eventDensity}
\alias{steadyDensity}
\alias{intervalDensity}
\title{
  Estimate steady-state density, density after trawling, average density between trawls
}

\description{
 \code{steadyDensity} calculates the steady-state density of benthic taxa, based on their carrying capacity K, on the fishing intensity (swept area ratio, sar), and on the species-specific parameters r and d (rate of increase and depletion fraction).

 \code{eventDensity} calculates the density of benthic taxa immediately after the n-th trawling event. 
 
 \code{intervalDensity} calculates the mean density of benthic taxa inbetween the n-th and n+1 th trawling event. 

}

\usage{
steadyDensity(K=1, sar=1, r=1, d=0.1)

eventDensity(K=1, sar=1, r=1, d=0.1, p=1-d, D0=K, n=1) 

intervalDensity(K=1, sar=1, r=1, d=0.1, p=1-d, D0=K, n=1)
}

\arguments{
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year]. One number, or, (for steadyDensity) a vector of the same length as \code{density}.
  }
  \item{r }{the rate of increase of each taxon, units e.g. [/year]. One number, or (for steadyDensity) a vector.
  }
  \item{d }{depletion fraction due to fishing. One number, or (for steadyDensity) a vector of same length as \code{r}.
  }
  \item{p }{fraction remaining after fishing. One number.
  }
  \item{K }{carrying capacity. One number, or (for steadyDensity) a vector of same length as \code{r}.
  }
  \item{D0 }{initial density. One number.
  }
  \item{n }{number of fishing events. One number.
  }
}

\value{
 returns the estimated densities, one number or a vector.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{manip}

\examples{

## ====================================================
## Density after n fishing events
## ====================================================

eventDensity(n=1,  K=100)     
eventDensity(n=10, K=100) 

eventDensity(n=10, d=0.2, K=80)      # density immediately after 10th trawling
eventDensity(n=11, d=0.2, K=80)/0.8  # density immediately before 11th trawling

intervalDensity(n=10, d=0.2, K=80)   # mean density in the interval

# apply the function for a sequence of events
tend   <- 20
S      <- 1
events <- seq(1/S, tend, by=1/S)

# Density immediately after all events
fishD  <- mapply(events, FUN=function(x) 
                         eventDensity(n=x, d=0.1, K=100, sar=S))
fishD2 <- mapply(events, FUN=function(x) 
                         eventDensity(n=x, d=0.2, K= 80, sar=S))

# Mean density inbetween events
fMean  <- mapply(events, FUN=function(x) 
                         intervalDensity(n=x, d=0.1, K=100, sar=S))
fMean2 <- mapply(events, FUN=function(x) 
                         intervalDensity(n=x, d=0.2, K= 80, sar=S))

# compare it with the dynamic model solution
times  <- seq(0, tend,length.out=1000)

Dens <- perturb(parms=list(r=c(1,1), K=c(100, 80), d=c(0.1, 0.2)), 
                times=times, events=events)

matplot.0D(Dens, lty=1, legend=FALSE)
points(events, fishD,  col=1, pch=18)
points(events, fishD2, col=2, pch=18)
lines(events+0.5/S, fMean, lty=2)
lines(events+0.5/S, fMean2, col=2, lty=2)
legend("right", col=1, lty=c(1,2,NA), pch=c(NA,NA,18), 
      legend=c("continuous", "interval average", "after perturbation"))

## ====================================================
## Steady-state
## ====================================================

# discrete parameter values and run

parms  <- c(r=0.5, d=0.2, K=1)
sar    <- 0.5

times  <- seq(from=0,      to=20, length.out=1000)
events <- seq(from=10/365, to=20, by=1/sar)

run.discrete   <- perturb(parms  = parms, 
                          times  = times, 
                          events = events)

# continuous simulation 

# Estimate continous mortality equivalent with event fishing
m              <- par_m(sar = sar, 
                        r   = parms["r"], 
                        d   = parms["d"])

# new parameter vector                        
parms.c        <- c(parms[c("r", "K")], m = m)

run.continuous <- logistic(parms = parms.c, 
                           times = times)

# show both simulations in one figure

plot(run.discrete, run.continuous)

# use steadyDensity to estimate the steady-state condition

std <- steadyDensity(r   = 0.5, 
                     d   = 0.2, 
                     K   = 1, 
                     sar = sar)
abline(h=std, lty=2, col="green")

}

\seealso{
\link{perturb} for how to run a disturbance model.

\link[Btrait]{Traits_nioz}, for trait databases in package Btrait.

\link[Btrait]{MWTL} for data sets on which fishing can be imposed.

\link[Btrait]{mapKey} for simple plotting functions.
}

\references{
Hiddink, JG, Jennings, S, Sciberras, M, et al., 2019. Assessing bottom trawling impacts based on the longevity of benthic invertebrates. J Appl Ecol 56: 1075-1084. https://doi.org/10.1111/1365-2664.13278 

Hiddink, Jan Geert, Simon Jennings, Marija Sciberras, Claire L. Szostek, Kathryn M. Hughes, Nick Ellis, Adriaan D. Rijnsdorp, Robert A. McConnaughey, Tessa Mazor, Ray Hilborn, Jeremy S. Collie, C. Roland Pitcher, Ricardo O. Amoroso, Ana M. Parma, Petri Suuronen, and Michel J. Kaiser, 2017.  Global analysis of depletion and recovery of seabed biota after bottom trawling disturbance. Proc. Nat. Aca. Sci, 114 (31) 8301-8306. https://doi.org/10.1073/pnas.161885811.

C.R. Pitcher, N. Ellis, S. Jennings, J.G. Hiddink, T. Mazor, M.J.Kaiser, M.I. Kangas, R.A. McConnaughey, A.M. Parma, A.D. Rijnsdorp, P. Suuronen, J.S. Collie, R. Amoroso, K.M. Hughes and R. Hilborn, 2017. Estimating the sustainability of towed fishing-gearimpacts on seabed habitats: a simple quantitative riskassessment method applicable to data-limited fisheries. Methods in Ecology and Evolution,8,472-480doi: 10.1111/2041-210X.12705

}
