\name{density }
\alias{density_perturb}
\alias{density_logistic}

\title{
  Estimates density at specific times for the logistic and perturbation model.
}

\description{
 \code{density_perturb} and \code{density_logistic} calculates the density at specific times for the perturbation and continuous model respectively. 

}

\usage{
density_perturb (K = 1, r = 1, d = 0.1, 
                parms = data.frame(K = K, r = r, d = d),
                sar = 1, D0 = parms[["K"]], times = 1, 
                tend_perturb = max(times), as.deSolve = FALSE)

density_logistic(K = 1, r = 1, m = 0.1, 
                parms = data.frame(K = K, r = r, m = m),
                D0 = parms[["K"]], times = 1,
                tend_perturb = max(times), as.deSolve = FALSE)
}

\arguments{
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year]. One number. Will be ignored if a column named "sar" is present in \code{parms}.
  }
  \item{r }{the rate of increase of each taxon, units e.g. [/year]. One number, or a vector.
  }
  \item{d }{depletion fraction due to fishing. One number, or a vector of same length as \code{r}.
  }
  \item{m }{mortality rate. One number, or a vector of same length as \code{r}.
  }
  \item{K }{carrying capacity. One number, or a vector of same length as \code{r}.
  }
  \item{parms }{a data.frame with all parameters. If it also includes a column labeled "times" or "sar", then the argument(s) with the same name will be ignored.
  }
  \item{D0 }{initial density. One number, or a vector of length = nrow(parms).
  }
  \item{times }{time at which to estimate the density. One number, or a vector. Ignored if \code{parms} contains a column labeled \code{"times"}.
  }
    \item{tend_perturb }{time at which the perturbation stops - from then on, no events will occur (\code{density_perturb}) or the mortality parameter will be set = 0 (\code{density_logistic}.
  }
  \item{as.deSolve }{if \code{TRUE}, will return an object of class deSolve, where the first column is time}
}

\value{
 
 If \code{as.deSolve} is \code{FALSE}, and \code{times} is part of the argument \code{parms}, 
 returns a \code{data.frame}, with the following columns:
 
 \itemize{
 
  \item for \code{density_logistic}:
  
  \itemize{
   \item \code{K, r, m}: the model parameters, 
   \item \code{D0}: the initial condition, 
   \item \code{times}, the times as in the argument \code{times}, 
   \item \code{density}, the requested density at \code{times}.
  }
 
 \item for \code{density_perturb}:
  
  \itemize{
   \item \code{K, r, d}: the model parameters, 
   \item \code{sar}: the swept area ratio, 
   \item \code{D0}: the initial condition, 
   \item \code{times}, the times as in the argument \code{times}, 
   \item \code{density}, the requested density at \code{times}.
  }
 }


 If \code{as.deSolve} is \code{FALSE}, and \code{times} is a vector: 
 returns a \code{data.frame}, with the following columns:
 
 \itemize{
 
  \item for \code{density_logistic}:
  
  \itemize{
   \item \code{K, r, m}: the model parameters, 
   \item \code{D0}: the initial condition, 
   \item \code{Dt_..}, columns with the density for each value in the argument \code{times} 
  }
 
 \item for \code{density_perturb}:
  
  \itemize{
   \item \code{K, r, d}: the model parameters, 
   \item \code{sar}: the swept area ratio, 
   \item \code{D0}: the initial condition, 
   \item \code{Dt_..}, columns with the density for each value in the argument \code{times} 
  }
 }

 If \code{as.deSolve} is \code{TRUE}, returns an object of class \code{deSolve}:
 
 \itemize{
 
  \itemize{
   \item \code{time}: the time as in the argument \code{times} 
   \item \code{..}, columns with the density for each row in the argument \code{parms}
  }
 
 }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{manip}

\examples{

## ====================================================
## Compare density_perturb with dynamic run_perturb
## ====================================================

times <- c(0.5, 1.5, 2.5, 8.5)

T1 <- density_perturb (K     = 100, 
                       times = times)      
T1

# dynamic run

rP <- run_perturb(parms = list(K = 100, r = 1, d = 0.1),
                  sar   = 1,
                  times = seq(0, 11, by = 0.01))

# compare

plot(rP, ylim = c(80, 100), las = 1, 
     main = "comparison density_perturb and run_perturb")
     
points(x   = T1$times, y = T1$density, 
       pch = 16, cex = 2)

legend("topright", 
        legend = c("run_perturb", 
                   "density_perturb"),
        lty = c(1, NA), pch = c(NA, 16))       
       
## ====================================================
## Compare density_logistic with dynamic run_logistic
## ====================================================

times <- c(0.5, 1.5, 2.5, 8.5)

T2 <- density_logistic(K     = 100, 
                       times = times)
T2

# dynamic run

rL <- run_logistic(parms = list(K = 100, r = 1, m = 0.1),
                   times = seq(0, 11, by = 0.01))

# compare

plot(rL, ylim = c(80, 100), las = 1, 
     main = "comparison density_logistic and run_logistic")
     
points(x   = T2$times, y = T2$density,   
       pch = 18, cex = 2)

legend("topright", 
        legend = c("run_logistic", 
                   "density_logistic"),
        lty = c(1, NA), pch = c(NA, 18))       

}

\seealso{
\link{density_trawl} for density before, after trawling and mean density inbetween trawling.

\link{steady_perturb} for steady-state calculations.

\link{run_perturb} for how to run a disturbance model.

\link{run_perturb} for how to run a disturbance model.

\link[Btrait]{Traits_nioz}, for trait databases in package Btrait.

\link[Btrait]{MWTL} for data sets on which fishing can be imposed.

\link[Btrait]{map_key} for simple plotting functions.
}

\references{
Hiddink, JG, Jennings, S, Sciberras, M, et al., 2019. Assessing bottom trawling impacts based on the longevity of benthic invertebrates. J Appl Ecol 56: 1075-1084. https://doi.org/10.1111/1365-2664.13278 

Hiddink, Jan Geert, Simon Jennings, Marija Sciberras, Claire L. Szostek, Kathryn M. Hughes, Nick Ellis, Adriaan D. Rijnsdorp, Robert A. McConnaughey, Tessa Mazor, Ray Hilborn, Jeremy S. Collie, C. Roland Pitcher, Ricardo O. Amoroso, Ana M. Parma, Petri Suuronen, and Michel J. Kaiser, 2017.  Global analysis of depletion and recovery of seabed biota after bottom trawling disturbance. Proc. Nat. Aca. Sci, 114 (31) 8301-8306. https://doi.org/10.1073/pnas.161885811.

C.R. Pitcher, N. Ellis, S. Jennings, J.G. Hiddink, T. Mazor, M.J.Kaiser, M.I. Kangas, R.A. McConnaughey, A.M. Parma, A.D. Rijnsdorp, P. Suuronen, J.S. Collie, R. Amoroso, K.M. Hughes and R. Hilborn, 2017. Estimating the sustainability of towed fishing-gearimpacts on seabed habitats: a simple quantitative riskassessment method applicable to data-limited fisheries. Methods in Ecology and Evolution,8,472-480doi: 10.1111/2041-210X.12705

}
