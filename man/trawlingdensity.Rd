\name{density during trawling events}
\alias{density_trawl}

\title{
  Estimates density before and after a trawling event, or the average density between trawls.
}

\description{
 \code{density_trawl} calculates the density of taxa immediately at the start and end after the n-th trawling event, and the mean density inbetween trawling events. 
 
}

\usage{
density_trawl   (K = 1, r = 1, d = 0.1, 
                parms = data.frame(K = K, r = r, d = d), 
                sar = 1, D0 = parms[["K"]], n = 1) 
}

\arguments{
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year]. One number.
  }
  \item{r }{the rate of increase of each taxon, units e.g. [/year]. One number, or a vector.
  }
  \item{d }{depletion fraction due to fishing. One number, or a vector of same length as \code{r}.
  }
  \item{K }{carrying capacity. One number, or a vector of same length as \code{r}.
  }
  \item{parms }{a data.frame with all parameters. 
  }
  \item{D0 }{initial density. One number, or a vector of length = nrow(parms).
  }
  \item{n }{number of fishing events. One number, or a vector.
  }
}

\value{
 Returns a \code{data.frame}, with the following columns:
 
 \itemize{
   \item \code{K, r, d, D0}: the model parameters, 
   \item \code{times}, the start time of the trawling events as in argument \code{n}, 
   \item \code{end_times}, the end time of the trawling events as in argument \code{n}, 
   \item \code{ntrawl_from}, the number of the trawl at the start of the trawling interval (\code{n}),
   \item \code{ntrawl_to}, the number of the trawl at the end of the trawling interval (\code{n+1}),
   \item \code{average_density}, the density averaged over the trawling interval,
   \item \code{start_density}, the density after trawling at the n-th trawling interval,
   \item \code{end_density}, the density at the end of the n-th trawling interval, before the next trawling.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{manip}

\examples{

## ====================================================
## Density for one fishing events
## ====================================================

D1 <- density_trawl(n = 1,  
                    K = 100)   
D1
             
## ====================================================
## Density for a sequence of fishing events
## ====================================================

D2 <- density_trawl(n = 1:10, 
                    K = 100) 
D2

## ====================================================
## Density for a sequence of fishing events 
## and multiple parameters
## ====================================================

D3 <- density_trawl(n = 1:5, 
                    K = c(10, 20))
D3

## ====================================================
## compare density_trawl with the dynamic model solution
## ====================================================

S      <- 1    # swept area ratio
tend   <- 20   

fishD  <- density_trawl(n   = 1:tend, 
                        d   = c(0.2, 0.1),
                        K   = c(40, 100),
                        sar = S)

head(fishD)

# dynamic run

times  <- seq(from = 0, to = tend, length.out = 1000)

parms  <- list(r = c(1,     1), 
               K = c(40,  100), 
               d = c(0.2, 0.1) )

Dens <- run_perturb(parms = parms, 
                    times = times, 
                    sar   = S)  # can only have one number for sar

# compare

matplot.0D(Dens, 
           lty = 1, legend = FALSE, 
           ylab = "density", 
           main = "comparison density_trawl and run_perturb")
           
with (fishD,
   points(times, fishD$start_density ,  
          col = factor(K), pch = 15))

with (fishD,
   points(end_times, fishD$end_density ,  
          col = factor(K), pch = 1))

with (fishD,
   points(x = 0.5*(times+end_times), y = fishD$average_density ,  
          col = factor(K), pch = 22))


legend("left", lty = c(1, NA, NA, NA), pch = c(NA, 15, 1, 22), 
      legend=c("run_perturb", "after perturbation", 
      "before perturbation", "interval average"))

}

\seealso{
\link{steady_perturb} for steady-state calculations.

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
