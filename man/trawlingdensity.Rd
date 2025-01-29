\name{density during trawling events}
\alias{trawl_perturb}
\alias{trawl_metier}

\title{
  Estimates density before and after a trawling event, or the average density between trawls.
}

\description{
 \code{trawl_perturb} and \code{trawl_metier} calculates the density of taxa immediately at the start and end after the n-th trawling event, and the mean density inbetween trawling events, for one metier and multiple metiers respectively. 
 
}

\usage{
trawl_perturb   (K = 1, r = 1, d = 0.1, 
                parms = data.frame(K = K, r = r, d = d), 
                sar = 1, D0 = parms[["K"]], n = 1) 

trawl_metier    (K = 1, r = 1, 
                 parms = data.frame(K = K, r = r), 
                 d = data.frame(0.1, 0.1), 
                 sar = data.frame(1, 2), 
                 D0 = parms[["K"]], n = 1) 
}

\arguments{
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year]. For \code{trawl_perturb}: one number. For \code{trawl_metier}, a data.frame with one value per metier. 
  }
  \item{r }{the rate of increase of each taxon, units e.g. [/year]. One number, or a vector.
  }
  \item{d }{depletion fraction due to fishing. For \code{trawl_perturb}: one number, or a vector of same length as \code{r}.  For \code{trawl_metier}, a data.frame with nrows = length(r), and ncols = number of metiers, containing one number per metier and per taxon.
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
   \item \code{times_end}, the end time of the trawling events as in argument \code{n}, 
   \item \code{ntrawl_from}, the number of the trawl at the start of the trawling interval (\code{n}),
   \item \code{ntrawl_to}, the number of the trawl at the end of the trawling interval (\code{n+1}),
   \item \code{density}, the density averaged over the trawling interval,
   \item \code{density_start}, the density after trawling at the n-th trawling interval,
   \item \code{density_end}, the density at the end of the n-th trawling interval, before the next trawling.
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

D1 <- trawl_perturb(n = 1,  
                    K = 100)   
D1
             
## ====================================================
## Density for a sequence of fishing events
## ====================================================

D2 <- trawl_perturb(n = 1:10, 
                    K = 100) 
D2

## ====================================================
## Density for a sequence of fishing events 
## and multiple parameters
## ====================================================

D3 <- trawl_perturb(n = 1:5, 
                    K = c(10, 20))
D3

## ====================================================
## compare trawl_perturb with the dynamic model solution
## ====================================================

S      <- 1    # swept area ratio
tend   <- 20   

fishD  <- trawl_perturb(n   = 1:tend, 
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
           main = "comparison trawl_perturb and run_perturb")
           
with (fishD,
   points(times, fishD$density_start ,  
          col = factor(K), pch = 15))

with (fishD,
   points(times_end, density_end ,  
          col = factor(K), pch = 1))

with (fishD,
   points(x = 0.5*(times+times_end), y = density ,  
          col = factor(K), pch = 22))


legend("left", lty = c(1, NA, NA, NA), pch = c(NA, 15, 1, 22), 
      legend=c("run_perturb", "after perturbation", 
      "before perturbation", "interval average"))

## ====================================================
## compare trawl_metier with the dynamic model solution
## ====================================================

S      <- data.frame(1, 3)    # swept area ratio


d   = data.frame(rbind(c(0.2, 0.1),
                       c(0.4, 0.2)))

tend   <- 20   

fishD  <- trawl_metier (n   = 1:tend, 
                        d   = d,
                        K   = c(40, 100),
                        sar = S)

head(fishD)

# dynamic run

times  <- seq(from = 0, to = 7, length.out = 100)

Dens <- density_metier(d   = d,
                       K   = c(40, 100),
                       sar = S,
                       times = times,
                       as.deSolve = TRUE)  
# compare

matplot.0D(Dens, 
           lty = 1, legend = FALSE, 
           ylab = "density", 
           main = "comparison trawl_metier and density_metier")
           
with (fishD,
   points(times, fishD$density_start ,  
          col = factor(K), pch = 15))

with (fishD,
   points(times_end, density_end ,  
          col = factor(K), pch = 1))

with (fishD,
   points(x = 0.5*(times+times_end), y = density ,  
          col = factor(K), pch = 22))


legend("left", lty = c(1, NA, NA, NA), pch = c(NA, 15, 1, 22), 
      legend=c("density_metier", "after perturbation", 
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
