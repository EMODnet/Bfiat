\name{critical values for fishing}
\alias{critical_sar}
\alias{critical_gpd}
\alias{critical_d}
\alias{critical_r}
\title{
 Critical values for fishing or taxon parameters. 
}

\description{
  \code{critical_sar} and \code{critical_gpd} estimates for one or more taxa the fishing intensity (\code{sar}) or gear penetration depth (\code{gpd}), above which the density/caryying capacity ratio decreases below a critical value (\code{crit_DK}). 
  
  \code{critical_r} and \code{critical_d} estimate for a combination of fishing intensity (\code{sar}) and gear penetration depth (\code{gpd}), the critical intrinsic rate of natural increase (\code{r}) or depletion fraction (d) at which the density, carrying capacity ratio reaches the critical value (\code{crit_DK}). 
}

\usage{
critical_sar (gpd = 1, fDepth = 1, uDepth = 0, d = NULL,
              r = 0.1, crit_DK = 0, ...) 
  
critical_gpd ( sar = 1, fDepth = 1, uDepth = 0, d = NULL,
               r = 0.1, crit_DK = 0, ...) 

critical_r (sar = 1, gpd = 1, fDepth = 1, uDepth = 0, d = NULL, 
            crit_DK = 0, ...) 

critical_d (sar = 1, r = 0.1, crit_DK = 0, ...) 
}

\arguments{
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year]. 
  One number or a vector.
  }
  \item{gpd }{gear penetration depth, units e.g. [cm].   
  One number or a vector.
  }
  \item{r }{the intrinsic rate of natural increase of a taxon, units e.g. [/year]. 
  One number or a vector.
  }
  \item{d }{depletion fraction due to fishing, one value or a vector. 
  }
  \item{crit_DK }{the steady-state density over carrying capacity ratio of the taxon that need to be matched. The default value of 0 indicates that the taxon is extinct.  A value of 0.01 indicates that the actual density is 1/100 of its potential density (carrying capacity). 
  }
  \item{fDepth }{fractional occurrence of species in sediment layers, dimensionless. 
  A vector of the same length as \code{uDepth}. The sum of \code{fDepth} should equal 1. Will be used to estimate the depletion fraction \code{d}. Will be ignored if \code{d} is given a value.
  }
  \item{uDepth }{depth of the upper position of the sediment layers, units e.g. [cm]. 
  A vector with length equal to the number of columns of \code{fDepth}. Will be ignored if \code{d} is given a value.
  }
  \item{... }{arguments passed to the \code{par_d} function
  }
}

\value{

These functions will either return one number, a vector or a matrix, depending on the input arguments
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{manip}

\examples{

## -----------------------------------------------------------------------
## Critical sar for beam trawling for all species from the MWTL dataset
## -----------------------------------------------------------------------

gpd_mud   <- 3.2 # beam trawl in muddy sediment
gpd_sand  <- 1.9 #               sandy

head(MWTL$fishing)
fDepth <- MWTL$fishing[,c("p0", "p0_5cm", "p5_15cm", "p15_30cm", "p30cm")]
uDepth <-               c( 0,    0,         5,        15,          30)

csar_sand <- critical_sar(gpd    = gpd_sand, 
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)

csar_mud  <- critical_sar(gpd    = gpd_mud, 
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)

csar <- data.frame(taxon     = MWTL$fishing$taxon, 
                   csar_sand = csar_sand, 
                   csar_mud  = csar_mud)
                   
summary(csar)


## -----------------------------------------------------------------------
## Critical gear penetration depth for species from the MWTL dataset
## -----------------------------------------------------------------------

fDepth <- MWTL$fishing[,c("p0", "p0_5cm", "p5_15cm", "p15_30cm", "p30cm")]
uDepth <-               c( 0,    0,         5,        15,          30)

cgpd_01   <- critical_gpd(sar   = 1,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
cgpd_10   <- critical_gpd(sar   = 10,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
cgpd_100  <- critical_gpd(sar   = 100,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
# sar of 100 and critical density = 10 percent
cgpd_100b <- critical_gpd(sar     = 100,
                          fDepth  = fDepth, uDepth = uDepth, 
                          r       = MWTL$fishing$r,
                          crit_DK = 0.1)

cgpd <- data.frame(taxon     = MWTL$fishing$taxon, 
                   cgpd_01   = cgpd_01, 
                   cgpd_10   = cgpd_10,
                   cgpd_100  = cgpd_100,
                   cgpd_100b = cgpd_100b)

# a cgpd of 100 means: no limits                   
summary(cgpd)


}

\seealso{
\link{run_perturb} for how to run a disturbance model.

\link[Btrait]{Traits_nioz}, for trait databases in package Btrait.

\link[Btrait]{MWTL} for data sets on which fishing can be imposed.

\link[Btrait]{map_key} for simple plotting functions.
}

\references{
to be added
}