\name{Impact of bottom disturbance}
\alias{run_perturb}
\alias{run_logistic}
\title{
  Functions to estimate the impact of bottom disturbances on benthic taxa. Allows for non-regularly spaced events.
}

\description{
 \code{run_perturb} applies a logistic model to estimate the depletion and recovery of benthic taxa from (a sequence of) bottom disturbing (e.g. trawling) events.  Fishing can be implemented by specifiying the perturbation events.

 \code{run_logistic} applies a logistic model with constant mortality.
}

\usage{
run_perturb (parms, times, sar = 1, 
             tstart_perturb = min(times) + 0.5/sar, tend_perturb = max(times),
             events = NULL, taxon_names = parms[["taxon"]], 
             D0 = parms[["K"]], addsum = FALSE, verbose = FALSE)
   
run_logistic(parms, times, 
             tstart_perturb = min(times), tend_perturb=max(times), 
             taxon_names = parms[["taxon"]], 
             D0 = parms[["K"]], addsum = FALSE, verbose = FALSE)
}

\arguments{
  \item{parms }{a vector, list, or data.frame with logistic parameters: for the function \code{run_perturb}, the parameters should be \code{r} (the rate of increase), \code{K} (carrying capacity), and \code{d}, the depletion factor during the event. For function \code{run_logistic}, it should contain either the parameters \code{r, K} and \code{m}, the constant mortality rate, or the parameters \code{K, r, sar}, and \code{d}, In the latter case, the parameter \code{m} will be estimated based on \code{sar, and d}. 
    }    
  \item{taxon_names }{a vector with names of the taxa whose characteristics are present in \code{parms}. Should be of the same length as (the number of rows of) \code{parms}. When not present or \code{NULL}, the taxa in the output will be labeled \code{"tax_1"}, \code{"tax_2"}, etc.
  }
  \item{sar }{fishing intensity, estimated as Swept Area Ratio, units e.g. [m2/m2/year], one number. Only used when \code{events} is \code{NULL}.
  }
  \item{times }{a vector with the output times of the model (sorted). The units of times should be consistent with the unit of parameter \code{r}; i.e. if \code{r} is expressed per year (1/y), then the time unit will be years.
  }
  \item{tstart_perturb, tend_perturb }{time at which the perturbation starts and stops - before                     tstart_perturb and after tend_perturb, no events will occur (\code{run_perturb}) or the mortality parameter will be set = 0 (\code{run_logistic}.
  }
  \item{events }{If not \code{NULL}, a vector with the times at which the disturbing event occurs. Event times that do not fall within \code{times} will be ignored. If \code{NULL}, then events will be estimated from \code{times}, \code{sar}, and \code{tend_perturb} as:  \code{seq(from = times[1], by = 1/sar, to = tend_perturb)}
  }
  \item{D0 }{initial conditions. If \code{D0} is not given, it will be set equal to the carrying capacity (= \code{parms["K"]}.
  }
  \item{addsum }{if \code{TRUE} an extra column with the summed densities will be added.
  }
  \item{verbose }{if \code{TRUE} will write warnings to the screen.
  }
}

\value{
  \code{run_perturb} and \code{run_logistic} return an object (matrix) of class \code{deSolve} (see \link[deSolve]{ode}), and contains, for each output time (rows), the density of all species (columns). If \code{addsum=TRUE}, the last column will be the summed density.
  
  \code{run_perturb} has the instantaneous depletion values in attribute \code{depletion}.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\details{
The function \code{run_perturb} implements the logistic function to describe density or biomass evolution of benthic taxa inbetween disturbance events.

The disturbance event itself creates an instantaneous density reduction for (some of) the taxa (as defined in parameter d).

Inbetween disturbances, the differential equation describing the change of \code{C} in time (\code{t}) is:

\ifelse{html}{\out{dC/dt = r * C * (1-C/K)}}{\eqn{\frac{dC}{dt}=r \cdot C \cdot (1-\frac{C}{K})}}

where \code{r} is the intrinsic rate of increase (units [/time]), K is the carrying capacity (same units as C).

This differential equation has as analytical solution:

\ifelse{html}{\out{C(t)=(C0*K) / (C0 + (K-C0) * exp(-r* t))}}{\eqn{C(t)=\frac{(C_0 \cdot K)}{C_0 + (K-C_0) \cdot exp^{(-r \cdot t)}}}}

where \code{C0} is the initial condition, and \code{t} is the time. 
 
During the disturbance at time T, the density or biomass C is instantaneously reduced with a factor \code{d}:

\ifelse{html}{\out{C{T+} = C{T-} * (1-d)}}{\eqn{C|_{T^+}=C|_{T^-} \cdot (1-d)}}

where (T-) and (T+) are the time immediately before (-) and after (+) the disturbance. 
The reduced density is then used as initial condition for the dynamics in the next section of time until the next disturbance event.
}
\keyword{manip}

\examples{

## ====================================================
## A small model - units in days
## ====================================================

# model parameters
 Trawl.pars <- data.frame(
   taxon = c("sp1", "sp2", "sp3", "sp4"),  # name of taxa
   r     = c(0.01,  0.005,  0.05, 0.001),  # rate of increase (/day)
   K     = c(   1,     10,   0.5,    10),  # carrying capacity (density or biomass)
   d     = c( 0.1,    0.1,   0.5,   0.8))  # instantaneous depletion fraction

# run the model for 5 years 
 times   <- 0:1825                          # time in days (consistent with pars)

# trawling frequency once or twice per year
 trawl1  <- seq(from=92, by=  365, to=365*5)
 trawl2  <- seq(from=92, by=182.5, to=365*5)

 T1 <- run_perturb (parms  = Trawl.pars, 
                    times  = times, 
                    events = trawl1)      
 T2 <- run_perturb (parms  = Trawl.pars, 
                    times  = times, 
                    events = trawl2)      
               
# show results
 plot(T1, T2, lty=1)

 head(T1)

# Depletion for each species at each trawling event
 attributes(T1)$depletion

## ----------------------------------------------------
## Randomly placed trawling events 
## ----------------------------------------------------

 trawl1r <- runif(n= 5, min=0, max=365*5)
 trawl2r <- runif(n=10, min=0, max=365*5)

 T1r <- run_perturb (parms  = Trawl.pars, 
                     times  = times, 
                     events = trawl1r)      
 T2r <- run_perturb (parms  = Trawl.pars, 
                     times  = times, 
                     events = trawl2r)      
 plot(T1r, T2r, lty=1)
 plot(T1 , T1r, lty=1)

 attributes(T1)$trawl

## ----------------------------------------------------
## Merging the small dataset with traits
## ----------------------------------------------------

# Note: no trait data for "sp.4"
# fuzzy-coded traits

 Btraits <- data.frame(
   taxon   = c("sp1","sp2","sp3","sp5","sp6"),
   T1_M1   = c(0    , 0   ,   0 , 0.2 ,    1),
   T1_M2   = c(1    , 0   , 0.5 , 0.3 ,    0),
   T1_M3   = c(0    , 1   , 0.5 , 0.5 ,    0),
   T2_M1   = c(0    , 0   ,   1 , 0.5 ,    1),
   T2_M2   = c(1    , 1   ,   0 , 0.5 ,    0)
)

# trait metadata (meaning of column names T1_M1 etc...)
 Btraits.lab <- data.frame(
   colname  = c("T1_M1","T1_M2","T1_M3","T2_M1","T2_M2"),
   trait    = c("T1"   ,"T1"   ,"T1"   ,"T2"   ,"T2"),
   modality = c("M1"   ,"M2"   ,"M3"   ,"M1"   ,"M2"), 
   score    = c(0      , 0.5   , 1     , 0.2   , 2)
)

# merge run with fuzzy traits, keep fuzzy code
 Trun <- get_trait_model(model = T1, 
                         trait = Btraits)
 plot(Trun)

# expand the run with crisp trait information
# (by passing trait_class and trait_score)
 Trun.crisp <- get_trait_model(model       = T1, 
                               trait       = Btraits, 
                               trait_class = Btraits.lab$trait, 
                               trait_score = Btraits.lab$score, 
                               scalewithvalue = FALSE,  # to have total numbers
                               verbose     = FALSE
                            )
 plot(Trun.crisp)

# categorical traits need to be converted to fuzyy (wide) format first

 Bcategory <- data.frame(
   taxon   = c("sp1","sp2","sp3","sp4"),
   C1      = c(  "A",  "B",  "A",  "C")
 )
 
 Bcat.fuzzy <- crisp2fuzzy(Bcategory)

 Trun.cat <- get_trait_model(model = T1, 
                             trait = Bcat.fuzzy)

 plot(Trun.cat)

## ====================================================
## Continuous mortality
## ====================================================

 Cont.pars <- data.frame(
   taxon = c("sp1",  "sp2",  "sp3",  "sp4"),  # name of taxa
   r     = c(0.01,   0.005,   0.05,  0.001),  # rate of increase (/day)
   K     = c(   1,      10,    0.5,     10),  # carrying capacity (density or biomass)
   m     = c( 0.005, 0.001,  0.025, 0.0004))  # daily mortality rate (/day)

# run the model for 5 years 
 times   <- 0:1825                          # time in days (consistent with pars)

# One run with fishing all times
 T1 <- run_logistic (parms  = Cont.pars, 
                     times  = times)

# Run where fishing tops after 3 years
 T2 <- run_logistic(parms        = Cont.pars,
                    tend_perturb = 3*365,
                    times        = times)
               
# show results
 plot(T1, T2, lty=1)

 tail(T1)

## ====================================================
## Trawling in Dutch part of the Northsea - 
## note: units are in years here!
## ====================================================

# ----------------------------------------------------
# 1. Select the data for one station; 
# ----------------------------------------------------
# use density data, averaged over all years.

 spB02 <- get_density( 
               data        = MWTL$density,
               subset      = (station == "BREEVTN02"), 
               averageOver = year, 
               descriptor  = station, 
               taxon       = taxon,
               value       = density)

# ----------------------------------------------------
# 2. find corresponding parameter values
# ----------------------------------------------------

# Select the traits that we will use to generate parameter values
# Generate unknown traits based on taxonomic closeness

  spTrait <- get_trait(taxon    = spB02$taxon, 
                       trait    = Traits_nioz, 
                       taxonomy = Taxonomy)

  Traits.ext <- na.omit(spTrait)   # remove the unknown traits

  Traits.lab <- attributes(Traits_nioz)$description

# The proportion of species occurrence shallower than 5 cm
# is used to estimate the depletion factor "d"

  head (Traits.lab)  # "ET1.M1" and ET1.M2: < 5cm DEEP
  prop5cm <- rowSums(Traits.ext[, c("ET1.M1", "ET1.M2")])

# convert traits to crisp values, to estimate mean lifespan
  sp.Trait.crisp   <- fuzzy2crisp(trait       = Traits.ext, 
                                  trait_class = Traits.lab$trait, 
                                  trait_score = Traits.lab$value)
  Life.span <- sp.Trait.crisp$Life.span

# Parameters r and d:
  MWTLparms <- data.frame(taxon = Traits.ext$taxon, 
                          d     = 0.8 * prop5cm,
                          r     = 5.31/Life.span)  # /year

# merge with sp, to add the carrying capacity (taken as the mean density)

  MWTLparms <- merge(MWTLparms, spB02, 
                     by = "taxon") 
                     
  colnames(MWTLparms)[which(colnames(MWTLparms)=="density")] <- "K"

# ----------------------------------------------------
# 3. run the disturbance, trawling 1 and 2 x per year
# ----------------------------------------------------

  times   <- seq(0, 5, by=1/365)                          # time in years
  trawl1r <- runif(n=5,  min=0, max=5)
  trawl2r <- runif(n=10, min=0, max=5)


  T1r <- run_perturb (parms  = MWTLparms, 
                      times  = times, 
                      events = trawl1r,
                      addsum = TRUE)      # to have "sum"   
  T2r <- run_perturb (parms  = MWTLparms, 
                      times  = times, 
                      events = trawl2r,
                      addsum = TRUE)    
                      
  plot(T1r, T2r, lty=1, which="Abra alba")
  plot(T1r, T2r, lty=1, which="sum")

# ----------------------------------------------------
# 4. The evolution of feeding types due to trawling
# ----------------------------------------------------
# from Traits.lab we see:
#   RT6.M1       Feeding type     Deposit feeder
#   RT6.M2       Feeding type     Suspension feeder
#   RT6.M3       Feeding type     Herbivore/Grazer
#   RT6.M4       Feeding type     Carnivore/Scavenger

# use this to select only the traits related to feeding type

FT <- Traits.ext[, c("taxon", "RT6.M1", "RT6.M2", "RT6.M3", "RT6.M4")]
colnames(FT)[-1] <- c("Deposit", "Suspension", "Herbivore", "Carnivore")

head(FT)

# merge perturbation output with the FT traits; use absolute numbers

FTrun <- get_trait_model(model          = T1r, 
                         trait          = FT, 
                         scalewithvalue = FALSE)  # keep absolute values
head(FTrun)

with(as.data.frame(FTrun), 
   matplot(x=times, y=cbind(Deposit, Suspension, Herbivore, Carnivore), 
           main="FEEDING TYPES station BREEVTN02", ylab="ind/m2", 
           xlab="time, years", type="l", lty=1, lwd=2))

legend("top", col=1:4, lty=1, lwd=2,
   legend = c("Deposit", "Suspension", "Herbivore", "carnivore"))   


}

\seealso{
\link{get_trait_model} for how to combine the output of a disturbance model with traits.

\link[Btrait]{Traits_nioz}, for trait databases in package Btrait..

\link[Btrait]{get_trait_density}, for the function on which \code{get_trait_model} is based.

\link[Btrait]{get_trait} for functions from package Btrait to extract traits

\link[Btrait]{MWTL} for data sets on which fishing can be imposed

\link[Btrait]{map_key} for simple plotting functions
}
