## ====================================================================
## ====================================================================
## Estimating trawling parameters
## ====================================================================
## ====================================================================


## ====================================================================
# a function to estimate depletion rate based on depth occurrence of species
## ====================================================================

par_d <- function(gpd   =    1 ,   # gear penetration depth, /cm
                  m_d   = 0.075,   # mortality parameter, /cm
                  m_max =   0.45,  # maximum depletion in any layer
                  fDepth,          # fractional occurrence of species in depth layers
                  uDepth           # upper position of depth layers, cm
                        ) {    
  if (!is.vector(uDepth)) 
    stop("'uDepth' should be a vector with the upper position of each depth layer")

  uL <- length(uDepth)

  if (is.vector(fDepth) & ! length(fDepth) == uL) 
    stop("'uDepth' and 'fDepth' should be of equal length if 'fDepth' is a vector")
  
  if (! is.numeric(fDepth))
    fDepth <- sapply(fDepth, as.numeric)

  if (is.vector(fDepth) ) 
    fDepth <- matrix(ncol=length(uDepth), data=fDepth)
  
  if (! is.data.frame(fDepth) & ! is.matrix(fDepth)) 
    stop("'fDepth' should be a data.frame or matrix with species x depth fraction")

  if (ncol(fDepth) != uL) 
    stop("number of columns of 'fDepth' should be = length of 'uDepth'")
    
  if (length(gpd) == 1 & nrow(fDepth) != 1){
    gpd <- rep(gpd, times=nrow(fDepth))
    
  } else if (nrow(fDepth) ==1 & length(gpd) > 1) 
    fDepth <- matrix(ncol=ncol(fDepth), nrow=length(gpd), data=fDepth, byrow=TRUE)
  
  if (length(gpd) != nrow(fDepth))
      stop ("'fDepth' and 'gpd' not compatible; length of gpd should be = 1 or = nrow(fDepth")

  if (is.null(m_max)) m_max <- 1
  if (length(gpd) == 1){
    cn <-  m_d*pmax(0, gpd - uDepth)
    if (length(cn) != ncol(fDepth)) 
      stop ("'fDepth' and 'uDepth' not compatible; length of uDepth should be = ncol(fDepth")
    AR <- sweep(fDepth, MARGIN=2, STATS=cn, FUN="*") 
  } else{
    cn <- outer(X=gpd, Y=uDepth, FUN = function(x,y) m_d*pmax(0, x - y))
    if (ncol(cn) != ncol(fDepth)) 
      stop ("'fDepth' and 'uDepth' not compatible; length of uDepth should be = ncol(fDepth")
    AR <- fDepth*cn
  }
  
  AR <- rowSums(AR)
  pmin(m_max, AR)
}

## ====================================================================
# a function to estimate the rate of increase based on longevity
## ====================================================================

par_r <- function(longevity = 1  # lifetime, years
                        ) {    
   ri <- 5.32/longevity
   ri[is.infinite(ri)] <- NA
   ri
  }

## ====================================================================
# a function to estimate carrying capacity, based on species abundances,
# swept area ratios  (units: year) and model parameters (units: days)
## ====================================================================

par_K <- function(density,  # current densities
                  sar,      # Swept area ratio, units m2/m2/year
                  r,        # /year, rate of increase
                  d)        # -, depletion fraction due to fishing
                         {    
  # To check if all have same length or are compatible
  DATA <- data.frame(density=density, sar=sar, r=r, d=d)
  isn  <- which(sar <= 0) 
  Ki <- DATA$density/steadyDensity(K=1, r=DATA$r, sar=DATA$sar, d=DATA$d)
  if (length(isn)) Ki[isn] <- DATA$density[isn]
  Ki[Ki < 0]          <- NA   # These densities are not sustainable 
  Ki[is.infinite(Ki)] <- NA   # Cannot be calculated (r = sar*d)
  Ki
}

## ====================================================================
# a function to estimate constant fishing mortality,
# based on the parameters of the discrete (event-driven) model
## ====================================================================

par_m <- function(sar, 
                  r, 
                  d){
  
  Fcd <- function(m, r, S, d){
    D0d <- 1-d      # initial condition in discrete
    D0c <- 1-d/2    # initial condition in continuous
    
    D1d <- D0d/(D0d+(1-D0d)*exp(-r/S))
    D1c <- D0c*(r-m)/(r*D0c+((r-m)-r*D0c)*exp(-(r-m)/S))
    
    D1d*(1-d/2)-D1c  # mean before &after fishing in discrete - continuous
  }
  
  D <- pmin(d, 1-1e-5)
  DD <- data.frame(sar=sar, r=r, d=D)
  
  # check for NA in inputs
  wisna <- apply(DD, MARGIN=1, FUN=function(x) any(is.na(x)))
  DD[wisna,] <- 1
  
  m <- sapply(1:nrow(DD), FUN=function(i) 
    uniroot(f=Fcd, interval=c(0, 1000), r=DD$r[i], S=DD$sar[i], d=DD$d[i])$root)
  m[d>=1] <- NA
  m[wisna] <- NA
  m
}


par_m_2 <- function(sar, 
                  r, 
                  d){

  n <- 1000
  
  Fcd <- function(m, r, S, d){
    
   D0d <- eventDensity(sar=S, r=r, d=d, n=n)      # initial condition in discrete
   D0c <- (r-m)/(r-m*exp(-(r-m)*n/S))    # initial condition in continuous

   # average density in interval after fishing - only part within log
   # (K/r*log(1-D0/K*(1-3^(r/S))))  
   D1d <- D0d*(1-exp(-r/S))                   # discrete
   D1c <- D0c*r/(r-m)*(1-exp(-(r-m)/S)) # continuous
   if (is.nan(D1c)) D1c <- 0
   D1d-D1c  # mean before &after fishing in discrete - continuous
  }
  
  D  <- pmin(d, 1-1e-5)                # depletion rate
  DD <- data.frame(sar=sar, r=r, d=D)  # so that all have the same length
  
  # check for NA in inputs
  wisna <- apply(DD, MARGIN=1, FUN=function(x) any(is.na(x)))
  DD[wisna,] <- 1
  
  # check for NA in inputs
  wis0 <- which(DD$d==0)
  DD[wis0,] <- 1
  
  m <- sapply(1:nrow(DD), FUN=function(i) 
        uniroot(f=Fcd, interval=c(0, 1000), r=DD$r[i], S=DD$sar[i], d=DD$d[i])$root)
  m[d>=1] <- NA
  m[wisna] <- NA
  m[wis0] <- 0
  m
}

