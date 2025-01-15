## ====================================================================
## ====================================================================
## Estimating trawling parameters
## ====================================================================
## ====================================================================

## ====================================================================
# a function to estimate depletion rate based on depth occurrence of species
## ====================================================================

par_d <- function(gpd   =     1 ,  # gear penetration depth, cm
                  m_d   =   0.055, # mortality parameter, /cm
                  m_max =   0.42,  # maximum depletion in any layer
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
    gpd <- rep(gpd, times = nrow(fDepth))
    
  } else if (nrow(fDepth) ==1 & length(gpd) > 1) 
    fDepth <- matrix(ncol = ncol(fDepth), nrow = length(gpd), 
                     data = fDepth, byrow = TRUE)
  
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

par_r <- function(age_at_maturity = NULL,  # age at maturity, years
                  longevity       = NULL)  {     # lifetime, years
  
  if (is.null(age_at_maturity) & ! is.null(longevity)){                       
     ri <- 5.31/longevity
  
  } else if (is.null(age_at_maturity)) {
    stop ("either longevity or age_at_maturity should have a value")
  
  } else {
     ri <- 5.31*0.7844/age_at_maturity
  }
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
                  d,        # -, depletion fraction due to fishing
                  t_density = NULL) # time at which density was estimated; NULL=steady-state
                         {    
  # To check if all have same length or are compatible
  DATA <- data.frame(density=density, sar=sar, r=r, d=d, K = 1)
  
  isn  <- which(sar <= 0) 
  
  if (is.null(t_density))     # steady state at time infinite, using K=1
     Ki   <- DATA$density/steady_perturb(K = 1, 
                                         r = DATA$r, d = DATA$d,  
                                         sar   = DATA$sar)$density
  
  else if (t_density > 0) {   # Density is estimated after t_density years of fishing  
     
    ########### NEED TO CHANGE THIS  ################
    # value at t_density, starting with K=1
    
    #SUBROUTINE eventdensity2(nspec, B0, sar,                         & 
    #                           K, r, d, eventnr, B)
    DATA$times = t_density
    D_t <- density_perturb2(parms = DATA[, c("r", "d", "K", "times")], 
                           sar   = DATA$sar)$density
     Ki <- DATA$density / D_t 
  } else Ki <- DATA$density
  
  if (length(isn)) Ki[isn] <- DATA$density[isn]
  Ki[Ki < 0]          <- NA   # These densities are not sustainable 
  Ki[is.infinite(Ki)] <- NA   # Cannot be calculated (extinct)
  Ki
}

## ====================================================================
# a function to estimate constant fishing mortality,
# based on the parameters of the discrete (event-driven) model
## ====================================================================

par_m <- function(sar, 
                  r,  
                  K, 
                  d, 
                  refD = K){
    
    # Assume Dens = 1
    Fcd <- function(m, r, S, d, K, refD){
      D0d <- 1-d      # initial condition in discrete
      D0c <- 1-d/2    # initial condition in continuous
      
      D1d <- D0d/(D0d + (K/refD-D0d)*exp(-r/S))
      D1c <- (r-m)/(r*D0c+((r-m)*K/refD-r*D0c)*exp(-(r-m)/S))
      
      D1d-D1c  # mean before &after fishing in discrete & continuous should be =
    }
    
    D  <- pmin(d, 1-1e-5)
    DD <- data.frame(sar=sar, r=r, d=D, K = K, refD = refD)  # all have equal length
    
    # check for NA in inputs
    wisna      <- apply(DD, MARGIN=1, FUN=function(x) any(is.na(x)))
    DD[wisna,] <- 1
    
    m <- sapply(1:nrow(DD), FUN=function(i) 
      uniroot(f=Fcd, interval=c(0, 1000), 
              r=DD$r[i], S=DD$sar[i], d=DD$d[i], K = DD$K[i], 
              refD = DD$refD[i], tol=1e-8)$root)
    m[d>=1]  <- NA
    m[wisna] <- NA
    m
  }

# older version with K = 1
par_m_1 <- function(sar,   
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
    uniroot(f=Fcd, interval=c(0, 1000), 
            r=DD$r[i], S=DD$sar[i], d=DD$d[i], tol=1e-8)$root)
  m[d>=1] <- NA
  m[wisna] <- NA
  m
}

