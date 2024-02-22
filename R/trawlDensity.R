## ====================================================================
## ====================================================================
## Estimating trawling density
## ====================================================================
## ====================================================================

steadyDensity <- function(K=1, sar=1, r=1, d=0.1){
  m <- par_m(sar=sar, r=r, d=d)
  ST <- K*(1-m/r)
  ST
}

eventDensity <- function(K=1, sar=1, r=1, d=0.1, p=1-d, D0=K, n=1){ # n=number of fishing events
  if (n==0) return(D0)
  if (abs(trunc(n)-n) > 1e-8) stop ("'n' should be an integer")
  a <- (1-exp(-r/sar))/K
  i <- 1:n
  p^n/(sum(a * p^(n-i)*exp(-(i-1)*r/sar))+1/D0*exp(-n*r/sar))
}

intervalDensity <- function(K=1, sar=1, r=1, d=0.1, p=1-d, D0=K, n=1){ # n=number of fishing events
  Di = eventDensity(K=K, sar=sar, r=r, d=d, p=p, D0=D0, n=n)  #density after trawling
  return(K/r*log(1-Di/K*(1-exp(r/sar))))
}

