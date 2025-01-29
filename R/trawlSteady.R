## ====================================================================
## ====================================================================
## Estimating steady-state trawling density in perturbation model
## ====================================================================
## ====================================================================

steady_metier  <- function (K = 1, r = 1,  
                parms = data.frame(K = K, r = r), 
                d = data.frame(0.1, 0.1), sar = data.frame(1, 2), 
                D0 = parms[["K"]], tol = 1e-06){
  
  # value of d: weighted mean; value of sar: sum of sar
  rS <- rowSums(sar)
  ss <- as.matrix(sar/rowSums(sar))
  dd <- rowSums(sweep(as.matrix(d), MARGIN = 2, STATS = ss, FUN= "*"))
  sp <- steady_perturb(K = parms$K, r = parms$r, 
                       d = dd, sar = rS, D0 = D0, tol = tol)
  names(d)   <- paste("d", 1:ncol(d), sep="_")
  names(sar) <- paste("sar", 1:ncol(sar), sep="_")
  sp <- sp[, c("K", "r", "density", "ntrawl", "time")]
  sp <- cbind(d, sar, sp)  
  sp
}

steady_perturb <- function(K     = 1, r = 1, d = 0.1,  
                          parms = data.frame(K = K, r = r, d = d),
                          sar   = 1, D0 = parms[["K"]], 
                          tol = 1e-6){
    
  parms <- as.data.frame(parms)
  pn <- names(parms)
  if (any(!c("K", "r", "d") %in% pn))
    stop ("'parms' should contain values for 'K', 'r', 'd' ")
  
    # So that all these parameters are the same length
  parms <- data.frame(parms[, c("K", "r", "d")], sar = sar, D0 = D0)
  
    NAs <- unique(unlist(apply(parms, 
                               MARGIN = 2, 
                               FUN = function(x) which(is.na(x)))))
    
    if (length(NAs)) 
      parms[NAs, ] <- 0
    
#    if (any (is.na(parms)))
#      stop ("cannot estimate steady-state density: input contains NAs")
    
    if (tol[1] < 1e-10)
      warning ("tol too small, taken as 1e-10")
    
    tol[1] <- max(tol[1], 1e-10)
    
    nspec <- as.integer(nrow(parms))
    
    DD <- .Fortran("perturb_steady", nspec = nspec, 
                   sar  = as.double(parms$sar), 
                   K    = as.double(parms$K), 
                   r    = as.double(parms$r), 
                   d    = as.double(parms$d), 
                   atol = as.double(1e-10), 
                   rtol = as.double(tol[[1]]),
                   
                   steadydens   = as.double (parms$D0),
                   steadytrawl  = as.integer(rep(0,  times = nspec)),
                   steadybefore = as.double(rep(0.,  times = nspec)),
                   steadyafter  = as.double(rep(0.,  times = nspec)),
                   steadymean   = as.double(rep(0.,  times = nspec)),
                   steadytimes  = as.double (rep(0., times = nspec)) 
                   )
    
    
    data.frame(parms,
               density_before = DD$steadybefore, # density before trawl
               density_after  = DD$steadyafter,  # density before trawl
               density        = DD$steadymean,   # density averaged
               ntrawl         = DD$steadytrawl,
               time           = DD$steadytimes)
    
}

## ====================================================================
## ====================================================================
## Estimating steady-state trawling density in perturbation model
## ====================================================================
## ====================================================================

steady_logistic <- function(K     = 1, r = 1, m = 0.1, 
                           parms = data.frame(K = K, r = r, m = m), 
                           D0    = parms[["K"]], 
                           tol   = 1e-6){
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  
  # So that all these parameters are the same length
  parms <- data.frame(parms, D0 = D0)
  
  if (any (is.na(parms)))
    stop ("cannot estimate steady-state density: input contains NAs")
  
  pn <- names(parms)
  if (any(!c("K", "r", "m") %in% pn))
    stop ("'parms' should contain values for 'K', 'r', 'm' ")
  
  # So that all these parameters are the same length
  parms <- data.frame(parms[, c("K", "r", "m")], D0 = D0)
  
  if (tol < 1e-10)
    warning ("tol too small, taken as 1e-10")
  tol <- max(tol, 1e-10)
  
  with(parms, {
    
    Ds <- K*(1-m/r)
    
    DD <- pmin(K, Ds + tol)
    DD <- pmax(DD, tol)
    Ds[Ds < 0] <- 0
    DD <- (K*(r-m)*D0 - DD*r*D0)/(DD*((r-m)*K-r*D0))
    ts <- -log(DD)/(r-m)
    ts[Ds == K] <- 0
    ts[ts < 0] <- NA
    ts[is.infinite(ts)] <- NA
    
    return(data.frame(parms, 
                      density = pmax(0, Ds), 
                      time    = ts))
  })
}


## ====================================================================
## General function (not exported)
## ====================================================================

steady_density <- function(K   = 1, r = 1, d = 0.1, m = NULL, 
                          parms = data.frame(K = K, r = r, d = d, m = m),
                          sar = 1, ...){
  if (is.null(m)){
    
   return(steady_perturb(parms = parms, 
                         sar = sar, 
                         D0 = parms["K"], ...)$density)
  }

  if (! is.data.frame(parms)) 
    parms <- as.data.frame(parms)  
  with(parms, return(K*(1-m/r)))
}

