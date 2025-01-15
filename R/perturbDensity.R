
## ====================================================================
## ====================================================================
## density at particular time for the perturbation model
## ====================================================================
## ====================================================================

density_perturb <- function(K   = 1, r = 1, d = 0.1, 
                          parms = data.frame(K = K, r = r, d = d),
                          sar   = 1, 
                          D0    = parms[["K"]], 
                          times = 1,  # time to estimate density
                          tend_perturb = max(times),
                          as.deSolve = FALSE){
  pn <- names(parms)
  if (any(!c("K", "r", "d") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r', 'd' in discrete models")
  }
  
  if ("times" %in% names(parms)){
    dens <- density_perturb2(parms = parms, 
                           sar   = sar, 
                           D0    = D0)$density
    Res <- data.frame(parms, sar = sar, D0 = D0, density = dens)
    
  } else {
    
    if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
    
    # So that all these parameters are the same length
    if (! "sar" %in% names(parms))
       parms  <- data.frame(parms, sar = sar, D0 = D0)

    pn <- names(parms)
    if (any(!c("K", "r", "d", "sar") %in% pn)){
      stop ("'parms' should contain values for 'K', 'r', 'd', 'sar' ")
       }
       
    # So that all these parameters are the same length
       parms  <- data.frame(parms[, c("K", "r", "d", "sar")], D0 = D0)
       
       nspec  <- nrow(parms)
       
       # check for NA in inputs
       wisna          <- apply(parms, MARGIN=1, FUN=function(x) any(is.na(x)))
       parms[wisna, ] <- 1
       ntimes <- length(times)
       
       DD <- .Fortran("logistictrawl2", nspec = nspec, ntimes = ntimes,
                      B0    = as.double(parms$D0),
                      K     = as.double(parms$K), 
                      r     = as.double(parms$r), 
                      d     = as.double(parms$d), 
                      sar   = as.double(parms$sar),
                      times = as.double(times), 
                      tend_perturb = as.double(tend_perturb[1]),
                      B = matrix(nrow = nspec, ncol = ntimes, data = -999.)
       )
       
       if (length(wisna))
         DD$B[wisna,] <- NA
       
       if (! as.deSolve){
         colnames(DD$B) <- paste("Dt", times, sep="")
         Res  <- data.frame(parms, DD$B)
         attributes(Res)$times <- times
       } else {
           Res <- cbind(times, t(DD$B))
           colnames(Res)[1] <- "times"
           colnames(Res)[-1] <- row.names(parms)
           class(Res) <- c("deSolve", class(Res))
       }
   }
  return(Res)
}  

## ====================================================================
## ====================================================================
## density at particular time - not very optimal
## ====================================================================
## ====================================================================

density_perturb2 <- function(K   = 1, r = 1, d = 0.1, times = 1,
                            parms = data.frame(K = K, r = r, d = d, times = times),
                            sar   = 1, 
                            D0    = parms[["K"]]
                            ){ # time to estimate density
  
  parms$n <- trunc(parms$times/sar)+1
  
  #density after trawling at start 
  Di = eventDensity2(parms = parms, 
                    sar   = sar, 
                    D0    = D0)  
  
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  
  # So that all these parameters are the same length
  parms  <- data.frame(parms, sar = sar, D0 = D0)
  
  tt   <- parms$times - (parms$n-1)/parms$sar
  dens <- with(parms, K*Di/(Di + (K-Di)*exp(-r*tt)))
  
  data.frame(parms, density = dens)
}  

