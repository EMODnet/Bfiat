
## ====================================================================
## ====================================================================
## density at particular time for the perturbation model
## ====================================================================
## ====================================================================

density_perturb <- function(K   = 1, r = 1, d = 0.1, 
                          parms = data.frame(K = K, r = r, d = d),
                          sar   = 1, 
                          D0    = parms[["K"]], 
                          times = 1,  # time(s) to estimate density
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
       
       DD <- .Fortran("perturb_times2", nspec = nspec, ntimes = ntimes,
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


## ====================================================================
## ====================================================================
## density at particular time for the perturbation model with various metiers
## ====================================================================
## ====================================================================

density_metier <- function(K   = 1, r = 1, 
                           parms = data.frame(K = K, r = r),
                           d     = data.frame(0.1, 0.1), 
                           sar   = data.frame(1, 2), 
                           D0    = parms[["K"]], 
                           times = 1,  # time to estimate density
                           tend_perturb = max(times),
                           as.deSolve = FALSE){
  pn <- names(parms)
  if (any(!c("K", "r") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r' in density_metier")
  }
  
  # So that all these parameters are the same length
  parms  <- data.frame(parms, D0 = D0)
  
  nspec <- nrow(parms)
  
  if (! (is.data.frame(d) | is.matrix(d)))
      stop ("'d' should be a data.frame or a matrix")
  if (! (is.data.frame(sar) | is.matrix(sar)))
    stop ("'sar' should be a data.frame or a matrix")
  
  nmetier <- ncol(sar)
  if (ncol(d) != nmetier) 
    stop ("'sar' and 'd' not compatible : should have the same number of columns")

  if (nrow(d) == 1 & nspec > 1)
    d <- matrix(nrow  = nspec, 
                ncol  = nmetier, 
                byrow = TRUE, 
                data  = unlist(d))
  else   if (nrow(d) != nspec)
    stop ("'d' should have same number of rows as 'parms'")
  
  else if (is.data.frame(d))
    d <- as.matrix(d)
  
  if (nrow(sar) == 1 & nspec > 1)
    sar <- matrix(nrow  = nspec, 
                  ncol  = nmetier, 
                  byrow = TRUE, 
                  data  = unlist(sar))
  else   if (nrow(sar) != nspec)
    stop ("'sar' should have same number of rows as 'parms'")
  
  if (is.data.frame(sar))
    sar <- as.matrix(sar)
  
  tend_perturb <- rep(tend_perturb, length.out = nmetier)

  # check for NA in inputs
    wisna          <- apply(parms, MARGIN=1, FUN=function(x) any(is.na(x)))
    parms[wisna, ] <- 1
  
    ntimes <- length(times)
    
    DD <- .Fortran("metier_time", 
                   nspec = nspec, nmetier = nmetier, ntimes = ntimes,
                   B0    = as.double(parms$D0),
                   K     = as.double(parms$K), 
                   r     = as.double(parms$r), 
                   d     = as.double(d),  
                   sar   = as.double(sar),  
                   times = as.double(times), 
                   tend_perturb = as.double(tend_perturb),
                   B     = matrix(nrow = nspec, ncol = ntimes, data = -999.)
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
   return(Res)
}  


