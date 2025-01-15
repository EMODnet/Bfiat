
## ====================================================================
## ====================================================================
## density at particular time for the logistic model
## ====================================================================
## ====================================================================

density_logistic <- function(K = 1, r = 1, m = 0.1, 
                            parms = data.frame(K = K, r = r, m = m),
                            D0    = parms[["K"]], 
                            times = 1, # time to estimate density
                            tend_perturb = max(times),
                            as.deSolve = FALSE){ 
  
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  pn <- names(parms)
  
  if (any(!c("K", "r", "m") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r', 'm' in logistic models")
  }
  
  TT <- parms$times
  
  if (! is.null(TT)){
    # So that all these parameters are the same length
    parms  <- data.frame(parms[, c("K", "r", "m", "times")], D0 = D0)
    
      rn <- parms$r - parms$m
      D  <- D0
      DD <- rn *parms$K *D/(parms$r*D + 
                        (rn *parms$K-parms$r*D)*exp(-rn*(parms$times)))
      
      return(data.frame(parms, density = DD))
      
  } else {
    
 
  
    nspec  <- nrow(parms)
  
  # check for NA in inputs
    wisna          <- apply(parms, MARGIN=1, FUN=function(x) any(is.na(x)))
    parms[wisna, ] <- 1
  
    ntimes <- length(times)
  
    DD <- .Fortran("logistic", nspec = nspec, ntimes = ntimes,
                 B0    = as.double(parms$D0),
                 K     = as.double(parms$K), 
                 r     = as.double(parms$r), 
                 m     = as.double(parms$m), 
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

