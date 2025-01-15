

## ====================================================================
## ====================================================================
## density after and between trawling events
## ====================================================================
## ====================================================================

density_trawl <- function(K   = 1, r = 1, d = 0.1, 
                          parms = data.frame(K = K, r = r, d = d),
                          sar = 1, 
                          D0  = parms[["K"]], 
                          n   = 1){ # n=number of fishing events
  
  #density after trawling at start 
  Di = eventDensity(parms = parms, 
                    sar   = sar, 
                    D0    = D0, 
                    n     = n)  
  
  #density after trawling at end 
  Di2 = eventDensity(parms = parms, 
                     sar   = sar, 
                     D0    = D0, 
                     n     = n+1) 
  
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  
  start_t <- (n-1)/sar
  end_t   <- n/sar
  
  # So that all these parameters are the same length
  parms  <- data.frame(parms[,c("K", "r", "d")], sar = sar, D0 = D0)
  
  #  return(K/r*log(1-Di/K*(1-exp(r/sar))))
  mean_dens <-   with(parms, K + sar*K/r*log(((K-Di)*exp(-r/sar)+Di)/K))
  
  X1 <- data.frame(parms, mean_dens)
  X1 <- reshape(X1, direction = "long", 
                varying = list(6:ncol(X1)), 
                v.names = "density")  # the mean
  
  X2 <- data.frame(parms, Di)
  X2 <- reshape(X2, direction = "long", 
                varying = list(6:ncol(X2)), 
                v.names = "density_start")
  
  X3 <- data.frame(parms, Di2/(1-parms$d))
  X3 <- reshape(X3, direction = "long", 
                varying = list(6:ncol(X3)), 
                v.names = "density_end")
  XX <- merge(merge(X1, X2), X3)
  
  XX$ntrawl_from   <- n[XX$time]
  XX$ntrawl_to     <- n[XX$time]+1
  XX$times         <- start_t[XX$time]
  XX$end_times     <- end_t[XX$time]
  
  XX <- XX[order(XX$id, XX$time),]
  
  return(XX[,c(names(parms), "times", "end_times", "ntrawl_from", "ntrawl_to",
               "density", "density_start", "density_end")])                   
}

## ====================================================================
## ====================================================================
## Density after a trawling event
## ====================================================================
## ====================================================================

eventDensity <- function(K   = 1, r = 1, d = 0.1, 
                         parms = data.frame(K = K, r = r, d = d),
                         sar = 1, 
                         D0  = parms[["K"]], 
                         n   = 1){ # n=number of fishing events
  
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  
  pn <- names(parms)
  if (any(!c("K", "r", "d") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r', 'd' in discrete models")
  }
  
  # So that all these parameters are the same length
  parms  <- data.frame(parms, sar = sar, D0 = D0)
  nspec  <- nrow(parms)
  nevent <- length(n)
  
  # check for NA in inputs
  wisna          <- apply(parms, MARGIN=1, FUN=function(x) any(is.na(x)))
  parms[wisna, ] <- 1
  
  n[which(is.na(n))] <- 0
  if (any(abs(trunc(n) - n) > 1e-8) )
      stop ("'n' should be an integer, one value")
  
  ii <- which(n == 0)
  if (length(ii)) n[ii] <- 1 

  if (length(n) > 1) n <- n[order(n)]

  DD <- .Fortran("eventdensity", nspec = nspec, nevent = nevent, 
                     B0 = as.double(parms$D0),
                     sar = as.double(parms$sar),
                     K  = as.double(parms$K), 
                     r  = as.double(parms$r), 
                     d  = as.double(parms$d), 
                     eventnr = as.integer(n), 
                     B = matrix(nrow=nspec, ncol=nevent, data=as.double(-999.))
                     )
      
    if (length(ii))
      DD$B[ii] <- parms$D0[ii]
    if (length(wisna))
      DD$B[wisna] <- parms$D0[wisna]
    DD$B  
}

## ====================================================================
## same but now for each species its own event!
## ====================================================================

eventDensity2 <- function(K   = 1, r = 1, d = 0.1, n   = 1,
                         parms = data.frame(K = K, r = r, d = d, n = n),
                         sar = 1, 
                         D0  = parms[["K"]]
                         ){ 
  
  if (! is.data.frame(parms)) parms <- as.data.frame(parms)  
  
  pn <- names(parms)
  if (any(!c("K", "r", "d", "n") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r', 'd', 'n' ")
  }
  # So that all these parameters are the same length
  parms  <- data.frame(parms, sar = sar, D0 = D0)
  nspec  <- nrow(parms)

  # check for NA in inputs
  wisna          <- apply(parms, MARGIN=1, FUN=function(x) any(is.na(x)))
  parms[wisna, ] <- 1
  
  ii <- which(parms$n == 0)
  if (length(ii)) parms$n[ii] <- 1 
  
# eventdensity2(nspec, B0, sar,                                & 
#                  K, r, d, eventnr, B)
  DD <- .Fortran("eventdensity2", nspec = nspec, 
                 B0 = as.double(parms$D0),
                 sar = as.double(parms$sar),
                 K  = as.double(parms$K), 
                 r  = as.double(parms$r), 
                 d  = as.double(parms$d), 
                 eventnr = as.integer(parms$n), 
                 B = rep(-999., times = nspec)
  )
  
  if (length(ii))
    DD$B[ii] <- parms$D0[ii]
  if (length(wisna))
    DD$B[wisna] <- parms$D0[wisna]
  DD$B  
}

