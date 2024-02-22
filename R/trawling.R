## ====================================================================
## ====================================================================
## Model of impact of bottom disturbances on density or biomass
## ====================================================================
## ====================================================================

perturb <- function(parms,                 # list/data.frame with r, k, d
                    times,                 # output times, in years
                    events=NULL,           # times of trawling events
                    taxon=parms["taxon"],  # names of taxa
                    Cini=parms["K"],       # initial condition
                    addsum=FALSE, 
                    verbose=FALSE){
  
  if (is.vector(parms)) 
    parms <- as.list(parms)
  else parms <- as.data.frame(parms)
  
  pn <- names(parms)
  if (any(!c("K", "r", "d") %in% pn))
    stop ("'parms' should contain values for 'K', 'r', 'd' ")
  r  <- parms[["r"]]
  K  <- parms[["K"]]
  d  <- parms[["d"]]
  
  include <- (1:length(K))
  omit    <- unique(c(which(is.na(r)), which(is.na(K)), which(is.na(d))))
  
  if (length(omit)){
    if (verbose) warning("parameters not given for",length(omit),
      " variables- check attributes()$omit")
    include <- include[-omit] 
  }
  
  C0  <- Cini
  nC0 <- 1  
  
  if (is.data.frame(Cini) | is.matrix(Cini)) {
    if (nrow(Cini) != nrow(parms))
      stop(" 'Cini' should have as many rows as 'parms'")
    nC0 <- ncol(Cini)
    
    # If more runs need to be done: 
    # model runs with relative density first, i.e. with C0=1
    # then these are multiplied with the actual initial densities (using sweep)
    if (nC0 > 1) {
      if (any(is.character(Cini[,1]))) 
        stop("'first column of Cini contains strings - should be numeric")
      C0 <- rep(1, times=length(r))
    }
  } else if (length(Cini) == 1) 
    C0 <- rep(Cini, times=length(r))

  C0 <- unlist(C0)
  r  <- r[include]
  K  <- K[include]
  d  <- d[include]
  C0 <- C0[include]

  # trawling events must lie within times  
  t0     <- times[1]
  tend   <- times[length(times)]
  trawl  <- events[events >= t0 & events <= tend]
  events <- sort(unique(c(trawl, tend)))
  
  # trawls that are not in times are added to times
  tnew   <- trawl[!trawl %in% times]
  times  <- sort(unique(c(times, trawl)))
  
  nspec  <- as.integer(length(C0))
  nevent <- as.integer(length(events))
  ntimes <- as.integer(length(times))
  B0     <- as.double(C0)
    
  DD <- .Fortran("logistictrawl", nspec=nspec, nevent=nevent, 
          ntimes=ntimes, B0=B0, K=as.double(K), r=as.double(r), 
          d=as.double(d), times=as.double(times), events=as.double(events), 
          B=matrix(nrow=nspec, ncol=ntimes, data=-0.999), 
          dTrawl=matrix(ncol=nspec, nrow=nevent, data=-0.999))
  
  Bt <- cbind(times, t(DD$B)) # add dynamic results to times

  if (any(is.nan(Bt))) {
    Bt[is.nan(Bt)] <- 0
    if (verbose) warning("carrying capacity of one species probably=0")
  }
  
  # remove trawls not in 'times'
  if (length(tnew))                             
    Bt <- Bt[-which(times%in% tnew), ]
  
  taxon <- unlist(taxon)
  if (is.null(taxon))   
    taxon <- paste("tax", 1:(ncol(Bt)-1), sep="_")
  else if (is.factor(taxon)) 
    taxon <- as.character(taxon)[include]
  else 
    taxon <- taxon[include]

  colnames(Bt) <- c("times", as.vector(unlist(taxon)))

  if (length(trawl)) {
    dTrawl           <- DD$dTrawl
    colnames(dTrawl) <- unlist(taxon)
    dTrawl           <- cbind(times=trawl[order(trawl)], 
                              dTrawl[1:length(trawl),]) 
  } else dTrawl <- NULL

  ADDsum <- function(run)  # rowSums does not work on a vector
    if (ncol(run) >2)  
      return(cbind(run, sum=rowSums(run[,-1])))
    else
      return(cbind(run, sum=run[,-1]))

  if (nC0 == 1){  # Model applied for only one 'station' - return a matrix
   
    if (addsum) Bt <- ADDsum(Bt)

   # instantaneous depletions 
    attributes(Bt)$depletion <- dTrawl
  
    class(Bt) <- c("deSolve", "matrix")
    return(Bt)
    
  } else { # Model applied for many 'station' - return a list of matrices

    RES <- list()  
    for (i in 1:nC0){ 
      C0 <- unlist(Cini[include,i])
      ii <- which(C0 > 0)      
      # Multiply columns (relative densities) with true initial values
      res <- sweep(x      = Bt[,(ii+1)], # relative species densities in columns
                   MARGIN = 2, 
                   STATS  = C0[ii], 
                   FUN    = "*")
      res <- cbind(Bt[,1], res)
      if (addsum) res <- ADDsum(res)

     # instantaneous depletions are also multiplied with initial densities
  
      if (length(trawl)) {
        dT <- sweep(x      = dTrawl[1:length(trawl), (ii+1)], 
                    MARGIN = 2, 
                    STATS  = C0[ii], 
                    FUN    = "*")
        attributes(res)$depletion <- cbind(times=dTrawl[,1], dT) 
      } else attributes(res)$depletion <- NULL
  
      class(res) <- c("deSolve", "matrix")
      RES[[i]] <- res
    }
    names(RES) <- colnames(Cini)
    return(RES)  
  }
}

## ====================================================================
## ====================================================================
## Model of impact of bottom disturbances on density or biomass
## ====================================================================
## ====================================================================

logistic <- function(parms,                  # list/data.frame with r, k, m
                     times,                  # output times, in years
                     tendPerturb=max(times), # time at which perturbation stops
                     taxon=parms["taxon"],   # names of taxa
                     Cini=parms["K"],        # initial condition
                     addsum=FALSE, 
                     verbose=FALSE){
  
  if (is.vector(parms)) 
    parms <- as.list(parms)
  else parms <- as.data.frame(parms)
  
  pn <- names(parms)
  if (any(!c("K", "r", "m") %in% pn))
    stop ("'parms' should contain values for 'K', 'r', 'm' in continuos logistic model ")
  r  <- parms[["r"]]
  K  <- parms[["K"]]
  m  <- parms[["m"]]
  
  include <- (1:length(K))
  omit    <- unique(c(which(is.na(r)), which(is.na(K)), which(is.na(m))))
  
  if (length(omit)){
    if (verbose) warning("parameters not given for",length(omit),
      " variables- check attributes()$omit")
    include <- include[-omit] 
  }
  
  C0  <- Cini
  nC0 <- 1  
  
  if (is.data.frame(Cini) | is.matrix(Cini)) {
    if (nrow(Cini) != nrow(parms))
      stop(" 'Cini' should have as many rows as 'parms'")
    nC0 <- ncol(Cini)
    
    # If more runs need to be done: 
    # model runs with relative density first, i.e. with C0=1
    # then these are multiplied with the actual initial densities (using sweep)
    if (nC0 > 1) {
      if (any(is.character(Cini[,1]))) 
        stop("'first column of Cini contains strings - should be numeric")
      C0 <- rep(1, times=length(r))
    }
  } else if (length(Cini) == 1) 
    C0 <- rep(Cini, times=length(r))

  C0 <- unlist(C0)
  r  <- r[include]
  K  <- K[include]
  m  <- m[include]
  C0 <- C0[include]

  nspec  <- as.integer(length(C0))
  ntimes <- as.integer(length(times))
  B0     <- as.double(C0)
 
  t1    <- times[times <= tendPerturb]  # times with perturbations
  lent1 <- length(t1)
  t1    <- unique(c(t1, tendPerturb))  # times with perturbations
  t2    <- times[times > tendPerturb]             # times without perturbations
  Logistic <- function(i, t){
    rn <- r[i]-m[i]
    D  <- C0[i]
    rn*K[i]*D/(r[i]*D + (rn*K[i]-r[i]*D)*exp(-rn*(t)))
  }
  
  B  <- outer(X=1:nspec, Y=t1, FUN=function(X,Y)Logistic(X,Y))
  if (length(t2) >= 1){  
   t2 <- t2-tendPerturb 
   m[] <- 0
   C0  <- B[, ncol(B)]
   B2 <- outer(X=1:nspec, Y=t2, FUN=function(X,Y)Logistic(X,Y))
   if (nrow(B2) > 1)
     B <- cbind(B[,1:lent1], B2)
   else B <- c(B[1:lent1], B2)
   if (is.vector(B)) B <- matrix(nrow=1, data=B)
  }
#  DD <- .Fortran("logistic", nspec=nspec, 
#          ntimes=ntimes, B0=B0, K=as.double(K), r=as.double(r), 
#          d=as.double(d), times=as.double(times), events=as.double(events), 
#          B=matrix(nrow=nspec, ncol=ntimes, data=-0.999), 
#          mort=matrix(ncol=nspec, nrow=times, data=-0.999))
  
  Bt <- cbind(times, t(B)) # add dynamic results to times

  if (any(is.nan(Bt))) {
    Bt[is.nan(Bt)] <- 0
    if (verbose) warning("carrying capacity of one species probably=0")
  }
  
  taxon <- unlist(taxon)
  if (is.null(taxon))   
    taxon <- paste("tax", 1:(ncol(Bt)-1), sep="_")
  else if (is.factor(taxon)) 
    taxon <- as.character(taxon)[include]
  else 
    taxon <- taxon[include]

  colnames(Bt) <- c("times", as.vector(unlist(taxon)))

  ADDsum <- function(run)  # rowSums does not work on a vector
    if (ncol(run) >2)  
      return(cbind(run, sum=rowSums(run[,-1])))
    else
      return(cbind(run, sum=run[,-1]))

  if (nC0 == 1){  # Model applied for only one 'station' - return a matrix
   
    if (addsum) Bt <- ADDsum(Bt)

    class(Bt) <- c("deSolve", "matrix")
    return(Bt)
    
  } else { # Model applied for many 'station' - return a list of matrices

    RES <- list()  
    for (i in 1:nC0){ 
      C0 <- unlist(Cini[include,i])
      ii <- which(C0 > 0)      
      # Multiply columns (relative densities) with true initial values
      res <- sweep(x      = Bt[,(ii+1)], # relative species densities in columns
                   MARGIN = 2, 
                   STATS  = C0[ii], 
                   FUN    = "*")
      res <- cbind(Bt[,1], res)
      if (addsum) res <- ADDsum(res)

      class(res) <- c("deSolve", "matrix")
      RES[[i]] <- res
    }
    names(RES) <- colnames(Cini)
    return(RES)  
  }
}



#####################
## old version of perturb, in R

# Logistic <- function(B0, times)   # analytical solution
#     (B0*K) / (B0 + (K-B0) * exp(-r*(times)))  # gives NaN for B and K=0
 
#   tout   <- NULL     # output times
#   Bt     <- NULL     # density versus time
#   dTrawl <- NULL     # depletion during each trawling

#  # if simulation starts with an event: reduce density
#   if (events[1] == t0) {
#    dTrawl <- C0*d           # depletion
#    C0     <- C0*(1-d)       # new density
#    events <- events[-1]     # remaining events
#   } 
  
#  # loop over events
#   for (t.event in events){
#    tt <- times[times >= t0 & times < t.event]  # time inbetween events
#    if (t.event == tend)                        # run till end
#      tt <- c(tt, tend)      
#    BN <- outer(C0, tt-t0, FUN=Logistic)        # simulate (tt-t0=time since event)
#    dTrawl <- rbind(dTrawl, BN[, ncol(BN)]*d)   # depletion during trawl
#    C0 <- BN[, ncol(BN)] * (1-d)                # new density after trawling
#    Bt <- rbind(Bt, cbind(tt, t(BN)))           # add results 
#    t0 <- t.event                               # new initial time
#   }
 
 