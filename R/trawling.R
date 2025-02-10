## ====================================================================
## ====================================================================
## Model of impact of bottom disturbances on density or biomass
## ====================================================================
## ====================================================================

run_perturb <- function(
                    parms,                     # list/data.frame with r, k, d
                    times,                     # output times, in years
                    sar = 1,
                    tstart_perturb = min(times) + 0.5/sar,
                    tend_perturb = max(times), # time at which perturbation stops
                    events = NULL,             # times of trawling events
                    taxon_names  = parms[["taxon"]], # names of taxa
                    D0     = parms[["K"]],       # initial condition
                    addsum = FALSE, 
                    verbose = FALSE){
  
  if (is.vector(parms)) 
    parms <- as.list(parms)
  else parms <- as.data.frame(parms)
  
  if (length(sar) > 1) stop ("sar should be one number")
  
  if (is.null(events)) 
    events <- seq(from = max(times[1], tstart_perturb), 
                  by   = 1/sar, 
                  to   = tend_perturb)
  
  pn <- names(parms)
  if (any(!c("K", "r", "d") %in% pn))
    stop ("'parms' should contain values for 'K', 'r', 'd' ")
  r  <- parms[["r"]]
  K  <- parms[["K"]]
  d  <- parms[["d"]]
  
  include <- (1:length(K))
  omit    <- unique(c(which(is.na(r)), 
                      which(is.na(K)), 
                      which(is.na(d))))
  
  if (length(omit)){
    if (verbose) warning("parameters not given for",length(omit),
      " variables- check attributes()$omit")
    include <- include[-omit] 
  }
  
  C0  <- D0
  nC0 <- 1  
  
  if (is.data.frame(D0) | is.matrix(D0)) {
    if (nrow(D0) != nrow(parms))
      stop(" 'D0' should have as many rows as 'parms'")
    nC0 <- ncol(D0)
    
    # If more runs need to be done: 
    # model runs with relative density first, i.e. with C0=1
    # then these are multiplied with the actual initial densities (using sweep)
    if (nC0 > 1) {
      if (any(is.character(D0[,1]))) 
        stop("'first column of D0 contains strings - should be numeric")
      C0 <- rep(1, times=length(r))
    }
  } else if (length(D0) == 1) 
    C0 <- rep(D0, times=length(r))

  C0 <- unlist(C0)
  r  <- r[include]
  K  <- K[include]
  d  <- d[include]
  C0 <- C0[include]

  # trawling events must lie within times  
  t0     <- times[1]
  tend   <- times[length(times)]
  trawl  <- events[events >= t0 & events <= tend]
  events <- sort(unique(c(trawl, tend+0.1)))
  
  # trawls that are not in times are added to times
  tnew   <- trawl[!trawl %in% times]
  times  <- sort(unique(c(times, trawl)))
  
  nspec  <- as.integer(length(C0))
  nevent <- as.integer(length(events))
  ntimes <- as.integer(length(times))
  B0     <- as.double(C0)
    
  DD <- .Fortran("perturb_times", nspec=nspec, nevent=nevent, 
          ntimes=ntimes, B0=B0, K=as.double(K), r=as.double(r), 
          d=as.double(d), times=as.double(times), events=as.double(events), 
          B=matrix(nrow=nspec, ncol=ntimes, data=-0.999), 
          dTrawl=matrix(nrow=nspec, ncol=nevent, data=-0.999))
  
  Bt <- cbind(times, t(DD$B)) # add dynamic results to times

  if (any(is.nan(Bt))) {
    Bt[is.nan(Bt)] <- 0
    if (verbose) warning("carrying capacity of one species probably=0")
  }
  
  # remove trawls not in 'times'
  if (length(tnew))                             
    Bt <- Bt[-which(times%in% tnew), ]
  
  taxon <- unlist(taxon_names)
  if (is.null(taxon))   
    taxon <- paste("tax", 1:(ncol(Bt)-1), sep="_")
  else if (is.factor(taxon)) 
    taxon <- as.character(taxon)[include]
  else 
    taxon <- taxon[include]

  colnames(Bt) <- c("time", as.vector(unlist(taxon)))

  if (length(trawl)) {
    dTrawl           <- t(DD$dTrawl)
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
      C0 <- unlist(D0[include,i])
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
    names(RES) <- colnames(D0)
    return(RES)  
  }
}

## ====================================================================
## ====================================================================
## Model of impact of bottom disturbances on density or biomass
## ====================================================================
## ====================================================================

run_logistic <- function(
                  parms,                      # list/data.frame with r, k, m
                  times,                      # output times, in years
                  tstart_perturb = min(times),
                  tend_perturb = max(times),  # time at which perturbation stops
                  taxon_names   = parms[["taxon"]], # names of taxa
                  D0      = parms[["K"]],     # initial condition
                  addsum  = FALSE, 
                  verbose = FALSE){
  
  if (is.vector(parms)) 
    parms <- as.list(parms)
  else parms <- as.data.frame(parms)
  
  pn <- names(parms)
  
  if (!"m"  %in% pn)  # try to estimate mortality from r, d, and sar
    if (all(c("sar", "r", "d") %in% pn) )
      parms$m <- with(parms, par_m(sar=sar, r=r, d = d, K=K)) ## ADDED K=K
  pn <- names(parms)
  
  if (any(!c("K", "r", "m") %in% pn)){
    stop ("'parms' should contain values for 'K', 'r', 'm' in continuos logistic model ")
  }
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
  
  C0  <- D0
  nC0 <- 1  
  
  if (is.data.frame(D0) | is.matrix(D0)) {
    if (nrow(D0) != nrow(parms))
      stop(" 'D0' should have as many rows as 'parms'")
    nC0 <- ncol(D0)
    
    # If more runs need to be done: 
    # model runs with relative density first, i.e. with C0=1
    # then these are multiplied with the actual initial densities (using sweep)
    if (nC0 > 1) {
      if (any(is.character(D0[,1]))) 
        stop("'first column of D0 contains strings - should be numeric")
      C0 <- rep(1, times=length(r))
    }
  } else if (length(D0) == 1) 
    C0 <- rep(D0, times=length(r))

  C0 <- unlist(C0)
  r  <- r[include]
  K  <- K[include]
  m  <- m[include]
  C0 <- C0[include]

  nspec  <- as.integer(length(C0))
  ntimes <- as.integer(length(times))

  Logistic <- function(i, t){
    rn <- r[i]-m[i]
    D  <- C0[i]
    rn*K[i]*D/(r[i]*D + (rn*K[i]-r[i]*D)*exp(-rn*(t)))
  }
  
  Logistic2 <- function(i, t){  # No mortality
    rn <- r[i]
    D  <- C0[i]
    rn*K[i]*D/(r[i]*D + (rn*K[i]-r[i]*D)*exp(-rn*(t)))
  }
  
  tb    <- times[1]
  t0    <- times[times < tstart_perturb]  # times with perturbations
  B0    <- NULL
  if (length(t0)) {
    t0 <- t0 - tb
    B0 <- outer(X   = 1:nspec, 
                Y   = t0, 
                FUN = function(X,Y) Logistic2(X,Y)) # no fish mortality
    C0 <- outer(X   = 1:nspec, 
                Y   = tstart_perturb, 
                FUN = function(X,Y) Logistic2(X,Y)) # no fish mortality
    tb <- tstart_perturb
  }
  
  t1    <- times[times >= tstart_perturb & 
                 times <= tend_perturb]  # times with perturbations
  lent1 <- length(t1)
  
  t1    <- t1 - tb
  
  B  <- outer(X   = 1:nspec, 
              Y   = t1, 
              FUN = function(X,Y)Logistic(X,Y))
  
  C0 <- outer(X   = 1:nspec, 
              Y   = tend_perturb, 
              FUN = function(X,Y) Logistic(X, Y)) # no fish mortality
  
  if (nrow(B) > 1)
    B <- cbind(B0, B)
  else 
    B <- c(B0, B)
  
  t2    <- times[times > tend_perturb]   # times without perturbations
  
  if (length(t2) >= 1){  
   t2 <- t2 - tend_perturb 
   B2  <- outer(X   = 1:nspec, 
                Y   = t2, 
                FUN = function(X,Y) Logistic2(X,Y))
   if (nrow(B2) > 1)
     B <- cbind(B, B2)
   else 
     B <- c(B, B2)
   
  }
  
  if (is.vector(B)) B <- matrix(nrow = 1, data = B)
  
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
  
  taxon <- unlist(taxon_names)
  if (is.null(taxon))   
    taxon <- paste("tax", 1:(ncol(Bt)-1), sep="_")
  else if (is.factor(taxon)) 
    taxon <- as.character(taxon)[include]
  else 
    taxon <- taxon[include]

  colnames(Bt) <- c("time", as.vector(unlist(taxon)))

  ADDsum <- function(run)  # rowSums does not work on a vector
    if (ncol(run) >2)  
      return(cbind(run, sum=rowSums(run[,-1])))
    else
      return(cbind(run, sum=run[,-1]))

  if (nC0 == 1){  # Model applied for only one 'station' - return a matrix
   
    if (addsum) Bt <- ADDsum(Bt)

    class(Bt) <- c("deSolve", "matrix")
    return(Bt)
    
  } else { # Model applied for many 'station's - return a list of matrices

    RES <- list()  
    for (i in 1:nC0){ 
      C0 <- unlist(D0[include,i])
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
    names(RES) <- colnames(D0)
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
 
 
