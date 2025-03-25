
# =============================================================================
# Find critical values for fishing parameters
# =============================================================================

critical_sar <- function(
    gpd    = 1,
    fDepth = 1, uDepth = 0, d = NULL, 
    r      = 0.1,
    crit_DK = 0,  #critical D/K
    ...) {  #arguments passed to par_d function
  
  # depletion for the gpd sequence
  if (is.null (d)){
    if (is.matrix(fDepth) | is.data.frame(fDepth)){
      if (any(abs(rowSums(na.omit(fDepth)) - 1) > 1e-6)) 
        stop (" sum of rows of fDepth should be = 1")
      d <- outer(Y = gpd, X = 1:nrow(fDepth), 
                 FUN = function(i, gpd) 
                   par_d(gpd, uDepth = uDepth, fDepth = fDepth[i,], ...))
    }  else {
      if (sum(fDepth) != 1) stop (" sum of fDepth should be = 1")
      d <-  par_d(gpd, uDepth = uDepth, fDepth = fDepth, ...)
    }
  }
  # Steady-state density can be estimated using this formula
  # (from Pritchard, 2017): D = 1+ sar/r * log(1-d)
  
  # For DK = crit_DK
  sar <- (crit_DK - 1.) / log(1-d)*r
  sar[is.infinite(sar)] <- NA
  sar
}

# =============================================================================

critical_r <- function(
    sar     = 1,
    gpd     = 1,
    fDepth  = 1, uDepth = 0, d = NULL,
    crit_DK = 0,
    ...) {  #arguments passed to par_d function
  
  # depletion for the gpd sequence
  if(is.null(d)){
    if (is.matrix(fDepth) | is.data.frame(fDepth)){
      if (any(abs(rowSums(na.omit(fDepth)) - 1) > 1e-6)) 
          stop (" sum of rows of fDepth should be = 1")
      d <- outer(Y = gpd, X = 1:nrow(fDepth), 
                 FUN = function(i, gpd) 
                        par_d(gpd, uDepth = uDepth, fDepth = fDepth[i,], ...))
   } else {
     if (sum(fDepth) != 1) stop (" sum of fDepth should be = 1")
     d <-  par_d(gpd, uDepth = uDepth, fDepth = fDepth, ...)
   }
  }
  
  # Steady-state density/K can be estimated using this formula
  # (from Pritchard, 2017): D = 1+ sar/r * log(1-d)
  # (D-1) / log(1-d)/sar
  # For DK = crit_DK
  r <- (crit_DK - 1.) / log(1-d)/sar
  r[is.infinite(r)] <- NA
  r
}

# =============================================================================

critical_d <- function(
    sar     = 1, r = 0.1,
    crit_DK = 0,
    ...) {  #arguments passed to par_d function
  
  # critical depletion for the sar
  # (from Pritchard, 2017): D = 1+ sar/r * log(1-d)
  
  # For D/K = crit_DK
  d <- 1 - exp((crit_DK - 1.) / sar * r)
  d[is.infinite(d)] <- NA
  d
}

# =============================================================================

critical_gpd <- function(
    sar     = 1, 
    fDepth  = 1, uDepth = 0, d = NULL,  
    r       = 0.1,
    crit_DK = 0,
    ...) {  #arguments passed to par_d function
  
  # critical depletion for the sar
  
  #  D = 1+ sar/r * log(1-d)
  
  if (length(r) > 1) {
    if (length(sar) > 1)
      stop("cannot solve for both r and sar a vector")
    if (is.vector(fDepth))
      fDepth <- matrix(nrow  = length(r), ncol = length(fDepth), 
                       byrow = TRUE, data = fDepth)
    else if (nrow(fDepth) != length(r))
      stop ("'fDepth'  and 'r' not compatible")
  }
    
  if (is.null(d)){
    # For DK = crit_DK
    d <- critical_d(sar = sar, r = r, crit_DK = crit_DK, ...)
    m_max <- 0.42
    d <- pmin(d, m_max)
  }
  
  zz <- which(is.na(d))
  if (length(zz)) 
    d[zz] <- 0.1
  fDepth[is.na(fDepth)] <- 1
  
  if (is.vector(fDepth)) 
    fDepth <- matrix(ncol = length(uDepth), nrow = length(d),
                     byrow = TRUE, data = fDepth)
  if (nrow(fDepth) != length(d))  
    fDepth <- matrix(ncol = length(uDepth), nrow = length(d),
                     byrow = TRUE, data = fDepth)
  froot <- function(gpd, dd){
    dd - par_d(gpd = gpd, fDepth = fDepth, uDepth = uDepth)
  }
  froot2 <- function(gpd, i){
    d[i] - par_d(gpd = gpd, fDepth = fDepth[i,], uDepth = uDepth)
  }
  if (length(d) == 1)
    res <- uniroot(f=froot, interval = c(0, 100), d = d)$root
  else
    res <- sapply(1:length(d), FUN = function(i) 
      uniroot(f=froot2, interval = c(0, 100), i = i)$root)
  
  if (length(zz)) res[zz] <- NA
  
  res
}
