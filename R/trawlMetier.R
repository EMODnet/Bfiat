

## ====================================================================
## ====================================================================
## density after and between trawling events - more than one metier
## ====================================================================
## ====================================================================

trawl_metier <- function(K   = 1, r = 1, 
                         parms = data.frame(K = K, r = r),
                         d = data.frame(0.1, 0.1), 
                         sar = data.frame(1, 2),
                         D0  = parms[["K"]], 
                         n   = 1){ # n=number of fishing events
  
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
  
  # check for NA in inputs
  wisna          <- apply(parms, 
                          MARGIN = 1, 
                          FUN    = function(x) any(is.na(x)))
  parms[wisna, ] <- 1
  
  n[which(is.na(n))] <- 0
  
  if (any(abs(trunc(n) - n) > 1e-8) )
      stop ("'n' should be an integer")
  
  ii <- which(n <= 0)
  if (length(ii)) n[ii] <- 1 
 
  if (length(n) > 1) n <- unique(n[order(n)])
  nevent <- length(n)
  
  DD <- .Fortran("metier_event", nspec = nspec, nmetier = nmetier ,
                     nevent = nevent, 
                     B0  = as.double(parms$D0),
                     sar = as.double(sar),
                     K   = as.double(parms$K), 
                     r   = as.double(parms$r), 
                     d   = as.double(d), 
                     eventnr = as.integer(n), 
                     B     = matrix(nrow = nspec, ncol = nevent, 
                                    data = as.double(-999.)),
                     Bend  = matrix(nrow = nspec, ncol = nevent, 
                                    data = as.double(-999.)),
                     Bmean = matrix(nrow = nspec, ncol = nevent, 
                                    data = as.double(-999.)),
                                                 
                     times = matrix(nrow = nspec, ncol = nevent, 
                                    data = as.double(-999.)),
                     tend  = matrix(nrow = nspec, ncol = nevent, 
                                    data = as.double(-999.))
  )
  if (length(ii)){
    DD$B[, ii]     <- parms$D0
    DD$Bend[, ii]  <- parms$D0
    DD$Bmean[, ii] <- parms$D0
    DD$times[, ii] <- NA
  }
  if (sum(wisna)){
    DD$B    [wisna, ] <- NA
    DD$Bend [wisna, ] <- NA
    DD$Bmean[wisna, ] <- NA
    DD$times[wisna, ] <- NA
  }
  
  # from wide to long format.
  nc <- ncol(parms) + 1 
    
  X1 <- data.frame(parms, DD$Bmean)
  X1 <- reshape(X1, direction = "long", 
                varying = list(nc:ncol(X1)), 
                v.names = "density")  # the mean
    
  X2 <- data.frame(parms, DD$B)
  X2 <- reshape(X2, direction = "long", 
                varying = list(nc:ncol(X2)), 
                v.names = "density_start")
  XX <- merge(X1, X2)
    
  X3 <- data.frame(parms, DD$Bend)
  X3 <- reshape(X3, direction = "long", 
                varying = list(nc:ncol(X3)), 
                v.names = "density_end")
  XX <- merge(XX, X3)
    
  X4 <- data.frame(parms, DD$times)
  X4 <- reshape(X4, direction = "long", 
                varying = list(nc:ncol(X4)), 
                v.names = "times")
  XX <- merge(XX, X4)
  
  X5 <- data.frame(parms, DD$tend)
  X5 <- reshape(X5, direction = "long", 
                varying = list(nc:ncol(X5)), 
                v.names = "times_end")
  XX <- merge(XX, X5)
  XX <- XX[order(XX$id, XX$time),]
  XX$ntrawl_from   <- n
  XX$ntrawl_to     <- n+1
  row.names(XX) <- NULL
  
  return(XX[,c(names(parms), "times", "times_end", "ntrawl_from", "ntrawl_to",
                 "density", "density_start", "density_end")])                   
    
    }


