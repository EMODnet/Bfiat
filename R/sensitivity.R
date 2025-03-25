
# =============================================================================
# Sensititivity of taxa to fishing/or areas
# =============================================================================

sensitivity_taxon <- function(
             fDepth = 1, uDepth = 0, r = 0.1,
             sar.seq = seq(0,    10, length.out = 200),
             gpd.seq = seq(0,     5, length.out = 200),
             ...) {  #arguments passed to par_d function
  
  # depletion for the gpd sequence
  if (sum(fDepth) != 1) stop (" sum of fDepth should be = 1")
  d.seq <- par_d(gpd = gpd.seq, fDepth = fDepth, uDepth = uDepth, ...)
  
  DK     <- outer(sar.seq, d.seq, 
                  FUN = function(s, d) 1/par_K(density = 1, sar = s, r = r, d = d))
  
#  rM_sar <- apply(DKK, MARGIN = 1, FUN = function(x) gpd.seq[which.min(x)])
  rM_gpd <- critical_gpd(sar     = sar.seq, 
                         fDepth  = fDepth, uDepth = uDepth, r = r, 
                         crit_DK = 0)
  rM_gpd[rM_gpd <= min(gpd.seq) ] <- NA
  rM_gpd[rM_gpd >= max(gpd.seq) ] <- NA
  
  #  rM_gpd <- apply(DKK, MARGIN = 2, FUN = function(x) sar.seq[which.min(x)])
  rM_sar <- critical_sar(gpd     = gpd.seq, 
                         d       = d.seq, r = r, 
                         crit_DK = 0)
  rM_sar[rM_sar <= min(sar.seq) ] <- NA
  rM_sar[rM_sar >= max(sar.seq) ] <- NA
  
  DK[is.infinite(DK)] <- NA
  
  RES <- list(sar = sar.seq, gpd = gpd.seq, DK = DK, 
              critical_sar = rM_sar, critical_gpd = rM_gpd, 
              r = r, Depth_mean = sum(fDepth*uDepth)/sum(fDepth), 
              d_mean = mean(d.seq))
  
  description <- data.frame(
    name = c("sar", "gpd", "DK", 
             "critical_sar", "critical_gpd", 
             "r", "Depth_mean", "d_mean"), 
    description = c("sequence of swept area ratios (fisheries)",
                    "sequence of gear penetration depths",
                    "density/carrying capacity, corresponding to sar (row) and gpd (column)",
                    "critical value of sar for each gpd, above which the taxon is extinct",
                    "critical value of gpd for each sar, above which the taxon is extinct",
                    "the intrinsic rate of natural increase of the taxon",
                    "the mean living depth of the taxon (sum(fDepth*uDepth))",
                    "the mean depletion fraction (for each gpd)"),
    units = c("/year", "cm", "-", "/yr", "cm", "/yr", "cm", "-"))
  attributes(RES)$description <- description
  class(RES) <- c("sensitivity_taxon", class(RES))
  RES
}


####

sensitivity_area <- function(sar = 1, gpd = 1,
    r.seq = seq(1e-1, 5, length.out = 200),
    d.seq = seq(0,   0.5, length.out = 200)
    ) { 
  
  # depletion for the gpd sequence
  DK     <- outer(r.seq, d.seq, 
                  FUN = function(r, d) 
                    1/par_K(density = 1, sar = sar, r = r, d = d))
  
#  rM_d <- apply(DKK, MARGIN = 1, FUN = function(x) d.seq[which.min(x)])
  rM_d <- critical_d(sar     = sar, 
                     r       = r.seq, 
                     crit_DK = 0)
  rM_d[rM_d <= min(d.seq) ] <- NA
  rM_d[rM_d >= max(d.seq) ] <- NA
  
#  rM_r <- apply(DKK, MARGIN = 2, FUN = function(x) r.seq[which.min(x)])
  rM_r <- critical_r(sar     = sar, 
                     d       = d.seq, 
                     crit_DK = 0)
  rM_r[rM_r <= min(r.seq) ] <- NA
  rM_r[rM_r >= max(r.seq) ] <- NA
  
  DK[is.infinite(DK)] <- NA
  
  RES <- list(r = r.seq, d = d.seq, DK = DK, 
              critical_r = rM_r, critical_d = rM_d, 
              sar = sar, gpd = gpd)
  
  description <- data.frame(
    name = c("r", "d", "DK", 
             "critical_r", "critical_d", 
              "sar", "gpd"), 
    description = c("sequence of r, the intrinsic rate of natural increase of the taxon",
                    "sequence of d, the depletion fraction",
                    "density/carrying capacity, corresponding to r (row) and d (column)",
                    "critical value of r for each d, below which the taxon is extinct",
                    "critical value of d for each r, above which the taxon is extinct",
                    "swept area ratio of the area",
                    "gear penetration depth of the area"),
    units = c("/year", "-", "-", "/yr", "-", "/yr", "cm"))
  attributes(RES)$description <- description
  class(RES) <- c("sensitivity_area", class(RES))
  RES
}

