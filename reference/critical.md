# Critical values for fishing or taxon parameters.

`critical_sar` and `critical_gpd` estimates for one or more taxa the
fishing intensity (`sar`) or gear penetration depth (`gpd`), above which
the density/caryying capacity ratio decreases below a critical value
(`crit_DK`).

`critical_r` and `critical_d` estimate for a combination of fishing
intensity (`sar`) and gear penetration depth (`gpd`), the critical
intrinsic rate of natural increase (`r`) or depletion fraction (d) at
which the density, carrying capacity ratio reaches the critical value
(`crit_DK`).

## Usage

``` r
critical_sar (gpd = 1, fDepth = 1, uDepth = 0, d = NULL,
              r = 0.1, crit_DK = 0, ...) 
  
critical_gpd ( sar = 1, fDepth = 1, uDepth = 0, d = NULL,
               r = 0.1, crit_DK = 0, ...) 

critical_r (sar = 1, gpd = 1, fDepth = 1, uDepth = 0, d = NULL, 
            crit_DK = 0, ...) 

critical_d (sar = 1, r = 0.1, crit_DK = 0, ...)
```

## Arguments

- sar :

  fishing intensity, estimated as Swept Area Ratio, units e.g.
  \[m2/m2/year\]. One number or a vector.

- gpd :

  gear penetration depth, units e.g. \[cm\]. One number or a vector.

- r :

  the intrinsic rate of natural increase of a taxon, units e.g.
  \[/year\]. One number or a vector.

- d :

  depletion fraction due to fishing, one value or a vector.

- crit_DK :

  the steady-state density over carrying capacity ratio of the taxon
  that need to be matched. The default value of 0 indicates that the
  taxon is extinct. A value of 0.01 indicates that the actual density is
  1/100 of its potential density (carrying capacity).

- fDepth :

  fractional occurrence of species in sediment layers, dimensionless. A
  vector of the same length as `uDepth`. The sum of `fDepth` should
  equal 1. Will be used to estimate the depletion fraction `d`. Will be
  ignored if `d` is given a value.

- uDepth :

  depth of the upper position of the sediment layers, units e.g. \[cm\].
  A vector with length equal to the number of columns of `fDepth`. Will
  be ignored if `d` is given a value.

- ... :

  arguments passed to the `par_d` function

## Value

These functions will either return one number, a vector or a matrix,
depending on the input arguments

## Author

Karline Soetaert \<karline.soetaert@nioz.nl\>

## See also

[run_perturb](https://EMODnet.github.io/Bfiat/reference/trawlingfun.md)
for how to run a disturbance model.

[Traits_nioz](https://EMODnet.github.io/Btrait/reference/Btraitdata.html),
for trait databases in package Btrait.

[MWTL](https://EMODnet.github.io/Btrait/reference/MWTLdensitydata.html)
for data sets on which fishing can be imposed.

[map_key](https://EMODnet.github.io/Btrait/reference/BtraitMap.html) for
simple plotting functions.

## References

to be added

## Examples

``` r

## -----------------------------------------------------------------------
## Critical sar for beam trawling for all species from the MWTL dataset
## -----------------------------------------------------------------------

gpd_mud   <- 3.2 # beam trawl in muddy sediment
gpd_sand  <- 1.9 #               sandy

head(MWTL$fishing)
#>                   taxon  p0 p0_5cm p5_15cm p15_30cm p30cm Age.at.maturity     r
#> 1 Abludomelita obtusata 0.5    0.5     0.0        0     0        0.500000 5.120
#> 2             Abra alba 0.0    0.5     0.5        0     0        0.500000 5.120
#> 3           Abra nitida 0.0    1.0     0.0        0     0        0.500000 5.120
#> 4       Abra prismatica 0.0    1.0     0.0        0     0        2.000000 1.280
#> 5           Abra tenuis 0.0    1.0     0.0        0     0        2.000000 1.280
#> 6 Abyssoninoe hibernica 0.0    0.5     0.5        0     0        3.333333 0.768
fDepth <- MWTL$fishing[,c("p0", "p0_5cm", "p5_15cm", "p15_30cm", "p30cm")]
uDepth <-               c( 0,    0,         5,        15,          30)

csar_sand <- critical_sar(gpd    = gpd_sand, 
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)

csar_mud  <- critical_sar(gpd    = gpd_mud, 
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)

csar <- data.frame(taxon     = MWTL$fishing$taxon, 
                   csar_sand = csar_sand, 
                   csar_mud  = csar_mud)
                   
summary(csar)
#>        taxon       csar_sand          csar_mud      
#>  Length   :400   Min.   :  5.799   Min.   :  3.306  
#>  N.unique :400   1st Qu.: 13.515   1st Qu.:  7.827  
#>  N.blank  :  0   Median : 23.852   Median : 13.896  
#>  Min.nchar:  5   Mean   : 37.407   Mean   : 21.606  
#>  Max.nchar: 34   3rd Qu.: 46.388   3rd Qu.: 26.448  
#>                  Max.   :438.392   Max.   :259.250  
#>                  NAs    :26        NAs    :26       


## -----------------------------------------------------------------------
## Critical gear penetration depth for species from the MWTL dataset
## -----------------------------------------------------------------------

fDepth <- MWTL$fishing[,c("p0", "p0_5cm", "p5_15cm", "p15_30cm", "p30cm")]
uDepth <-               c( 0,    0,         5,        15,          30)

cgpd_01   <- critical_gpd(sar   = 1,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
cgpd_10   <- critical_gpd(sar   = 10,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
cgpd_100  <- critical_gpd(sar   = 100,
                          fDepth = fDepth, uDepth = uDepth, 
                          r      = MWTL$fishing$r)
# sar of 100 and critical density = 10 percent
cgpd_100b <- critical_gpd(sar     = 100,
                          fDepth  = fDepth, uDepth = uDepth, 
                          r       = MWTL$fishing$r,
                          crit_DK = 0.1)

cgpd <- data.frame(taxon     = MWTL$fishing$taxon, 
                   cgpd_01   = cgpd_01, 
                   cgpd_10   = cgpd_10,
                   cgpd_100  = cgpd_100,
                   cgpd_100b = cgpd_100b)

# a cgpd of 100 means: no limits                   
summary(cgpd)
#>        taxon        cgpd_01       cgpd_10          cgpd_100      
#>  Length   :400   Min.   :100   Min.   : 1.127   Min.   : 0.1160  
#>  N.unique :400   1st Qu.:100   1st Qu.: 2.549   1st Qu.: 0.2698  
#>  N.blank  :  0   Median :100   Median : 4.778   Median : 0.5319  
#>  Min.nchar:  5   Mean   :100   Mean   : 5.465   Mean   : 1.0229  
#>  Max.nchar: 34   3rd Qu.:100   3rd Qu.: 7.286   3rd Qu.: 0.9075  
#>                  Max.   :100   Max.   :20.095   Max.   :15.3686  
#>                  NAs    :8     NAs    :8        NAs    :8        
#>    cgpd_100b      
#>  Min.   : 0.1044  
#>  1st Qu.: 0.2430  
#>  Median : 0.4792  
#>  Mean   : 0.9566  
#>  3rd Qu.: 0.8188  
#>  Max.   :15.3321  
#>  NAs    :8        

```
