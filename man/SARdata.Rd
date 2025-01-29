\name{Fishing intensity data}
\alias{SAR}
\title{
  Swept-area ratios in the Northsea, OSPAR data.
}
\description{
 \code{SAR} has average swept area ratios (2009-2020) for the Northsea (long [-2, 8]; lat [50, 60]), with resolution of 0.05 dg, (source: OSPAR).
}

\format{

 \itemize{
 \item{**\code{SAR}**} is a list with 5 items.  
  \itemize{
   \item{lon, x}, a \code{vector} (201 elements) with the longitudes (WGS84), dgE.  
   \item{lat, y}, a \code{vector} (201 elements) with the latitudes (WGS84), dgN.  
   \item{sar}, a \code{matrix} of dimensions (201 x 201) with annual swept area ratios for the surface (0-2cm).  
   \item{subsar}, a \code{matrix} of dimensions (201 x 201) with annual swept area ratios for the subsurface (>2cm).  
   \item{gpd}, a \code{matrix} of dimensions (201 x 201) with estimated gear penetration depths ([cm]).  
  }
  Data source: OSPAR (ICES, 2022)
 }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}

\examples{

##-----------------------------------------------------
## Show contents of the fishing data set
##-----------------------------------------------------

metadata(SAR)
head(SAR$sar)

##-----------------------------------------------------
## The fishing data - asp sets the good aspect ratio
##-----------------------------------------------------

# Aspect ratio of the graph
asp <- 1/cos((mean(SAR$lat) * pi)/180)

par(mfrow=c(2,2))
with (SAR, image2D(x=lon, y=lat, z=sar, 
      xlab = expression(""^o~E), 
      ylab = expression(""^o~N),
      clab = "#/y", 
      main = "Swept Area Ratio (GPD: 0-2cm)", 
      log="c", clim=c(0.0001,50), 
      colkey=list(length=0.5, width=0.5),
      asp=asp, las=1)
)      
  
with (SAR, image2D(x=lon, y=lat, z=subsar, 
      xlab = expression(""^o~E), 
      ylab = expression(""^o~N), 
      clab = "#/y",     
      main = "Swept Area Ratio (GPD: >2cm)",  
      log="c", clim=c(0.0001,50), 
      colkey=list(length=0.5, width=0.5),  
      asp=asp, las=1) 

)
  
with (SAR, image2D(x=lon, y=lat, z=gpd, 
      xlab = expression(""^o~E), 
      ylab = expression(""^o~N),
      clab = "cm",
      main = "Gear penetration depth", 
      colkey=list(length=0.5, width=0.5),
      asp=asp, las=1)
)      
  
hist (log10(SAR$sar), n=100, axes=FALSE, 
      xlab = "swept area ratio, /yr", 
      main = "histogram")
axis(side=1, at=-4:2, labels=10^(-4:2))
axis(side=2)

}


\seealso{
\link[Btrait]{Traits_nioz} for trait data from package Btrait.

\link[Btrait]{MWTL} for density data from package Btrait.
}

\details{

 Original data in \code{SAR} were compiled by ICES upon request by OSPAR. 
 
 These data are based on relevant VMS and logbook data. 
 Several Quality Checks were undertaken on the full VMS dataset to produce an overview QC report. 
 All R scripts and SQL code used to access and process the VMS data are available on GitHub 
 (https://github.com/ices-eg/wg_WGSFD).

 The full data set comprised sar, and subsar data from 2009 till 2020, and from
 [-43.675, 43.7250 dg E] and [36.02, 81.72 dg N]. 
 
 From these data, yearly averages were estimated, and the range restricted to be within 
 [-2, 8 dgE] & [50, 60 dg N].
 
 Mean gear penetration depth was estimated based on the metiers:

The metiers are Aggregated into beam trawl (TBB), dredge (DRB), demersal seine (SN), and otter trawl (OT), based on the 
metier layers: OT_CRU, OT_DMF, OT_MIX, OT_MIX_CRU, OT_MIX_DMF_BEN, OT_MIX_DMF_PEL, OT_MIX_CRU_DMF,
OT_SPF, TBB_CRU, TBB_DMF, TBB_MOL, DRB_MOL, SDN_DMF, SSC_DMF.

The gear penetration depths used were the mean values over sand and muddy sediments: 

\itemize{
 \item in sand: 3.5, 1.1, 1.1, 1.9 cm for DRB, OT, SN, TBB respectively
 \item in mud : 5.4, 2.0, 2.0, 3.2 cm for DRB, OT, SN, TBB respectively
 \item average (used here): 4.45, 1.55, 1.55, 2.55 cm for DRB, OT, SN, TBB respectively
}


}
\keyword{datasets}

\references{
(older version:)

ICES Technical Service, Greater North Sea and Celtic Seas Ecoregions, 29 August 2018
sr.2018.14 Version 2: 22 January 2019 https://doi.org/10.17895/ices.pub.4508
OSPAR request on the production of spatial data layers of fishing intensity/pressure.
}