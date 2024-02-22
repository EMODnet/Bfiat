\name{Fishing intensity data}
\alias{SAR}
\title{
  Swept-area ratios in the Northsea, OSPAR data.
}
\description{
 \code{SAR} has average swept area ratios (2009-2020) for the Northsea (long [-1.025, 8.025]; lat [50.025 59.025]), with resolution of 0.05 dg, (source: OSPAR).
}

\format{

 \itemize{
 \item{**\code{SAR}**} is a list with 5 items.  
  \itemize{
   \item{lon, x}, a \code{vector} (182 elements) with the longitudes (WGS84), dgE.  
   \item{lat, y}, a \code{vector} (181 elements) with the latitudes (WGS84), dgN.  
   \item{sar}, a \code{matrix} of dimensions (182 x 181) with annual swept area ratios for the surface (0-2cm).  
   \item{subsar}, a \code{matrix} of dimensions (182 x 181) with annual swept area ratios for the subsurface (>2cm).  
   \item{gpd}, a \code{matrix} of dimensions (182 x 181) with estimated gear penetration depths ([cm]).  
  }
  Data source: OSPAR (ICES, 2018)
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
      log="c", clim=c(0.001,100), 
      colkey=list(length=0.5, width=0.5),
      asp=asp, las=1)
)      
  
with (SAR, image2D(x=lon, y=lat, z=subsar, 
      xlab = expression(""^o~E), 
      ylab = expression(""^o~N), 
      clab = "#/y",     
      main = "Swept Area Ratio (GPD: >2cm)",  
      log="c", clim=c(0.001,100), 
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
 [-1.025, 8.025 dgE] & [50.025, 59.025 dg N].
 
 Mean gear penetration depth was estimated as:
 (1.*(sar-subsar) + 3*subsar)/sar
 
 where we assumed a mean penetration of 1 cm for the (0-2 cm) sar values, and 
 a mean penetration of 3 cm for the (>3 cm) subsar values.
}
\keyword{datasets}

\references{
ICES Technical Service, Greater North Sea and Celtic Seas Ecoregions, 29 August 2018
sr.2018.14 Version 2: 22 January 2019 https://doi.org/10.17895/ices.pub.4508
OSPAR request on the production of spatial data layers of fishing intensity/pressure.
}