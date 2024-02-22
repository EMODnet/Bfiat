\name{Datasets for trawling assessment}
\alias{NIOZdata}
\title{
  Biological data for running the BFIAT tool.
}

\description{
 \code{NIOZdata} provides macrofauna data at 4 sites for running the BFIAT model tools. 
}

\format{
 \itemize{
 \item{**\code{NIOZdata}**} is a list with 5 tables .  
  \itemize{
   \item{density} a \code{data.frame} with the station names, taxon, density, biomass and date.  
   \item{stations} a \code{data.frame} with the station names, x and y positions, depth and the number of samples on which the data were based.
   \item{taxonomy} a \code{data.frame} with column names \code{taxon}, \code{genus}, \code{family}, \code{order}, \code{class}, \code{phylum} and \code{AphiaID}.  
   \item{fishing} a \code{data.frame} with the species traits necessary for estimating the fishing parameters; columns are \code{taxon, p0, p0_5, p5_15, p15_30, p30, swim, Life.span, r}.  
   \item{contours} a \code{list} with the x, y, z values of the contours (for use with \link[Btrait]{mapKey} or  \link[Btrait]{mapLegend}).  
   \item{abiotics} a \code{data.frame} with abiotic conditions (depth, D50, silt percentage).  
   \item{sar} a \code{data.frame} with fishing intensities for the stations, per year; columns are \code{station, year, sar, subsar, gpd}.  
  }
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}

\examples{

##-----------------------------------------------------
## Show contents of the NIOZ data sets
##-----------------------------------------------------

metadata(NIOZdata$density)
head(NIOZdata$density)

metadata(NIOZdata$fishing)
head(NIOZdata$fishing)

metadata(NIOZdata$sar)
NIOZdata$sar
}

\references{

Dataset from the Vlakte van de Raan:

J. Tiano, J. Depestele, G. Van Hoey, J. Fernandes, P. van Rijswijk, and K. Soetaert, 2022. Trawling effects on biogeochemical processes are mediated by fauna in high energy biogenic reef-inhabited coastal sediments. Biogeosciences, 19, 2583-2598, 2022

Dataset from the Frisian Front:

Tiano, J. C., R. Witbaard, M.J.N. Bergman, P. van Rijswijk, A. Tramper, D. van Oevelen, K. Soetaert, 2019. Acute impacts of bottom trawl gears on benthic metabolism and nutrient cycling. ICES journal of Marine Science, fsz060, https://doi.org/10.1093/icesjms/fsz060. 

Tiano, J., van der Reijden, K, O'Flynn, S, Beauchard, O, van der Ree, S,  van der Wees, J, Ysebaert, T, Soetaert, K., 2020. Experimental assessment of bottom trawling finds significant effects on epifauna and juvenile infauna. Marine Environmental Research, 159, 11 p., 104964. https://doi.org/10.1016/j.marenvres.2020.104964

Dataset from the Dogger Bank and Fladen Ground:

Emil De Borger, Ulrike Braeckman, Karline Soetaert, 2021. Rapid organic matter cycling in North Sea sediments. Continental shelf research, vol. 214, 2021, doi:10.1016/j.csr.2020.104327.
}

\seealso{
\link[Btrait]{Traits_nioz} for trait data from package Btrait.

\link[Btrait]{MWTL} for density data from package Btrait.
}

\details{

 Data from four North Sea stations was used in \code{NIOZdata}. 
 \itemize{

  \item {Vlakte_van_de_Raan} The sandy, high metabolism site uses averaged data from 59 boxcores from the Vlakte van de Raan (Tiano et al., 2022). Faunal samples for this site were processed at the Flanders Research Institute for Agriculture, Fisheries and Food (ILVO). 
  \item{FrieseFront} The high metabolism, muddy site, uses information from 25 boxcore samples in the Frisian Front (Tiano et al., 2019, 2020) with faunal samples processed at the Netherlands Institute for Sea Research (NIOZ) and Fugro. 
  \item{Doggerbank}Data for the low metabolism sandy site were collected from three boxcore samples from the Dogger Bank (De Borger et al., 2021) and were processed at the NIOZ. 
  \item{FladenGround} Data for the low metablism muddy sites were collected from three boxcore samples from the Fladen Grounds (De Borger et al., 2021) and were processed at the NIOZ. 
 }
 
}
\keyword{datasets}
