\name{Bioturbation and bioirrigation indices from Model output}
\alias{get_Db_model}
\alias{get_irr_model}
\title{
  From model output (time x taxon density) to time x bioturbation or time x bioirrigation potentials.
}

\description{
 \code{get_Db_model} and \code{get_irr_model} combine output of a biological model that describes how taxa evolve with time with a trait database to obtain how bioturbation and bioirrigation potentials evolve over time.  They estimate the bioturbation potential index (BPc) sensu Querios et al., 2013 and Solan et al., 2004, and the bioirrigation potential index (IPc) sensu Wrede et al., 2018, respectively.
}

\usage{
get_Db_model(model, trait=Btrait::Traits_Db, taxonomy = NULL,  
             weight, verbose = FALSE, na.rm=FALSE)  

get_irr_model(model, trait=Btrait::Traits_irr, taxonomy = NULL,  
            weight, verbose = FALSE, na.rm=FALSE)  
}

\arguments{
  \item{model}{ dynamic model output that describes how taxa change over time. This is a matrix with (time x taxon) information, such as generated from the function \link{run_perturb}. The first column should contain the time, and the other columns the taxon density over time, with the column names the actual taxon names. (see note for the consequences of this.)
  }
  \item{trait}{ (taxon x trait) data. For \code{get_Db_model}, \code{traits} should contain the columns \code{Ri} and \code{Mi} (with the reworking and mobility traits). For use in \code{get_irr_model}, \code{traits} should contain columns \code{BT} and \code{FT}, and \code{ID}, the burrowing type, feeding type, and injection depth respectively. The trait matrix should have the names of the taxa in its first column. A good choice of a trait database for estimating bioturbation potentials is \link[Btrait]{Traits_Db}; a good choice for estimating bioirrigation potentials is \link[Btrait]{Traits_irr}.
  }
  \item{verbose }{when TRUE, will write warnings to the screen.
  }
  \item{na.rm }{when TRUE, will remove species for which trait could not be estimated (as not present in the trait database).
  }
  \item{taxonomy }{taxonomic information; first column will be matched with \code{taxon}, regardless of its name.
  }
  \item{weight }{a \code{data.frame} with taxon (first column) versus weight (2nd column) characteristics for each taxon in the \code{model}.
  }
}

\value{
\code{get_Db_model} returns the time x Db-index matrix of class \code{deSolve}. Both the bioturbation index for each taxon ("BPC_taxon.."), and the summed index is returned ("BPc").

\code{get_irr_model} returns the time x Irr-index density matrix of class \code{deSolve}. Both the irrigation index for each taxon ("IPC_taxon.."), and the summed index is returned ("IPc").

The result has in its attributes a vector called "Factor", which is the constant factor with which the model output has been multiplied; for \code{get_Db_model}, this is estimated as sqrt(weight)*Mi*Ri; for \code{get_irr_model}, this is estimated as (weight)^0.75*FT*BT*ID.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\note{
The names of the columns in \code{model} should have the taxon names, that are also found in the trait database. This can create a problem if the model output is a \code{data.frame} and the taxon names contain blancs. In this case, the names will have "." rather than " ". This can be solved either (1) by making the names in the trait database consistent with the model output (use R-function \code{make.names()} to convert names with spaces to good names, as used in data.frames), or (2) by adding, to model an attribute called taxon.names that has the actual taxon names (in the order of the columns); (attributes(model)$taxon.names  <- taxon.names).

\itemize{
 \item{get_Db_model}: The formula for estimating the bioturbation Index for taxon i (as in Querios et al., 2013) is:

BPc_i = sqrt(Wi) * density_i * Ri*Wi

 \item{get_irr_model}: The formula for estimating the bioirrigation Index for taxon i (as in Wrede et al., 2018) is:

IPc_i = (Wi)^(0.75) * density_i * BTi*FTi*IDi

 \item{The stations Index} is the sum of all species indices.
}

}

\seealso{
\link{run_perturb} for the disturbance model.

\link{get_trait_model} for deriving functional traits from a perturbation or logistic model.

\link[Btrait]{Traits_Db}, for trait databases in packge Btrait..

\link[Btrait]{get_trait_density}, for the function on which \code{get_trait_model} is based.

\link[Btrait]{get_trait} for functions from package Btrait to extract traits
}

\keyword{manip}

\examples{

## ====================================================
## A small model
## ====================================================

 Trawl.pars <- data.frame(
    taxon = c("sp.1","sp.2","sp.3","sp.4"),  # name of taxa
    r     = c(   3,      1,    20,    0.1),  # rate of increase (/year)
    K     = c(   1,     10,   0.5,    10 ),  # carrying capacity (density or biomass)
    d     = c( 0.1,    0.1,   0.5,   0.8))  # instantaneous depletion

# run the model for 5 years 
 times   <- seq(0, 5, by=1/365)             # time in years (consistent with pars)

# trawling frequency once / twice per year
 trawl1  <- seq(from=1/4, by=1,   to=5)

 trawl2  <- seq(from=1/4, by=0.5, to=5)

# run the model
 trawl1run <- run_perturb (parms  = Trawl.pars, 
                           times  = times, 
                           events = trawl1)     

 trawl2run <- run_perturb (parms  = Trawl.pars, 
                           times  = times, 
                           events = trawl2)     

 plot(trawl1run, trawl2run)
 
##-----------------------------------------------------
## Small dataset: taxonomy
##-----------------------------------------------------

Btaxonomy <- data.frame(
  species = c("sp.1","sp.2","sp.3","sp.4","sp.5","sp.6"),
  genus   = c( "g.1", "g.2", "g.2", "g.2", "g.3", "g.4"),
  family  = c( "f.1", "f.1", "f.1", "f.1", "f.2", "f.3"),
  order   = c( "o.1", "o.1", "o.1", "o.1", "o.2", "o.2"),
  class   = c( "c.1", "c.1", "c.1", "c.1", "c.1", "c.1")
  )

##-----------------------------------------------------
## Db and Irr traits
##-----------------------------------------------------

# Note: no data for "sp.4"
DbTraits <- data.frame(
  taxon   = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  Ri      = c(1     , 1    ,   3  ,   2  ,     4),
  Mi      = c(1     , 2    ,   2  ,   3  ,     3)
)

IrrTraits <- data.frame(
  taxon   = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  BT      = c(1     , 1    ,   3  ,   2  ,     3),
  FT      = c(3     , 2.5  ,   2  ,   3  ,     2),
  ID      = c(1     , 2    ,   1  ,   3  ,     2)
)

# Weight of the species - it is assumed this does not change over time
Weight <- data.frame(taxon=c("sp.1","sp.2","sp.3","sp.4","sp.5","sp.6"),
                     weight =c(0.1,   0.2,     3,  0.04,   0.5,     6))
                     
DbRun  <- get_Db_model (model       = trawl1run,  
                        trait       = DbTraits, 
                        weight      = Weight)
                      
head(DbRun)   # contains NA for sp.4 that is not in DbTraits
attributes(DbRun)$notrait

# use taxonomy to also estimate Db for sp.4
DbRun1 <- get_Db_model (model       = trawl1run,  
                        trait       = DbTraits, 
                        weight      = Weight,
                        taxonomy    = Btaxonomy)
head(DbRun1)  # sp.4 estimated based on taxonomic closeness
 
# Show results - note: no sp4 for DbRun, so only one line
plot(DbRun1, DbRun)

# Irrigation potential index
IrrRun1 <- get_irr_model (model       = trawl1run,  
                          trait       = IrrTraits, 
                          weight      = Weight,
                          taxonomy    = Btaxonomy)

plot(IrrRun1)

}
\references{
Queiros, Ana M., Silvana N. R. Birchenough, Julie Bremner, Jasmin A. Godbold, Ruth E. Parker, Alicia Romero-Ramirez, Henning Reiss, Martin Solan, Paul J. Somerfield, Carl Van Colen, Gert Van Hoey, Stephen Widdicombe, 2013. A bioturbation classification of European marine infaunal invertebrates. Ecology and Evolution 3 (11), 3958-3985

Solan M, Cardinale BJ, Downing AL, Engelhardt KAM, Ruesink JL, Srivastava DS. 2004. Extinction and ecosystem function in the marine benthos. Science 306:1177-80.

A. Wrede, J.Beermann, J.Dannheim, L.Gutow, T.Brey, 2018. Organism functional traits and ecosystem supporting services - A novel approach to predict bioirrigation. Ecological indicators, 91, 737-743.
}
