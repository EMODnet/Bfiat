## ====================================================================
## ====================================================================
## Merging model output with traits
## ====================================================================
## ====================================================================

getTraitModel <- function(model,         
                          trait,
                          trait.class = NULL, trait.score = NULL, 
                          taxonomy = NULL,  
                          t.column = 1, 
                          scalewithvalue = TRUE, 
                          verbose = FALSE){
  traitModel <- getTraitDensity(
            wide        = model, 
            trait       = trait, 
            trait.class = trait.class, 
            trait.score = trait.score, 
            taxonomy    = taxonomy, 
            t.column    = t.column, 
            scalewithvalue = scalewithvalue, 
            verbose     = verbose) 
  
  Atts <- attributes(traitModel)
  traitModel <- as.matrix(traitModel)
  attributes(traitModel)$notrait <- Atts$notrait
  class(traitModel) <- c("deSolve", "matrix")
  return(traitModel)
}

## ====================================================================
## Estimating bioturbation/bioirrigation from model output
## ====================================================================

getDbModel <- function(model,         
                       trait=Btrait::Traits_Db,
                       taxonomy = NULL,  
                       weight, # two-columned data.frame (taxon, weight)
                       verbose = FALSE, na.rm=FALSE)  
getIndexModel(model=model, trait=trait, taxonomy=taxonomy, 
              weight=weight, verbose=verbose, na.rm=na.rm, type="BPc")

## ====================================================================

getIrrModel <- function(model,         
                       trait=Btrait::Traits_irr,
                       taxonomy = NULL,  
                       weight, # two-columned data.frame (taxon, weight)
                       verbose = FALSE, na.rm=FALSE)  
getIndexModel(model=model, trait=trait, taxonomy=taxonomy, 
              weight=weight, verbose=verbose, na.rm=na.rm, type="IPc")

## ====================================================================

getIndexModel <- function(model,         
                       trait,
                       taxonomy = NULL,  
                       weight, # two-columned data.frame (taxon, weight)
                       verbose = FALSE, na.rm=FALSE,
                       type = "BPc") { 

# Get the traits for the taxa - use taxonomic closeness to extend traits
  taxcname <- colnames(trait)[1]  # name of column with taxa in trait
  taxa     <- attributes(model)$taxon.names                                   
  if (is.null(taxa))                                                     
    taxa     <- colnames(model)[-1] # taxa in the model run
  if (type == "BPc") 
    traitnames <- c(taxcname, "Ri", "Mi")
  else if (type == "IPc")
    traitnames <- c(taxcname, "BT", "FT", "ID")
  S <- try(Trait <- trait[,traitnames], silent = TRUE)
  if (inherits(S, "try-error"))
    stop ("Trait does not have all columns, should contain", paste(traitnames, collapse = ", "))
  
  Data_Db  <- getTrait(
        taxon    = taxa, 
        trait    = Trait, 
        taxonomy = taxonomy,
        verbose  = verbose)

  # taxon for which trait could not be derived
  notrait <- attributes(Data_Db)$notrait 
  
  Data_Db <- cbind(Data_Db, 1:nrow(Data_Db))  # add order
  nc      <- ncol(Data_Db)
  
  if (is.null(ncol(weight)))
    stop ("weight should be a matrix or dataframe with taxon in first column and weight in second column") 
  colnames(weight)[2] <- "Weight"
  Data_Db <- merge(Data_Db, weight, by=1, all.x=TRUE)
  Noweight <- NULL
  if (any(is.na(Data_Db$Weight))){
    if (verbose) warning("Weight not known for ", length(which(is.na(Data_Db$Weight)))," taxa - set to 0")
    Data_Db$Weight[is.na(Data_Db$Weight)] <- 0
    Noweight <- Data_Db$taxon[is.na(Data_Db$Weight)]
   }
  Data_Db <- Data_Db[order(Data_Db[,nc]), ]   # restore order

  if (type == "BPc") 
    Data_FAC <- with(Data_Db, sqrt(Weight)*Mi*Ri)  # THIS WILL REMAIN CONSTANT
  else if (type == "IPc")
    Data_FAC <- with(Data_Db, Weight^(0.75)*BT*FT*ID)  
  
  # multiply all rows with fac
  model_BPC <- sweep(model[,-1], MARGIN=2, STATS=Data_FAC, FUN="*")
  colnames(model_BPC) <- paste(type, taxa, sep="_")
  if (na.rm) {
    isna <- apply(model_BPC, MARGIN=2, FUN=function(x)!any(is.na(x)))
    model_BPC <- model_BPC[,isna]
  }
  if (is.matrix(model_BPC))
    model_BPC <- cbind(model_BPC, rowSums(model_BPC, na.rm=TRUE))
  else
    model_BPC <- cbind(model_BPC, model_BPC)
  colnames(model_BPC)[ncol(model_BPC)] <- type
  model_BPC <- cbind(model[, 1], model_BPC)
  model_BPC <- as.matrix(model_BPC)
  colnames(model_BPC)[1] <- "times"
  class(model_BPC) <- c("deSolve", "matrix")
  attributes(model_BPC)$Factor <- Data_FAC 
  attributes(model_BPC)$Noweight <- Noweight 
  attributes(model_BPC)$notrait <- notrait
  return(model_BPC)
}
