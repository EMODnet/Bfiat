## ====================================================================
## ====================================================================
## Merging model output with traits
## ====================================================================
## ====================================================================

get_trait_model <- function(model,         
                            trait,
                            trait_class = NULL, trait_score = NULL, 
                            taxon_names = colnames(model)[-1],
                            taxonomy = NULL,  
                            taxon_column = 1, 
                            scalewithvalue = TRUE, 
                            verbose = FALSE){
  wide <- model
  colnames(wide)[-1] <- taxon_names 
  traitModel <- get_trait_density(
            wide        = model, 
            trait       = trait, 
            trait_class = trait_class, 
            trait_score = trait_score, 
            taxonomy    = taxonomy, 
            taxon_column    = taxon_column, 
            scalewithvalue = scalewithvalue, 
            verbose     = verbose) 
  
  Atts              <- attributes(traitModel)
  traitModel        <- as.matrix(traitModel)
  attributes(traitModel)$notrait <- Atts$notrait
  class(traitModel) <- c("deSolve", "matrix")
  return(traitModel)
}

## ====================================================================
## Estimating bioturbation/bioirrigation from model output
## ====================================================================

get_Db_model <- function(model,         
                       trait=Btrait::Traits_Db,
                       taxon_names = colnames(model)[-1],
                       taxonomy = NULL,  
                       weight, # two-column data.frame (taxon, weight) or vector
                       
                       verbose = FALSE, na.rm=FALSE)  
getIndexModel(model=model, trait=trait, taxonomy=taxonomy, 
              weight=weight, taxon_names=taxon_names, 
              verbose=verbose, na.rm=na.rm, type="BPc")

## ====================================================================

get_irr_model <- function(model,         
                       trait=Btrait::Traits_irr,
                       taxon_names = colnames(model)[-1],
                       taxonomy = NULL,  
                       weight, # two-column data.frame (taxon, weight) or vector
                       verbose = FALSE, na.rm=FALSE)  
getIndexModel(model=model, trait=trait, taxonomy=taxonomy, 
              weight=weight, taxon_names=taxon_names, 
              verbose=verbose, na.rm=na.rm, type="IPc")

## ====================================================================

getIndexModel <- function(model,         
                       trait,
                       taxonomy = NULL,  
                       weight, # two-columned data.frame (taxon, weight)
                       taxon_names = colnames(model)[-1],
                       verbose = FALSE, na.rm=FALSE,
                       type = "BPc") { 

# number of trajectories in model (first column is time)  
  nTrajectory <- ncol(model) - 1

# taxon names, can also be passed as attribute to model
  
  taxcname <- colnames(trait)[1]  # name of column with taxa in trait
  taxa     <- taxon_names                                   
  if (is.null(taxa))                                                     
    taxa     <- colnames(model)[-1] # taxa in the model run
  if (length(taxa) != nTrajectory)
    stop ("the length of taxon_names should be = ncol(model) - 1")
  
# Get the traits for the taxa - use taxonomic closeness to extend traits
  if (type == "BPc") 
    traitnames <- c(taxcname, "Ri", "Mi")
  else if (type == "IPc")
    traitnames <- c(taxcname, "BT", "FT", "ID")
  
  S <- try(Trait <- trait[,traitnames], silent = TRUE)
  if (inherits(S, "try-error"))
    stop ("Trait does not have all columns, should contain", paste(traitnames, collapse = ", "))
  
  Data_Db  <- get_trait(
        taxon    = taxa, 
        trait    = Trait, 
        taxonomy = taxonomy,
        verbose  = verbose)
  # one taxon can be present several times...
  if (nrow(Data_Db) != length(taxa)){  # taxa may have been present multiple times
    row.names(Data_Db) <- Data_Db[,1]
    Data_Db <- Data_Db[taxa, ]
  }
  
  # taxon for which trait could not be derived -  will be added to attributes
  notrait <- attributes(Data_Db)$notrait 
  
  Data_Db <- cbind(Data_Db, 1:nrow(Data_Db))  # add order
  nc      <- ncol(Data_Db)
  
  # weight data: either a vector, or a two-columned matrix
  if (is.vector(weight)){
    Data_Db$Weight <- weight
  } else if (ncol(weight) != 2){
    stop ("weight should either be a vector or a matrix or dataframe with taxon in 1st column and weight in 2nd column") 
  } else {colnames(weight)[2] <- "Weight"
  
    Data_Db <- merge(Data_Db, weight, 
                     by = 1, all.x = TRUE)
    Data_Db <- Data_Db[order(Data_Db[ ,nc]), ]   # restore order
  }  
  
  Noweight <- NULL
  
  if (any(is.na(Data_Db$Weight))){
    if (verbose) 
      warning("Weight not known for ", length(which(is.na(Data_Db$Weight)))," taxa - set to 0")
    Data_Db$Weight[is.na(Data_Db$Weight)] <- 0
    Noweight <- Data_Db$taxon[is.na(Data_Db$Weight)]
   }

# Part in formula that remains constant 
  if (type == "BPc") 
    Data_FAC <- with(Data_Db, sqrt(Weight)*Mi*Ri)  
  else if (type == "IPc")
    Data_FAC <- with(Data_Db, Weight^(0.75)*BT*FT*ID)  
  
# multiply all rows in model with the constant factor
  model_BPC <- sweep(model[,-1], 
                     MARGIN = 2, 
                     STATS  = Data_FAC, 
                     FUN    = "*")
  colnames(model_BPC) <- paste(type, taxa, sep="_")

# remove NAs    
  if (na.rm) {
    isna <- apply(model_BPC, 
                  MARGIN = 2, 
                  FUN    = function(x)!any(is.na(x)))
    model_BPC <- model_BPC[,isna]
  }
  
# add summed values  
  if (is.matrix(model_BPC))
    model_BPC <- cbind(model_BPC, 
                       rowSums(model_BPC, na.rm = TRUE))
  else
    model_BPC <- cbind(model_BPC, model_BPC)
  
  colnames(model_BPC)[ncol(model_BPC)] <- type
  model_BPC <- cbind(model[, 1], # time column
                     model_BPC)
  
  model_BPC <- as.matrix(model_BPC)
  
  colnames(model_BPC)[1] <- "time"
  
  class(model_BPC) <- c("deSolve", "matrix")
  
  attributes(model_BPC)$Factor <- Data_FAC 
  attributes(model_BPC)$Noweight <- Noweight 
  attributes(model_BPC)$notrait <- notrait
  return(model_BPC)
}
