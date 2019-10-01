#' Annotate the CpG ID's with genetic information
#'
#' This funtion annotates the results of an EWAS with the "IlluminaHumanMethylation450kanno.ilmn12.hg19" package.
#' The function calls the annotation data frame from the IlluminaHumanMethylation450kanno.ilmn12.hg19 bioconductor package
#' and returns a new data frame, with all the CpG, gene and location information as well as the EWAS model results
#'
#' IMPORTANT: The CpG identifier column needs to be called "ID"!
#' 
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value, standard error and beta coefficient.
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19 
#' @importFrom minfi getAnnotation
#' @return \code{annotateCpG} as an data table.
#' @export

annotateCpG <- function(EWAS){
  
  # Get the annotation for the Illumina 450k
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  Annot = minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  Annot$ID         <- rownames(Annot)
  rownames(Annot)  <- NULL
  
  # Make sure that the CpG columns is called "ID" for both data frames
  #names(EWAS)[1] <- "ID"
  EWAS_Model <- merge(Annot, EWAS, by="ID")
  rm(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  EWAS_Model
}
