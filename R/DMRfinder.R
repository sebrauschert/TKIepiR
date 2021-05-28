#' Differentially methylated regions (DMR) finder
#'
#' This function is a wrapper for the \code{DMRcate} bioconductor package to identify differentially methylated regions.
#' The input needs to be a result file from an EWAS, containing p-values, beta-coefficients and standard errors.
#' 
#' IMPORTANT: The CpG identifier column needs to be called "ID"!
#' 
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value, standard error and beta coefficient.
#' @param annotate The default is that the
#' @param p.column.name The name of the p-value column
#' @param beta.column.name The name of the beta coefficient column
#' @param se.column.name The name of the standard error column
#' @param C This is as per \code{DMRcate} package: Scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma.
#' Empirical testing shows that, for 450k data when lambda=1000, near- optimal prediction of sequencing-derived DMRs is obtained when C is approxi- mately 2, i.e.
#' 1 standard deviation of Gaussian kernel = 500 base pairs. Should be a lot larger for sequencing data - suggest C=50. Cannot be < 0.2
#' @param lambda This is according to the \code{DMRcate} package: Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR bookend definition;
#' gaps >= lambda between significant CpG sites will be in separate DMRs. Support is truncated at 5*lambda. Default is 1000 nucleotides. See details for further info.
#' @param pcutoff P-value threshold for selecting significan DMRs. Defaults to the DMRcate packages default "fdr"
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import Hmisc
#' @importFrom DMRcate dmrcate extractRanges
#' @return \code{DMRfinder} as a data table.
#' @export


DMRfinder <- function(EWAS, annotate = TRUE, p.column.name = p.column.name, beta.column.name = beta.column.name, se.column.name = se.column.name, C = 2, lambda = 100, pcutoff = "fdr"){
  
  ## In the same source file (to remind you that you did it) add:
  #if(getRversion() >= "2.15.1")  utils::globalVariables("naresid.omit")
  #if(getRversion() >= "3.1.0") utils::suppressForeignCheck("localvariable")
  
  if (annotate==TRUE){
    EWAS <- as.data.frame(EWAS)
    EWAS <- annotateCpG(EWAS)
  }
  if (annotate==FALSE){
    EWAS <- as.data.frame(EWAS[,c("ID", p.column.name, beta.column.name, se.column.name)])
    EWAS <- annotateCpG(EWAS)
  }
  # Prepare a data set to match the format of the dmracte algorithm
  DMRset        <- EWAS[,c("ID", "chr","pos",beta.column.name, p.column.name, se.column.name)]
  DMRset$stat   <- as.numeric(as.character(DMRset[,beta.column.name]))/as.numeric(as.character(DMRset[,se.column.name]))
  DMRset$infdr  <- p.adjust(as.numeric(as.character(DMRset[, p.column.name])), method="fdr")
  DMRset$betafc <- as.numeric(as.character(DMRset[,beta.column.name]))
  DMRset$is.sig <- DMRset$infdr <= 0.05

  # Create a list object, as this is required for DMRcate
  model         <- list(DMRset$ID, DMRset$stat, DMRset$chr, DMRset$pos, DMRset$betafc, DMRset$infdr, DMRset$is.sig)
  names(model)  <- c("ID", "stat", "chr", "pos", "betafc", "infdr", "is.sig")


  annotated <- data.frame(ID = model$ID, stat = model$stat,
                          CHR = model$chr, pos = model$pos, betafc = model$betafc,
                          indfdr = model$infdr, is.sig=model$infdr < 0.05)

  annotated <- annotated[order(annotated$CHR, annotated$pos),]
  #class(annotated) <- "annot"
  
  # The newest version of dmrcate calls the object CpGannotated rather than annot.
  class(annotated) <- "CpGannotated"

  annotated$dmrcate                  <- DMRcate::dmrcate(annotated, lambda=lambda, C=C, pcutoff=pcutoff)
  annotated$dmrcate$results.ranges   <- DMRcate::extractRanges(annotated$dmrcate, genome="hg19")
  as.data.frame(annotated$dmrcate$results.ranges)
}
