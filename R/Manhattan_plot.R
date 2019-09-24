#' Annotate the CpG ID's with genetic information
#'
#' This funtion annotates the results of an EWAS with the "IlluminaHumanMethylation450kanno.ilmn12.hg19" package.
#' The function calls the annotation data frame from the IlluminaHumanMethylation450kanno.ilmn12.hg19 bioconductor package
#' and returns a new data frame, with all the CpG, gene and location information as well as the EWAS model results
#'
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value, standard error and beta coefficient.
#' @return \code{annotateCpG} as an data table.
#' @examples
#'  annotateCpG(data)
#' @export

annotateCpG <- function(EWAS){

  # Get the annotation for the Illumina 450k
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  Annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  Annot$ID         <- rownames(Annot)
  rownames(Annot)  <- NULL

  # Make sure that the CpG columns is called "ID" for both data frames
  names(EWAS)[1] <- "ID"
  EWAS_Model <- merge(Annot, EWAS, by="ID")
  EWAS_Model
}



#' Create a Manhattan plot
#'
#'This function utilizes both the \code{annotateCpG()} function of this package, as well as the \code{qqman} CRAN package. It
#'creates a Manhattan plot with annotated gene names based on the \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} bioconductor package.
#'
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value,
#' @param annotate If TRUE (default) this will annotatate the EWAS results file with the gene and position info, using the \code{annotateCpG()} function of
#' this package
#' @param p.column.name This specifies the column for the model p-value in the EWAS result file
#' @param title Specify a title for the Manhattan plot
#' @param col.scheme Specify two colors for the alternating color scheme of the Manhattan plot
#' @return \code{ManhattanPlot} as an image.
#' @examples
#'  ManhattanPlot(data)
#' @export

ManhattanPlot <- function(EWAS, annotate=TRUE, p.column.name=p.column.name, title=title, col.scheme=c("turquoise4", "springgreen4")){

  if (annotate==TRUE){
    EWAS <- as.data.frame(EWAS)
    #EWAS1 <- EWAS[-which(EWAS[,1] %in% ""),]
    EWAS <- annotateCpG(EWAS)
  }
  # We need to extract the following columns for the Manhattan plot, as they are required for the qqman package:
  # a) CpG name
  # b) Nearest Gene
  # c) Chromosome
  # d) Position
  # e) P value
  # f) SNP
  # The SNP column orignially is intended to contain the rs SNP number in a GWAS analysis, but here we put the Gene name in this position, so we can
  # Annotate it in the Manhattan plot

  Mod.Man        <- EWAS[, c("ID","CpG_rs", "chr", "pos",p.column.name, "UCSC_RefGene_Name")]
  names(Mod.Man) <- c("ID","Gene..nearest.Gene.", "CHR", "BP", "P", "SNP" ) # SNP and Gene name changed for annotation

  # Making sure the columns have the right file type (numeric)
  Mod.Man$CHR <- as.numeric(as.character(substring(Mod.Man$CHR, 4)))
  Mod.Man$P   <- as.numeric(as.character(Mod.Man$P))

  # Remove X and Y chromosomes
  Mod.Man     <- data.frame(subset(Mod.Man, Mod.Man$CHR %nin% NA))

  # qqman Manhattan plot function
  qqman::manhattan(Mod.Man, main=title, annotatePval = (0.05/462925),annotateTop=F, col=col.scheme)
}
