#' Create a Manhattan plot
#'
#'This function utilizes both the \code{annotateCpG()} function of this package, as well as the \code{qqman} CRAN package. It
#'creates a Manhattan plot with annotated gene names based on the \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} bioconductor package.
#' 
#' IMPORTANT: The CpG identifier column needs to be called "ID"!
#' 
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value,
#' @param annotate If TRUE (default) this will annotatate the EWAS results file with the gene and position info, using the \code{annotateCpG()} function of
#' this package
#' @param p.column.name This specifies the column for the model p-value in the EWAS result file
#' @param title Specify a title for the Manhattan plot
#' @param col.scheme Specify two colors for the alternating color scheme of the Manhattan plot
#' @importFrom qqman manhattan
#' @return \code{ManhattanPlot} as an image.
#' @export

ManhattanPlot <- function(EWAS, annotate=TRUE, p.column.name=p.column.name, title="title", col.scheme=c("turquoise4", "springgreen4")){

  if (annotate==TRUE){
    EWAS <- as.data.frame(EWAS)
    #EWAS1 <- EWAS[-which(EWAS[,1] %in% ""),]
    EWAS <- annotateCpG(EWAS)
  }
  if (annotate==FALSE){
    EWAS <- as.data.frame(EWAS[,c("ID", p.column.name)])
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
  qqman::manhattan(Mod.Man, main=title, annotatePval = (0.05/length(Mod.Man$P)),annotateTop=F, col=col.scheme)
}
