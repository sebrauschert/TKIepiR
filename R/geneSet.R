#' Gene set enrichment analysis
#'
#' This function is a wrapper for the \code{methylGSA} bioconductor package, which is a very easy to use method
#' for gene set enrichment analysis
#'
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value, standard error and beta coefficient.
#' @param p.val.col Columna name for the p-value in the EWAS file
#' @param minsize Mimimum amount of genes in GO pathway
#' @param maxsize MAximum amount of genes in GO pathway  
#' @param sig.gut Significance cut off for the gene set enrichment: Only CpGs with a value equal or smaller than this will 
#' be included
#' @param plottitle Title for the heatmap-barplot, if \code{plot.it} is TRUE
#' @param plot.it If true, plot will be returned rather than the table with GO terms
#' @return \code{geneSet} as a data table or heatmap-type barplot.
#' @examples
#'  geneSet(EWAS)
#' @export


geneSet <- function(EWAS, p.val.col = p.val.col, minsize = 100, maxsize = 1000, sig.cut = 0.001, plottitle="title", plot.it=FALSE){
  
  EWAS            <- as.data.frame(EWAS)
  EWAS            <- subset(EWAS, as.numeric(EWAS[,p.val.col]) %nin% NA)
  cpg.pval        <- as.numeric(EWAS[, p.val.col])
  
  names(EWAS)[1] <- "ID"
  names(cpg.pval) <- EWAS$ID
  
  # Perform the gene set enrichment analysis
  results = methylGSA::methylgometh(cpg.pval = cpg.pval, sig.cut = sig.cut, 
                         minsize = minsize, maxsize = maxsize, topDE=10)
  if (plot.it == TRUE){
  methylGSA::barplot(results, colorby="padj", num= 10) + theme_minimal() + ggtitle(plottitle)
  } 
  else{ 
    results
  }
}
