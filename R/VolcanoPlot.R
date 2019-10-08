#' Create a Volcano plot
#'
#'This function utilizes the \code{annotateCpG()} function of this package. It
#'creates a Volcano plot with annotated gene names based on the \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} bioconductor package.
#' 
#' IMPORTANT: The CpG identifier column needs to be called "ID"!
#' 
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value,
#' @param annotate If TRUE (default) this will annotatate the EWAS results file with the gene and position info, using the \code{annotateCpG()} function of
#' this package
#' @param p.column.name This specifies the column for the model p-value in the EWAS result file
#' @param beta.column.name The name of the beta coefficient column
#' @param title Specify a title for the Volcano plot
#' @return \code{VolcanoPlot} as an image.
#' @export


VolcanoPlot <- function(EWAS, annotate=TRUE, p.column.name = p.column.name, beta.column.name = beta.column.name, title=title){
  
  if (annotate==TRUE){
    EWAS <- as.data.frame(EWAS)
    #EWAS1 <- EWAS[-which(EWAS[,1] %in% ""),]
    EWAS <- annotateCpG(EWAS)
  }
  if (annotate==FALSE){
    EWAS <- as.data.frame(EWAS[,c("ID", p.column.name)])
    EWAS <- annotateCpG(EWAS)
  }
  
  EWAS[,p.column.name] <- as.numeric(as.character(EWAS[,p.column.name]))
  EWAS[,beta.column.name] <- as.numeric(as.character(EWAS[,beta.column.name]))
  
  # Create a column indicating Bonferroni significance so it can be color coded in the volcano plot
  EWAS$SIGNI <- ifelse((EWAS[,p.column.name] < 0.05/length(EWAS[,p.column.name])),"significant","not significant")
  
  EWAS$GENE <- ifelse(EWAS$SIGNI %in% "significant", EWAS$UCSC_RefGene_Name, NA)
  
  ggplot() + 
    geom_point(data = EWAS, mapping = aes(y=EWAS[,beta.column.name], x=-log(EWAS[,p.column.name], base=10),colour=EWAS$SIGNI),size=2) +
    
    geom_vline(aes(xintercept=-log(0.05/475429,base=10)),col="red")+
    geom_vline(aes(xintercept=-log(0.05,base=10)),col="darkgrey")+
    geom_hline(aes(yintercept=0),linetype="dotted")+
    
    # Annotate the significant CpGs with the nearest gene name
    geom_text_repel(
      aes(y=EWAS[, beta.column.name], x=-log(EWAS[,p.column.name], base=10), label =as.character(EWAS$GENE)),segment.color = "black",
      force=10, size=2) +
    labs(y= expression(beta*-Coefficient), x=expression(-log[10](italic(p))), color="Bonferroni corrected:") +
    
    coord_flip() +
    ggtitle(title)
  
}
