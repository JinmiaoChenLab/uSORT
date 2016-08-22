#' An expression averaging function for duplicated genes
#' 
#' An expression averaging function which takes the average value of duplicated genes as 
#' Monocle's CellDataSet object cannot handle duplicated genes.
#' 
#' @param exprs a matrix of normalized expression data (TPM or FPKM) containing n-rows of cells 
#' and m-cols of genes 
#' @return an averaged expression data with unique gene names
#' @author MaiChan Lau
#' @export

avg_same_genes <- function(exprs = NULL){
  
  dup_flag<-duplicated(rownames(exprs))
  unique_names<-unique(rownames(exprs)[dup_flag==TRUE])
  if(length(unique_names)>0){
    for(i in 1:length(unique_names)){
      exprs<-avg_function(unique_names[i],exprs)
    }
  }
  
  return (exprs)
}


#' An expression averaging function for an unique gene
#' 
#' An expression averaging function which collapses multiple expression with same gene name to the 
#' average value
#' 
#' @param gene_name a character of gene name
#' @param exprs a matrix of normalized expression data (TPM or FPKM) containing n-rows of cells 
#' and m-cols of genes 
#' @return an expression data with average expression value computed for \code{gene_name}
#' @author MaiChan Lau
#' @export

avg_function <- function(gene_name, exprs){
  
    if(is.null(ncol(exprs))) dim(exprs)<-c(length(exprs),1)
    
    id<-which(rownames(exprs) %in% gene_name)
    avg_val<-colMeans(exprs[id,])
    exprs[id[1],]<-avg_val
    exprs<-exprs[-id[2:length(id)],]
    return(exprs)
}

