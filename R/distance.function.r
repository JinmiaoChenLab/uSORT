#' A distance function

#' This function allows you to compute cell-to-cell distances
#' 
#' @param expr an expresssion matrix containing n-rows of cells and m-cols of genes
#' @param flag an integer indicating the type of distance metric: (1) euclidean distance,
#' (2) correlations, (3) jaccard, (4) input as distance
#' @param method a function, a registry entry, or a mnemonic string referencing the 
#' proximity measure, parameter for \code{proxy::dist} method
#' @return a n-by-n pairwise distance matrix for n cells
#' @author MaiChan Lau
#' @importFrom proxy dist
#' @export


distance.function<-function(expr, flag=1, method= 'eJaccard'){
  ## Define parameters
  n<-nrow(expr)
  
  ## Compute Distance
  # euclidean distance
  if(flag==1){
    vecs<-as.matrix(expr)
    g2=rowSums(vecs*vecs)
    colRepeat<-do.call("cbind", replicate(n, g2, simplify = FALSE))
    vv<-2*as.matrix(vecs)%*%(t(as.matrix(vecs)))
    rowRepeat<-do.call("rbind", replicate(n, g2, simplify = FALSE))
    D<-sqrt(abs(colRepeat-vv+rowRepeat))
    colnames(D)<-row.names(D)  
  }
  # correlations
  else if (flag==2){
    D <- as.dist(1 - cor(t(expr)))
  }
  # jaccard
  else if (flag==3){
    D <- dist(expr, method = method)
  }
  # input as distance
  else {
    D =  expr
  }
  
  
  return(D)
}
