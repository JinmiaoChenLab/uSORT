#' @author MaiChan Lau
#' @export
#' 
STS_sortingcost<-function(expr = NULL, dist_metric_flag=NULL){
  ## Define parameters & constant
  n<-nrow(expr)
  I<-diag(n)
  
  # a strictly increasing column vector 
  X<-matrix(ncol=1,nrow=n)
  for (i in 1:n){
    X[i]<-i-(n+1)/2
  }
  
  d <- distance.function(expr)
  
  cost<-sum(diag(I%*%d%*%t(I)%*%X%*%t(X)))
  
  return(cost)
}



