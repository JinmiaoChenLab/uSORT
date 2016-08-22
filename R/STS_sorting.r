#' A sorting function using the Side-to-side (STS) algorithm
#' 
#' This function allows you to order cells according to the Side-to-side (STS) 
#' algorithm proposed in SPIN paper Tsafrir, D. et al. (2005). Bioinformatics, 21(10), page 2301 to 2308.
#' 
#' @param d an expresssion matrix containing n-by-n pairwise distances computed for n cells
#' @param max_iter maximum number of iteration
#' @return a list containinig a vector \code{ordering} of the permutated cell ordering, and a cost value \code{cost}
#' @author MaiChan Lau
#' @export
#' 
STS_sorting<-function(d, max_iter=10){
  ## Define parameters & constant
  n<-nrow(d)
  I<-diag(n)
  # a strictly increasing column vector 
  X<-matrix(ncol=1,nrow=n)
  for (i in 1:n){
    X[i]<-i-(n+1)/2
  }
  
  ## STS algorithm
  # Step 1
  t<-0
  Dt<-d
  global_permutation<-I
  prev_permutation<-I
  cost<-sum(diag(I%*%Dt%*%t(I)%*%X%*%t(X)))
  
  flag<-'on'
  
  while(flag=='on' && t<max_iter){
    
    # Step 2
    St<-Dt%*%X
    
    # Step 3
    descending.order<-order(St, decreasing=T)
    Pt<-I[descending.order,]
    global_permutation<-global_permutation[descending.order,]
    cost<-sum(diag(global_permutation%*%d%*%t(global_permutation)%*%X%*%t(X)))
    
    ## Step 4
    if(isTRUE(all.equal((Pt%*%St),(prev_permutation%*%St)))){
      flag<-'off'
    }
    else{
      Dt<-Dt[descending.order,descending.order]
      t<-t+1
      prev_permutation<-Pt
    }
  } #end while-loop
  
  
  ordering<-apply(global_permutation,1,which.max)
  return(list('ordering'=ordering, 'cost'=cost))
}



