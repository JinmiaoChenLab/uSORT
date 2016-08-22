#' A neighborhood sorting wrapper function
#' 
#' A sorting wrapper function seeks for an (locally) optimal cell permutation by 
#' repeating the neighborhood sorting operation (i.e. neighborhood_sorting) on randomly chosen starting ordering
#' 
#' @param expr an expresssion matrix containing n-rows of cells and m-cols of genes
#' @param repeats number of randomly chosen starting ordering to be tested
#' @return a list containinig a matrix \code{permutated.expr} of the permutated expression matrix
#' and the best cost value \code{best.cost}
#' @author MaiChan Lau
#' @export


neighborhood_sorting_wrapper<-function(expr, sigma_width=1, no_randomization=10){
    
    ## Initialization
    best.cost<-1e20
    set.seed(123)
    mat_size <- nrow(expr)
    s <- seq(mat_size)
    m <- rep(s,mat_size)
    i <- matrix(m,ncol=mat_size,nrow = mat_size); i <- t(i)
    j <- matrix(m,ncol=mat_size,nrow = mat_size)
    
    G = exp(-(i-j)^2/sigma_width/mat_size)
    
    for (i in 1:10){
        G <- G/(matrix(rep(colSums(G),mat_size), ncol = mat_size, nrow = mat_size, byrow = T))
        G <- G/(matrix(rep(rowSums(G),mat_size), ncol = mat_size, nrow = mat_size, byrow = F))
    }
    G <- (G + t(G))/2
    weights_mat <- G
    
    for (i in 1:no_randomization){
      
        if(i==1) temp.expr<-expr
        else{
          rand_id<-sample.int(mat_size)
          temp.expr<-expr[rand_id,,drop=F]
        }  
        temp.dist<-distance.function(temp.expr)
        temp.SPIN.res<-neighborhood_sorting(temp.dist, weights_mat)
        
        if(temp.SPIN.res$cost<best.cost){
            temp.best.ordering<-temp.expr[temp.SPIN.res$ordering,,drop=F]
            best.ordering<-temp.best.ordering
            best.cost<-temp.SPIN.res$cost
        }
      
    }#end for-loop
    
    return(list('permutated.expr'=best.ordering,'best.cost'=best.cost))
}