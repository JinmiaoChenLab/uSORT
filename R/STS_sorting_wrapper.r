#' A STS sorting wrapper function
#' 
#' A sorting wrapper function seeks for an (locally) optimal cell permutation by 
#' repeating the STS sorting operation (i.e. STS_sorting) on randomly chosen starting ordering
#' 
#' @param expr an expresssion matrix containing n-rows of cells and m-cols of genes
#' @param repeats number of randomly chosen starting ordering to be tested
#' @return a list containinig a matrix \code{permutated.expr} of the permutated expression matrix
#' and the best cost value \code{best.cost}
#' @author MaiChan Lau
#' @export


STS_sorting_wrapper<-function(expr, no_randomization=10){
    ## Define parameters
    n<-nrow(expr)
    ## Initialization
    best.cost<-1e20
    set.seed(123)
    for (i in 1:no_randomization){
      
        if(i==1) temp.expr<-expr
        else{
          rand_id<-sample.int(n)
          temp.expr<-expr[rand_id,,drop=F]
        }  
        temp.dist<-distance.function(temp.expr)
        temp.SPIN.res<-STS_sorting(temp.dist)
        
        if(temp.SPIN.res$cost<best.cost){
            temp.best.ordering<-temp.expr[temp.SPIN.res$ordering,,drop=F]
            best.ordering<-temp.best.ordering
            best.cost<-temp.SPIN.res$cost
        }
      
    }#end for-loop
    
    return(list('permutated.expr'=best.ordering,'best.cost'=best.cost))
}