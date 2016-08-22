#' A sorting function using the neighborhood algorithm
#' 
#' This function allows you to order cells according to the neighborhood
#' algorithm proposed in SPIN paper Tsafrir, D. et al. (2005). Bioinformatics, 21(10), page 2301 to 2308.
#' 
#' @param d an expresssion matrix containing n-by-n pairwise distances computed for n cells
#' @param max_iter maximum number of iteration
#' @return a list containinig a vector \code{ordering} of the permutated cell ordering, and a cost value \code{cost}
#' @author MaiChan Lau
#' @export
#' 
neighborhood_sorting<-function(d,  weights_mat=NULL, max_iter=100){
    
    if(is.null(weights_mat))  stop("weights_mat is not available for neighborhood sorting!")
    
    #
    mat_size <- nrow(d)
    I<-diag(mat_size)
    t<-0
    Dt<-d
    global_permutation<-I
    prev_permutation<-I
    cost<-sum(diag(I%*%Dt%*%t(I)%*%weights_mat))
    
    flag<-'on'
    
    while(flag=='on' && t<max_iter){
        
        ## ============
        mismatch = Dt %*% weights_mat;
        val <- Biobase::rowMin(mismatch)
        mn <- apply(mismatch,1,which.min)
        main_diag <- diag(mismatch)   
        sort_score <- seq(mat_size)
        mx <- max(val)
        sort_score <- mn - 0.1*sign((mat_size/2-mn))*val/(mx)
        
        # Sorting the matrix
        sorted_ind <- order(sort_score)
        val <- sort_score[sorted_ind]
        
        # update of ordering
        Pt<-I[sorted_ind,]
        global_permutation<-global_permutation[sorted_ind,]
        cost<-sum(diag(global_permutation%*%d%*%t(global_permutation)%*%weights_mat))
        #cat('@t:', t, ' cost=', cost, '\n')
        ## ============
        
        ## 
        if(isTRUE(all.equal((Pt%*%mismatch),(prev_permutation%*%mismatch)))){
            flag<-'off'
        }
        else{
            Dt<-Dt[sorted_ind,sorted_ind]
            t<-t+1
            prev_permutation<-Pt
        }
    } #end while-loop
    
    #cat(val)
    ordering<-apply(global_permutation,1,which.max)
    return(list('ordering'=ordering, 'cost'=cost))
}



