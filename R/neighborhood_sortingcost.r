#' @author MaiChan Lau
#' @export
#' 
neighborhood_sortingcost<-function(expr = NULL, dist_metric_flag=NULL, sigma_width=1){
    ## Define parameters & constant
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
    
    I<-diag(mat_size)
    
    d <- distance.function(expr)
    
    cost<-sum(diag(I%*%d%*%t(I)%*%weights_mat%*%t(weights_mat)))
    
    return(cost)
}



