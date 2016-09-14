#' A wrapper function for autoSPIN sorting method
#'
#' A wrapper function for autoSPIN method which implements optimized local refinement using the selected
#' SPIN sorting method, i.e. STS or Neighborhood.
#'
#' @param data An log2 transformed expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param data_type A character string indicating the type of progression, i.e. 'linear' (strictly linear) or
#' 'cyclical' (cyclically linear).
#' @param sorting_method A character string indicating the choice of SPIN sorting method, i.e. 'STS' (Side-to-Side)
#' or 'Neighborhood'.
#' @param alpha A fraction value denoting the size of locality used for calculating the summed local variance.
#' @param sigma_width An integer number denoting the degree of spread of the gaussian distribution which
#' is used for computing weight matrix for Neighborhood sorting method.
#' @param no_randomization An integer number indicating the number of repeated sorting, each of which
#' uses randomly selected initial cell position.
#' @param window_perc_range A fraction value indicating the range of window size to be examined during local refinement.
#' @param window_size_incre_perct A fraction value indicating the step size at each iteration
#' for incrementing window size.
#' @return A data frame containing single column of ordered sample IDs.
#' @export
#'
#' @examples
#' set.seed(15)
#' da <- iris[sample(150, 150, replace = FALSE), ]
#' rownames(da) <- paste0('spl_',seq(1,nrow(da)))
#' d <- da[,1:4]
#' dl <- da[,5,drop=FALSE]
#' res <- autoSPIN(data = d)
#' dl <- dl[match(res$SampleID,rownames(dl)),]
#' annot <- data.frame(id = seq(1,nrow(res)), label=dl, stringsAsFactors = FALSE)
#' #ggplot(annot, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
autoSPIN <-function(data,
                    data_type = c('linear', 'cyclical'),
                    sorting_method = c("STS", "neighborhood"),
                    alpha = 0.2,
                    sigma_width = 1,
                    no_randomization = 20,
                    window_perc_range = c(0.1, 0.9),
                    window_size_incre_perct = 0.05) {

    n <- nrow(data)
    sorting_method <- match.arg(sorting_method)
    data_type <- match.arg(data_type)

    ## =====Global sorting=========================
    if(sorting_method == 'STS') {
        global_res <- STS_sorting_wrapper(data,
                                          no_randomization= no_randomization)
    }else if(sorting_method == 'neighborhood') {
        global_res <- neighborhood_sorting_wrapper(data,
                                                   no_randomization = no_randomization,
                                                   sigma_width = sigma_width)
    }
    else{stop('Sorting type can onlyl be STS or neighborhood')}
    global_ordering <- global_res$permutated.expr
    global_distVar <- summed_local_variance(global_ordering,
                                            alpha=alpha,
                                            data_type=data_type)

    ## =====Local sorting=========================
    cat('  Local sorting...\n')
    lowest_distVar<-global_distVar
    incremental_window_size<-ceiling(window_perc_range[1]*n)
    list_distVar<-c()
    list_windowSize<-c()

    ## Check if starting window size is larger than limit
    if(incremental_window_size>=ceiling(n*window_perc_range[2])){
        cat('  Warning: starting window size > allowable size, \n')
        cat('  No local sorting refinement performed\n')
        global_temp<-global_ordering;break
    }

    ## Iterate through incrementally enlarging local window sizes
    while(incremental_window_size < ceiling(n*window_perc_range[2])){

        if(length(list_windowSize)==0){
            cat(paste0('  Testing for window size = ',
                       incremental_window_size, ","))
        }else{
            cat(paste0(incremental_window_size, ","))
        }
        stop_flag=FALSE # set to TRUE when the size of unhandled (local sorting)
        list_windowSize<-c(list_windowSize,incremental_window_size)

        # Apply STS_sorting to local data window
        # Iterate through from left to right every windows with size <= incremental_window_size
        for (i in 1:floor(n/incremental_window_size))
        {
            startID<-(i-1)*incremental_window_size+1
            endID<-i*incremental_window_size

            if(endID>=(n-1) | ((i+1)*incremental_window_size)>n)
            {endID<-n; stop_flag=TRUE}

            # Get local data window for sorting
            data_window<-global_ordering[startID:endID,,drop=FALSE]

            if(sorting_method == 'STS'){
                local_res <- STS_sorting_wrapper(data_window,
                                                 no_randomization=no_randomization)
            }else if(sorting_method == 'neighborhood') {
                local_res<-neighborhood_sorting_wrapper(data_window,
                                                        no_randomization= no_randomization,
                                                        sigma_width = sigma_width)
                }
            local_ordering<-local_res$permutated.expr

            if(i==1) temp_combined_ordering<-local_ordering
            else temp_combined_ordering<-rbind(temp_combined_ordering,
                                               local_ordering)

            if(stop_flag==TRUE) break
        }

        temp_combined_dist <- distance.function(temp_combined_ordering)

        # Check if reversal operation is needed
        one_sided_length<-ceiling(incremental_window_size)*0.5
        iter_no<-floor(n/incremental_window_size)

        # Iterate through from left to right every windows with size <= incremental_window_size
        for (i in 1:iter_no)
        {
            if(iter_no<=1) break

            startID<-(i-1)*incremental_window_size+1
            endID<-i*incremental_window_size

            if(endID>n  | ((i+1)*incremental_window_size)>n) endID<-n
            idx<-seq(1,n)
            idx[startID:endID]<-rev(idx[startID:endID])
            partially_rev_dist<-temp_combined_dist[idx,idx]

            # 1st window is compared to the following RHS window
            if (i==1){
                # Reversed local window reversed: sum of distance across neighboring windows
                between_dist_rev<-sum(partially_rev_dist[(endID-one_sided_length):endID,
                                                         (endID+1):(endID+one_sided_length)])

                # Original local window: sum of distance across neighboring windows
                between_dist<-sum(temp_combined_dist[(endID-one_sided_length):endID,
                                                     (endID+1):(endID+one_sided_length)])

                if(between_dist>between_dist_rev) {
                    temp_combined_ordering<-temp_combined_ordering[idx,,drop=FALSE]
                    temp_combined_dist[idx,idx]
                }
            }else{
                # Reversed local window reversed: sum of distance across neighboring windows
                between_dist_rev<-sum(partially_rev_dist[(startID-one_sided_length):startID,
                                                         startID:(startID+one_sided_length)])
                # Original local window: sum of distance across neighboring windows
                between_dist<-sum(temp_combined_dist[(startID-one_sided_length):startID,
                                                     startID:(startID+one_sided_length)])

                if(between_dist>between_dist_rev) {
                    temp_combined_ordering<-temp_combined_ordering[idx,,drop=FALSE]
                    temp_combined_dist[idx,idx]
                }
            }

        }


        temp_combined_distVar<-summed_local_variance(temp_combined_ordering,
                                                     alpha=alpha,
                                                     data_type=data_type)


        list_distVar<-c(list_distVar,temp_combined_distVar)

        # Update global ordering if local sorting achieve better fine-grained smoothness
        if(temp_combined_distVar<lowest_distVar){
            best_local_ordering<-temp_combined_ordering
            lowest_distVar<-temp_combined_distVar
        }

        incremental_window_size<-incremental_window_size +
            ceiling(n*window_size_incre_perct)
    }
    cat('\n')

    if(exists("best_local_ordering")){
        global_ordering<-best_local_ordering
        }
    ordering <- data.frame('SampleID'=rownames(global_ordering))

    return(ordering)
}





#' A distance function

#' A distance function computes cell-to-cell distance matrix.
#'
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param method A character string indicating the distance function.
#'
#' @return A matrix containing n-by-n cell distance.
distance.function<-function(expr,
                            method= c('Euclidean','Correlation','eJaccard','none')){
    ## Define parameters
    n<-nrow(expr)
    method_name <- match.arg(method)

    switch(method_name,
           Euclidean = {
               vecs<-as.matrix(expr)
               g2=rowSums(vecs*vecs)
               colRepeat<-do.call("cbind", replicate(n, g2,
                                                     simplify = FALSE))
               vv<-2*as.matrix(vecs)%*%(t(as.matrix(vecs)))
               rowRepeat<-do.call("rbind", replicate(n, g2,
                                                     simplify = FALSE))
               D<-sqrt(abs(colRepeat-vv+rowRepeat))
               colnames(D)<-row.names(D)
           },
           Correlation = {
               D <- as.dist(1 - cor(t(expr)))
           },
           eJaccard = {
               D <- dist(expr, method = method)
           },
           none = {
               D =  expr
           })

    return(D)
}


#' A summed local variance function
#'
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param alpha A fraction value indicating the size of window for local variance measurement.
#' @param data_type A character string indicating the type of progression, i.e. 'linear' (strictly linear) or
#' 'cyclical' (cyclically linear).
#'
#' @return A numeric value of the summed local variance.
summed_local_variance <- function(expr = NULL,
                                       alpha = NULL, data_type = 'linear'){

    d <-distance.function(expr)
    if(data_type == 'linear') {
        d_var <-summed_local_variance_linear(d, alpha=alpha)
    }else if(data_type == 'cyclical'){
        d_var <- summed_local_variance_cyclical(d, alpha=alpha)
    }
    return (d_var)

}



#' A summed local variance function for strictly linear data type
#'
#' @param d A cell-to-cell distance matrix.
#' @param alpha A fraction value indicating the size of window for local variance measurement.
#'
#' @return A numeric value of the summed local variance.
summed_local_variance_linear<-function(d, alpha=0.3){

    ## Initialization
    num_cell<-nrow(d)
    total.var<-0
    window_size<-ceiling(num_cell*alpha)

    for(i in (as.integer(window_size/2)):(num_cell-(as.integer(window_size/2)))){
        window_around_cell<-d[i,]
        startID <- i - (as.integer(window_size/2))
        endID <- i + (as.integer(window_size/2))
        if(startID<1)startID<-1
        if(endID>num_cell) endID <- num_cell
        #cat(i,'\t',startID,'\t',endID,'\n')

        window_around_cell<-window_around_cell[startID:endID]
        window_around_cell <- window_around_cell[-i]
        ind.cell.var<-var(window_around_cell)
        total.var<-total.var+ind.cell.var
    }

    return(total.var)
}



#' A summed local variance function for cyclical linear data type
#'
#' @param d A cell-to-cell distance matrix.
#' @param alpha A fraction value indicating the size of window for local variance measurement.
#'
#' @return A numeric value of the summed local variance.
summed_local_variance_cyclical<-function(d, alpha=0.3){

    ## Initialization
    num_cell<-nrow(d)
    total.var<-0
    window_size<-ceiling(num_cell*alpha)

    for(i in 1:num_cell){
        window_around_cell<-d[i,]
        window_around_cell_extended <- c(window_around_cell,
                                         window_around_cell,
                                         window_around_cell)

        startID <- i - (as.integer(window_size/2))
        endID <- i + (as.integer(window_size/2))

        if(startID<1){
            extended_startID <- num_cell-i
        }else{
            extended_startID <- num_cell+i}

        extended_endID <- num_cell + endID

        window_around_cell_extended<-window_around_cell_extended[extended_startID:extended_endID]
        window_around_cell_extended <- window_around_cell_extended[-i]
        ind.cell.var<-var(window_around_cell_extended)
        total.var<-total.var+ind.cell.var
        #cat(i, '\t', 'ind.cell.var = ', ind.cell.var, '\n')
    }

    return(total.var)
}

