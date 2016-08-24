#' A sorting wrapper function
#'
#' A sorting wrapper function includes 2 additional features for the original SPIN method,
#' (1) optimization with multiple randomly chosen starting ordering in
#' \code{STS_sorting_wrapper}, and (2) local sorting refinement through
#' \code{c2c_dist_uneveness}
#'
#' @param expr an expresssion matrix containing n-rows of cells and m-cols of genes
#' @param alpha window size which is given in  percentage of total no. of cells,
#' parameter for \code{c2c_dist_uneveness}
#' @param window_perc_range range of window size given in percentages of total no. of cells
#' indicating the smallest and largest window sizes for local sorting refinement
#' @param window_size_incre_perct size increment in percentage for local sorting refinement
#' @param sorting_method character specifying the type of sorthing algorithm, i.e. STS or neighborhood
#' @param data_type
#' @param sigma_width
#' @param no_randomization
#' @param baseNm base name used for saving output plot
#'
#' @return a matrix of the permutated expression matrix
#' @author MaiChan Lau
#' @export
autoSPIN <-function(expr,
                    data_type = c('linear', 'cyclical'),
                    sorting_method = c("STS", "neighborhood"),
                    baseNm = 'uSPIN',
                    alpha = 0.2,
                    sigma_width = 1,
                    no_randomization = 10,
                    window_perc_range = c(0.1, 0.9),
                    window_size_incre_perct = 0.05) {

    n <- nrow(expr)
    cat('alpha = ', alpha, '\n')
    sorting_method <- match.arg(sorting_method)
    data_type <- match.arg(data_type)

    ## =====Global sorting=========================
    cat('Global sorting ...... in progress\n')
    if(sorting_method == 'STS') {
        global_res<-STS_sorting_wrapper(expr, no_randomization= no_randomization)
    }else if(sorting_method == 'neighborhood') {
        global_res<-neighborhood_sorting_wrapper(expr, no_randomization= no_randomization, sigma_width = sigma_width)
    }
    else stop('Sorting type can onlyl be STS or neighborhood')
    global_ordering<-global_res$permutated.expr
    global_distVar<-c2c_dist_uneveness_wrapper(global_ordering, alpha=alpha,
                                               data_type=data_type)

    ## =====Local sorting=========================

    cat('Local sorting ...... in progress\n')
    lowest_distVar<-global_distVar
    incremental_window_size<-ceiling(window_perc_range[1]*n)
    list_distVar<-c()
    list_windowSize<-c()

    ## Check if starting window size is larger than limit
    if(incremental_window_size>=ceiling(n*window_perc_range[2])){
        cat('Warning: starting window size > allowable size, \n')
        cat('no local sorting refinement performed\n')
        global_temp<-global_ordering;break
    }

    ## Iterate through incrementally enlarging local window sizes
    while(incremental_window_size < ceiling(n*window_perc_range[2])){

        if(length(list_windowSize)==0) cat('Testing for window size=',incremental_window_size,'...\t')
        else cat(incremental_window_size,'...\t')
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
            data_window<-global_ordering[startID:endID,,drop=F]

            if(sorting_method == 'STS') local_res<-STS_sorting_wrapper(data_window, no_randomization=no_randomization)
            else if(sorting_method == 'neighborhood') local_res<-neighborhood_sorting_wrapper(data_window, no_randomization= no_randomization, sigma_width = sigma_width)
            local_ordering<-local_res$permutated.expr

            if(i==1) temp_combined_ordering<-local_ordering
            else temp_combined_ordering<-rbind(temp_combined_ordering,local_ordering)

            if(stop_flag==TRUE)break
        }


        temp_combined_dist<-distance.function(temp_combined_ordering)

        # Check if reversal operation is needed
        one_sided_length<-ceiling(incremental_window_size)
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
                between_dist_rev<-sum(partially_rev_dist[(endID-one_sided_length):endID,(endID+1):(endID+one_sided_length)])

                # Original local window: sum of distance across neighboring windows
                between_dist<-sum(temp_combined_dist[(endID-one_sided_length):endID,(endID+1):(endID+one_sided_length)])

                if(between_dist>between_dist_rev) {
                    temp_combined_ordering<-temp_combined_ordering[idx,,drop=F]
                    temp_combined_dist[idx,idx]
                }
            }

            # All other windows are compared to the preceeding LHS window
            else{
                # Reversed local window reversed: sum of distance across neighboring windows
                between_dist_rev<-sum(partially_rev_dist[(startID-one_sided_length):startID,startID:(startID+one_sided_length)])
                # Original local window: sum of distance across neighboring windows
                between_dist<-sum(temp_combined_dist[(startID-one_sided_length):startID,startID:(startID+one_sided_length)])

                if(between_dist>between_dist_rev) {
                    temp_combined_ordering<-temp_combined_ordering[idx,,drop=F]
                    temp_combined_dist[idx,idx]
                }
            }

        }#end i


        temp_combined_distVar<-c2c_dist_uneveness_wrapper(temp_combined_ordering, alpha=alpha,
                                                          data_type=data_type)


        list_distVar<-c(list_distVar,temp_combined_distVar)

        # Update global ordering if local sorting achieve better fine-grained smoothness
        if(temp_combined_distVar<lowest_distVar){
            best_local_ordering<-temp_combined_ordering
            lowest_distVar<-temp_combined_distVar
        }

        incremental_window_size<-incremental_window_size + ceiling(n*window_size_incre_perct)
    }# end while loop
    cat('\n')

    pdf(paste0('window.size.optimization.',baseNm, ".pdf"))
    plot (list_windowSize, list_distVar, xlab='Window size (no. of cells)',
          ylab='Smoothness of local C2C distance', type='l', main='Finding optimal window size')
    points(list_windowSize, list_distVar, pch=20, col=ifelse(list_distVar==min(list_distVar), 'red','black'))
    dev.off()
    cat(list_distVar, '\n')
    if(exists("best_local_ordering")){global_ordering<-best_local_ordering}
    ordering <- data.frame('SampleID'=rownames(global_ordering), 'GroupID' = rep('untitled',n))

    ## Printout cost value using preliminary genes
    exp.cost <- fluidigmSC::readExpObject('preprocessed_exp.fso')
    exp.cost <- fluidigmSC::updateGeneListFromFile(exp.cost, gene_list_file='pca.genes.txt')
    exp.cost <- fluidigmSC::updateSampleListFromList(exp.cost, sample_list = ordering)
    if(sorting_method == 'STS') {
        cost <- STS_sortingcost(expr = t(exp.cost$log2ex_data))
    }else{
        cost <- neighborhood_sortingcost(expr = t(exp.cost$log2ex_data), sigma_width=1)
    }
    cat('autoSPIN sorting cost (PCA genes) @',baseNm,' = ', cost, '\n')

    return(ordering)

}


