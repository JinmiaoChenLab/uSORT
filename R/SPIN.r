#' A SPIN sorting wrapper function
#'
#' A SPIN sorting wrapper function provides a R version of SPIN [Tsafrir et al. 2005].
#'
#' @param expr an expresssion matrix containing n-rows of cells and m-cols of genes
#' @param baseNm a base name used for saving output plot
#' @param sorting_method character specifying the type of sorting algorithm, i.e. STS or neighborhood
#' @param sigma_width width size used for computing weight matrix in Neighborhood sorting method
#' @return a matrix of the permutated expression matrix
#' @author SPIN
#' @export
SPIN <- function(expr, baseNm = 'SPIN', sorting_method = NULL, sigma_width = 1){

    n <- nrow(expr)

    if(sorting_method == 'STS'){
        global_res <- STS_sorting_wrapper(expr, no_randomization= 1)
    }else if(sorting_method == 'neighborhood') {
        global_res<-neighborhood_sorting_wrapper(expr, no_randomization= 1, sigma_width = sigma_width)
    }else {
        stop('Sorting type can onlyl be STS or neighborhood')
    }
    ordering<-global_res$permutated.expr

    return(data.frame('SampleID'=rownames(ordering), 'GroupID' = rep('untitled',n)))
}


