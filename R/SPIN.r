#' A wrapper function for SPIN sorting method
#'
#' A wrapper function for SPIN method  provides a R version of SPIN [Tsafrir et al. 2005].
#'
#' @param data An log2 transformed expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param sorting_method A character string indicating the choice of sorting method, i.e. 'STS' (Side-to-Side)
#' or 'Neighborhood'.
#' @param sigma_width An integer number determining the degree of spread of the gaussian distribution which
#' is used for computing weight matrix for Neighborhood sorting method.
#'
#' @return A data frame containing single column of ordered sample IDs.
#' @export
#'
#' @examples
#' set.seed(15)
#' da <- iris[sample(150, 150, replace = FALSE), ]
#' rownames(da) <- paste0('spl_',seq(1,nrow(da)))
#' d <- da[,1:4]
#' dl <- da[,5,drop=FALSE]
#' res <- SPIN(data = d)
#' dl <- dl[match(res$SampleID,rownames(dl)),]
#' annot <- data.frame(id = seq(1,nrow(res)), label=dl, stringsAsFactors = FALSE)
#' #ggplot(annot, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
SPIN <- function(data, sorting_method = c("STS", "neighborhood"),
    sigma_width = 1) {

    n <- nrow(data)
    sorting_method <- match.arg(sorting_method)
    switch(sorting_method, STS = {
        global_res <- STS_sorting_wrapper(data, no_randomization = 1)
    }, neighborhood = {
        global_res <- neighborhood_sorting_wrapper(data,
                                                   no_randomization = 1,
                                                   sigma_width = sigma_width)
    })

    ordering <- global_res$permutated.expr
    return(data.frame(SampleID = rownames(ordering),
                      GroupID = rep("untitled",
        n)))
}



#' A wrapper function for Side-to-Side (STS) sorting.
#'
#' A wrapper function for Side-to-Side (STS) sorting as proposed in [Tsafrir et al. 2005].
#'
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param no_randomization An integer number indicating the number of repeated sorting, each
#' of which uses a randomaly selected initial cell ordering.
#'
#' @return A list containing \code{permutated.expr}(data frame) and \code{best.cost}(a numeric value).
STS_sorting_wrapper <- function(expr, no_randomization = 10) {
    ## Define parameters
    n <- nrow(expr)
    ## Initialization
    best.cost <- 1e+20
    set.seed(123)
    for (i in 1:no_randomization) {

        if (i == 1)
            temp.expr <- expr else {
            rand_id <- sample.int(n)
            temp.expr <- expr[rand_id, , drop = FALSE]
        }
        temp.dist <- distance.function(temp.expr)
        temp.SPIN.res <- STS_sorting(temp.dist)

        if (temp.SPIN.res$cost < best.cost) {
            temp.best.ordering <- temp.expr[temp.SPIN.res$ordering,
                , drop = FALSE]
            best.ordering <- temp.best.ordering
            best.cost <- temp.SPIN.res$cost
        }

    }  #end for-loop

    return(list(permutated.expr = best.ordering, best.cost = best.cost))
}



#' A sorting function using the Side-to-Side (STS) algorithm
#'
#' @param d A matrix containing n-by-n cell distance.
#' @param max_iter An integer number indicating the maximum number of iteration if sorting
#' does not converge.
#'
#' @return A list containing \code{ordering}(a vector of re-ordered sequence) and
#' \code{cost}(a numeric value).
STS_sorting <- function(d, max_iter = 10) {
    ## Define parameters & constant
    n <- nrow(d)
    I <- diag(n)
    # a strictly increasing column vector
    X <- matrix(ncol = 1, nrow = n)
    for (i in 1:n) {
        X[i] <- i - (n + 1)/2
    }

    ## STS algorithm Step 1
    t <- 0
    Dt <- d
    global_permutation <- I
    prev_permutation <- I
    cost <- sum(diag(I %*% Dt %*% t(I) %*% X %*% t(X)))

    flag <- "on"

    while (flag == "on" && t < max_iter) {

        # Step 2
        St <- Dt %*% X

        # Step 3
        descending.order <- order(St, decreasing = TRUE)
        Pt <- I[descending.order, ]
        global_permutation <- global_permutation[descending.order,
            ]
        cost <- sum(diag(global_permutation %*% d %*%
                             t(global_permutation) %*%
                             X %*% t(X)))

        ## Step 4
        if (isTRUE(all.equal((Pt %*% St), (prev_permutation %*%
            St)))) {
            flag <- "off"
        } else {
            Dt <- Dt[descending.order, descending.order]
            t <- t + 1
            prev_permutation <- Pt
        }
    }  #end while-loop


    ordering <- apply(global_permutation, 1, which.max)
    return(list(ordering = ordering, cost = cost))
}



#' A cost computation function for Side-to-Side (STS) algorithm
#'
#' @param method A character string indicating the distance function.
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#'
#' @return A numeric value of sorting cost.
#' @export
#'
#' @examples
#' set.seed(15)
#' da <- iris[sample(150, 150, replace = FALSE), ]
#' d <- da[,1:4]
#' randomOrdering_cost <- STS_sortingcost(d, method= 'Euclidean')
#' randomOrdering_cost
#'
#' da <- iris
#' d <- da[,1:4]
#' properOrdering_cost <- STS_sortingcost(d, method= 'Euclidean')
#' properOrdering_cost
STS_sortingcost <- function(expr = NULL, method = c("Euclidean",
    "Correlation", "eJaccard", "none")) {
    ## Define parameters & constant
    n <- nrow(expr)
    I <- diag(n)

    # a strictly increasing column vector
    X <- matrix(ncol = 1, nrow = n)
    for (i in 1:n) {
        X[i] <- i - (n + 1)/2
    }

    d <- distance.function(expr, method = method)

    cost <- sum(diag(I %*% d %*% t(I) %*% X %*% t(X)))

    return(cost)
}



#' A wrapper function for Neighborhood sorting.
#'
#' A wrapper function for Neighborhood sorting as proposed in [Tsafrir et al. 2005].
#'
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param sigma_width An integer number determining the degree of spread of the gaussian distribution which is used for computing weight matrix for Neighborhood sorting method.
#' @param no_randomization An integer number indicating the number of repeated sorting, each of which uses a randomaly selected initial cell ordering.
#'
#' @return A list containing \code{permutated.expr}(data frame) and \code{best.cost}(a numeric value).
neighborhood_sorting_wrapper <- function(expr, sigma_width = 1,
    no_randomization = 10) {

    ## Initialization
    best.cost <- 1e+20
    set.seed(123)
    mat_size <- nrow(expr)
    s <- seq(mat_size)
    m <- rep(s, mat_size)
    i <- matrix(m, ncol = mat_size, nrow = mat_size)
    i <- t(i)
    j <- matrix(m, ncol = mat_size, nrow = mat_size)

    G = exp(-(i - j)^2/sigma_width/mat_size)

    for (i in 1:10) {
        G <- G/(matrix(rep(colSums(G), mat_size), ncol = mat_size,
            nrow = mat_size, byrow = TRUE))
        G <- G/(matrix(rep(rowSums(G), mat_size), ncol = mat_size,
            nrow = mat_size, byrow = FALSE))
    }
    G <- (G + t(G))/2
    weights_mat <- G

    for (i in 1:no_randomization) {

        if (i == 1)
            temp.expr <- expr else {
            rand_id <- sample.int(mat_size)
            temp.expr <- expr[rand_id, , drop = FALSE]
        }
        temp.dist <- distance.function(temp.expr)
        temp.SPIN.res <- neighborhood_sorting(temp.dist, weights_mat)

        if (temp.SPIN.res$cost < best.cost) {
            temp.best.ordering <- temp.expr[temp.SPIN.res$ordering,
                , drop = FALSE]
            best.ordering <- temp.best.ordering
            best.cost <- temp.SPIN.res$cost
        }

    }  #end for-loop

    return(list(permutated.expr = best.ordering, best.cost = best.cost))
}



#' A sorting function using the Neighborhood algorithm
#'
#' @param d A matrix containing n-by-n cell distance.
#' @param weights_mat A weight matrix of size n-by-n.
#' @param max_iter An integer number indicating the maximum number of iteration if sorting
#' does not converge.
#'
#' @return A list containing \code{ordering}(a vector of re-ordered sequence) and
#' \code{cost}(a numeric value).
neighborhood_sorting <- function(d, weights_mat = NULL, max_iter = 100) {

    #
    mat_size <- nrow(d)
    I <- diag(mat_size)
    t <- 0
    Dt <- d
    global_permutation <- I
    prev_permutation <- I
    cost <- sum(diag(I %*% Dt %*% t(I) %*% weights_mat))

    flag <- "on"

    while (flag == "on" && t < max_iter) {

        ## ============
        mismatch = Dt %*% weights_mat
        val <- Biobase::rowMin(mismatch)
        mn <- apply(mismatch, 1, which.min)
        main_diag <- diag(mismatch)
        sort_score <- seq(mat_size)
        mx <- max(val)
        sort_score <- mn - 0.1 * sign((mat_size/2 - mn)) * val/(mx)

        # Sorting the matrix
        sorted_ind <- order(sort_score)
        val <- sort_score[sorted_ind]

        # update of ordering
        Pt <- I[sorted_ind, ]
        global_permutation <- global_permutation[sorted_ind, ]
        cost <- sum(diag(global_permutation %*% d %*%
                             t(global_permutation) %*%
                             weights_mat))
        # cat('@t:', t, ' cost=', cost, '\n') ============

        ##
        if (isTRUE(all.equal((Pt %*% mismatch), (prev_permutation %*%
            mismatch)))) {
            flag <- "off"
        } else {
            Dt <- Dt[sorted_ind, sorted_ind]
            t <- t + 1
            prev_permutation <- Pt
        }
    }  #end while-loop

    # cat(val)
    ordering <- apply(global_permutation, 1, which.max)
    return(list(ordering = ordering, cost = cost))
}




#' A cost computation function for Neighborhood algorithm
#'
#' @param expr An expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param sigma_width An integer number determining the degree of spread of the gaussian distribution which
#' is used for computing weight matrix for Neighborhood sorting method.
#' @param method A character string indicating the distance function.
#'
#' @return A numeric value of sorting cost.
#' @export
#'
#' @examples
#' set.seed(15)
#' da <- iris[sample(150, 150, replace = FALSE), ]
#' d <- da[,1:4]
#' randomOrdering_cost <- neighborhood_sortingcost(d, method= 'Euclidean')
#' randomOrdering_cost
#'
#' da <- iris
#' d <- da[,1:4]
#' properOrdering_cost <- neighborhood_sortingcost(d, method= 'Euclidean')
#' properOrdering_cost
neighborhood_sortingcost <- function(expr = NULL, sigma_width = 1,
    method = c("Euclidean", "Correlation", "eJaccard", "none")) {
    ## Define parameters & constant
    mat_size <- nrow(expr)
    s <- seq(mat_size)
    m <- rep(s, mat_size)
    i <- matrix(m, ncol = mat_size, nrow = mat_size)
    i <- t(i)
    j <- matrix(m, ncol = mat_size, nrow = mat_size)

    G = exp(-(i - j)^2/sigma_width/mat_size)

    for (i in 1:10) {
        G <- G/(matrix(rep(colSums(G), mat_size), ncol = mat_size,
            nrow = mat_size, byrow = TRUE))
        G <- G/(matrix(rep(rowSums(G), mat_size), ncol = mat_size,
            nrow = mat_size, byrow = FALSE))
    }
    G <- (G + t(G))/2
    weights_mat <- G

    I <- diag(mat_size)

    d <- distance.function(expr)

    cost <- sum(diag(I %*% d %*% t(I) %*% weights_mat %*% t(weights_mat)))

    return(cost)
}
