#' sWanderlust
#'
#' autoSPIN guided wanderlust. Specifically, we use autoSPIN to help find the starting point for wanderlust.
#'
#' @param data data Input data matrix.
#' @param data_type The data type which guides the autoSPIN sorting, including \code{linear}, \code{cyclical}.
#' @param SPIN_option SPIN contains two options including \code{STS}(default), \code{neighborhood}.
#' @param alpha alpha parameter for autoSPIN, default is 0.2.
#' @param sigma_width Sigma width parameter for SPIN, default is 1.
#' @param diffusionmap_components Number of components from diffusion map used for wanderlust analysis, default is 4.
#' @param l Number of nearest neighbors, default is 15.
#' @param k Number of nearest neighbors for repeating graphs, default is 15, should be less than or equal to l.
#' @param num_waypoints Number of waypoint used for wanderlust, default is 150.
#' @param flock_waypoints The number of times for flocking the waypoints, default is 2.
#' @param waypoints_seed The seed for reproducing the results.
#'
#' @export
#'
#' @author Hao Chen
#' @return a vector of the sorted oder.
#' @examples
#' set.seed(15)
#' shuffled_iris <- iris[sample(150, 150, replace = FALSE), ]
#' data <- shuffled_iris[,1:4]
#' data_label <- shuffled_iris[,5]
#' wishbone <- sWanderlust(data = data, num_waypoints = 100)
sWanderlust <- function(data, data_type = c("linear", "cyclical"),
    SPIN_option = c("STS", "neighborhood"), alpha = 0.2, sigma_width = 1,
    diffusionmap_components = 4, l = 15, k = 15, num_waypoints = 150,
    flock_waypoints = 2, waypoints_seed = 2711) {

    ## autoSPIN to find the starting point
    data_type <- match.arg(data_type)
    SPIN_option <- match.arg(SPIN_option)
    autoSPIN_order <- autoSPIN(data, data_type = data_type,
                               sorting_method = SPIN_option,
                               sigma_width = sigma_width,
                               alpha = alpha)

    ## try both ends as starting point with wanderlust and find the
    ## best wanderlust preprocess
    medColSums <- median(rowSums(data))
    data <- t(apply(data, 1, function(x) x/sum(x) * medColSums))
    nonzero_genes <- colSums(data) > 0
    data <- data[, nonzero_genes]
    df <- diffusionmap(data, data_type = "scSeq")
    dm <- df$diffusion_eigenvectors

    # Forward first cell
    startingCell1 <- as.character(autoSPIN_order$SampleID[1])
    cat(paste0("  Try Starting cell: ", startingCell1))
    wishbone1 <- Rwanderlust(data = dm[, 2:(diffusionmap_components +
        1)], s = startingCell1, l = l, k = k, num_waypoints = num_waypoints,
        flock_waypoints = flock_waypoints, waypoints_seed = waypoints_seed)

    wanderlust_ordering1 <- wishbone1$Order
    wanderlust_exp1 <- data[wanderlust_ordering1, ]
    c1 <- STS_sortingcost(expr = wanderlust_exp1)

    # Reverse first cell
    startingCell2 <- as.character(autoSPIN_order$SampleID[nrow(autoSPIN_order)])
    cat(paste0("  Try Starting cell: ", startingCell2))
    wishbone2 <- Rwanderlust(data = dm[, 2:(diffusionmap_components +
        1)], s = startingCell2, l = l, k = k, num_waypoints = num_waypoints,
        flock_waypoints = flock_waypoints, waypoints_seed = waypoints_seed)

    wanderlust_ordering2 <- wishbone2$Order
    wanderlust_exp2 <- data[wanderlust_ordering2, ]
    c2 <- STS_sortingcost(expr = wanderlust_exp2)

    # Pick an optimal ordering
    if (c1 < c2) {
        sWanderlust_order <- wanderlust_ordering1
        cat("Starting cell ", startingCell1, ", cost=",
            c1, "(driver genes)\n")
    } else {
        sWanderlust_order <- wanderlust_ordering2
        cat("Starting cell ", startingCell2, ", cost=",
            c2, "(driver genes)\n")
    }

    return(sWanderlust_order)
}

#' a wrapper of wanderlust for sWanderlust
#'
#' @param data Input data matrix.
#' @param s The ID of starting point.
#' @param diffusionmap_components Number of components from diffusion map used for wanderlust analysis, default is 4.
#' @param l Number of nearest neighbors, default is 15.
#' @param k Number of nearest neighbors for repeating graphs, default is 15, should be less than or equal to l.
#' @param num_graphs Number of repreated graphs.
#' @param num_waypoints Number of waypoint used for wanderlust, default is 150.
#' @param waypoints_seed The seed for reproducing the results.
#' @param flock_waypoints The number of times for flocking the waypoints, default is 2.
#'
#' @return sorted order.
#'
#' @author Hao Chen
wanderlust_wrapper <- function(data, s, diffusionmap_components = 4,
    l = 15, k = 15, num_graphs = 1, num_waypoints = 150,
    waypoints_seed = 123,
    flock_waypoints = 2) {
    ## wanderlust preprocess
    if (missing(s) || is.null(s)) {
        message("No starting point specified, the first row will be used!")
        s <- row.names(data)[1]
    }

    if (is.numeric(s)) {
        s <- row.names(data)[s]
    }
    cat(paste0("  Starting cell set to be: ", s))
    medColSums <- median(rowSums(data))
    data <- t(apply(data, 1, function(x) x/sum(x) * medColSums))
    nonzero_genes <- colSums(data) > 0
    data <- data[, nonzero_genes]
    df <- diffusionmap(data, data_type = "scSeq")
    dm <- df$diffusion_eigenvectors

    res <- Rwanderlust(data = dm[, 2:(diffusionmap_components + 1)],
                       s = s, l = l, k = k, num_graphs = num_graphs,
                       num_waypoints = num_waypoints,
                       waypoints_seed = waypoints_seed,
                       flock_waypoints = flock_waypoints)

    return(res)
}


#' R inplementation of wanderlust
#'
#' @param data Input data matrix.
#' @param s Starting point ID.
#' @param l l nearest neighbours.
#' @param k k nearest neighbours, k < l.
#' @param num_waypoints Number of waypoints to guide the trajectory detection.
#' @param flock_waypoints The number of times for flocking the waypoints, default is 2.
#' @param num_graphs Number of repreated graphs.
#' @param waypoints_seed The seed for reproducing the results.
#' @param metric Distance calculation metric for nearest neighbour detection.
#' @param voting_scheme The scheme of voting.
#' @param band_sample Boolean, if band the sample
#' @param partial_order default NULL
#' @param verbose Boolean, if print the details
#'
#' @return a list containing Trajectory, Order, Waypoints
#' @author Hao Chen
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' set.seed(15)
#' shuffled_iris <- iris[sample(150, 150, replace = FALSE), ]
#' data <- shuffled_iris[,1:4]
#' data_label <- shuffled_iris[,5]
#' wishbone <- Rwanderlust(data = data, num_waypoints = 100, waypoints_seed = 2)
#' pd1 <- data.frame(id = wishbone$Trajectory, label=data_label, stringsAsFactors = FALSE)
#' pd2 <- data.frame(id = seq_along(row.names(data)), label=data_label, stringsAsFactors = FALSE)
#' #ggplot(pd1, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
#' #ggplot(pd2, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
Rwanderlust <- function(data, s, l = 15, k = 15, num_graphs = 1,
    num_waypoints = 250, waypoints_seed = 123, flock_waypoints = 2,
    metric = "euclidean", voting_scheme = "exponential",
    band_sample = FALSE,
    partial_order = NULL, verbose = TRUE) {
    if (missing(s) || is.null(s)) {
        message("No starting point specified, the first row will be used!")
        s <- row.names(data)[1]
    }

    if (is.numeric(s)) {
        s <- row.names(data)[s]
    }

    if (verbose)
        cat("  Building lNN graph...\n")

    ## Construct nearest neighbour graph
    lnn_time <- system.time({
        lnn <- nn2(data, data, l + 1, searchtype = "standard")
    })

    if (verbose)
        cat(paste0("    lNN computed in: ", round(lnn_time[3],
            2), " seconds\n"))

    ## generate klNN graphs and iteratively refine a trajectory in
    ## each
    trajectory <- NULL
    for (graph_iter in seq_len(num_graphs)) {
        if (k != l) {
            klnn <- spdists_klnn(lnn, k)  # randomly select k neighbours
        } else {
            # remove self neighbour
            klnn <- list()
            klnn$nn.idx <- lnn$nn.idx[, -1]
            klnn$nn.dists <- lnn$nn.dists[, -1]
        }

        ## transform klnn to a sparse matrix object
        sm_i <- rep(1:nrow(klnn$nn.idx), each = ncol(klnn$nn.idx))
        sm_j <- as.vector(t(klnn$nn.idx))
        klnn <- sparseMatrix(i = sm_i, j = sm_j,
                             x = as.vector(t(klnn$nn.dists)))
        ## Make the graph undirected
        klnn[cbind(sm_j, sm_i)] <- klnn[cbind(sm_i, sm_j)]

        ## Detect landmarks and trajectory
        traj_l <- trajectory_landmarks(klnn, data, s,
                                       partial_order = partial_order,
            waypoints = num_waypoints, waypoints_seed = waypoints_seed,
            metric = metric, flock_waypoints = flock_waypoints,
            band_sample = band_sample)

        traj <- traj_l[[1]]
        dist <- traj_l[[2]]
        landmarks <- traj_l[[3]]

        ## calculate weighed trajectory
        W = weighting_scheme(dist, voting_scheme)

        ## save initial solution - start point's shortest path distances
        t <- list(traj[1, ])
        t <- c(t, list(colSums(traj * W)))

        ## iteratively realign trajectory (because landmarks moved)
        if (verbose)
            cat("  Iteratively realign trajectory: \n")
        converged <- FALSE
        user_break <- FALSE
        realign_iter <- 1
        while (converged == FALSE && user_break == FALSE) {
            if (verbose)
                cat(paste0("    Running iterations... ", realign_iter,
                  "\n"))
            realign_iter <- realign_iter + 1
            traj <- do.call(rbind, dist)
            traj <- realign_trajectory(t, dist, landmarks, traj,
                1, length(dist), realign_iter)
            ## calculate weighed trajectory
            t <- c(t, list(colSums(traj * W)))
            ## check for convergence
            fpoint_corr <- cor(t[[realign_iter + 1]], t[[realign_iter]],
                method = "pearson")
            if (verbose)
                cat(paste0("      Correlation with previous iteration: ",
                  round(fpoint_corr, 6), "\n"))
            converged <- fpoint_corr > 0.9999

            if (realign_iter > 20) {
                # break after too many alignments - something is wrong
                user_break <- TRUE
                warning(paste("Force exit after", realign_iter,
                  "iterations.\n"))
            }

        }

        iter_traj <- t[[realign_iter + 1]]  ## +1 here
        iter_traj <- (iter_traj - min(iter_traj))/(max(iter_traj) -
            min(iter_traj))
        trajectory <- rbind(trajectory, iter_traj)
    }

    # Normalize the iter_trajectory
    if (verbose)
        cat("  Wanderlust Sorting Done!\n")

    ## average the graph num_graphs
    trajectory <- colSums(trajectory)/nrow(trajectory)
    order <- row.names(data)[order(trajectory)]
    waypoints <- row.names(data)[landmarks]

    return(list(Trajectory = trajectory, Order = order,
                Waypoints = waypoints))
}


# Randomly removing l-k edges for each iteration of graph_num
spdists_klnn <- function(lnn, k) {
    l <- ncol(lnn$nn.dists)
    n <- nrow(lnn$nn.dists)
    if (k > l - 1)
        stop("k is too big, more than the number of nearest neighbors!")

    neighbourDice <- seq_len(l)[-1]  # remove self neighbour
    sampleID <- do.call(rbind, lapply(1:n, function(i) {
        cbind(i, sample(neighbourDice, k, replace = FALSE))
    }))
    klnn <- list()
    klnn$nn.idx <- matrix(lnn$nn.idx[sampleID], byrow = TRUE, ncol = k)
    klnn$nn.dists <- matrix(lnn$nn.dists[sampleID], byrow = TRUE,
        ncol = k)
    return(klnn)
}


#' determining initial trajectory and landmarks
#'
#' @param knn A sparse matrix of knn.
#' @param data data.
#' @param s The ID of starting point.
#' @param partial_order A vector of IDs specified as recommended waypoints, NULL to ignore.
#' @param waypoints Either the number of waypoints, or specify the waypoint IDs.
#' @param waypoints_seed Random sampling seed, for reproducible results.
#' @param metric Distance calculation metric for nearest neighbour detection.
#' @param flock_waypoints Iteration of using nearest points around waypoint to adjust its position.
#' @param band_sample if give more chance to nearest neighbours of starting point in randomly waypoints selection.
#'
#' @return a list
#' @importFrom plyr laply
#' @importFrom igraph graph_from_adjacency_matrix shortest.paths E
trajectory_landmarks <- function(knn, data, s, partial_order = NULL,
    waypoints = 250, waypoints_seed = 123, metric = "euclidean",
    flock_waypoints = 2, band_sample = FALSE) {

    ## randomly select one starting point
    if (length(s) > 1)
        s <- sample(s, 1)
    if (is.character(s)) {
        s <- match(s, row.names(data))
    } else {
        stop("wrong starting cell!")
    }

    dijkstra <- knn[s, ]  ##** different from python code

    if (length(waypoints) == 1 && is.numeric(waypoints)) {
        cat("    Randomly Select waypoints...\n")
        n_opts <- seq_len(nrow(data))
        ## ?? band_sample doesn't make sense for me
        if (band_sample) {
            cat("    Band samples...\n")
            n_opts <- NULL
            window_size <- 0.1
            prc <- 0.998
            max_dist <- max(dijkstra)
            while (prc > 0.08) {
                band <- (dijkstra >= ((prc - window_size) * max_dist)) &
                  (dijkstra <= prc * max_dist)
                n_opts <- c(n_opts, sample(which[band], min(sum(band),
                  waypoints - length(partial_order)), replace = FALSE))
                prc <- prc - window_size
            }
        }

        ## use FlowSOM to aid the waypoints selection
        if (is.numeric(waypoints_seed))
            set.seed(waypoints_seed)
        waypoints <- sample(n_opts, waypoints - length(partial_order),
            replace = FALSE)
    }

    # Flock wayoints
    if (flock_waypoints) {
        cat("    Flock waypoints...\n")
        for (f in seq_len(flock_waypoints)) {
            n_neighbors <- 20
            idx <- nn2(data, data[waypoints, , drop = FALSE], n_neighbors,
                searchtype = "standard")$nn.idx
            for (i in seq_along(waypoints)) {
                med_data = matrix(apply(data[idx[i, ], ], 2, median),
                  nrow = 1)
                waypoints[i] <- nn2(data, med_data, 1)$nn.idx[1,
                  1]
            }
        }
    }

    if (!(s %in% partial_order)) {
        ## partial_order includes start point
        partial_order <- c(s, partial_order)
    }

    ## add extra landmarks if user specified l <-
    ## unique(c(partial_order, waypoints)) ## start with starting
    ## point
    l <- c(partial_order, waypoints)
    # print(table(l))

    ## calculate all shortest paths
    g <- graph_from_adjacency_matrix(knn, mode = "undirected",
        weighted = TRUE)
    dist <- list()
    for (i in seq_len(length(l))) {
        sp_dist <- shortest.paths(g, v = l[i], weights = E(g)$weight)[1,
            ]
        unreachable <- is.infinite(sp_dist)
        sp_dist[unreachable] <- max(c(sp_dist[!unreachable], laply(dist,
            max)))
        dist[[i]] <- sp_dist
    }

    ## adjust paths according to partial order by redirecting leave
    ## out here with partial_order=NULL
    traj <- do.call(rbind, dist)
    if (length(l) > length(partial_order)) {
        traj <- realign_trajectory(dist, dist, l, traj, length(partial_order) +
            1, length(l), 1)
    }

    return(list(traj, dist, l))
}


realign_trajectory <- function(t, dist, landmarks, traj, start_range,
    end_range, realign_iter) {
    traj_distVector <- t[[realign_iter]]
    for (i in c(start_range:end_range)) {
        ## find position of landmark in previous iteration
        dist_to_waypoint_i <- traj_distVector[landmarks[i]]
        before_indices <- which(traj_distVector < dist_to_waypoint_i)
        if (length(before_indices) > 0) {
            ## convert all cells before starting point to the negative
            traj[i, before_indices] = -dist[[i]][before_indices]
            # set zero to position of starting point
            traj[i, ] = traj[i, ] + dist_to_waypoint_i
        }
    }

    traj <- traj - min(traj)
    return(traj)
}



weighting_scheme <- function(dist, voting_scheme = c("uniform",
    "exponential", "linear")) {
    voting_scheme <- match.arg(voting_scheme)
    W_full <- do.call(rbind, dist)
    switch(voting_scheme, uniform = {
        W_full[] <- 1
    }, exponential = {
        std <- apply(W_full, 2, function(x) {
            sqrt(sum((x - mean(x))^2)/length(x))
        })
        sdv <- mean(std) * 3
        W_full <- exp(-0.5 * (W_full/sdv)^2)
    }, linear = {
        W_full <- max(W_full) - W_full
    })

    ## The weighing matrix must be a column stochastic operator
    W_full <- apply(W_full, 2, function(x) {
        x/sum(x)
    })
    return(W_full)
}


#' @importFrom RSpectra eigs_sym
diffusionmap <- function(data, data_type = c("scSeq", "massCyt"),
    knn = 10, epsilon = 1, n_diffusion_components = 10,
    n_pca_components = 15,
    markers = NULL) {
    cat("  Runing Diffusion Map...")
    data_type <- match.arg(data_type)
    if (n_pca_components > min(dim(data))) {
        n_pca_components <- min(dim(data))
    }
    switch(data_type, scSeq = {
        data_pca <- pca(data, n_components = n_pca_components)$loadings
        ## center and scale looks weird here but just follow the python
        ## implementation
        data <- data - min(data)
        data <- data/max(data)
        data <- data %*% data_pca[, 1:n_pca_components]
    }, massCyt = {
        if (is.null(markers)) markers <- colnames(data)
        data <- data[, markers]
    })
    ## Nearest neighbors
    N <- nrow(data)
    nbrs <- nn2(data, data, knn, searchtype = "standard")

    ## transform klnn to a sparse adjacency matrix
    sm_i <- rep(1:nrow(nbrs$nn.idx), each = ncol(nbrs$nn.idx))
    sm_j <- as.vector(t(nbrs$nn.idx))
    dist <- as.vector(t(nbrs$nn.dists))
    # Symmetrize W and convert to affinity (with selfloops)
    W <- sparseMatrix(i = sm_i, j = sm_j, x = dist)
    n_i <- NULL
    n_j <- NULL
    n_x = NULL
    for (i in 1:nrow(W)) {
        for (j in 1:ncol(W)) {
            if (i != j && (W[i, j] + W[j, i]) > 0) {
                n_i <- c(n_i, i)
                n_j <- c(n_j, j)
                n_x <- c(n_x, W[i, j] + W[j, i])
            } else if (i == j) {
                n_i <- c(n_i, i)
                n_j <- c(n_j, j)
                n_x <- c(n_x, 0)
            }
        }
    }
    W <- sparseMatrix(i = n_i, j = n_j, x = exp(-n_x/(epsilon^2)))

    # create D
    D <- Matrix::colSums(W)
    D[D != 0] = 1/D[D != 0]
    # Symmetric markov normalization
    D <- sparseMatrix(i = 1:N, j = 1:N, x = sqrt(D))
    P <- D
    TT <- D %*% W %*% D
    TT <- (TT + Matrix::t(TT))/2

    # Eigen value decomposition
    DV <- eigs_sym(TT, k = n_diffusion_components, opts = list(tol = 1e-04,
        maxiter = 1000, retvec = TRUE))

    D <- DV$values
    V <- DV$vectors
    inds <- order(D, decreasing = TRUE)
    D <- D[inds]
    V <- V[, inds]
    V = P %*% V

    # Normalize, The Frobenius norm
    V <- apply(V, 2, function(x) x/sqrt(sum(x^2)))

    row.names(V) <- row.names(data)
    return(list(diffusion_eigenvectors = V, diffusion_eigenvalues = D))
}


pca <- function(x, n_components = 100) {
    # Make sure data is zero mean
    x <- (x - min(x))/max(x)
    if (ncol(x) < nrow(x)) {
        c <- cov(x)
    } else {
        # if #gene > #cell, we better use this matrix for the eigen
        # decomposition
        c <- x %*% t(x)/nrow(x)
    }
    r <- eigen(c)
    l <- r$values
    m <- -r$vectors

    # Sort eigenvectors in descending order
    ind <- order(l, decreasing = TRUE)
    if (n_components < 1) {
        n_components <- which(cumsum(l/sum(l)) >= n_components)[1]
        cat(paste0("Embedding into ", n_components, " dimensions.\n"))
    }
    if (n_components > ncol(m)) {
        n_components <- ncol(m)
        cat(paste0("Target dimensionality reduced to ", n_components,
            ".\n"))
    }

    m <- m[, ind[1:n_components], drop = FALSE]
    l <- l[ind[1:n_components]]

    # Apply mapping on the data
    if (ncol(x) >= nrow(x)) {
        m <- t(x) %*% m/matrix(sqrt(nrow(x) * l), byrow = TRUE,
            nrow = nrow(t(x)), ncol = ncol(m))
    }

    return(list(loadings = m, eigenvalues = l))
}
















