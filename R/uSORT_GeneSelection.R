#' Gene selection using PCA technique
#'
#' @param data A matrix of data.frame with row.name of cells, and col.name of genes
#'
#' @return a vector of the names of selected genes.
#' @export
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' exprs <- uSORT_preProcess(exprs_file = file)
#' exp_trimmed <- t(exprs$exprs_log_trimed)
#' PCA_selected_genes <- pca_gene_selection(exp_trimmed)
pca_gene_selection <- function(data) {
    data <- as.matrix(data)
    pca <- prcomp(data)
    ## select PCs
    eigen_values <- pca$sdev^2
    select_PCs <- elbow_detection(scores = eigen_values)
    cat("  No. of PC selected = ", length(select_PCs), "\n")
    ## select genes
    rotation_matrix <- as.data.frame(abs(pca$rotation[, select_PCs,
        drop = FALSE]))
    gene_scores <- rowSums(rotation_matrix)
    PCA_genes <- elbow_detection(scores = gene_scores)
    cat("  No. of PCA genes = ", length(PCA_genes), "\n")

    return(as.character(colnames(data)[PCA_genes]))
}


#' A elbow detection function
#'
#' A elbow detection function detects the elbow/knee of a given vector of values.
#' Values will be sorted descendingly before detection, and the ID of those values
#' above the elbow will be returned.
#'
#' @param scores A vector of numeric scores.
#' @param if_plot Boolean determine if plot the results.
#' @return a vector of selected elements IDs
#'
#' @export
#' @examples
#' scores <- c(10, 9 ,8, 6, 3, 2, 1, 0.1)
#' elbow_detection(scores, if_plot = TRUE)
elbow_detection <- function(scores, if_plot = FALSE) {

    num_scores <- length(scores)
    if (num_scores < 2) {
        stop("Input scores must be a vector with length more than 1!")
    }
    scores <- data.frame(id = seq_len(num_scores), value = scores)
    sorted_scores <- scores[order(scores$value, decreasing = TRUE),
        ]

    ## use distance to diagonal line to determine the elbow point
    xy_coordinates <- cbind(x = seq_len(num_scores), y = sorted_scores$value)
    start_point <- xy_coordinates[1, ]
    end_point <- xy_coordinates[num_scores, ]
    x1 <- start_point[1]
    x2 <- end_point[1]
    y1 <- start_point[2]
    y2 <- end_point[2]
    a <- y1 - y2
    b <- x2 - x1
    c <- x1 * y2 - x2 * y1
    dist_to_line <- abs(a * xy_coordinates[, "x"] + b * xy_coordinates[,
        "y"] + c)/sqrt(a^2 + b^2)
    best_point_id <- which.max(dist_to_line)
    score_cutoff <- xy_coordinates[best_point_id, "y"]
    select_ID <- scores$id[which(scores$value >= score_cutoff)]

    if (if_plot) {
        plot(seq(nrow(scores)), sorted_scores$value,
             col = ifelse(sorted_scores$value >=
            score_cutoff, "red", "black"), xlab = "ID", ylab = "Score",
            main = paste0("Optimal number = ", length(select_ID),
                " with cutoff value = ", round(score_cutoff, digits = 4)))
    }

    return(select_ID)
}



#' A feature/ gene selection function
#'
#' A feature/ gene selection function (1) removes sparsely expressed genes, (2) identifies
#' differentially expressed genes based on preliminary cell ordering, (3) removes highly
#' dispersed genes from the identified DEGs, (4) further picks genes which are expected
#' to have large expression difference on the 2 extreme ends of preliminary cell ordering
#'
#' @param cds a Monocle's CellDataSet object
#' @param min_expr the minimum expression value
#' @param scattering.cutoff.prob probability used for removing largely dispersed genes
#' @param driving.force.cutoff a value used for removing genes which do not change much along cell progress along cell progress path
#' @param qval_cutoff a user-defined adjusted p-value below which genes are retained
#' @param data_type a character indicating the type of underlying cell progression, i.e. linear or cyclical.
#' @param nCores Number of cores to use.
#' @return integer
#' @author MaiChan Lau
#' @importFrom VGAM sm.ns
#' @importFrom parallel detectCores
#' @importFrom monocle detectGenes
#' @export
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' #exprs <- uSORT_preProcess(exprs_file = file)
#' #exp_raw <- t(exprs$exprs_raw)
#' #exp_trimmed <- t(exprs$exprs_log_trimed)
#' #cds <- uSORT:::EXP_to_CellDataSet(exp_trimmed, exp_raw)
#' #driver_genes <- driving_force_gene_selection(cds = cds)
driving_force_gene_selection <- function(cds, scattering.cutoff.prob = 0.75,
    driving.force.cutoff = NULL, qval_cutoff = 0.05, min_expr = 0.1,
    data_type = c("linear", "cyclical"), nCores = 1) {

    if (nCores > detectCores()) {
        nCores <- max(detectCores - 1, 1)
    }

    ## =====Gene trimming: sparse/ dropout genes==========
    data_type <- match.arg(data_type)
    cds <- detectGenes(cds, min_expr = min_expr)
    trimmed_cds <- cds
    ## =====Gene trimming: scattered genes==========
    cat("  Removing scattered genes...\n")
    degree.scattering <- scattering_quantification_per_gene(trimmed_cds)
    scattering.cutoff <- quantile(degree.scattering, scattering.cutoff.prob)
    degree.scattering <- degree.scattering[degree.scattering <
        scattering.cutoff]
    gene.retained.2 <- names(degree.scattering)
    degree.scattering.output <- data.frame(GeneID = names(degree.scattering),
        Scattering = degree.scattering)
    cat("  Scattering.cutoff @prob", scattering.cutoff.prob, " =",
        scattering.cutoff, "\n")
    trimmed_cds <- trimmed_cds[gene.retained.2, ]

    no.genes <- nrow(trimmed_cds)
    cat("  No. of non scattered genes = ", no.genes, "\n")

    ## =====Identifying driver genes using Monocle's VGAM
    ## functions==========
    cat("  Identifying driver genes...\n")
    diff_test_res <- differentialGeneTest1(trimmed_cds, cores = nCores,
        fullModelFormulaStr = "expression~VGAM::sm.ns(Pseudotime, df=3)")
    # write.table(diff_test_res,file='driver.genes.pval.txt',row.names
    # = T, col.names = NA, quote=T, sep='\t')
    no.missing.models <- nrow(diff_test_res[diff_test_res$status ==
        "FAIL", ])
    cat("  No. of genes failed to fit vgam model = ", no.missing.models,
        "\n")
    driver.genes <- rownames(diff_test_res)[diff_test_res$qval <
        qval_cutoff]
    trimmed_cds <- trimmed_cds[as.character(driver.genes), ]
    no.genes <- nrow(trimmed_cds)
    cat(paste0("  No. of driver genes @qval(", qval_cutoff, ") = ",
        no.genes, "\n"))

    ## =====Filter driver genes for large driving force==========
    change <- ifelse(data_type == "cyclical", "cyclic.driving.force",
        "linear.driving.force")
    trimmed_diff_test_res <- diff_test_res[rownames(diff_test_res) %in%
        rownames(trimmed_cds@assayData$exprs), ]
    driving_force <- data.frame(GeneID = rownames(trimmed_diff_test_res),
        Score = trimmed_diff_test_res[, change])
    driving_force$Score[is.infinite(driving_force$Score)] <- 0
    driving_force <- driving_force[order(driving_force$Score,
                                         decreasing = TRUE),
        , drop = FALSE]

    if (nrow(driving_force) <= 1)
        stop("No. of filtered driver genes < 2!")

    if (!is.null(driving.force.cutoff)) {
        ## Based on user-input cutoff value
        selected_driver_genes <- driving_force[driving_force$Score >
            driving.force.cutoff, "GeneID"]
        cat("  Driving force cutoff = ", driving.force.cutoff,
            "\n")
    } else {
        ## Based on elbow detection method
        selected_driver_genes <- elbow_detection(scores = driving_force$Score)
        selected_driver_genes <- driving_force$GeneID[selected_driver_genes]
        if (length(selected_driver_genes) == 0)
            stop("Unable to get turning point while identifying driver genes!\n")
    }

    cat("  No. of genes selected for final sorting = ",
        length(selected_driver_genes),
        "\n")
    return(as.character(selected_driver_genes))
}


#' differential gene test
#'
#' modified from FludigmSC pacakge
#'
#' @param cds Input object.
#' @param fullModelFormulaStr Full model formula.
#' @param reducedModelFormulaStr Reduced model formula.
#' @param cores Number of cores will be used.
#'
#' @importFrom Biobase esApply
#' @return test results
differentialGeneTest1 <- function(cds,
                                  fullModelFormulaStr = "expression~sm.ns(Pseudotime, df=3)",
    reducedModelFormulaStr = "expression~1", cores = 1) {
    if (cores > 1) {
        diff_test_res <- mcesApply1(cds, 1, diff_test_helper1,
            cores = cores, fullModelFormulaStr = fullModelFormulaStr,
            reducedModelFormulaStr = reducedModelFormulaStr,
            expressionFamily = cds@expressionFamily,
            lowerDetectionLimit = cds@lowerDetectionLimit)
        diff_test_res
    } else {
        diff_test_res <- esApply(cds, 1, diff_test_helper1,
                                 fullModelFormulaStr = fullModelFormulaStr,
            reducedModelFormulaStr = reducedModelFormulaStr,
            expressionFamily = cds@expressionFamily,
            lowerDetectionLimit = cds@lowerDetectionLimit)
        diff_test_res
    }
    diff_test_res <- do.call(rbind.data.frame, diff_test_res)
    diff_test_res$qval <- p.adjust(diff_test_res$pval, method = "BH")
    diff_test_res
}


#' @importFrom Biobase multiassign exprs
#' @importFrom parallel makeCluster stopCluster
#' @importMethodsFrom BiocGenerics clusterEvalQ parRapply parCapply
mcesApply1 <- function(X, MARGIN, FUN, cores = 1, ...) {
    parent <- environment(FUN)
    if (is.null(parent))
        parent <- emptyenv()
    e1 <- new.env(parent = parent)
    multiassign(names(pData(X)), pData(X), envir = e1)
    environment(FUN) <- e1
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {
        require(VGAM)
    })
    if (MARGIN == 1) {
        res <- parRapply(cl, exprs(X), FUN, ...)
    } else {
        res <- parCapply(cl, exprs(X), FUN, ...)
    }
    stopCluster(cl)
    res
}


#' A modified monocle's function
#'
#' A modified monocle's function for 'compareModels' which identifies and removes genes
#' whose reduced_models is better than full_models in term of likelihood
#'
#' @param expr_matrix Expression matrix.
#' @param krange krange.
#' @param method method function.
#' @param ... Other parameters.
#'
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM lrtest
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom methods new
#' @importFrom fpc pamk
#' @importFrom cluster pam
#' @importFrom stats as.dist as.formula cor cov dist logLik median p.adjust prcomp predict quantile var
#' @importFrom utils read.table
clusterGenes1 <- function(expr_matrix, krange, method = function(x) {
    as.dist((1 - cor(t(x)))/2)
}, ...) {
    expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0,
        ]
    expr_matrix <- t(scale(t(log10(expr_matrix))))
    expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) ==
        FALSE, ]
    expr_matrix[is.na(expr_matrix)] <- 0
    n <- method(expr_matrix)
    pamk.best <- pamk(n, krange = krange)

    cat("optimal k =", pamk.best$nc, "\n")
    clusters <- pam(n, pamk.best$nc, ...)
    class(clusters) <- "list"
    clusters$exprs <- expr_matrix
    clusters
}


#' A modified monocle's function
#'
#' A modified monocle's function for 'compareModels' which identifies and removes genes
#' whose reduced_models is better than full_models in term of likelihood
#'
#' @param full_models a Monocle's vgam full model
#' @param reduced_models a Monocle's vgam reduced/ null model
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM lrtest
compareModels1 <- function(full_models, reduced_models) {
    stopifnot(length(full_models) == length(reduced_models))
    test_res <- mapply(function(x, y) {
        if (is.null(x) == FALSE && is.null(y) == FALSE) {
            if (logLik(x) < logLik(y)) {
                data.frame(status = "FAIL", pval = 1)
            } else {
                lrt <- lrtest(x, y)
                pval = lrt@Body["Pr(>Chisq)"][2, ]
                data.frame(status = "OK", pval = pval)
            }
        } else {
            data.frame(status = "FAIL", pval = 1)
        }
    }, full_models, reduced_models, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    test_res <- do.call(rbind.data.frame, test_res)
    test_res$qval <- p.adjust(test_res$pval, method = "BH")
    test_res
}


#' A modified monocle's helper function
#'
#' A modified monocle's function for 'diff_test_helper1' which includes more attempts
#' on finding models and also compute max. magnitude change in expression values predicted
#' by GLM model
#'
#' @param x an expression data
#' @param fullModelFormulaStr a Monocle's model structure
#' @param reducedModelFormulaStr a Monocle's model structure
#' @param expressionFamily a Monocle's family character
#' @param lowerDetectionLimit a threshold value
#' @param type_ordering a character indicating the type of underlying cell progression,
#' i.e. linear or circular
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM vgam
diff_test_helper1 <- function(x, fullModelFormulaStr, reducedModelFormulaStr,
    expressionFamily, lowerDetectionLimit = 0.1, type_ordering = "linear") {
    leftcensored <- x < lowerDetectionLimit
    x[x < lowerDetectionLimit] <- lowerDetectionLimit
    if (expressionFamily@vfamily %in% c("zanegbinomialff", "negbinomial",
        "poissonff", "quasipoissonff")) {
        expression <- round(x)
        integer_expression <- TRUE
    } else if (expressionFamily@vfamily %in% c("gaussianff")) {
        expression <- x
    } else {
        expression <- log10(x)
        integer_expression <- FALSE
    }


    p <- parent.env(environment())
    # print(ls(p))

    for (i in 1:20) {
        test_res <- tryCatch({

            full_model_fit <- suppressWarnings(
                vgam(as.formula(fullModelFormulaStr),
                family = expressionFamily))
            reduced_model_fit <- suppressWarnings(
                vgam(as.formula(reducedModelFormulaStr),
                family = expressionFamily))


            if (integer_expression) {
                pred_res <- predict(full_model_fit, type = "response")
                pred_res[pred_res < lowerDetectionLimit] <- lowerDetectionLimit
            } else {
                pred_res <- 10^(predict(full_model_fit, type = "response"))
                pred_res[pred_res < log10(lowerDetectionLimit)] <- log10(lowerDetectionLimit)
            }

            pred <- data.frame(Pseudotime = p$Pseudotime,
                               expectation = log10(pred_res))
            pred <- pred[order(pred$Pseudotime), ]
            cyclic.driving.force <- diff(range(pred$expectation))
            linear.driving.force <- abs(pred$expectation[1] -
                                            pred$expectation[nrow(pred)])
            res <- compareModels1(list(full_model_fit),
                                  list(reduced_model_fit))
            res$cyclic.driving.force <- cyclic.driving.force
            res$linear.driving.force <- linear.driving.force
            if (is.na(cyclic.driving.force) | is.na(linear.driving.force)) {
                res$pval = 1
                res$qval = 1
                res$status = "FAIL"
            } else res$status = "OK"
            res

        }, error = function(e) {
            print(e)
            data.frame(status = "FAIL", pval = 1, qval = 1,
                       cyclic.driving.force = 0,
                linear.driving.force = 0)
        })
        if (test_res$status != "FAIL")
            break
    }

    test_res
}



#' An expression scattering measurement function
#'
#' An expression scattering measurement function computes the level of scattering
#' for individual genes along the cell ordering
#'
#' @param CDS a Monocle's CellDataSet object
#' @return integer
#' @author MaiChan Lau
scattering_quantification_per_gene <- function(CDS = NULL) {

    x <- CDS@assayData$exprs
    if (CDS@expressionFamily@vfamily %in% c("zanegbinomialff",
        "negbinomial", "poissonff", "quasipoissonff")) {
        expression <- round(x)
        integer_expression <- TRUE
    } else if (CDS@expressionFamily@vfamily %in% c("gaussianff")) {
        expression <- x
    } else {
        expression <- log10(x)
        integer_expression <- FALSE
    }

    log_exprs <- log10(x)
    log_exprs[!is.finite(log_exprs)] <- log10(CDS@lowerDetectionLimit)

    total.dispersion <- apply(log_exprs, 1, variability_per_gene)
    return(total.dispersion)
}

#' A utility function for scattering_quantification_per_gene
#'
#' A utility function for scattering_quantification_per_gene which computes the degree
#' of scattering for single gene, whereby the value is computed by summing over the
#' local values of smaller local windows
#'
#' @param logExp a log-scale expression vector of a gene
#' @param min_expr a minimum expression value
#' @param window_size_perct a window size (in % of total no. of cells) used for computing
#' dispersion level
#' @param nonZeroExpr_perct a minimum amount of cells (in % of window size) with non-zero
#' expression, otherwise the associated window will be assigned to 0 disperson value
#' @return integer
#' @author MaiChan Lau
variability_per_gene <- function(logExp = NULL, min_expr = 0.1,
    window_size_perct = 0.1, nonZeroExpr_perct = 0.1) {

    num_cell <- length(logExp)
    window_size <- ceiling(num_cell * window_size_perct)
    ind.cell.var <- 0
    total.dispersion <- 0

    for (i in 1:num_cell) {

        data_window <- logExp
        startID <- (i - 1) * window_size + 1
        endID <- startID + (window_size)
        if (endID > num_cell)
            endID <- num_cell

        data_window <- data_window[startID:endID]
        if (sum(data_window > log10(min_expr)) <= nonZeroExpr_perct *
            window_size)
            ind.cell.var <- 0 else {
            window_nonZero <- data_window[data_window > log10(min_expr)]
            ind.cell.var <- diff(range(window_nonZero))/diff(range(logExp))
        }
        total.dispersion <- total.dispersion + ind.cell.var
        if (endID == num_cell)
            break
    }

    return(total.dispersion)
}
