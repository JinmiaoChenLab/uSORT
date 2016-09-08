#' A data loading and pre-processing function
#'
#' A data loading and pre-processing function which firstly identifies outlier cells and scarcely expressed genes.
#'
#' @param exprs_file Input file name in txt format, with rownames of cells and colnames of genes.
#' @param log_transform Boolean, if TRUE log transform the data.
#' @param remove_outliers Boolean, if TRUE remove the outliers.
#' @param lod A value of limit of detection in the unit of TPM/CPM/RPKM. It will be used as the starting value
#' for outlier cell detection and the basis for removing scarce genes.
#'
#' @return A list containing \code{exprs_raw}(data frame) and \code{exprs_log_trimed}(data.frame).
#' @export
#'
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' exprs <- uSORT_preProcess(exprs_file = file)
uSORT_preProcess <- function(exprs_file, log_transform = TRUE, 
    remove_outliers = TRUE, lod = 1) {
    
    exp <- fluidigmSC_readLinearExp(exprs_file)
    
    ## Outlier detection
    threshold <- lod
    step <- 4
    fine_step <- 1
    num_fine_test <- 4
    pct_goodsample_threshold <- 0.5
    quantile_threshold <- 0.95
    low_quantile_threshold <- 0.15
    min_gene_number <- 20
    outlier <- fluidigmSC_identifyExpOutliers(exp$log2ex_data, 
        exp$expression_data_raw, threshold, step = step, fine_step = fine_step, 
        num_fine_test = num_fine_test, pct_goodsample_threshold = pct_goodsample_threshold, 
        quantile_threshold = quantile_threshold, low_quantile_threshold = low_quantile_threshold, 
        min_gene_number = min_gene_number, lod = lod)
    # cat(outlier,'\n')
    cat(paste0("  ", length(outlier), " outliers were identified\n", 
        sep = ""))
    
    ## update exp
    log2ex_data <- exp$log2ex_data
    log2ex_data <- log2ex_data[, !colnames(log2ex_data) %in% outlier]
    
    ## update log2ex_avg
    expr_raw <- exp$expression_data_raw
    expression_data <- expr_raw[, !colnames(expr_raw) %in% outlier]
    expression_data_avg <- apply(expression_data, 1, mean)
    log2ex_avg_data <- as.data.frame(log2(expression_data_avg))
    
    ## remove lowly expressed genes
    cutoff <- lod * 2
    retained_genes <- fluidigmSC_removeGenesByLinearExpForAllType(log2ex_data, 
        log2ex_avg_data, threshold = cutoff)
    
    ## update exp
    log2ex_data <- log2ex_data[rownames(log2ex_data) %in% retained_genes, 
        ]
    return(list(exprs_raw = expr_raw, exprs_log_trimed = log2ex_data))
}



#' A gene finding function
#'
#' A gene finding function looking for genes in the target set x from the source set y,
#' reproduced from FluidigmSC package.
#'
#' @param x A vector of characters representing gene names (target genes).
#' @param y A vector of characters representing gene names (source genes).
#' @param ignore_case Boolean, if TRUE ignores letter case.
#'
#' @return A vector of characters representing gene names.
fluidigmSC_isElementIgnoreCase <- function(x, y, ignore_case = TRUE) {
    x <- as.vector(x)
    y <- as.vector(y)
    if (ignore_case == TRUE) {
        x <- tolower(x)
        y <- tolower(y)
    }
    return(is.element(x, y))
}


#' An expression reading function
#'
#' An expression reading function which imports expression data from .txt file, and then
#' computes log2 transformed data, reproduced from FluidigmSC package.
#'
#' @param exp_file Input file name in txt format, with rownames of cells and colnames of genes.
#' @param lod A value of limit of detection in the unit of TPM/CPM/RPKM. It will be used as the starting value
#' for outlier cell detection and the basis for removing scarce genes.
#'
#' @return A list containing \code{expression_data_raw}(data frame), \code{log2ex_data}(data frame),
#' and \code{log2ex_avg_data}(data frame).
fluidigmSC_readLinearExp <- function(exp_file = TRUE, lod = 1) {
    expression_data <- read.table(exp_file, row.names = 1, header = TRUE, 
        sep = "\t")
    expression_data <- as.matrix(expression_data)
    gene_ids <- as.vector(rownames(expression_data))
    expression_data <- apply(expression_data, 2, as.numeric)
    expression_data[is.na(expression_data)] <- -1
    expression_data_raw <- expression_data
    
    expression_data[expression_data < lod] <- lod
    
    log2ex_data <- as.data.frame(log2(expression_data))
    
    expression_data_avg <- apply(expression_data, 1, mean)
    log2ex_avg_data <- as.data.frame(log2(expression_data_avg))
    
    rownames(expression_data_raw) <- gene_ids
    rownames(log2ex_data) <- gene_ids
    rownames(log2ex_avg_data) <- gene_ids
    
    return(list(expression_data_raw = expression_data_raw, log2ex_data = log2ex_data, 
        log2ex_avg_data = log2ex_avg_data))
}


#' An outlier detection function
#'
#' An outlier detection function identifies cells with median expression below that of the bulk,
#' reproduced from FluidigmSC package.
#'
#' @param log2ex_data A data frame containing log2 tranformed expression values, with rownames of genes
#' and colnames of cells.
#' @param expression_data_raw A data frame containing raw expression values, with rownames of genes
#' and colnames of cells.
#' @param threshold A value in raw expression used as the starting threshold value.
#' @param step An integer number indicating the increment of threshold value at each iteration.
#' @param fine_step An integer number indicating the increment of threshold value at each iteration,
#' at the refining stage.
#' @param num_fine_test An integer number indicating the number of iteration of the refining stage.
#' @param pct_goodsample_threshold A fraction value indicating the minimum percentage of samples
#' on which the representative genes are detectable.
#' @param quantile_threshold A probability of gene detection rate above which a sample is considered
#' as good sample.
#' @param low_quantile_threshold A probability of average gene expression value below which a sample
#' is taken as an outlier.
#' @param min_gene_number An integer indicating the minimum size of representative genes.
#' @param lod  A value of limit of detection in the unit of TPM/CPM/RPKM.
#'
#' @return A vector of character stating the IDs of outlier cells.
fluidigmSC_identifyExpOutliers <- function(log2ex_data, expression_data_raw, 
    threshold, step, fine_step, num_fine_test, pct_goodsample_threshold = 0.5, 
    quantile_threshold = 0.95, low_quantile_threshold = 0.25, min_gene_number = 25, 
    lod) {
    last_threshold <- threshold
    good_sample_found = FALSE
    org_log2ex_data <- log2ex_data
    for (i in 0:100) {
        log2ex_data <- fluidigmSC_removeGenesByLinearExpForAllType_log2(log2ex_data, 
            threshold)
        
        expression_data_trimmed <- expression_data_raw[rownames(expression_data_raw) %in% 
            rownames(log2ex_data), ]
        
        gene_num <- length(rownames(log2ex_data))
        total_samples <- as.vector(colnames(log2ex_data))
        total_sample_num <- length(total_samples)
        gene_detection <- fluidigmSC_analyzeGeneDetection(expression_data_trimmed, 
            lod)
        good_samples <- rownames(gene_detection)[as.vector(gene_detection[, 
            2]) >= quantile_threshold]
        bad_samples <- rownames(gene_detection)[as.vector(gene_detection[, 
            2]) < quantile_threshold]
        num_good_samples <- length(good_samples)
        pct_good_samples <- num_good_samples/total_sample_num
        if (pct_good_samples >= pct_goodsample_threshold || gene_num < 
            min_gene_number || i == 100) {
            good_sample_found <- pct_good_samples >= pct_goodsample_threshold
            if (threshold != last_threshold) 
                last_threshold <- threshold - step
            break
        }
        threshold <- threshold + step
    }
    if (good_sample_found) {
        threshold <- last_threshold + fine_step
        log2ex_data <- org_log2ex_data
        for (i in 1:num_fine_test) {
            log2ex_data <- fluidigmSC_removeGenesByLinearExpForAllType_log2(log2ex_data, 
                threshold)
            expression_data_trimmed <- expression_data_raw[rownames(expression_data_raw) %in% 
                rownames(log2ex_data), ]
            
            gene_num <- length(rownames(log2ex_data))
            total_samples <- as.vector(colnames(log2ex_data))
            total_sample_num <- length(total_samples)
            gene_detection <- fluidigmSC_analyzeGeneDetection(expression_data_trimmed, 
                lod)
            good_samples <- rownames(gene_detection)[as.vector(gene_detection[, 
                2]) >= quantile_threshold]
            bad_samples <- rownames(gene_detection)[as.vector(gene_detection[, 
                2]) < quantile_threshold]
            num_good_samples <- length(good_samples)
            pct_good_samples <- num_good_samples/total_sample_num
            if (pct_good_samples >= pct_goodsample_threshold || 
                i == 100) {
                good_gene_exp <- as.matrix(log2ex_data[, is.element(colnames(log2ex_data), 
                  good_samples)])
                median_gene_exp <- apply(log2ex_data, 2, median)
                exp_quantile_cutoff <- quantile(good_gene_exp, 
                  low_quantile_threshold)
                outliers <- colnames(log2ex_data)[median_gene_exp <= 
                  exp_quantile_cutoff]
                break
            }
            threshold <- threshold + fine_step
        }
    } else {
        cat(paste("\nCould not find more than ", min_gene_number, 
            " commonly expressed genes.  Please manually select the outlier candidates in the boxplot.\n", 
            sep = ""))
        exp_quantile_cutoff = 0
        outliers <- character(0)
    }
    return(outliers)
}



#' A gene trimming function
#'
#' A gene trimming function removes genes whose average expression value is below the log2(threshold),
#' and also present in at least 10% of total cells;  reproduced from FluidigmSC package.
#'
#' @param log2ex_data A data frame containing log2 tranformed expression values, with rownames of genes
#' and colnames of cells.
#' @param log2ex_avg_data A data frame containing log2 tranformed average expression values for individual
#' gene.
#' @param threshold A limit of detection in the unit of TPM/CPM/RPMK.
#'
#' @return A vector of character containing gene names of those passed the filtering.
fluidigmSC_removeGenesByLinearExpForAllType <- function(log2ex_data, 
    log2ex_avg_data, threshold) {
    pct_detection <- 0.1
    
    removed_genes <- rownames(log2ex_avg_data)[log2ex_avg_data <= 
        log2(threshold)]
    
    gene_names <- rownames(log2ex_data)
    sample_names <- colnames(log2ex_data)
    detected_genes <- c()
    
    
    if (length(sample_names) > 1) {
        log2ex_data_numeric <- apply(log2ex_data, 1, as.numeric)
        log2ex_data_numeric <- log2ex_data_numeric >= log2(threshold)
        log2ex_data_numeric <- apply(log2ex_data_numeric, 2, sum)
    } else {
        log2ex_data_numeric <- log2ex_data
        log2ex_data_numeric <- log2ex_data_numeric >= log2(threshold)
    }
    gene_detection <- log2ex_data_numeric/length(sample_names)
    detected_genes <- gene_names[gene_detection >= pct_detection]
    
    undetected_genes_by_percent <- gene_names[!is.element(gene_names, 
        detected_genes)]
    removed_genes <- c(removed_genes, undetected_genes_by_percent)
    detected_genes <- gene_names[!fluidigmSC_isElementIgnoreCase(gene_names, 
        removed_genes)]
    
    return(detected_genes)
}

#' A gene trimming function
#'
#' A gene trimming function removes genes whose average expression value is below the log2(threshold);
#' reproduced from FluidigmSC package.
#'
#' @param log2ex_data A data frame containing log2 tranformed expression values, with rownames of genes
#' and colnames of cells.
#' @param threshold A limit of detection in the unit of TPM/CPM/RPMK.
#'
#' @return A vector of character containing gene names of those passed the filtering.
fluidigmSC_removeGenesByLinearExpForAllType_log2 <- function(log2ex_data, 
    threshold) {
    gene_list <- data.frame(GeneID = rownames(log2ex_data))
    log2ex_max_avg_data <- apply(log2ex_data, 1, mean)
    log2ex_max_avg_data <- as.data.frame(log2ex_max_avg_data)
    
    removed_genes <- rownames(log2ex_max_avg_data)[log2ex_max_avg_data <= 
        log2(threshold)]
    
    
    gene_list <- gene_list[!fluidigmSC_isElementIgnoreCase(gene_list[, 
        1], removed_genes), , drop = FALSE]
    log2ex_data <- log2ex_data[rownames(log2ex_data) %in% gene_list[, 
        1], ]
    return(log2ex_data)
}


#' A gene detection function
#'
#' A gene detection function computes the fraction of genes detected in each cell,
#' reproduced from FluidigmSC package.
#'
#' @param expression_data A data frame containing raw expression values, with rownames of genes
#' and colnames of cells.
#' @param threshold A limit of detection in the unit of TPM/CPM/RPMK.
#'
#' @return A data frame containing a column of number of genes detected, and a column of the
#' corresponding percentage of gene detection, rownames of cells.
fluidigmSC_analyzeGeneDetection <- function(expression_data, threshold = 1) {
    
    samples <- as.vector(colnames(expression_data))
    genes <- as.vector(rownames(expression_data))
    
    expression_data <- as.matrix(expression_data)
    expression_data <- apply(expression_data, 2, as.numeric)
    cutoff <- c()
    cutoff <- expression_data >= threshold
    
    gene_num <- apply(cutoff, 2, sum, na.rm = TRUE)
    total_genes <- length(genes)
    pct_gene_detection <- gene_num/total_genes
    all <- cbind(gene_num, pct_gene_detection, samples)
    all <- all[order(as.numeric(all[, 1])), ]
    
    gene_num <- as.numeric(all[, 1])
    samples <- all[, 3]
    
    gene_detection <- as.data.frame(all[, 1:2])
    rownames(gene_detection) <- rownames(all)
    colnames(gene_detection) <- c("number_of_genes", "pct_of_genes")
    invisible(gene_detection)
}

