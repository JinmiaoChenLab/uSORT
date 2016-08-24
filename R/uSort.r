#' uSort: A self-refining ordering pipeline for gene selection
#'
#' This package is designed to uncover the intrinsic cell progression path
#' from single-cell RNA-seq data.
#'
#' This package incorporates data pre-processing, preliminary PCA gene selection, preliminary
#' cell ordering, feature selection, refined cell ordering, and post-analysis interpretation
#' and visualization. The uSort workflow can be implemented through calling the main function
#' \code{\link{uSort}}.
#'
#' @docType package
#' @name uSort
#'
NULL


#' The main function of uSort pacakge
#'
#' The main function of \code{uSort-pacakge} which provides a workflow of sorting scRNA-seq data.
#'
#' @param result_directory the directory to which output results will be written
#' @param data_type a character indicating the type of underlying cell progression, i.e. linear or cyclical
#' @param exprs_data
#' @param if_preprocess
#' @param preliminary_sort_method
#' @param refine_sort_method
#' @param scattering_cutoff_prob
#' @param driving_force_cutoff
#' @param sorting_method
#' @param alpha
#' @param sigma_width
#' @param no_randomization
#' @param diffusionmap_components
#' @param l
#' @param k
#' @param num_waypoints
#' @param waypoints_seed
#' @param qval_cutoff_featureSelection a user-defined adjusted p-value below which genes are selected for driver gene selection
#'
#' @return a matrix \code{res} of the permutated norm_exprs. If choose 'writeResults = TRUE', results will be saved into files under \code{result_directory}
#'
#' @author MaiChan Lau
#'
#' @import fluidigmSC
#' @import monocle
#' @importFrom gplots heatmap.2
#' @importFrom gplots colorpanel
#' @importFrom ggplot2 ggsave
#' @importFrom VGAM sm.ns
#' @export
uSort <- function(exprs_data,
                  if_preprocess = TRUE,
                  preliminary_sort_method = c("autoSPIN, sWanderlust, monocle, SPIN"),
                  refine_sort_method = c("autoSPIN, sWanderlust, monocle, SPIN"),
                  result_directory = file.path(exprs_data),
                  # gene selection parameters
                  scattering_cutoff_prob = 0.75,
                  driving_force_cutoff = NULL,
                  qval_cutoff_featureSelection = 0.05,
                  # autoSPIN parameters
                  data_type = c("linear", "cyclical"),
                  sorting_method = c("STS", "neighborhood"),
                  alpha = 0.2,
                  sigma_width = 1,
                  no_randomization = 10,
                  # Wanderlust parameters
                  diffusionmap_components = 4,
                  l=15, k=15,
                  num_waypoints = 150,
                  waypoints_seed = 2711){

    ## Parameter check
    preliminary_sort_method <- match.arg(preliminary_sort_method)
    refine_sort_method <- match.arg(refine_sort_method)
    data_type <- match.arg(data_type)
    sorting_method <- match.arg(sorting_method)
    fluidigmSC_functions_update()
    if(!dir.exists(result_directory))
        dir.create(result_directory)
    setwd(result_directory)

    ## Data loading and pre-processing
    exp_trimmed <- uSort_preProcess(exprsData, preProcess = TRUE)
    saved_filename <- itn_saveData(exp_trimmed,
                                   save_title = 'preprocessed_exp',
                                   init_filename = 'preprocessed_exp')


    ## Preliminary gene selection using PCA
    pca_res <- prcomp(as.matrix(t(exp_trimmed$log2ex_data)))
    eigenVal <- pca_res$sdev ^ 2
    eigenVal <- data.frame('PC' = seq(length(eigenVal)), 'Score' = eigenVal)
    opt_PC <- elbow_detection(scores = eigenVal,
                              xlab.text = 'PC#',
                              ylab.text = 'Variance',
                              baseNm = 'optimal_number_of_PC')
    pca_conmpnent_number <- nrow(opt_PC)
    cat('No. of PC selected = ', no_component, '\n')
    rotation_matrix <- as.data.frame(abs(pca_res$rotation[, 1:no_component, drop = FALSE]))
    rotation_pc_sums <- rowSums(rotation_matrix)
    sum_order <- order(rotation_pc_sums, decreasing = TRUE)
    gene_score <- data.frame('GeneID' = names(rotation_pc_sums)[sum_order],
                             'Score'  = rotation_pc_sums[sum_order])
    PCA_genes <- elbow_detection(scores = gene_score,
                                 baseNm = 'PCA.gene.selection',
                                 ylab.text = paste0('Sum of (abs)PC1-', no_component, ' loadings'),
                                 xlab.text = 'Gene#')
    colnames(PCA_genes) <- 'GeneID'
    cat('no. of PCA genes = ', nrow(PCA_genes), '\n')
    write.table(PCA_genes, file = 'pca_genes.txt', row.names = F,
                col.names = T, quote = T, sep = '\t')

    ## Preliminary uSort sorting
    exp_PCA_genes <- updateGeneListFromList(exp_trimmed, PCA_genes)
    cat('\nPreliminary sorting:\n')
    preliminary_ordering <- sorting_wraper(exp_PCA_genes,
                                           data_type = data_type,
                                           method = preliminary_sort_method,
                                           local_sort = T,
                                           sorting_method = sorting_method,
                                           sigma_width = sigma_width,
                                           no_randomization = no_randomization,
                                           data_type = data_type,
                                           alpha = alpha)
    exp_preliminary_ordering <- updateSampleListFromList(exp_trimmed, preliminary_ordering)
    write.table(exp_preliminary_ordering$sample_list,
                file = 'preliminary_ordering.txt',
                row.names = F,
                col.names = T,
                quote = T,
                sep = '\t')


    ## Plot cell-to-cell distance heatmap
    exp_PCA_genes <- updateSampleListFromList(exp_PCA_genes, preliminary_ordering)
    c2c.dist <- distance.function(t(exp_PCA_genes$log2ex_data))
    pdf('C2C.distHeatmap_preliminary.pdf')
    col <- gplots::colorpanel(10, "red", "yellow", "blue")
    gplots::heatmap.2(c2c.dist,
                      dendrogram = 'none',
                      trace = 'none',
                      col = col,
                      Rowv = F,
                      Colv = F,
                      labCol = NA ,
                      labRow = NA ,
                      keysize = 1.2,
                      cexRow = 1.1,
                      main = paste('Cell to cell distance \ncomputed on ',
                                   nrow(PCA_genes), ' genes'))
    dev.off()

    ## Refined gene selection by driver genes selection
    cds <- EXP_to_CellDataSet(EXP = exp_preliminary_ordering)
    driver_genes_refined_order <- feature_selection(cds = cds,
                                                    scattering.cutoff.prob = scattering_cutoff_prob,
                                                    driving.force.cutoff = driving_force_cutoff,
                                                    data_type = data_type,
                                                    qval_cutoff = qval_cutoff_featureSelection)

    ## Refined uSort sorting
    exp_driver_genes_refined_ordering <- updateGeneListFromList(exp_preliminary_ordering,
                                                                driver_genes_refined_order)
    cat('Refined sorting:\n')
    refined_ordering <- sorting_wraper(exp_driver_genes_refined_ordering,
                                       data_type = data_type,
                                       method = preliminary_sort_method,
                                       local_sort = T,
                                       sorting_method = sorting_method,
                                       sigma_width = sigma_width,
                                       no_randomization = no_randomization,
                                       data_type = data_type,
                                       alpha = alpha)

    exp_refined_ordering <- updateSampleListFromList(exp_driver_genes_refined_ordering, refined_ordering)
    write.table(exp_refined_ordering$sample_list,
                file = 'refined.cell.ordering.txt',
                row.names = F,
                col.names = T,
                quote = T,
                sep = '\t')

    ## Plot driver gene profiles on final ordering
    pdf('final.driver.genes.profiles.pdf', width = 8)
    heatmap.2(as.matrix( exp_refined_ordering$log2ex_data),
              dendrogram = 'row',
              trace = 'none',
              col = bluered,
              Rowv = T,
              Colv = F,
              scale = 'row',
              cexRow = 0.5,
              margins = c(8, 8))
    dev.off()

    ## Plot cell-to-cell distance heatmap
    c2c.dist <- distance.function(t(exp_refined_ordering$log2ex_data))
    pdf('C2C.distHeatmap.pdf')
    col <- gplots::colorpanel(10, "red", "yellow", "blue")
    gplots::heatmap.2(c2c.dist,
                      dendrogram = 'none',
                      trace = 'none',
                      col = col,
                      Rowv = F,
                      Colv = F,
                      labCol = NA ,
                      labRow = NA ,
                      keysize = 1.2,
                      cexRow = 1.1,
                      main = paste('Cell to cell distance \ncomputed on ',
                                   nrow(driver_genes_refined_order), ' genes'))
    dev.off()

    return(exp_refined_ordering$sample_list)
}


#' Data loading and pre-processign for uSort
#'
#' @param exprsData
#' @param preProcess
uSort_preProcess <- function(exprsData, preProcess = TRUE){

    if(file.exists(exprs_data)){
        exp_noOutlier <- fluidigmSC::identifyOutliers(exprsData)  ## load data from file
    }else if(is.object(exprs_data)){
        exp_noOutlier <- exprsData                                ## data is a object
    }else{
        stop("Cannot loading data, please check your exprsData parameter!")
    }

    cat('\nOriginal data: no. of genes = ',
        nrow(exprsData$org_data),
        '\t no. of samples = ',
        ncol(exprsData$org_data),
        '\n'
        )
    if(preProcess){
        exp_trimmed <- pre_processing(exp_noOutlier)
        cat(
            'After SINGuLAR pre-processing: no. of genes = ',
            nrow(exp_trimmed$gene_list),
            '\t no. of samples = ',
            nrow(exp_trimmed$sample_list),
            '\n'
        )
        return(exp_trimmed)
    }else{
        exp_trimmed <- exp_noOutlier
    }

    return(exp_trimmed)
}


#' wrapper of all avaliable sorting methods in uSort
#'
#'
#' @param data
#' @param method
#' @param data_type
#' @param sorting_method
#' @param alpha
#' @param sigma_width
#' @param no_randomization
#' @param diffusionmap_components
#' @param num_waypoints
#' @param waypoints_seed
sorting_wraper <- function(data,
                           method = c("autoSPIN, sWanderlust, monocle, SPIN"),
                           # autoSPIN parameters
                           data_type = c("linear", "cyclical"),
                           sorting_method = c("STS", "neighborhood"),
                           alpha = 0.2,
                           sigma_width = 1,
                           no_randomization = 10,
                           # Wanderlust parameters
                           diffusionmap_components = 4,
                           l=15, k=15,
                           num_waypoints = 150,
                           waypoints_seed = 2711){

    data_type <- match.arg(data_type)
    method <- match.arg(method)
    data_type <- match.arg(data_type)
    sorting_method <- match.arg(sorting_method)
    switch(method,
           autoSPIN = {
               order <- autoSPIN(t(data$log2ex_data),
                                 sorting_method = sorting_method,
                                 baseNm = 'preliminary',
                                 sigma_width = sigma_width,
                                 no_randomization = no_randomization,
                                 data_type = data_type,
                                 alpha = alpha)
           },
           SPIN = {
               order <- SPIN(t(data$log2ex_data),
                             baseNm = 'SPIN',
                             sorting_method = sorting_method,
                             sigma_width = sigma_width)

           },
           monocle = {
               g <- data$gene_list
               cds <- EXP_to_CellDataSet(EXP = data)
               f <- fData(cds)
               f$gene_short_name <- rownames(f)
               fData(cds) <- f
               cds <- setOrderingFilter(cds, ordering_genes = g[,1])
               cds <- reduceDimension(cds, use_irlba=FALSE)
               cds <- orderCells(cds, num_paths=num_path, reverse=TRUE)
               pseudotime <- pData(cds)
               pseudotime <- pseudotime[order(pseudotime$Pseudotime), ,drop=FALSE]
               order <- rownames(pseudotime)
           },
           sWanderlust = {
               order <- sWanderlust(data,
                                    alpha = alpha,
                                    sorting_method = sorting_method,
                                    baseNm = 'uSPIN',
                                    sigma_width = sigma_width,
                                    data_type = data_type,
                                    diffusionmap_components = diffusionmap_components,
                                    l=l, k=k,
                                    num_waypoints = num_waypoints,
                                    waypoints_seed = waypoints_seed)
           })

    return(order)
}



