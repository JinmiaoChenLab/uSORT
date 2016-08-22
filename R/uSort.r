#' uSort: A self-refining ordering pipeline for gene selection
#'
#' This package is designed to uncover the intrinsic cell progression path
#' from single-cell RNA-seq data. 
#'
#' This package incorporates data pre-processing, preliminary PCA gene selection, preliminary
#' cell ordering, feature selection, refined cell ordering, and post-analysis interpretation
#' and visualization. The uSort workflow can be implemented through calling the main function
#' \code{\link{uSort_main}}.
#'
#' @docType package
#' @name uSort
#'
NULL

#' A main function of uSort
#'
#' A main function of uSort \code{uSort_main} which provides a workflow of sorting scRNA-seq data,
#' including data preprocess with merging methods of multiple fcs file, logicle transformation,
#' dimension reduction with PCA, isomap or tsne(default), and a kernal-based local maxima
#' clustering combined with SVM for subpopulation detection. The intermediate results can be saved
#' into seperate files and the cluster results can be visualized in heatmaps and scatter plots.
#'
#' @param resDir the directory to which output results will be written
#' @param exprsFile a character of input expression filename
# @param no_component number of principal components used for preliminary gene selection
#' @param exp.outlier.removed a EXP object of preprocessed expression data
#' @param driving.force.cutoff a value used for removing genes which do not
#' change much along cell progress
#' @param data_type a character indicating the type of underlying cell progression, 
#' i.e. linear or cyclical
#' @param qval_cutoff_featureselection a user-defined adjusted p-value below which
#' genes are selected for driver gene selection
#' @return a matrix \code{res} of the permutated norm_exprs. If choose 'writeResults = TRUE', results
#' will be saved into files under \code{resDir}
#' @author MaiChan Lau
#' @import fluidigmSC
#' @import monocle
#' @importFrom gplots heatmap.2
#' @importFrom gplots colorpanel
#' @importFrom ggplot2 ggsave
#' @importFrom VGAM sm.ns
#' @export

uSort_main<-function(resDir = NULL, exprsFile = NULL,scattering.cutoff.prob=0.75, 
                     exp.outlier.removed = NULL, driving.force.cutoff = NULL,
                     data_type = NULL, 
                     qval_cutoff_featureselection = 0.05, sigma_width = 1,
                     sorting_method = NULL, no_randomization = NULL,
                     alpha=NULL){
    ## =====Inputs check==========
    preprocess_flag=FALSE
    if(is.null(exp.outlier.removed)){
        
        if(is.null(exprsFile)) stop("No exprs file selected!")
        else
            #exprs <- read.table(exprsFile, row.names=1, header = TRUE)
            preprocess_flag <- TRUE
        
    }
    else exp_noOutlier <- exp.outlier.removed
    
    if(is.null(resDir)) stop("No output directory selected")
    
    ## =====Modify some fluidigmSC's functions==========
    assignInNamespace("itn_isSubSetIgnoreCase",itn_isSubSetIgnoreCase,ns="fluidigmSC",, envir=as.environment("package:fluidigmSC"))
    assignInNamespace("itn_getMatchingOriginalName",itn_getMatchingOriginalName,ns="fluidigmSC")
    assignInNamespace("itn_anova_by_factor",itn_anova_by_factor,ns="fluidigmSC")
    assignInNamespace("identifyOutliers",identifyOutliers,ns="fluidigmSC")
    assignInNamespace("itn_displayOutlierGUI",itn_displayOutlierGUI,ns="fluidigmSC")
    assignInNamespace("itn_ScatterPlot",itn_ScatterPlot,ns="fluidigmSC")
    assignInNamespace("itn_saveData",itn_saveData,ns="fluidigmSC")
    
    
    ## =====Pre-processing==========
    setwd(resDir)
    if(preprocess_flag == TRUE) exp_noOutlier <- fluidigmSC::identifyOutliers(exprsFile)
    exp_trimmed <- pre_processing(exp_noOutlier)
    cat('\nOriginal data: no. of genes = ', nrow(exp_trimmed$org_data),
        '\t no. of samples = ', ncol(exp_trimmed$org_data),'\n')
    cat('After SINGuLAR pre-processing: no. of genes = ', nrow(exp_trimmed$gene_list),
        '\t no. of samples = ', nrow(exp_trimmed$sample_list),'\n')
    save_title <- 'preprocessed_exp'
    saved_filename <- itn_saveData(exp_trimmed, 
                                   save_title = save_title, init_filename = save_title)

    ## =====Preliminary gene selection using PCA==========
    pca_res<-prcomp(as.matrix(t(exp_trimmed$log2ex_data)))
    eigenVal <- pca_res$sdev^2
    eigenVal <- data.frame('PC'=seq(length(eigenVal)),'Score'=eigenVal)
    opt_PC <- elbow_detection(scores = eigenVal, xlab.text='PC#',
                                 ylab.text = 'Variance',
                                 baseNm = 'optimal_number_of_PC')
    
    no_component <- nrow(opt_PC)   
    cat('No. of PC selected = ', no_component,'\n')
    PC1.to.X<-as.data.frame(abs(pca_res$rotation[,1:no_component,drop=F]))
    # sum of PCs
    PC1.to.X$total<-rowSums(PC1.to.X[,1:no_component,drop=F])

    PC1.to.X<-PC1.to.X[,-c(1:no_component),drop=F]
    sorted_PC1.to.X<-PC1.to.X[order(PC1.to.X[,1],decreasing=T),,drop=F]
    pca <- data.frame('GeneID'=rownames(sorted_PC1.to.X),'Score'=sorted_PC1.to.X$total)
    PCA_genes <- elbow_detection(scores = pca, baseNm = 'PCA.gene.selection',
                                 ylab.text = paste0('Sum of (abs)PC1-',no_component,' loadings'), 
                                 xlab.text= 'Gene#')

    # The top 100 pca genes
    #PCA_genes <- pca[1:100,'GeneID',drop=FALSE]
    colnames(PCA_genes) <- 'GeneID'
    #pca <- PCA(exp_trimmed, display_plots = FALSE)
    #PCA_genes <- fldm_pca$pca_ranked_genes[seq(100),'GeneID',drop=FALSE]
    
    cat('no. of pca genes = ', nrow(PCA_genes),'\n')  
    write.table(PCA_genes,file='pca.genes.txt',row.names = F, col.names = T, quote=T, sep='\t')

    #PCA_genes<-read.table('../commonPCAgenes_run1n16n17n18.txt', header = T)
    ## =====Preliminary uSort sorting==========
    exp.PCA.genes <- fluidigmSC::updateGeneListFromList(exp_trimmed, PCA_genes)
    cat('\nPreliminary sorting:\n')
    preliminary.ordering <- sorting_wrapper(t(exp.PCA.genes$log2ex_data), local_sort=T, sorting_method = sorting_method,
                                            baseNm = 'preliminary', sigma_width=sigma_width,
                                            no_randomization=no_randomization, data_type=data_type,
                                            alpha=alpha)
    exp.preliminary.ordering <- updateSampleListFromList(exp_trimmed, preliminary.ordering)
    write.table(exp.preliminary.ordering$sample_list,file='preliminary.ordering.txt',row.names = F, col.names = T, quote=T, sep='\t')
    
    ## Plot cell-to-cell distance heatmap
    exp.PCA.genes <- updateSampleListFromList(exp.PCA.genes,preliminary.ordering)
    c2c.dist<-distance.function(t(exp.PCA.genes$log2ex_data))
    pdf('C2C.distHeatmap_preliminary.pdf')
    col<-gplots::colorpanel( 10, "red", "yellow", "blue" )
    gplots::heatmap.2(c2c.dist,dendrogram='none',trace='none',col=col,
                      Rowv=F,Colv=F,labCol=NA ,labRow=NA ,keysize=1.2,
                      cexRow =1.1,
                      main=paste('Cell to cell distance \ncomputed on ',nrow(PCA_genes), ' genes'))
    dev.off()
    
    ## =====Select driver genes for refined sorting==========
    cds <- EXP_to_CellDataSet(EXP = exp.preliminary.ordering)
    driver.genes.refined.sorting <- feature_selection(cds = cds, scattering.cutoff.prob = scattering.cutoff.prob, driving.force.cutoff=driving.force.cutoff,
                                                      data_type=data_type,
                                                      qval_cutoff=qval_cutoff_featureselection)
    
    
    ## =====Refined uSort sorting==========
    exp.refined.ordering <- fluidigmSC::updateGeneListFromList(exp.preliminary.ordering, driver.genes.refined.sorting)
    cat('Refined sorting:\n')
    refined.ordering <- sorting_wrapper(t(exp.refined.ordering$log2ex_data), local_sort=T, sorting_method = sorting_method,
                                      baseNm = 'refined', sigma_width=sigma_width,
                                      no_randomization=no_randomization, data_type=data_type,
                                      alpha=alpha)
    exp.refined.ordering <- updateSampleListFromList(exp.refined.ordering, refined.ordering)
    write.table(exp.refined.ordering$sample_list, file='refined.cell.ordering.txt',row.names = F,col.names = T, quote=T, sep='\t')

    ## ===== For STS sorting: Refining with Wanderlust==========
    if(sorting_method=='STS' && data_type == 'linear'){
        ## command
        command = "python3.4"
        path2script="../wanderlust.py"
        ## Prepare data file
        tdata <- t(exp.refined.ordering$log2ex_data)
        write.csv(tdata, file = "log2ex_drivergenes_refinedOrdering_logTrans.csv")
        
        # Forward first cell
        cat('Refined sorting using Wanderlust:\n')
        startingCell <- as.character(refined.ordering$SampleID[1])
        args = startingCell
        allArgs = c(path2script, args)
        output = system2(command, args=allArgs, stdout=TRUE)
        wanderlust_trajectory <- read.csv("wanderlust_trajectory.csv",
                                          sep = "\t", header = FALSE, row.names = 1)
        colnames(wanderlust_trajectory) <- 'pseudotime'
        sorted_wanderlust_trajectory <- wanderlust_trajectory[order(wanderlust_trajectory$pseudotime),,drop=F]
        wanderlust_ordering1 <- data.frame('SampleID'= rownames(sorted_wanderlust_trajectory), 'GroupID'='untitled')
        wanderlust.exp1 <- updateSampleListFromList(exp.refined.ordering, wanderlust_ordering1)
        c1 <- STS_sortingcost(expr = t(wanderlust.exp1$log2ex_data))
        #cat('Starting cell ', startingCell,', cost=', c1 ,'\n')
        
        # Reverse first cell
        startingCell <- as.character(refined.ordering$SampleID[nrow(refined.ordering)])
        args = startingCell
        allArgs = c(path2script, args)
        output = system2(command, args=allArgs, stdout=TRUE)
        wanderlust_trajectory <- read.csv("wanderlust_trajectory.csv",
                                          sep = "\t", header = FALSE, row.names = 1)
        colnames(wanderlust_trajectory) <- 'pseudotime'
        sorted_wanderlust_trajectory <- wanderlust_trajectory[order(wanderlust_trajectory$pseudotime),,drop=F]
        wanderlust_ordering2 <- data.frame('SampleID'= rownames(sorted_wanderlust_trajectory), 'GroupID'='untitled')
        wanderlust.exp2 <- updateSampleListFromList(exp.refined.ordering, wanderlust_ordering2)
        c2 <- STS_sortingcost(expr = t(wanderlust.exp2$log2ex_data))
        #cat('Starting cell ', startingCell,', cost=', c2 ,'\n')
        

        # Pick an optimal ordering
        if(c1 < c2){
            write.table(wanderlust_ordering1, file='wanderlust_ordering.txt', row.names = F, col.names = T,
                        quote = F, sep = '\t')      
            exp.refined.ordering <- updateSampleListFromList(exp.refined.ordering, wanderlust_ordering1)
            cat('Starting cell ', startingCell,', cost=', c1 ,'(driver genes)\n')
        }else{
            write.table(wanderlust_ordering2, file='wanderlust_ordering.txt', row.names = F, col.names = T,
                        quote = F, sep = '\t')      
            exp.refined.ordering <- updateSampleListFromList(exp.refined.ordering, wanderlust_ordering2)
            cat('Starting cell ', startingCell,', cost=', c2 ,'(driver genes)\n')
        
        }
        exp.cost <- fluidigmSC::updateGeneListFromFile(exp.refined.ordering, gene_list_file='pca.genes.txt')
        STScost <- STS_sortingcost(expr = t(exp.cost$log2ex_data))
        cat('STS sorting cost of final ordering (PCA genes) = ', STScost, '\n')

    }else if(sorting_method=='neighborhood' && data_type == 'cyclical'){
        exp.cost <- fluidigmSC::updateGeneListFromFile(exp.refined.ordering, gene_list_file='pca.genes.txt')
        Neighborhoodcost <- neighborhood_sortingcost(expr = t(exp.cost$log2ex_data), sigma_width=1)
        cat('Neighborhood sorting cost of final ordering (PCA genes) = ', Neighborhoodcost, '\n')
        STScost <- STS_sortingcost(expr = t(exp.cost$log2ex_data))
        cat('STS sorting cost of final ordering (PCA genes) = ', STScost, '\n')
    }

    ##==========================================
    ## Plot driver gene profiles on final ordering
    pdf('final.driver.genes.profiles.pdf',width=8)     
    heatmap.2(as.matrix(exp.refined.ordering$log2ex_data),dendrogram='row',trace='none',col = bluered,
              Rowv=T,Colv=F,scale = 'row',cexRow=0.5,  margins = c(8, 8))
    dev.off()
    
    ## Plot cell-to-cell distance heatmap
    c2c.dist<-distance.function(t(exp.refined.ordering$log2ex_data))
    pdf('C2C.distHeatmap.pdf')
    col<-gplots::colorpanel( 10, "red", "yellow", "blue" )
    gplots::heatmap.2(c2c.dist,dendrogram='none',trace='none',col=col,
                      Rowv=F,Colv=F,labCol=NA ,labRow=NA ,keysize=1.2,
                      cexRow =1.1,
                      main=paste('Cell to cell distance \ncomputed on ',nrow(driver.genes.refined.sorting), ' genes'))
    dev.off()
    
    
    return(exp.refined.ordering$sample_list)
}


