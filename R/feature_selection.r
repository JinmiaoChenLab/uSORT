#' A feature/ gene selection function
#'
#' A feature/ gene selection function (1) removes sparsely expressed genes, (2) identifies
#' differentially expressed genes based on preliminary cell ordering, (3) removes highly
#' dispersed genes from the identified DEGs, (4) further picks genes which are expected
#' to have large expression difference on the 2 extreme ends of preliminary cell ordering
#'
#' @param cds a Monocle's CellDataSet object
#' @param min_expr the minimum expression value
#' @param scattering.cutoff.prob probability used for removing largely dispersed
#' genes
#' @param driving.force.cutoff a value used for removing genes which do not
#' change much along cell progress along cell progress path
#' @param qval_cutoff a user-defined adjusted p-value below which
#' genes are retained
#' @param data_type a character indicating the type of underlying cell progression, 
#' i.e. linear or cyclical
#' @return integer
#' @author MaiChan Lau
#' @importFrom VGAM sm.ns
#' @importFrom parallel detectCores
#' @export

feature_selection <- function(cds = NULL, scattering.cutoff.prob = 0.75,
                              driving.force.cutoff = NULL, qval_cutoff = NULL,
                              min_expr = 0.1, data_type = NULL){
    
    ## =====Gene trimming: sparse/ dropout genes==========
    cds <- monocle::detectGenes(cds, min_expr = min_expr)

#     cat('\nRemoving genes with significant dropout@dropout', max_dropout_percent,'......\n')
#     gene.retained.1 <- rownames(cds)[cds@featureData@data$num_cells_expressed >= ncol(cds)*max_dropout_percent]
#     trimmed_cds <- cds[gene.retained.1,]
#     no.genes <- nrow(trimmed_cds)
#     cat('No. of genes remained = ', no.genes,'\n')
#     write.table(data.frame('GeneID'=gene.retained.1),
#                 file='nondropout.genes.txt',row.names = F,col.names = T,sep='\t')
    
    trimmed_cds <- cds
    ## =====Modify some Monocle's functions==========
    assignInNamespace("compareModels",compareModels,ns="monocle")
    assignInNamespace("diff_test_helper",diff_test_helper,ns="monocle")
    assignInNamespace("differentialGeneTest",differentialGeneTest,ns="monocle")

    
    ## =====Gene trimming: scattered genes==========
    cat('\nRemoving scattered genes ......\n')
    
    degree.scattering <- scattering_quantification_per_gene(trimmed_cds)
    write.table(degree.scattering,file='full.degree.scattering.txt',
                row.names = F, col.names = T, quote=T, sep='\t')
    
    scattering.cutoff <- quantile(degree.scattering,scattering.cutoff.prob)
    degree.scattering <- degree.scattering[degree.scattering < scattering.cutoff]
    gene.retained.2 <- names(degree.scattering)
    degree.scattering.output <- data.frame('GeneID'=names(degree.scattering), 'Scattering'=degree.scattering)
    cat('scattering.cutoff @prob', scattering.cutoff.prob, ' =' , scattering.cutoff,'\n')
    
    write.table(degree.scattering.output,file='degree.scattering.txt',
                row.names = F, col.names = T, quote=T, sep='\t')
    
    trimmed_cds <- trimmed_cds[gene.retained.2,]
    #trimmed_cds <- trimmed_cds['H2-Aa',]
     
    no.genes <- nrow(trimmed_cds); cat('No. of non scattered genes = ', no.genes,'\n')
    write.table(data.frame('GeneID'=gene.retained.2),file='nonscattered.genes.txt',row.names = F, col.names = T, quote=T, sep='\t')
    
    
    ## =====Identifying driver genes using Monocle's VGAM functions==========
    cat('Identifying driver genes ......\n')
    #nCores <- parallel::detectCores()
    nCores <- 5
    diff_test_res <- monocle::differentialGeneTest(trimmed_cds, cores=nCores,
                                                   fullModelFormulaStr = "expression~VGAM::sm.ns(Pseudotime, df=3)")
    write.table(diff_test_res,file='driver.genes.pval.txt',row.names = T, col.names = NA, quote=T, sep='\t')
    no.missing.models <- nrow(diff_test_res[diff_test_res$status == 'FAIL',])
    cat('No. of genes failed to fit vgam model = ', no.missing.models,'\n')
    driver.genes <- rownames(diff_test_res)[diff_test_res$qval < qval_cutoff]
    trimmed_cds<-trimmed_cds[as.character(driver.genes),]
    no.genes<-nrow(trimmed_cds); cat('No. of driver genes @qval(',qval_cutoff, ')= ', no.genes,'\n')
    write.table(data.frame('GeneID'=driver.genes),file='driver.genes.txt',
                row.names = F, col.names = T, quote=T, sep='\t')
    
    
    ## =====Filter driver genes for large driving force==========
    if(data_type == 'cyclical') change <- 'cyclic.driving.force'
    else change <- 'linear.driving.force'
    
    trimmed_diff_test_res <- diff_test_res[rownames(diff_test_res)%in% rownames(trimmed_cds@assayData$exprs),]
    driving.force <- data.frame('GeneID'=rownames(trimmed_diff_test_res),'Score'=trimmed_diff_test_res[,change])
    driving.force[is.infinite(driving.force$Score),'Score'] <- 0
    driving.force<-driving.force[order(driving.force$Score,decreasing=T),,drop=F]
    driving.force.output <- driving.force
    colnames(driving.force.output) <- c('GeneID','Driving.force')
    write.table(driving.force.output,file='driving.force.txt',row.names = F, col.names = T, quote=T, sep='\t')
    
    if (nrow(driving.force)<=1) stop("No. of filtered driver genes < 2!")
    
    ## Based on user-input cutoff value
    if(!is.null(driving.force.cutoff)){
        selected.driver.genes <- as.character(driving.force[driving.force$Score > driving.force.cutoff,'GeneID'])
        dim(selected.driver.genes) <- c(length(selected.driver.genes),1)
        cat('Driving force cutoff = ',driving.force.cutoff, '\n')
        pdf('driver_genes_selection.pdf')
        plot(seq(nrow(driving.force)),driving.force$Score, xlab='Gene#', 
             ylab=expression(paste('driving force (', Delta, 'log2exprs)',sep=' ')),
             col=ifelse((driving.force$Score) >driving.force.cutoff, 'red','black'),
             main=paste0('No. of genes selected = ', nrow(selected.driver.genes)))
        dev.off()
        
    }
    
    ## Based on elbow detection method
    else{
        selected.driver.genes <- elbow_detection(scores = driving.force, baseNm = 'driver.genes.selection',
                                                 ylab.text = 'Driving force', xlab.text = 'Gene#')
        colnames(selected.driver.genes)
        if(nrow(selected.driver.genes)==0) stop('Unable to get turning point while identifying driver genes!\n')

    }
    
    cat('No. of genes selected for final sorting = ', nrow(selected.driver.genes),'\n')
    write.table(data.frame('GeneID'=selected.driver.genes, 'GroupID'='untitled'),file='selected.driver.genes.txt',row.names = F, col.names = T, quote=T, sep='\t')
    return(data.frame('GeneID'=selected.driver.genes))
}
