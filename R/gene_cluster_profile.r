#' A function generating gene subsets
#'
#' A function generating gene subsets
#'
#' @param exp an EXP object 
#' @param geneCluster a data.frame containing clusterID of a set of interesting genes
#' @param cell.ordering a data.frame of cell ordering; in the format of 
#' sample_list of EXP object
#' @param resDir an output directory
#' @author MaiChan Lau
#' @importFrom gplots heatmap.2
#' @export
gene_cluster_profile <- function(exp = NULL, cell.ordering = NULL, geneCluster = NULL,
                                 resDir=NULL){

    setwd(resDir)
    exp <- updateSampleListFromList(exp, cell.ordering)
    colnames(geneCluster) <- c('GeneID', 'ClusterID')
    
    k <- unique(geneCluster$ClusterID)
    
    for (i in 1:length(k)){
        
        cat('Plotting gene cluster# ',k[i],'\n')
        geneSet <- as.character(geneCluster[geneCluster$ClusterID==k[i],'GeneID'])
        log2ex<- exp$log2ex_data
        geneSet_log2ex <- log2ex[rownames(log2ex) %in% geneSet,]
        N <- ceiling(nrow(geneSet_log2ex)*0.1)
        geneSet_log2ex <-geneSet_log2ex[sample(nrow(geneSet_log2ex), N) ,]
        fileNm <- paste0('geneSet',k[i],'.Profile.pdf')        
        pdf(fileNm)
        heatmap.2(as.matrix(geneSet_log2ex),dendrogram='row',trace='none',
                  Rowv=T,Colv=F,scale = 'row',cexRow=1.0, )
        dev.off()
    }
    
    
    
}

