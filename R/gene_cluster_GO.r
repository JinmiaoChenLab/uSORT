#' A gene ontology function
#'
#' A gene ontology function performs gene set enrichment analysis using topGO package
#'
#' @param resDir the directory to which output results will be written
#' @param GO_controls parameters for gene enrichment study
#' @param gene.clusterID a dataframe consists of a column 'GeneID' for a list of
#' genes & a column of 'ClusterID' for the corresponding cluster assignment of
#' each gene
#' @param universal_geneSet a list of universal or backgroupd genes
#' @author MaiChan Lau
#' @export
gene_cluster_GO <- function(resDir=NULL, GO_controls = NULL, universal_geneSet = NULL,
                            geneClusterID_fname = NULL, baseNm = 'GO'){
    if(is.null(GO_controls)) stop("No GO specifications!")
    if(is.null(universal_geneSet)){
       setwd(resDir)
        preprocessed_exp<-fluidigmSC::readExpObject('preprocessed_exp.fso')
        universal_geneSet <- data.frame('GeneID'=rownames(preprocessed_exp$org_data))
        
    }    
    gene.clusterID <- read.table(geneClusterID_fname,header=T)
    
    
    bground_gene <- data.frame('GeneID'= universal_geneSet$GeneID[!universal_geneSet$GeneID %in% gene.clusterID$GeneID],
                               'ClusterID'=0)
    universe_genes <- rbind(bground_gene, gene.clusterID)
    cat('Universe genes = ', nrow(universe_genes), '\n')
    cat('Gene enrichment ......\n')
    topGO_geneCluster(geneList = universe_genes, GO_controls=GO_controls, resDir=getwd(), baseNm=baseNm)
    

       
}

