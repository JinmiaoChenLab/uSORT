#' A gene enrichment function
#'
#' A gene enrichment function which is based on topGO package. It identifies enriched function 
#' of each gene clusters.
#' 
#' @param geneList a data.frame of gene names and cluster assignment
#' @param GO_controls a data.frame containing the choice of the ontology type & mapping
#' @param resDir the directory to which output results will be written
#' @author MaiChan Lau
#' @import topGO
#' @export
topGO_geneCluster <- 
    function(geneList = NULL, GO_controls = NULL, resDir=NULL, baseNm = 'GO')
    {
        library(topGO)
        geneList_input <-geneList$ClusterID
        names(geneList_input)<-geneList$GeneID
        k <- unique(geneList_input)
        
        for (i in 1:length(k)){
            if(k[i]==0) next;
            cat('GO @cluster# ',k[i],'\n')
            
            geneSelectionFun <-function(clusterID){return(clusterID == k[i])}
            GOdata <- new("topGOdata",
                          ontology = as.character(GO_controls$onto),
                          allGenes = geneList_input,geneSelectionFun=geneSelectionFun,
                          nodeSize = 5,
                          annot = annFUN.org,
                          mapping = as.character(GO_controls$mapping),
                          ID = as.character(GO_controls$ID))   
            
            #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
            #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
            resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
            
#             allRes <- GenTable(GOdata, resultFis=resultFis, classicKS = resultKS, elimKS = resultKS.elim,
#                                orderBy = "resultFis", topNodes = 10)
            allRes <- GenTable(GOdata, resultFis=resultFis, 
                               orderBy = "resultFis", topNodes = 10)
            outputFn <- paste0(baseNm,'_cluster',k[i],'.GO.txt')
            setwd(resDir)
            write.table(allRes, file=outputFn, row.names = F, col.names = T, sep='\t')
        }
        
        
    }

