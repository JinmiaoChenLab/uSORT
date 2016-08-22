#' @export
find_developing_genes <- function(res_dir = NULL, ordering_fname = NULL, baseNm = 'postanalysis',
                                  qval_cutoff = 0.05, gene_clustering = FALSE, rev_plot_gene_cluster = FALSE,
                                  krange = 2:10){
    setwd(res_dir)
    cat('\nFinding genes whose expression changed with the uSPIN ordering ......\n')
    all.genes <- read.table('nondropout.genes.txt',header=T)
    preprocessed_exp<-fluidigmSC::readExpObject('preprocessed_exp.fso')
    sorted_exp <- fluidigmSC::updateSampleListFromFile(preprocessed_exp, ordering_fname)
    sorted_exp <- fluidigmSC::updateGeneListFromList(sorted_exp, all.genes)
    cds <- EXP_to_CellDataSet(EXP = sorted_exp)
    cds <- monocle::detectGenes(cds, min_expr = 0.1)
    
    max_dropout_percent = 0.2
    trimmed.genes <- rownames(cds)[cds@featureData@data$num_cells_expressed >= ncol(cds)*max_dropout_percent]
    trimmed.cds <- cds[trimmed.genes,]
    noGene <- nrow(trimmed.cds) 
    #cat('\nRemoving sparsely expressed genes ...... #gene remained = ', noGene,'\n')
    assignInNamespace("compareModels",compareModels,ns="monocle")
    assignInNamespace("diff_test_helper",diff_test_helper,ns="monocle")
    
    #nCores <- parallel::detectCores()
    nCores <- 5
    full.models <- monocle::fitModel(trimmed.cds,  cores=nCores,
                                     modelFormulaStr="expression~VGAM::sm.ns(Pseudotime, df=3)")
    null.models <- monocle::fitModel(trimmed.cds, modelFormulaStr="expression~1")
    diff_test_res <- compareModels(full.models, null.models)
    write.table(diff_test_res,file=paste0(baseNm,'.developing.genes.pval.txt'),
                row.names = T, col.names = NA, quote=T, sep='\t')
    
    passQval_diff_test_res <- rownames(diff_test_res)[diff_test_res$qval < qval_cutoff]
    fileNm <- paste0(baseNm,'.developing.genes.qval',qval_cutoff,'.txt')
    write.table(data.frame('GeneID'=passQval_diff_test_res),
                file=fileNm,row.names = F, col.names = F, quote=T, sep='\t')   
    
    if (gene_clustering ==  TRUE){
        library(fpc)
        library(cluster)
        cat('Clustering driver genes @q-val <', qval_cutoff,'......\n')
        passQval_fullmodels<-full.models[names(full.models)%in% as.character(passQval_diff_test_res)]
        expression_curve_matrix <- monocle::responseMatrix(passQval_fullmodels)
        #gene.clusters <- monocle::clusterGenes(expression_curve_matrix, k=k)
        assignInNamespace("clusterGenes",clusterGenes,ns="monocle")
        gene.clusters <- clusterGenes(expr_matrix = expression_curve_matrix, 
                                      krange = krange)
        plot_width <- length(unique(gene.clusters$clustering))*4
        
        gene.clusters.profiles <- monocle::plot_clusters(cds, gene.clusters)
        fileNm <-paste0(baseNm,'.regression.profiles.qval',qval_cutoff,'.pdf')
        ggplot2::ggsave(fileNm, plot=gene.clusters.profiles, width=plot_width)
        fileNm <-paste0(baseNm,'.clusterID.qval',qval_cutoff,'.txt')
        gene.clustersID <-data.frame('GeneID'=names(gene.clusters$clustering),'ClusterID'=gene.clusters$clustering)
        write.table(gene.clustersID,file=fileNm,row.names=F,col.names=T,sep='\t')
        
        if(rev_plot_gene_cluster == TRUE){
            cat('Plotting for flipped regression profiles\n')
            
            rev_ordering <- apply(sorted_exp$sample_list,2,rev)
            
            rev_sorted_exp <- fluidigmSC::updateSampleListFromList(sorted_exp,rev_ordering)
            rev_cds <- EXP_to_CellDataSet(EXP = rev_sorted_exp)
            
            gene.clusters.profiles <- monocle::plot_clusters(rev_cds, gene.clusters)
            fileNm <-paste0(baseNm,'.regression.profiles.qval',qval_cutoff,'.rev.pdf')
            ggplot2::ggsave(fileNm, plot=gene.clusters.profiles, width = plot_width)
        }
    }
}
