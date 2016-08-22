#' A heatmap plotting function
#'
#' A heatmap plotting function 
#'
#' @param resDir the directory to which output plot will be stored
#' @param sorted_exp a EXP object with cells sorted 
#' @param gene.clusterID a dataframe consists of a column 'GeneID' for a list of
#' genes & a column of 'ClusterID' for the corresponding cluster assignment of
#' each gene
#' @param flipping a logical value stating if the final cell ordering is flipped
#' for heatmap plotting
#' @author MaiChan Lau
#' @export
gene_cluster_heatmap<- function (resDir=NULL, sorted_exp = NULL, gene.clusterID = NULL,
                                 flipping = FALSE)
{
    log2ex_data<- sorted_exp$log2ex_data
    gene.clusterID <- gene.clusterID[order(gene.clusterID$ClusterID),]
    subset_log2ex_data<-log2ex_data[ match(gene.clusterID$GeneID, rownames(log2ex_data)),]
    
    ## Reverse cell ordering
    if(flipping== TRUE){
        cat('cell ordering is flipped for heatmap\n')
        rev_id <- rev(seq(ncol(subset_log2ex_data)))
        subset_log2ex_data <- subset_log2ex_data[,rev_id]
    }
    
    ## Define colors & legend
    num_colors <- 256
    colors <- c("blue",'white', "red" )
    colors <- colorRampPalette(colors)(num_colors)
    k<-length(unique(gene.clusterID$ClusterID))
    colorPal <- cm.colors(k)
    
    geneClus_color<-vector(mode='character',length=nrow(gene.clusterID))
    for(i in 1:length(geneClus_color)){
        for (j in 1:k){
            if(gene.clusterID$ClusterID[i]==j) geneClus_color[i]<-colorPal[j]
        }
    }
    
    legendNm <- c()
    for (i in 1:k){legendNm[i]<- paste0("GeneGp#", i)}

    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 40)
    log2exprs_z <- apply(subset_log2ex_data,1,scale)
    rownames(log2exprs_z) <- colnames(subset_log2ex_data)
    log2exprs_z <- t(log2exprs_z)
    
    if(flipping== TRUE)plotNm <- paste0('driver.gene.clusters.heatmap.flipped.pdf')
    else plotNm <- paste0('driver.gene.clusters.heatmap.pdf')
    pdf(plotNm)
    heatmap.2(as.matrix(log2exprs_z),col=my_palette,dendrogram='none',trace='none',
              Rowv=F,Colv=F,scale = 'none',cexRow=1.0,  RowSideColors=geneClus_color,
              breaks=seq(-2,2,0.1),density.info="none")
    legend("topright",legend=legendNm,
           col=colorPal,pch=20,horiz=T, inset=c(0,-0.01),
           pt.cex=1.5, bty = 'n')
    dev.off()
    
}

