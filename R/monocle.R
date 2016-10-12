#' A wrapper function for Monocle sorting method
#'
#' @param log2_exp An log2 transformed expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param expression_data_raw A data frame containing raw expression values, with rownames of cells
#' and colnames of genes.
#' @param lod A value of limit of detection in the unit of TPM/CPM/RPKM.
#'
#' @importFrom Biobase fData pData fData<-
#' @importFrom monocle setOrderingFilter orderCells
#' @importFrom monocle reduceDimension
#' @return A data frame containing single column of ordered sample IDs.
#' @export
#'
#' @examples
#' set.seed(15)
#' da <- iris[sample(150, 150, replace = FALSE), ]
#' rownames(da) <- paste0('spl_',seq(1,nrow(da)))
#' d <- da[,1:4]
#' dl <- da[,5,drop=FALSE]
#' #res <- monocle_wrapper(log2_exp = d, expression_data_raw = d)
#' #dl <- dl[match(res,rownames(dl)),]
#' #annot <- data.frame(id = seq(1,length(res)), label=dl, stringsAsFactors = FALSE)
#' #ggplot(annot, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
monocle_wrapper <- function(log2_exp, expression_data_raw, lod = 1) {
    g <- as.character(colnames(log2_exp))
    cds <- EXP_to_CellDataSet(log2_exp = log2_exp,
                              expression_data_raw = expression_data_raw,
                              lod = lod)
    f <- fData(cds)
    f$gene_short_name <- rownames(f)
    fData(cds) <- f
    cds <- setOrderingFilter(cds, ordering_genes = g)
    cds <- reduceDimension(cds, use_irlba = FALSE)
    cds <- orderCells(cds, num_paths = 1, reverse = FALSE)
    pseudotime <- pData(cds)
    pseudotime <- pseudotime[order(pseudotime$Pseudotime), , drop = FALSE]
    order <- rownames(pseudotime)
}


#' A function for constructing a Monocle's CellDataSet object from an expression matrix
#'
#' @param log2_exp An log2 transformed expresssion matrix containing n-rows of cells and m-cols of genes.
#' @param expression_data_raw A data frame containing raw expression values, with rownames of cells
#' and colnames of genes.
#' @param lod A value of limit of detection in the unit of TPM/CPM/RPKM.
#'
#' @importFrom monocle newCellDataSet
#' @return A CellDataSet object.
EXP_to_CellDataSet <- function(log2_exp = NULL, expression_data_raw = NULL,
    lod = 1) {
    log2_exp <- t(log2_exp)
    expression_data_raw <- t(expression_data_raw)
    samples <- data.frame(SampleID = as.character(colnames(log2_exp)))
    genes <- data.frame(GeneID = as.character(rownames(log2_exp)))

    exprs <- expression_data_raw
    exprs[exprs < lod] <- lod
    exprs <- exprs[, as.character(samples$SampleID)]
    exprs <- exprs[as.character(genes$GeneID), ]

    Sample_sheet <- data.frame(samples, Pseudotime = seq(1, ncol(exprs)))
    rownames(Sample_sheet) <- Sample_sheet$SampleID
    # Sample_sheet$SampleID <- NULL
    pd <- new("AnnotatedDataFrame", data = Sample_sheet)

    Feature_sheet <- genes
    rownames(Feature_sheet) <- Feature_sheet$GeneID
    # Feature_sheet$GeneID <- NULL
    fd <- new("AnnotatedDataFrame", data = Feature_sheet)

    res <- newCellDataSet(as.matrix(exprs),
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.1)
    return(res)


}
