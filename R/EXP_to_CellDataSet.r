#' A data type conversion function
#' 
#' A data type conversion function transforms fluidigmSC's EXP object to Monocle's 
#' CellDataSet object.
#' 
#' @param EXP a fluidigmSC's EXP object 
#' @return a Monocle's CellDataSet object 
#' @author MaiChan Lau
#' @export

EXP_to_CellDataSet <- function(EXP = NULL, lod=1){
  
  ## Inputs check
  if (is.null(EXP)) stop("NO input!")
  
  ## Exprs data
  exprs <- EXP$org_data 
  exprs[exprs<lod] <- lod
  #exprs<-exprs[row.names(exprs)%in%selectedGenes[,1],]
  
  
  # Remove duplicated gene names
  exprs <- avg_same_genes(exprs)
  
  # Remove outlier cells
  exprs<-exprs[,as.character(EXP$sample_list$SampleID)]
  # Remove QC-failed genes
  exprs<-exprs[as.character(EXP$gene_list$GeneID),]
  
  
  ## PhenoData
  #Sample_sheet<-EXP$sample_list
  #Pseudotime <- data.frame('SampleID'=Sample_sheet$SampleID,'Pseudotime'=seq(1,ncol(exprs)))
  
  #Sample_sheet<-merge(Sample_sheet,Pseudotime,by='SampleID',sort=F)
  Sample_sheet<-data.frame(EXP$sample_list,'Pseudotime'=seq(1,ncol(exprs)))
  rownames(Sample_sheet) <- Sample_sheet$SampleID
  Sample_sheet$SampleID <- NULL
  pd<-new("AnnotatedDataFrame",data=Sample_sheet)
  
  
  ## FeatureData
  Feature_sheet <- EXP$gene_list
  rownames(Feature_sheet) <- Feature_sheet$GeneID
  Feature_sheet$GeneID <- NULL
  fd<-new("AnnotatedDataFrame",data=Feature_sheet)
  
#   print(dim(pd))
#   print(dim(fd))
#   print(dim(exprs))
#   tryCatch({
#     #res<-newCellDataSet(as.matrix(exprs),phenoData=pd,featureData=fd)  
#     res <- new("CellDataSet", exprs = as.matrix(exprs), phenoData = pd, featureData = fd)
#     res
#     #newCellDataSet(as.matrix(exprs),phenoData=pd,featureData=fd)  
#   },error = function(e) {
#     print(e)
#   })
  #res <- new("CellDataSet", exprs = as.matrix(exprs), phenoData = pd, featureData = fd)
  res<-monocle::newCellDataSet(as.matrix(exprs),phenoData=pd,featureData=fd,lowerDetectionLimit = 0.1)  
  return(res)    
  
  
}

