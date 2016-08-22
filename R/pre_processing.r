#' A data pre-processing function
#' 
#' A data pre-processing function removes genes which are lowly expressed and performs PCA analysis
#' based on FluidigmSC's autoAnalysis function.
#' @param exp a fluidigmSC's EXP object output from identifyOutliers function
#' @return a fluidigmSC's EXP object which is ready for downstream analysis
#' @author MaiChan Lau
#' @export

pre_processing<-function(exp = NULL){
  
  cat ('Data pre-processing (using SINGuLAR tool)...... \n')
  data_type <- as.character(exp$data_type[1, 1])
  lod <- as.numeric(as.character(exp$lod[1, 1]))
  
  ## Remove lowly expressed genes
  cutoff <- lod * 2
  exp <- fluidigmSC::removeGenesByLinearExp(exp, cutoff)
  
  return(exp)
}
