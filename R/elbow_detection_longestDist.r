#' A elbow detection function
#' 
#' A elbow detection function detects the elbow/ knee for gene selection, either
#' using PCA loading values or driving forces. If both elbow & knee are detected
#' elbow will be returned.
#' @param scores a data.frame constaining gene names in 'GeneID' column and scores (loading values or
#' driving forces) in the 'Score' column.
#' @param boundary_cut a fraction (at both sides) of the boundary to be omitted for 
#' differentiating elbow from knee   
#' the assumed test statistic / distribution of the data. Currently only "Normal" and "CSS" supported.
#' @param baseNm a character for naming output files
#' @return a list of selected genes
#' @author Mai Chan
#' @export

elbow_detection <- function(scores = NULL, boundary_cut = 0.1, baseNm = 'auto.elbow.detection',
                            ylab.text = 'Scores', xlab.text = 'Index'){
    
    original_colnames <- colnames(scores)
    colnames(scores) <- c('X','Y')
    sorted_scores <- scores[order(scores$Y, decreasing = T),]
    
    nPoints = nrow(sorted_scores)
    allCoord <- data.frame('no'=seq(nPoints), 'scores'=sorted_scores$Y)
    
    firstPoint <- allCoord[1,]
    lineVec <- allCoord[nPoints,] - firstPoint
    lineVecN <- lineVec / sqrt(sum(lineVec^2))
    
    minus <- function(x,y){x-y}
    vecFromFirst <- t(apply(as.matrix(allCoord), MARGIN = 1, minus, y = as.matrix(firstPoint)))
    scalarProduct <- apply(vecFromFirst, MARGIN = 1, pracma::dot, y=as.numeric(lineVecN))
    dim(scalarProduct) <- c(length(scalarProduct),1)
    vecFromFirstParallel = t(apply(scalarProduct, MARGIN = 1, function(x){x*as.numeric(lineVecN)}))
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine <- apply(vecToLine, MARGIN = 1, function(x){sqrt(sum(x^2))})
    
#     s <- sign(vecToLine)[,2]
#     m<-changepoint::cpt.mean(s)
#     plot(s)
#     sign_cpt <- changepoint::cpts(m)
#     if(length(sign_cpt) ==0 )idxOfBestPoint <-which.max(distToLine)
#     else{
#         if((sign_cpt/ nPoints > 0.1) | (sign_cpt/ nPoints < 0.9)){
#             if(mean(s[1:sign_cpt]) < 0) idxOfBestPoint = which.max(distToLine[1:sign_cpt])
#             else if (mean (s[sign_cpt: nPoints]) <0) idxOfBestPoint = which.max(distToLine[sign_cpt: nPoints])
#         }
#         else idxOfBestPoint <-which.max(distToLine)
#     }

    idxOfBestPoint <-which.max(distToLine)
    score_cutoff <- sorted_scores$Y[idxOfBestPoint]
    opt_no <- data.frame('X' = sorted_scores$X[sorted_scores$Y >= score_cutoff])
    colnames(opt_no) <- original_colnames[1]

#     if(length(grep("driver", baseNm))>0) {
#         ylab.text <- expression(paste('driving force (', Delta, 'log2exprs)',sep=' '))
#     } else if(length(grep("PCA", baseNm))>0) {
#         ylab.text <- 'sum of PC1-3 loadings'
#     } else ylab.text <- 'scores'
    
    pdf(paste0(baseNm,'_longestDist.pdf'))
    plot(seq(nrow(scores)),sorted_scores$Y, col=ifelse(sorted_scores$Y >= score_cutoff,'red','black'),
         xlab=xlab.text, ylab=ylab.text, pch="*", 
         main=paste0("Optimal no. selected = ", nrow(opt_no), "@cutoff ", round(score_cutoff, digits=4))) 
    dev.off()
    
    return(opt_no)
    
}
