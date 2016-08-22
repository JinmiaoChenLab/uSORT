#' A modified monocle's function
#' 
#' A modified monocle's function for 'compareModels' which identifies and removes genes 
#' whose reduced_models is better than full_models in term of likelihood
#' 
#' @param full_models a Monocle's vgam full model 
#' @param reduced_models a Monocle's vgam reduced/ null model
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM lrtest
#' @export

clusterGenes <- function (expr_matrix, krange, method = function(x) 
{
    as.dist((1 - cor(t(x)))/2)
}, ...) 
{
    expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 
                                   0, ]
    expr_matrix <- t(scale(t(log10(expr_matrix))))
    expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) == 
                                   FALSE, ]
    expr_matrix[is.na(expr_matrix)] <- 0
    n <- method(expr_matrix)
    pamk.best <- fpc::pamk(n, krange=krange)
    
    #     pdf(paste0('pam_optimalk_qval',qval_cutoff,'_avg_silwidth.pdf'))
    #     barplot(pamk.best$crit, xlab = 'No. of cluster', ylab = 'Average silhouette width',
    #             names.arg = c(1,krange))
    #     dev.off()
    
    cat('optimal k =', pamk.best$nc, '\n')
    clusters <- pam(n, pamk.best$nc, ...)
    #plot(clusters)
    class(clusters) <- "list"
    clusters$exprs <- expr_matrix
    clusters
}


#' A modified monocle's function
#' 
#' A modified monocle's function for 'compareModels' which identifies and removes genes 
#' whose reduced_models is better than full_models in term of likelihood
#' 
#' @param full_models a Monocle's vgam full model 
#' @param reduced_models a Monocle's vgam reduced/ null model
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM lrtest
#' @export

compareModels <- 
    function (full_models, reduced_models) 
    {
        stopifnot(length(full_models) == length(reduced_models))
        test_res <- mapply(function(x, y) {
            if (is.null(x) == FALSE && is.null(y) == FALSE) {
                if(logLik(x) < logLik(y)){
                    data.frame(status = "FAIL", pval = 1)
                }
                else {
                    lrt <- VGAM::lrtest(x, y)
                    pval = lrt@Body["Pr(>Chisq)"][2, ]
                    data.frame(status = "OK", pval = pval)
                }
            }
            else {
                data.frame(status = "FAIL", pval = 1)
            }
        }, full_models, reduced_models, SIMPLIFY = FALSE, USE.NAMES = TRUE)
        test_res <- do.call(rbind.data.frame, test_res)
        test_res$qval <- p.adjust(test_res$pval, method = "BH")
        test_res
    }


#' A modified monocle's helper function 
#' 
#' A modified monocle's function for 'diff_test_helper' which includes more attempts 
#' on finding models and also compute max. magnitude change in expression values predicted
#' by GLM model
#' 
#' @param x an expression data
#' @param fullModelFormulaStr a Monocle's model structure
#' @param reducedModelFormulaStr a Monocle's model structure
#' @param expressionFamily a Monocle's family character
#' @param lowerDetectionLimit a threshold value  
#' @param type_ordering a character indicating the type of underlying cell progression, 
#' i.e. linear or circular
#' @return test_res a dataframe containing status of modeling and adjusted p-value
#' @author MaiChan Lau
#' @importFrom VGAM vgam
#' @export
diff_test_helper <-
    function (x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, 
              lowerDetectionLimit=0.1, type_ordering = 'linear') 
    {
        leftcensored <- x < lowerDetectionLimit
        x[x<lowerDetectionLimit] <- lowerDetectionLimit
        if (expressionFamily@vfamily %in% c("zanegbinomialff", "negbinomial", 
                                            "poissonff", "quasipoissonff")) {
            expression <- round(x)
            integer_expression <- TRUE
        }
        else if (expressionFamily@vfamily %in% c("gaussianff")) {
            expression <- x
        }
        else {
            expression <- log10(x)
            integer_expression <- FALSE
        }
        
        
        p <- parent.env(environment())
        #print(ls(p))
        
        for (i in 1:20){
            test_res <- tryCatch({
                
                full_model_fit <- suppressWarnings(VGAM::vgam(as.formula(fullModelFormulaStr), 
                                                              family = expressionFamily))
                reduced_model_fit <- suppressWarnings(VGAM::vgam(as.formula(reducedModelFormulaStr), 
                                                                 family = expressionFamily))
                
                
                if (integer_expression) {
                    pred_res <- predict(full_model_fit, type = "response")
                    pred_res[pred_res < lowerDetectionLimit] <- lowerDetectionLimit
                }
                else {
                    pred_res <- 10^(predict(full_model_fit, type = "response"))
                    pred_res[pred_res < log10(lowerDetectionLimit)] <- log10(lowerDetectionLimit)
                }

                pred <- data.frame(Pseudotime = p$Pseudotime, expectation = log10(pred_res))
                pred <- pred[order(pred$Pseudotime),]
                cyclic.driving.force<-diff(range(pred$expectation))
                linear.driving.force<-abs(pred$expectation[1] -  pred$expectation[nrow(pred)])
#                 pred$scaled_expectation <- ((pred$expectation) - min(pred$expectation)) / ((max(pred$expectation) - min(pred$expectation)))
#                 pred$scaled_expectation <- pred$scaled_expectation  * diff(range(expression))               
#                 cyclic.driving.force<-diff(range(pred$scaled_expectation))
#                 linear.driving.force<-abs(pred$scaled_expectation[1] -  pred$scaled_expectation[nrow(pred)])

# pdf('test.pdf')
# plot(p$Pseudotime, log10(pred_res))
# dev.off()
# pdf('scaled.pdf')
# plot(p$Pseudotime, pred$scaled_expectation)
# dev.off()
# save(list=ls(),file = 'variables.RData')                
                res<-compareModels(list(full_model_fit), list(reduced_model_fit))
                res$cyclic.driving.force <- cyclic.driving.force
                res$linear.driving.force <- linear.driving.force
                if(is.na(cyclic.driving.force) | is.na(linear.driving.force)){
                    res$pval = 1; res$qval = 1; res$status = 'FAIL'}
                else res$status = 'OK'
                res
                
            }, error = function(e) {
                print(e)
                data.frame(status = "FAIL", pval = 1, qval = 1, cyclic.driving.force = 0, linear.driving.force=0)
            })
            if(test_res$status!= 'FAIL') break;
        }
        
        test_res
    }


