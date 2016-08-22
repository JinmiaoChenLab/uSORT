#' A modified function of fluidigmSC's identifyOutliers
#' 
#' A modified version of fluidigmSC's identifyOutliers function which removes the saving
#' feature. 
#' 
#' @param exp_file a character of input expression filename
#' @param lod limit of detection used for gene trimming
#' @return The updated EXP object 
#' @author fluidigmSC
#' @export

identifyOutliers <- function (exp_file = NULL, lod=1) 
{
    exp <- readLinearExp(exp_file, lod = lod)
    
    data_type <- as.character(exp$data_type[1, 1])
    threshold <- 0
    step <- 4
    fine_step <- 1
    num_fine_test <- 4
    pct_goodsample_threshold <- 0.5
    quantile_threshold <- 0.95
    low_quantile_threshold <- 0.15
    min_gene_number <- 20
    if (data_type == "LinearExp") 
        threshold <- as.numeric(exp$lod[1, 1])
    else threshold <- 1
    sample_list <- exp$sample_list
    sample_groups <- levels(sample_list[, 2])
    current_group_no <- 1
    cat("Identifying outliers ...\n")
    flush.console()
    exp_trimmed_all <- list()
    exp_threshold_all <- list()
    exp_outliers_all <- list()
    for (i in 1:length(sample_groups)) {
        sample_group_id <- sample_groups[i]
        exp_group <- retainSampleGroup(exp, sample_group_id)
        exp_group_parameter <- itn_identifyExpOutliers(exp_group, 
                                                       exp_file, threshold, step, fine_step, num_fine_test, 
                                                       pct_goodsample_threshold, quantile_threshold, low_quantile_threshold, 
                                                       min_gene_number)
        exp_trimmed_all[[i]] <- exp_group_parameter$exp_trimmed
        exp_threshold_all[[i]] <- exp_group_parameter$threshold
        exp_outliers_all[[i]] <- as.data.frame(exp_group_parameter$outliers)
        if (length(exp_outliers_all[[i]]) == 0) 
            exp_outliers_all[[i]] <- NA
        if (nrow(exp_group$outlier_list) > 0) {
            exp_trimmed_all[[i]] <- restoreAllOutlier(exp_trimmed_all[[i]])
            exp_outliers_all[[i]] <- exp_group$outlier_list
        }
        names(exp_outliers_all)[i] <- sample_groups[i]
    }
    
    while (current_group_no > 0) {
        sample_group_id <- sample_groups[current_group_no]
        exp_group <- retainSampleGroup(exp, sample_group_id)
        outlier_list <- itn_displayOutlierGUI(exp_trimmed_all, 
                                              exp_file, exp_outliers_all, exp_threshold_all, current_group_no)
        if (length(outlier_list) > 1) 
            exp_outliers_all <- outlier_list[[2]]
        outlier_list <- outlier_list[[1]]
        if (class(outlier_list) == "integer") 
            current_group_no <- outlier_list
        else current_group_no <- 0
    }
    
    if (class(outlier_list) != "data.frame") {
        cat("\nCancelled\n")
        flush.console()
        invisible(exp_with_outliers)
    }
    else {
        if (nrow(outlier_list) > 0) {
            cat(paste("\n", "Total ", nrow(outlier_list), " outliers were identified.\n", 
                      sep = ""))
            sample_groups <- unique(as.vector(outlier_list[, 
                                                           2]))
            if (length(sample_groups) > 1) {
                for (sample_group in sample_groups) {
                    outlier_num <- nrow(outlier_list[is.element(outlier_list[, 
                                                                             2], sample_group), ])
                    cat(paste("    ", outlier_num, " outliers were identified in sample group of ", 
                              sample_group, "\n", sep = ""))
                }
            }
            flush.console()
            if (nrow(exp$outlier_list) > 0) 
                exp <- restoreOutlierFromList(exp, exp$outlier_list)
            exp_with_outliers <- itn_updateExpByOutliers(exp, 
                                                         outlier_list = outlier_list, action = "add")
            save_title <- "Save Outlier-Removed Expression Data ..."
            save_msg <- "\nExpression data with identified outliers is saved.\n"
        }
        else {
            cat("\nNo outlier was identified.\n")
            flush.console()
            exp_with_outliers <- exp
            save_title = "Save Expression Data ..."
            save_msg <- "\nExpression data is saved.\n"
        }
        
    }
    invisible(exp_with_outliers)
}
    
