itn_anova_by_factor<-function (exp, aov_factor = "", group_difference = TRUE, p_FDR = TRUE) 
{
    expression_data <- exp$log2ex_data
    pValues <- apply(expression_data, 1, itn_anova_x, aovFactor = aov_factor, 
                     Tukey_HSD = TRUE, colName = FALSE)
    if (sum(is.na(pValues)) > 0) 
        pValues[is.na(pValues)] <- 1
    pValues <- t(pValues)
    if (group_difference) {
        diff <- apply(expression_data, 1, itn_anova_x, aovFactor = aov_factor, 
                      Tukey_HSD = TRUE, colName = FALSE, Difference = TRUE)
        diff <- t(diff)
    }
    x1 <- as.numeric(expression_data[1, ])
    colNames <- itn_anova_x(x1, aovFactor = aov_factor, Tukey_HSD = TRUE, 
                            colName = TRUE)
    colnames(pValues) <- colNames
    if (!p_FDR) {
        colnames(pValues) <- paste(colnames(pValues), "pValue", 
                                   sep = ".")
        pFDR <- pValues
    }
    else {
        pFDR <- c()
        for (j in ncol(pValues):1) {
            FDR <- p.adjust(pValues[, j], method = "BH")
            pFDR <- cbind(pValues[, j], FDR, pFDR)
            colnames(pFDR)[c(1, 2)] <- paste(colnames(pValues)[j], 
                                             c("pValue", "FDR"), sep = ".")
            if (j > 1 && group_difference) {
                pFDR <- cbind(diff[, j - 1], pFDR)
                colnames(pFDR)[1] <- colnames(pValues)[j]
            }
        }
    }
    Gene <- rownames(expression_data)
    pFDR <- cbind(Gene, pFDR)
    rownames(pFDR) <- 1:nrow(pFDR)
    return(pFDR)
}
