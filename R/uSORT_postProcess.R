#' Resluts parsing for uSORT
#'
#' Save result object into a RData file. Save cell to cell distance heatmap for
#' both preliminary and refined results. Creat plot of driver gene profiles on
#' final ordering using heatmap.
#'
#' @param uSORT_results Result object from \code{uSort} function, a list.
#' @param project_name A prefix for the saving files.
#' @param result_directory The path where to save the results.
#'
#' @importFrom gplots heatmap.2 colorpanel bluered
#' @return save the results.
#' @export
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' #remove the # symbol of the following codes to test
#' #uSORT_results <- uSORT(exprs_file = file,
#' # project_name = 'test',
#' # preliminary_sorting_method = 'autoSPIN',
#' # refine_sorting_method = 'sWanderlust',
#' # save_results = FALSE)
#' #uSORT_write_results(uSORT_results,
#' # project_name = 'test',
#' # result_directory = getwd())
uSORT_write_results <- function(uSORT_results, project_name, result_directory) {
    if (!dir.exists(result_directory)) {
        dir.create(result_directory)
    }
    setwd(result_directory)
    save(uSORT_results, file = paste0(project_name, ".RData"))

    ## preliminary cell to cell distance heat map
    pre_data <- uSORT_results$trimmed_log2exp[uSORT_results$preliminary_sorting_order,
        uSORT_results$preliminary_sorting_genes]
    priliminary_c2c_dist <- distance.function(pre_data)
    pdf(paste0(project_name, "_distance_heatmap_preliminary.pdf"))
    heatmap.2(priliminary_c2c_dist, dendrogram = "none", trace = "none",
        col = colorpanel(10, "red", "yellow", "blue"), Rowv = FALSE,
        Colv = FALSE, labCol = NA, labRow = NA, keysize = 1.2, cexRow = 1.1,
        main = "Cell to cell distance")

    dev.off()

    ## refined cell to cell distance heat map
    ref_data <- uSORT_results$driverGene_refinedOrder_log2exp
    refined_c2c_dist <- distance.function(ref_data)
    pdf(paste0(project_name, "_distance_heatmap_refined.pdf"))
    heatmap.2(refined_c2c_dist, dendrogram = "none", trace = "none",
        col = colorpanel(10, "red", "yellow", "blue"), Rowv = FALSE,
        Colv = FALSE, labCol = NA, labRow = NA, keysize = 1.2, cexRow = 1.1,
        main = "Cell to cell distance")

    dev.off()

    ## Plot driver gene profiles on final ordering
    pdf(paste0(project_name, "_final_driver_genes_profiles.pdf"),
        width = 8)
    heatmap.2(as.matrix(t(uSORT_results$driverGene_refinedOrder_log2exp)),
        dendrogram = "row", trace = "none", col = bluered, Rowv = TRUE,
        Colv = FALSE, scale = "row", cexRow = 0.5, margins = c(8,
            8))
    dev.off()
}
