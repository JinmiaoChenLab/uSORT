#' wrapper of all avaliable sorting methods in uSORT
#'
#' Sorting methods include \code{autoSPIN}, \code{sWanderlust}, \code{monocle},
#' \code{Wanderlust}, \code{SPIN}. Any of the sorting method can be called directly
#' using this funciton.
#'
#' @param data Input preprocessed data matrix with row.name of cells and col.name of genes.
#' @param data_raw Input raw data matrix with row.name of cells and col.name of genes, for monocle method.
#' @param method The name of the sorting method to use, including \code{autoSPIN}, \code{sWanderlust}, \code{monocle}, \code{Wanderlust}, \code{SPIN} and \code{none}.
#' @param data_type The type of the data, either \code{linear} or \code{cyclical}.
#' @param SPIN_option The runing option of SPIN, \code{STS} or \code{neighborhood}.
#' @param SPIN_sigma_width Sigma width for SPIN.
#' @param autoSPIN_alpha alpha for autoSPIN.
#' @param autoSPIN_randomization Number of randomization for autoSPIN.
#' @param wanderlust_start_cell The id of the starting cell for wanderlust.
#' @param wanderlust_dfmap_components The number of components from diffusionmap for wanderlust.
#' @param wanderlust_l The number of nearest neighbors used for wanderlust.
#' @param wanderlust_num_waypoints The number of waypoints for wanderlust.
#' @param wanderlust_waypoints_seed The seed for reproducible analysis.
#' @param wanderlust_flock_waypoints The bumber of flock times for wanderlust.
#'
#' @return return the order of sorting results.
#' @export
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' exprs <- uSORT_preProcess(exprs_file = file)
#' exp_trimmed <- t(exprs$exprs_log_trimed)
#' PCA_selected_genes <- pca_gene_selection(exp_trimmed)
#' exp_PCA_genes <- exp_trimmed[, PCA_selected_genes]
#' #order <- uSORT_sorting_wrapper(data = exp_PCA_genes, method = 'autoSPIN')
uSORT_sorting_wrapper <- function(data, data_raw,
    method = c("autoSPIN", "sWanderlust", "monocle", "Wanderlust",
               "SPIN", "none"),
    data_type = c("linear", "cyclical"),
    SPIN_option = c("STS", "neighborhood"),
    SPIN_sigma_width = 1,
    autoSPIN_alpha = 0.2,
    autoSPIN_randomization = 20,
    wanderlust_start_cell = NULL,
    wanderlust_dfmap_components = 4,
    wanderlust_l = 15,
    wanderlust_num_waypoints = 150,
    wanderlust_waypoints_seed = 2711,
    wanderlust_flock_waypoints = 2) {

    data_type <- match.arg(data_type)
    method <- match.arg(method)
    data_type <- match.arg(data_type)
    SPIN_option <- match.arg(SPIN_option)
    switch(method, autoSPIN = {
        res <- autoSPIN(data, data_type = data_type,
                        sorting_method = SPIN_option,
            sigma_width = SPIN_sigma_width, alpha = autoSPIN_alpha,
            no_randomization = autoSPIN_randomization)
        order <- res$SampleID
    }, SPIN = {
        res <- SPIN(data = data, sorting_method = SPIN_option,
            sigma_width = SPIN_sigma_width)
        order <- res$SampleID

    }, monocle = {
        order <- monocle_wrapper(log2_exp = data,
                                 expression_data_raw = data_raw)
    }, Wanderlust = {
        res <- wanderlust_wrapper(data, s = wanderlust_start_cell,
            diffusionmap_components = wanderlust_dfmap_components,
            l = wanderlust_l, k = wanderlust_l,
            num_waypoints = wanderlust_num_waypoints,
            waypoints_seed = wanderlust_waypoints_seed,
            flock_waypoints = wanderlust_flock_waypoints)
        order <- res$Order
    }, sWanderlust = {
        res <- sWanderlust(data, data_type = data_type,
                           alpha = autoSPIN_alpha,
            SPIN_option = SPIN_option, sigma_width = SPIN_sigma_width,
            diffusionmap_components = wanderlust_dfmap_components,
            l = wanderlust_l, k = wanderlust_l,
            num_waypoints = wanderlust_num_waypoints,
            waypoints_seed = wanderlust_waypoints_seed,
            flock_waypoints = wanderlust_flock_waypoints)
        order <- res
    }, none = {
        order <- NULL
    })

    return(as.character(order))
}



