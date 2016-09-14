#' uSORT: A self-refining ordering pipeline for gene selection
#'
#' This package is designed to uncover the intrinsic cell progression path
#' from single-cell RNA-seq data.
#'
#' This package incorporates data pre-processing, preliminary PCA gene selection, preliminary
#' cell ordering, feature selection, refined cell ordering, and post-analysis interpretation
#' and visualization. The uSORT workflow can be implemented through calling the main function or the GUI.
#' \code{\link{uSORT}}.
#'
#' @docType package
#' @name uSORT
#'
NULL



#' The main function of uSORT pacakge
#'
#' The main function of \code{uSORT-pacakge} which provides a workflow of sorting scRNA-seq data.
#'
#'
#' @param exprs_file Input file name in txt format, with rownames of cells and colnames of genes.
#' @param log_transform Boolean, if log transform the data.
#' @param remove_outliers Boolean, if remove the outliers.
#' @param preliminary_sorting_method Method name for preliminary sorting, including \code{autoSPIN}, \code{sWanderlust}, \code{monocle}, \code{Wanderlust}, \code{SPIN}, or \code{none}.
#' @param refine_sorting_method Method name for refined sorting, including \code{autoSPIN}, \code{sWanderlust}, \code{monocle}, \code{Wanderlust}, \code{SPIN}, or \code{none}.
#' @param project_name A character name as the prefix of the saved result file.
#' @param result_directory The directory indicating where to save the results.
#' @param nCores Number of cores that will be employed for drive gene selection (parallel computing), default is 1.
#' @param save_results Boolean determining if save the results.
#' @param reproduce_seed A seed used for reproducing the result.
#' @param scattering_cutoff_prob Scattering cutoff value probability for gene selection, default 0.75.
#' @param driving_force_cutoff Driving force cutoff value for gene selection, default NULL(automatically).
#' @param qval_cutoff_featureSelection Q value cutoff for gene selection, default 0.05.
#' @param pre_data_type The data type which guides the autoSPIN sorting, including \code{linear}, \code{cyclical}.
#' @param pre_SPIN_option SPIN contains two options including \code{STS}(default), \code{neighborhood}.
#' @param pre_SPIN_sigma_width Sigma width parameter for SPIN, default is 1.
#' @param pre_autoSPIN_alpha alpha parameter for autoSPIN, default is 0.2.
#' @param pre_autoSPIN_randomization Number of randomizations for autoSPIN, default is 20.
#' @param ref_data_type The data type which guides the autoSPIN sorting, including \code{linear}, \code{cyclical}.
#' @param ref_SPIN_option SPIN contains two options including \code{STS}(default), \code{neighborhood}.
#' @param ref_SPIN_sigma_width Sigma width parameter for SPIN, default is 1.
#' @param ref_autoSPIN_alpha alpha parameter for autoSPIN, default is 0.2.
#' @param ref_autoSPIN_randomization Number of randomizations for autoSPIN, default is 20.
#' @param pre_wanderlust_start_cell The name of starting cell for wanderlust, default is the first cell from the data.
#' @param pre_wanderlust_dfmap_components Number of components from diffusion map used for wanderlust analysis, default is 4.
#' @param pre_wanderlust_l Number of nearest neighbors used for wanderlust, default is 15.
#' @param pre_wanderlust_num_waypoints Number of waypoint used for wanderlust, default is 150.
#' @param pre_wanderlust_waypoints_seed The seed for reproducing the wanderlust results.
#' @param pre_wanderlust_flock_waypoints The number of times for flocking the waypoints, default is 2.
#' @param ref_wanderlust_start_cell The name of starting cell for wanderlust, default is the first cell from the data.
#' @param ref_wanderlust_dfmap_components Number of components from diffusion map used for wanderlust analysis, default is 4.
#' @param ref_wanderlust_l Number of nearest neighbors used for wanderlust, default is 15.
#' @param ref_wanderlust_waypoints_seed The seed for reproducing the wanderlust results.
#' @param ref_wanderlust_flock_waypoints The number of times for flocking the waypoints, default is 2.
#' @param ref_wanderlust_num_waypoints Number of waypoint used for wanderlust, default is 150
#'
#' @export
#' @return results object (a list)
#' @seealso \code{\link{uSORT-package}}, \code{\link{uSORT_GUI}}
#'
#' @examples
#' dir <- system.file('extdata', package='uSORT')
#' file <- list.files(dir, pattern='.txt$', full=TRUE)
#' #remove the # symbol of the following codes to test
#' #uSORT_results <- uSORT(exprs_file = file, project_name = "test",
#' # preliminary_sorting_method = "autoSPIN",
#' # refine_sorting_method = "sWanderlust",
#' # save_results = FALSE)
uSORT <- function(exprs_file,
                  log_transform = TRUE,
                  remove_outliers = TRUE,
                  preliminary_sorting_method = c("autoSPIN", "sWanderlust", "monocle", "Wanderlust", "SPIN", "none"),
                  refine_sorting_method = c("autoSPIN", "sWanderlust", "monocle", "Wanderlust", "SPIN", "none"),
                  project_name = "uSORT",
                  result_directory = getwd(),
                  nCores = 1,
                  save_results = TRUE,
                  reproduce_seed = 1234,
                  # gene selection parameters
                  scattering_cutoff_prob = 0.75,
                  driving_force_cutoff = NULL,
                  qval_cutoff_featureSelection = 0.05,
                  # autoSPIN parameters for preliminary sorting
                  pre_data_type = c("linear", "cyclical"),
                  pre_SPIN_option = c("STS", "neighborhood"),
                  pre_SPIN_sigma_width = 1,
                  pre_autoSPIN_alpha = 0.2,
                  pre_autoSPIN_randomization = 20,
                  # Wanderlust parameters for preliminary sorting
                  pre_wanderlust_start_cell = NULL,
                  pre_wanderlust_dfmap_components = 4,
                  pre_wanderlust_l = 15,
                  pre_wanderlust_num_waypoints = 150,
                  pre_wanderlust_waypoints_seed = 2711,
                  pre_wanderlust_flock_waypoints = 2,
                  # autoSPIN parameters for refine sorting
                  ref_data_type = c("linear", "cyclical"),
                  ref_SPIN_option = c("STS", "neighborhood"),
                  ref_SPIN_sigma_width = 1,
                  ref_autoSPIN_alpha = 0.2,
                  ref_autoSPIN_randomization = 20,
                  # Wanderlust parameters for refine sorting
                  ref_wanderlust_start_cell = NULL,
                  ref_wanderlust_dfmap_components = 4,
                  ref_wanderlust_l = 15,
                  ref_wanderlust_num_waypoints = 150,
                  ref_wanderlust_flock_waypoints = 2,
                  ref_wanderlust_waypoints_seed = 2711){

    ## parameter check
    set.seed(reproduce_seed)
    preliminary_sorting_method <- match.arg(preliminary_sorting_method)
    refine_sorting_method <- match.arg(refine_sorting_method)
    pre_data_type <- match.arg(pre_data_type)
    ref_data_type <- match.arg(ref_data_type)
    pre_SPIN_option <- match.arg(pre_SPIN_option)
    ref_SPIN_option <- match.arg(ref_SPIN_option)

    ## data loading and pre-processing
    cat("-Data loading and pre-processing...\n")
    exprs <- uSORT_preProcess(exprs_file = exprs_file,
                              log_transform = log_transform,
                              remove_outliers = remove_outliers)
    exp_raw <- t(exprs$exprs_raw)
    exp_trimmed <- t(exprs$exprs_log_trimed)

    ## preliminary gene selection using PCA
    cat("-Preliminary gene selection using PCA...\n")
    PCA_selected_genes <- pca_gene_selection(exp_trimmed)
    exp_PCA_genes <- exp_trimmed[, PCA_selected_genes]

    ## preliminary uSORT sorting
    cat(paste0("-Preliminary sorting using ", preliminary_sorting_method, "\n"))
    preliminary_sorting_order <-
        uSORT_sorting_wrapper(data = exp_PCA_genes,
                              data_raw = exp_raw,
                              method = preliminary_sorting_method,
                              data_type = pre_data_type,
                              SPIN_option = pre_SPIN_option,
                              SPIN_sigma_width = pre_SPIN_sigma_width,
                              autoSPIN_alpha = pre_autoSPIN_alpha,
                              autoSPIN_randomization = pre_autoSPIN_randomization,
                              wanderlust_start_cell = pre_wanderlust_start_cell,
                              wanderlust_dfmap_components = pre_wanderlust_dfmap_components,
                              wanderlust_l = pre_wanderlust_l,
                              wanderlust_num_waypoints = pre_wanderlust_num_waypoints,
                              wanderlust_waypoints_seed = pre_wanderlust_waypoints_seed,
                              wanderlust_flock_waypoints = pre_wanderlust_flock_waypoints)

    if(is.null(preliminary_sorting_order)){
        exp_preliminary_order <- exp_PCA_genes
        cds <- EXP_to_CellDataSet(exp_trimmed, exp_raw)
    }else{
        exp_preliminary_order <- exp_PCA_genes[preliminary_sorting_order, ]
        cds <- EXP_to_CellDataSet(exp_trimmed[preliminary_sorting_order, ],
                                  exp_raw)
    }

    ## refined gene selection by driver genes selection
    cat("-Refined gene selection by driver genes detection...\n")
    driver_genes <-
        driving_force_gene_selection(cds = cds,
                                     scattering.cutoff.prob = scattering_cutoff_prob,
                                     driving.force.cutoff = driving_force_cutoff,
                                     data_type = ref_data_type,
                                     qval_cutoff = qval_cutoff_featureSelection,
                                     nCores = nCores)


    ## refined uSORT sorting
    exp_driver_genes <- exp_trimmed[preliminary_sorting_order, driver_genes]
    cat(paste0("-Refined sorting using ", refine_sorting_method, "\n"))
    refined_sorting_order <-
        uSORT_sorting_wrapper(data = exp_driver_genes,
                              data_raw = exp_trimmed,
                              method = refine_sorting_method,
                              data_type = ref_data_type,
                              SPIN_option = ref_SPIN_option,
                              SPIN_sigma_width = ref_SPIN_sigma_width,
                              autoSPIN_alpha = ref_autoSPIN_alpha,
                              autoSPIN_randomization = ref_autoSPIN_randomization,
                              wanderlust_start_cell = ref_wanderlust_start_cell,
                              wanderlust_dfmap_components = ref_wanderlust_dfmap_components,
                              wanderlust_l = ref_wanderlust_l,
                              wanderlust_num_waypoints = ref_wanderlust_num_waypoints,
                              wanderlust_waypoints_seed = ref_wanderlust_waypoints_seed,
                              wanderlust_flock_waypoints = ref_wanderlust_flock_waypoints)

    if(is.null(refined_sorting_order)){
        exp_refined_order <- exp_driver_genes
    }else{
        exp_refined_order <- exp_driver_genes[refined_sorting_order, ]
    }

    ## organize result object and write results
    uSORT_results <-
        list(exp_raw = exp_raw,
             trimmed_log2exp = exp_trimmed,
             preliminary_sorting_genes = PCA_selected_genes,
             preliminary_sorting_order = preliminary_sorting_order,
             refined_sorting_genes = driver_genes,
             refined_sorting_order = refined_sorting_order,
             driverGene_refinedOrder_log2exp = exp_refined_order)

    cat("-Sorting was done...\n")
    if(save_results){
        uSORT_write_results(uSORT_results = uSORT_results,
                            project_name = project_name,
                            result_directory = result_directory)
    }
    invisible(uSORT_results)
}




