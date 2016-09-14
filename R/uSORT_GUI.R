#' The user friendly GUI for \code{uSORT-package}
#'
#' This GUI provides an easy way for applying the uSORT package.
#'
#' @author Hao Chen
#' @return the GUI for \code{uSORT-package}
#' @import tcltk
#' @export
#' @seealso \code{\link{uSORT-package}}, \code{\link{uSORT}}
#' @references \url{http://JinmiaoChenLab.github.io/uSORT/}
#' @examples
#' interactive()
#' #if(interactive()) uSORT_GUI()  # remove the hash symbol to run
uSORT_GUI <- function(){
    paras <- uSORT_parameters_GUI()

    if(paras$return_status == "YES"){
        cat("-Parameters for uSORT: \n")
        for(i in seq_len(length(paras))){
            cat(paste0("  ", names(paras)[i], " = ", paras[[i]], "\n"))
        }

        ## run uSORT
        uSORT(exprs_file = paras$exprs_file,
              log_transform = TRUE,
              remove_outliers = TRUE,
              preliminary_sorting_method = paras$preliminary_sorting_method,
              refine_sorting_method = paras$refine_sorting_method,
              project_name = paras$project_name,
              result_directory = paras$result_directory,
              # gene selection parameters
              scattering_cutoff_prob = 0.75,
              driving_force_cutoff = NULL,
              qval_cutoff_featureSelection = 0.05,
              # autoSPIN parameters for preliminary sorting
              pre_data_type = paras$pre_data_type,
              pre_SPIN_option = paras$pre_SPIN_option,
              pre_SPIN_sigma_width = as.numeric(
                  paras$pre_SPIN_sigma_width),
              pre_autoSPIN_alpha = as.numeric(
                  paras$pre_autoSPIN_alpha),
              pre_autoSPIN_randomization = as.numeric(
                  paras$pre_autoSPIN_randomization),
              # autoSPIN parameters for refine sorting
              ref_data_type = paras$ref_data_type,
              ref_SPIN_option = paras$ref_SPIN_option,
              ref_SPIN_sigma_width = as.numeric(
                  paras$ref_SPIN_sigma_width),
              ref_autoSPIN_alpha = as.numeric(
                  paras$ref_autoSPIN_alpha),
              ref_autoSPIN_randomization = as.numeric(
                  paras$ref_autoSPIN_randomization),
              # Wanderlust parameters for preliminary sorting
              pre_wanderlust_start_cell = as.numeric(
                  paras$pre_wanderlust_start_cell),
              pre_wanderlust_dfmap_components = as.numeric(
                  paras$pre_wanderlust_dfmap_components),
              pre_wanderlust_l = as.numeric(
                  paras$pre_wanderlust_l),
              pre_wanderlust_num_waypoints = as.numeric(
                  paras$pre_wanderlust_num_waypoints),
              pre_wanderlust_waypoints_seed = as.numeric(
                  paras$pre_wanderlust_waypoints_seed),
              pre_wanderlust_flock_waypoints = as.numeric(
                  paras$pre_wanderlust_flock_waypoints),
              # Wanderlust parameters for refine sorting
              ref_wanderlust_start_cell = as.numeric(
                  paras$ref_wanderlust_start_cell),
              ref_wanderlust_dfmap_components = as.numeric(
                  paras$ref_wanderlust_dfmap_components),
              ref_wanderlust_l = as.numeric(paras$ref_wanderlust_l),
              ref_wanderlust_num_waypoints = as.numeric(
                  paras$ref_wanderlust_num_waypoints),
              ref_wanderlust_waypoints_seed = as.numeric(
                  paras$ref_wanderlust_waypoints_seed),
              ref_wanderlust_flock_waypoints = as.numeric(
                  paras$ref_wanderlust_flock_waypoints)
        )
    }else{
        message("uSORT analysis cancelled!")
    }

}

#' The GUI for inputting paramters for uSORT
#'
#' This is a function for generating the GUI for uSORT,
#' it's called by \code{\link{uSORT_GUI}}. For internal use only.
#'
#' @return a list of parameters.
#'
#' @author Hao Chen
uSORT_parameters_GUI <- function(){
    ## global options
    sort_options <- c("autoSPIN", "sWanderlust", "monocle",
                      "Wanderlust", "SPIN", "none")
    data_type_options <- c("linear", "cyclical")
    SPIN_option <- c("STS", "neighborhood")

    ## tcl varaibles
    input_data_file <- tclVar("")
    project_name <- tclVar("uSORT")
    preliminary_sort_method <- tclVar(sort_options[1])
    refine_sort_method <- tclVar(sort_options[1])
    result_directory <- tclVar(getwd())

    pre_data_type <- tclVar(data_type_options[1])
    pre_SPIN_option <- tclVar(SPIN_option[1])
    pre_SPIN_sigma_width <- tclVar(1)
    pre_autoSPIN_alpha <- tclVar(0.2)
    pre_autoSPIN_randomization <- tclVar(20)
    # autoSPIN parameters for refine sorting
    ref_data_type <- tclVar(data_type_options[1])
    ref_SPIN_option <- tclVar(SPIN_option[1])
    ref_SPIN_sigma_width <- tclVar(1)
    ref_autoSPIN_alpha <- tclVar(0.2)
    ref_autoSPIN_randomization <- tclVar(20)
    # Wanderlust parameters for preliminary sorting
    pre_wanderlust_start_cell <- tclVar(1)
    pre_wanderlust_dfmap_components <- tclVar(4)
    pre_wanderlust_l <- tclVar(15)
    pre_wanderlust_num_waypoints <- tclVar(150)
    pre_wanderlust_waypoints_seed <- tclVar(2711)
    pre_wanderlust_flock_waypoints <- tclVar(2)
    # Wanderlust parameters for refine sorting
    ref_wanderlust_start_cell <- tclVar(1)
    ref_wanderlust_dfmap_components <- tclVar(4)
    ref_wanderlust_l <- tclVar(15)
    ref_wanderlust_num_waypoints <- tclVar(150)
    ref_wanderlust_waypoints_seed <- tclVar(2711)
    ref_wanderlust_flock_waypoints <- tclVar(2)
    # return status label
    return_status <- tclVar("YES")

    ## widget control functions (commonds)
    reset_input_data_file <- function(){
        InputFile <- NULL
        InputFile <- tk_choose.files(
            default = paste(tclvalue(result_directory), "fso",
                            sep = .Platform$file.sep),
            caption = "Select your data file", multi = FALSE,
            filters = matrix(c("{txt files}", "{.txt}"),1, 2), index = 1)
        if (!is.null(InputFile)){
            tclvalue(input_data_file) <- InputFile
            tclvalue(result_directory) <- dirname(InputFile)
        }

    }

    reset_preliminary_sort_parameter <- function(){
        method <- tclvalue(preliminary_sort_method)
        paras <- sorting_method_parameter_GUI(method)

        tclvalue(pre_data_type) <- paras$data_type
        tclvalue(pre_SPIN_option) <- paras$SPIN_option
        tclvalue(pre_SPIN_sigma_width) <- paras$SPIN_sigma_width
        tclvalue(pre_autoSPIN_alpha) <- paras$autoSPIN_alpha
        tclvalue(pre_autoSPIN_randomization) <- paras$autoSPIN_randomization
        tclvalue(pre_wanderlust_start_cell) <- paras$wanderlust_start_cell
        tclvalue(pre_wanderlust_dfmap_components) <- paras$wanderlust_dfmap_components
        tclvalue(pre_wanderlust_l) <- paras$wanderlust_l
        tclvalue(pre_wanderlust_num_waypoints) <- paras$wanderlust_num_waypoints
        tclvalue(pre_wanderlust_waypoints_seed) <- paras$wanderlust_waypoints_seed
        tclvalue(pre_wanderlust_flock_waypoints) <- paras$wanderlust_flock_waypoints

    }

    reset_refine_sort_parameter <- function(){
        method <- tclvalue(refine_sort_method)
        paras <- sorting_method_parameter_GUI(method)

        tclvalue(ref_data_type) <- paras$data_type
        tclvalue(ref_SPIN_option) <- paras$SPIN_option
        tclvalue(ref_SPIN_sigma_width) <- paras$SPIN_sigma_width
        tclvalue(ref_autoSPIN_alpha) <- paras$autoSPIN_alpha
        tclvalue(ref_autoSPIN_randomization) <- paras$autoSPIN_randomization
        tclvalue(ref_wanderlust_start_cell) <- paras$wanderlust_start_cell
        tclvalue(ref_wanderlust_dfmap_components) <- paras$wanderlust_dfmap_components
        tclvalue(ref_wanderlust_l) <- paras$wanderlust_l
        tclvalue(ref_wanderlust_num_waypoints) <- paras$wanderlust_num_waypoints
        tclvalue(ref_wanderlust_waypoints_seed) <- paras$wanderlust_waypoints_seed
        tclvalue(ref_wanderlust_flock_waypoints) <- paras$wanderlust_flock_waypoints

    }

    reset_result_directory <- function(){
        ResultDir <- NULL
        ResultDir <- tclvalue(tkchooseDirectory(
            title = "Choose your data dircetory ..."))
        if (!is.null(ResultDir))
            tclvalue(result_directory) <- ResultDir
    }

    quit <- function() {
        tclvalue(return_status) <- "NO"
        tkdestroy(tt)
        }

    reset <- function() {
        tclvalue(input_data_file) <- ""
        tclvalue(preliminary_sort_method) <- sort_options[1]
        tclvalue(refine_sort_method) <- sort_options[1]
        tclvalue(result_directory) <- getwd()
    }

    submit <- function() {
        ErrorFlag <- FALSE
        if (!file.exists(tclvalue(input_data_file))) {
            ErrorFlag <- TRUE
            tkmessageBox(title = "uSORT Package",
                         message = "Data file doesn't exit,
                         please provide correct data file.",
                         icon = "info", type = "ok")
        }

        if (!ErrorFlag) tkdestroy(tt)
    }

    ## GUI parameters and widgets
    box_length <- 40
    cell_width <- 4
    bt_width <- 8

    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, "uSORT Analysis")

    input_data_file_label <- tklabel(tt, text = "Input Data File :")
    input_data_file_entry <- tkentry(tt, textvariable = input_data_file,
                                     width = box_length)
    input_data_file_button <- tkbutton(tt, text = " Select... ",
                                       width = bt_width,
                                       command = reset_input_data_file)

    project_name_label <- tklabel(tt, text = "Project Name :")
    project_name_entry <- tkentry(tt, textvariable = project_name,
                                  width = box_length)

    preliminary_sort_method_label <- tklabel(tt,
                                             text = "Preliminary Sorting :")
    preliminary_sort_method_combo <- ttkcombobox(tt, values = sort_options,
                                                 textvariable = preliminary_sort_method,
                                                 state="readonly",
                                                 width = box_length-2)
    preliminary_sort_method_button <- tkbutton(tt, text = " Parameter ",
                                       width = bt_width,
                                       command = reset_preliminary_sort_parameter)

    refine_sort_method_label <- tklabel(tt, text = "Refine Sorting :")
    refine_sort_method_combo <- ttkcombobox(tt, values = sort_options,
                                            textvariable = refine_sort_method,
                                            state="readonly",
                                            width = box_length-2)
    refine_sort_method_button <- tkbutton(tt, text = " Parameter ",
                                          width = bt_width,
                                          command = reset_refine_sort_parameter)

    result_directory_label <- tklabel(tt, text = "Result Directory :")
    result_directory_entry <- tkentry(tt, textvariable = result_directory,
                                      width = box_length)
    result_directory_button <- tkbutton(tt, text = " Choose... ",
                                        width = bt_width,
                                        command = reset_result_directory)

    quit_button <- tkbutton(tt, text = "Quit", command = quit)
    reset_button <- tkbutton(tt, text = "Reset", command = reset)
    submit_button <- tkbutton(tt, text = "Submit", command = submit)

    ## construct GUI
    tkgrid(input_data_file_label, input_data_file_entry,
           input_data_file_button, padx = cell_width)
    tkgrid.configure(input_data_file_label, input_data_file_entry,
                     input_data_file_button, sticky = "e")

    tkgrid(preliminary_sort_method_label, preliminary_sort_method_combo,
           preliminary_sort_method_button, padx = cell_width)
    tkgrid.configure(preliminary_sort_method_label,
                     preliminary_sort_method_combo,
                     preliminary_sort_method_button, sticky = "e")

    tkgrid(refine_sort_method_label, refine_sort_method_combo,
           refine_sort_method_button, padx = cell_width)
    tkgrid.configure(refine_sort_method_label, refine_sort_method_combo,
                     refine_sort_method_button, sticky = "e")

    tkgrid(result_directory_label, result_directory_entry,
           result_directory_button, padx = cell_width)
    tkgrid.configure(result_directory_label, result_directory_entry,
                     result_directory_button, sticky = "e")

    tkgrid(project_name_label, project_name_entry, padx = cell_width)
    tkgrid.configure(project_name_label, project_name_entry, sticky = "e")

    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)

    tkgrid(quit_button, reset_button, submit_button, padx = cell_width)
    tkgrid.configure(quit_button, sticky = "e")
    tkgrid.configure(submit_button, sticky = "w")

    tkwait.window(tt)

    return(list(exprs_file = tclvalue(input_data_file),
                preliminary_sorting_method = tclvalue(preliminary_sort_method),
                refine_sorting_method = tclvalue(refine_sort_method),
                result_directory = tclvalue(result_directory),
                project_name = tclvalue(project_name),

                pre_data_type = tclvalue(pre_data_type),
                pre_SPIN_option = tclvalue(pre_SPIN_option),
                pre_SPIN_sigma_width = tclvalue(pre_SPIN_sigma_width),
                pre_autoSPIN_alpha = tclvalue(pre_autoSPIN_alpha),
                pre_autoSPIN_randomization = tclvalue(pre_autoSPIN_randomization),
                pre_wanderlust_start_cell = tclvalue(pre_wanderlust_start_cell),
                pre_wanderlust_dfmap_components = tclvalue(pre_wanderlust_dfmap_components),
                pre_wanderlust_l = tclvalue(pre_wanderlust_l),
                pre_wanderlust_num_waypoints = tclvalue(pre_wanderlust_num_waypoints),
                pre_wanderlust_waypoints_seed = tclvalue(pre_wanderlust_waypoints_seed),
                pre_wanderlust_flock_waypoints = tclvalue(pre_wanderlust_flock_waypoints),

                ref_data_type = tclvalue(ref_data_type),
                ref_SPIN_option = tclvalue(ref_SPIN_option),
                ref_SPIN_sigma_width = tclvalue(ref_SPIN_sigma_width),
                ref_autoSPIN_alpha = tclvalue(ref_autoSPIN_alpha),
                ref_autoSPIN_randomization = tclvalue(ref_autoSPIN_randomization),
                ref_wanderlust_start_cell = tclvalue(ref_wanderlust_start_cell),
                ref_wanderlust_dfmap_components = tclvalue(ref_wanderlust_dfmap_components),
                ref_wanderlust_l = tclvalue(ref_wanderlust_l),
                ref_wanderlust_num_waypoints = tclvalue(ref_wanderlust_num_waypoints),
                ref_wanderlust_waypoints_seed = tclvalue(ref_wanderlust_waypoints_seed),
                ref_wanderlust_flock_waypoints = tclvalue(ref_wanderlust_flock_waypoints),

                return_status = tclvalue(return_status)
                )
           )
}


#' GUI for sorting method paramters
#'
#' The parameters appeared on GUI are based on input method, this function is called by
#' \code{\link{uSORT_parameters_GUI}}. For internal use only.
#'
#' @param method method name.
#' @return a list of parameters.
#'
#' @author Hao Chen
sorting_method_parameter_GUI <- function(
    method = c("autoSPIN", "sWanderlust",
               "monocle", "Wanderlust",
               "SPIN", "none")){
    method <- match.arg(method)
    data_type_options <- c("linear", "cyclical")
    SPIN_options <- c("STS", "neighborhood")

    data_type <- tclVar(data_type_options[1])
    SPIN_option <- tclVar(SPIN_options[1])
    SPIN_sigma_width <- tclVar(1)
    autoSPIN_alpha <- tclVar(0.2)
    autoSPIN_randomization <- tclVar(20)
    wanderlust_start_cell <- tclVar(1)
    wanderlust_dfmap_components <- tclVar(4)
    wanderlust_l <- tclVar(15)
    wanderlust_num_waypoints <- tclVar(150)
    wanderlust_waypoints_seed <- tclVar(2711)
    wanderlust_flock_waypoints <- tclVar(2)

    ## GUI parameters and widgets
    box_length <- 40
    cell_width <- 4
    bt_width <- 8

    switch(method,
           autoSPIN = {
               reset <- function() {
                   tclvalue(data_type) <- data_type_options[1]
                   tclvalue(SPIN_option) <- SPIN_options[1]
                   tclvalue(SPIN_sigma_width) <- 1
                   tclvalue(autoSPIN_alpha) <- 0.2
                   tclvalue(autoSPIN_randomization) <- 10
               }

               ok <- function() {
                   tkdestroy(tt)
               }

               tt <- tktoplevel(borderwidth = 20)
               tkwm.title(tt, "uSORT: parameters for autoSPIN")

               data_type_label <- tklabel(tt, text = "data type :")
               data_type_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[1],
                                    variable = data_type,
                                    value = data_type_options[1]),
                      side = "left")
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[2],
                                    variable = data_type,
                                    value = data_type_options[2]),
                      side = "left")

               method_label <- tklabel(tt, text = "sorting method :")
               method_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(method_rbuts, text = SPIN_options[1],
                                    variable = SPIN_option,
                                    value = SPIN_options[1]),
                      side = "left")
               tkpack(tkradiobutton(method_rbuts, text = SPIN_options[2],
                                    variable = SPIN_option,
                                    value = SPIN_options[2]),
                      side = "left")

               alpha_label <- tklabel(tt, text = "alpha :")
               alpha_entry <- tkentry(tt, textvariable = autoSPIN_alpha,
                                      width = box_length)

               sigma_width_label <- tklabel(tt, text = "sigma width :")
               sigma_width_entry <- tkentry(tt,
                                            textvariable = SPIN_sigma_width,
                                            width = box_length)

               no_randomization_label <- tklabel(tt,
                                                 text = "no. randomization :")
               no_randomization_entry <- tkentry(tt,
                                                 textvariable=autoSPIN_randomization,
                                                 width = box_length)

               reset_button <- tkbutton(tt, text = "Reset", command = reset)
               ok_button <- tkbutton(tt, text = "OK", command = ok)

               tkgrid(data_type_label, data_type_rbuts, padx = cell_width)
               tkgrid.configure(data_type_label, sticky = "e")
               tkgrid.configure(data_type_rbuts, sticky = "w")

               tkgrid(method_label, method_rbuts, padx = cell_width)
               tkgrid.configure(method_label, sticky = "e")
               tkgrid.configure(method_rbuts, sticky = "w")

               tkgrid(alpha_label, alpha_entry, padx = cell_width)
               tkgrid.configure(alpha_label, alpha_entry, sticky = "e")

               tkgrid(sigma_width_label, sigma_width_entry,
                      padx = cell_width)
               tkgrid.configure(sigma_width_label, sigma_width_entry,
                                sticky = "e")

               tkgrid(no_randomization_label, no_randomization_entry,
                      padx = cell_width)
               tkgrid.configure(no_randomization_label,
                                no_randomization_entry, sticky = "e")

               tkgrid(reset_button, ok_button, padx = cell_width)
               tkgrid.configure(reset_button, sticky = "e")

               tkwait.window(tt)
           },
           sWanderlust = {
               reset <- function() {
                   tclvalue(data_type) <- data_type_options[1]
                   tclvalue(SPIN_option) <- SPIN_options[1]
                   tclvalue(SPIN_sigma_width) <- 1
                   tclvalue(autoSPIN_alpha) <- 0.2
                   tclvalue(wanderlust_dfmap_components) <- 4
                   tclvalue(wanderlust_l) <- 15
                   tclvalue(wanderlust_num_waypoints) <- 150
                   tclvalue(wanderlust_waypoints_seed) <- 2711
                   tclvalue(wanderlust_flock_waypoints) <- 2
               }

               ok <- function() {
                   tkdestroy(tt)
               }

               tt <- tktoplevel(borderwidth = 20)
               tkwm.title(tt, "uSORT: parameters for sWanderlust")

               data_type_label <- tklabel(tt, text = "data type :")
               data_type_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[1],
                                    variable = data_type,
                                    value = data_type_options[1]),
                      side = "left")
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[2],
                                    variable = data_type,
                                    value = data_type_options[2]),
                      side = "left")

               method_label <- tklabel(tt, text = "sorting method :")
               method_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(method_rbuts,
                                    text = SPIN_options[1],
                                    variable = SPIN_option,
                                    value = SPIN_options[1]),
                      side = "left")
               tkpack(tkradiobutton(method_rbuts,
                                    text = SPIN_options[2],
                                    variable = SPIN_option,
                                    value = SPIN_options[2]),
                      side = "left")

               alpha_label <- tklabel(tt, text = "alpha :")
               alpha_entry <- tkentry(tt, textvariable = autoSPIN_alpha,
                                      width = box_length)

               sigma_width_label <- tklabel(tt, text = "sigma width :")
               sigma_width_entry <- tkentry(tt,
                                            textvariable = SPIN_sigma_width,
                                            width = box_length)

               dfmap_components_label <- tklabel(tt,
                                                 text = "diffusionmap components :")
               dfmap_components_entry <- tkentry(tt,
                                                 textvariable = wanderlust_dfmap_components,
                                                 width = box_length)

               l_label <- tklabel(tt, text = "no. nearest neighbors :")
               l_entry <- tkentry(tt, textvariable = wanderlust_l,
                                  width = box_length)

               num_waypoints_label <- tklabel(tt, text = "no. waypoints :")
               num_waypoints_entry <- tkentry(tt,
                                              textvariable = wanderlust_num_waypoints,
                                              width = box_length)

               waypoints_seed_label <- tklabel(tt, text = "waypoints seed :")
               waypoints_seed_entry <- tkentry(tt,
                                               textvariable = wanderlust_waypoints_seed,
                                               width = box_length)

               flock_waypoints_label <- tklabel(tt,
                                                text = "flock waypoints :")
               flock_waypoints_entry <- tkentry(tt,
                                                textvariable = wanderlust_flock_waypoints,
                                                width = box_length)

               reset_button <- tkbutton(tt, text = "Reset",
                                        command = reset)
               ok_button <- tkbutton(tt, text = "OK",
                                     command = ok)

               tkgrid(data_type_label, data_type_rbuts, padx = cell_width)
               tkgrid.configure(data_type_label, sticky = "e")
               tkgrid.configure(data_type_rbuts, sticky = "w")

               tkgrid(method_label, method_rbuts, padx = cell_width)
               tkgrid.configure(method_label, sticky = "e")
               tkgrid.configure(method_rbuts, sticky = "w")

               tkgrid(alpha_label, alpha_entry, padx = cell_width)
               tkgrid.configure(alpha_label, alpha_entry, sticky = "e")

               tkgrid(sigma_width_label, sigma_width_entry,
                      padx = cell_width)
               tkgrid.configure(sigma_width_label, sigma_width_entry,
                                sticky = "e")

               tkgrid(dfmap_components_label, dfmap_components_entry,
                      padx = cell_width)
               tkgrid.configure(dfmap_components_label,
                                dfmap_components_entry, sticky = "e")

               tkgrid(l_label, l_entry, padx = cell_width)
               tkgrid.configure(l_label, l_entry, sticky = "e")

               tkgrid(num_waypoints_label, num_waypoints_entry,
                      padx = cell_width)
               tkgrid.configure(num_waypoints_label, num_waypoints_entry,
                                sticky = "e")

               tkgrid(waypoints_seed_label, waypoints_seed_entry,
                      padx = cell_width)
               tkgrid.configure(waypoints_seed_label, waypoints_seed_entry,
                                sticky = "e")

               tkgrid(flock_waypoints_label, flock_waypoints_entry,
                      padx = cell_width)
               tkgrid.configure(flock_waypoints_label,
                                flock_waypoints_entry, sticky = "e")

               tkgrid(reset_button, ok_button, padx = cell_width)
               tkgrid.configure(reset_button, sticky = "e")

               tkwait.window(tt)

           },
           monocle = {
               return(NULL)
           },
           Wanderlust = {
               reset <- function() {
                   tclvalue(data_type) <- data_type_options[1]
                   tclvalue(wanderlust_start_cell) <- ""
                   tclvalue(wanderlust_dfmap_components) <- 4
                   tclvalue(wanderlust_l) <- 15
                   tclvalue(wanderlust_num_waypoints) <- 150
                   tclvalue(wanderlust_waypoints_seed) <- 2711
                   tclvalue(wanderlust_flock_waypoints) <- 2
               }

               ok <- function() {
                   tkdestroy(tt)
               }

               tt <- tktoplevel(borderwidth = 20)
               tkwm.title(tt, "uSORT: parameters for sWanderlust")

               data_type_label <- tklabel(tt, text = "data type :")
               data_type_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[1],
                                    variable = data_type,
                                    value = data_type_options[1]),
                      side = "left")
               tkpack(tkradiobutton(data_type_rbuts,
                                    text = data_type_options[2],
                                    variable = data_type,
                                    value = data_type_options[2]),
                      side = "left")

               start_cell_label <- tklabel(tt, text = "start cell :")
               start_cell_entry <- tkentry(tt,
                                           textvariable = wanderlust_start_cell,
                                           width = box_length)

               dfmap_components_label <- tklabel(tt,
                                                 text = "diffusionmap components :")
               dfmap_components_entry <- tkentry(tt,
                                                 textvariable = wanderlust_dfmap_components,
                                                 width = box_length)

               l_label <- tklabel(tt, text = "no. nearest neighbors :")
               l_entry <- tkentry(tt, textvariable = wanderlust_l,
                                  width = box_length)

               num_waypoints_label <- tklabel(tt, text = "no. waypoints :")
               num_waypoints_entry <- tkentry(tt,
                                              textvariable = wanderlust_num_waypoints,
                                              width = box_length)

               waypoints_seed_label <- tklabel(tt,
                                               text = "waypoints seed :")
               waypoints_seed_entry <- tkentry(tt,
                                               textvariable = wanderlust_waypoints_seed,
                                               width = box_length)

               flock_waypoints_label <- tklabel(tt,
                                                text = "flock waypoints :")
               flock_waypoints_entry <- tkentry(tt,
                                                textvariable = wanderlust_flock_waypoints,
                                                width = box_length)

               reset_button <- tkbutton(tt, text = "Reset",
                                        command = reset)
               ok_button <- tkbutton(tt, text = "OK",
                                     command = ok)

               tkgrid(data_type_label, data_type_rbuts,
                      padx = cell_width)
               tkgrid.configure(data_type_label, sticky = "e")
               tkgrid.configure(data_type_rbuts, sticky = "w")

               tkgrid(start_cell_label, start_cell_entry, padx = cell_width)
               tkgrid.configure(start_cell_label, start_cell_entry,
                                sticky = "e")

               tkgrid(dfmap_components_label, dfmap_components_entry,
                      padx = cell_width)
               tkgrid.configure(dfmap_components_label,
                                dfmap_components_entry, sticky = "e")

               tkgrid(l_label, l_entry, padx = cell_width)
               tkgrid.configure(l_label, l_entry, sticky = "e")

               tkgrid(num_waypoints_label, num_waypoints_entry,
                      padx = cell_width)
               tkgrid.configure(num_waypoints_label,
                                num_waypoints_entry, sticky = "e")

               tkgrid(waypoints_seed_label, waypoints_seed_entry,
                      padx = cell_width)
               tkgrid.configure(waypoints_seed_label, waypoints_seed_entry,
                                sticky = "e")

               tkgrid(flock_waypoints_label, flock_waypoints_entry,
                      padx = cell_width)
               tkgrid.configure(flock_waypoints_label,
                                flock_waypoints_entry, sticky = "e")

               tkgrid(reset_button, ok_button, padx = cell_width)
               tkgrid.configure(reset_button, sticky = "e")

               tkwait.window(tt)

           },
           SPIN = {
               reset <- function() {
                   tclvalue(SPIN_option) <- SPIN_options[1]
                   tclvalue(SPIN_sigma_width) <- 1
               }

               ok <- function() {
                   tkdestroy(tt)
               }

               tt <- tktoplevel(borderwidth = 20)
               tkwm.title(tt, "uSORT: parameters for SPIN")

               method_label <- tklabel(tt, text = "sorting method :")
               method_rbuts <- tkframe(tt)
               tkpack(tkradiobutton(method_rbuts, text = SPIN_options[1],
                                    variable = SPIN_option,
                                    value = SPIN_options[1]),
                      side = "left")
               tkpack(tkradiobutton(method_rbuts, text = SPIN_options[2],
                                    variable = SPIN_option,
                                    value = SPIN_options[2]),
                      side = "left")

               sigma_width_label <- tklabel(tt, text = "sigma width :")
               sigma_width_entry <- tkentry(tt,
                                            textvariable = SPIN_sigma_width,
                                            width = box_length)

               reset_button <- tkbutton(tt, text = "Reset", command = reset)
               ok_button <- tkbutton(tt, text = "OK", command = ok)

               tkgrid(method_label, method_rbuts, padx = cell_width)
               tkgrid.configure(method_label, sticky = "e")
               tkgrid.configure(method_rbuts, sticky = "w")

               tkgrid(sigma_width_label, sigma_width_entry,
                      padx = cell_width)
               tkgrid.configure(sigma_width_label, sigma_width_entry,
                                sticky = "e")

               tkgrid(reset_button, ok_button, padx = cell_width)
               tkgrid.configure(reset_button, sticky = "e")

               tkwait.window(tt)

           },
           none = {
               return(NULL)
           })

    paras <- list(data_type = tclvalue(data_type),
                  SPIN_option = tclvalue(SPIN_option),
                  SPIN_sigma_width = tclvalue(SPIN_sigma_width),
                  autoSPIN_alpha = tclvalue(autoSPIN_alpha),
                  autoSPIN_randomization = tclvalue(autoSPIN_randomization),
                  wanderlust_start_cell = tclvalue(wanderlust_start_cell),
                  wanderlust_dfmap_components = tclvalue(wanderlust_dfmap_components),
                  wanderlust_l = tclvalue(wanderlust_l),
                  wanderlust_num_waypoints = tclvalue(wanderlust_num_waypoints),
                  wanderlust_waypoints_seed = tclvalue(wanderlust_waypoints_seed),
                  wanderlust_flock_waypoints = tclvalue(wanderlust_flock_waypoints))

    return(paras)

}





