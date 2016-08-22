itn_saveData <-
function (data, file = TRUE, save_title = "", init_filename = "", 
          back_to_r_console = FALSE) 
{
    model <- data
    #require("tcltk", quiet = TRUE)
    is_sc_object <- FALSE
    if (save_title == "") {
        save_title = "Save Object..."
    }
    if (is.list(data) && sum(match(names(data), "obj_type", nomatch = 0)) > 
        0) {
        is_sc_object <- TRUE
    }
    if (file == TRUE) file <- paste0(save_title,'.fso')
#     if (file == TRUE) {
#         if (is_sc_object) {
#             file <- tclvalue(tkgetSaveFile(title = save_title, 
#                                            filetypes = "{{FluidigmSC Object File} {.fso}}", 
#                                            initialfile = init_filename))
#             if (back_to_r_console == TRUE) 
#                 itn_BringToTop(-1)
#             if (file != "") {
#                 file_ext = tools::file_ext(file)
#                 if (file_ext == "") 
#                     file = paste(file, "fso", sep = ".")
#             }
#         }
#         else {
#             file <- tclvalue(tkgetSaveFile(title = save_title, 
#                                            filetypes = "{{Tab Delimited File} {.txt}}", 
#                                            initialfile = init_filename))
#             if (back_to_r_console == TRUE) 
#                 itn_BringToTop(-1)
#             if (file != "") {
#                 file_ext = tools::file_ext(file)
#                 if (file_ext == "") 
#                     file = paste(file, "txt", sep = ".")
#             }
#         }
#     }
    if (file == "") {
        stop("No file name is given.")
    }
    else {
        ops <- options()
        options(warn = -1)
        if (is_sc_object == FALSE) {
            rowname_flag <- TRUE
            colname_flag <- TRUE
            if (ncol(data) == 0 || is.na(as.numeric(rownames(data)[1])) == 
                FALSE) {
                rowname_flag <- FALSE
            }
            if (ncol(data) == 0 || (colnames(data))[1] == "V1") {
                colname_flag <- FALSE
            }
            if (rowname_flag == TRUE && colname_flag == TRUE) {
                new_col_names <- t(as.matrix(c("ID", colnames(data))))
                write.table(new_col_names, file = file, sep = "\t", 
                            quote = FALSE, row.names = FALSE, col.names = FALSE, 
                            append = FALSE)
                write.table(data, file = file, sep = "\t", quote = FALSE, 
                            row.names = TRUE, col.names = FALSE, append = TRUE)
            }
            else {
                write.table(data, file = file, sep = "\t", quote = FALSE, 
                            row.names = rowname_flag, col.names = colname_flag, 
                            append = FALSE)
            }
        }
        else {
            write.table("fluidigmSC Object File Version 1.0", 
                        file, sep = "\t", quote = FALSE, row.names = FALSE, 
                        col.names = FALSE, append = FALSE)
            for (i_name in 1:length(model)) {
                write.table("", file, sep = "\t", quote = FALSE, 
                            row.names = FALSE, col.names = FALSE, append = TRUE)
                begin_line <- paste("<start::", names(model)[i_name], 
                                    ">", sep = "")
                write.table(begin_line, file, sep = "\t", quote = FALSE, 
                            row.names = FALSE, col.names = FALSE, append = TRUE)
                i_model <- model[[i_name]]
                rowname_flag <- TRUE
                colname_flag <- TRUE
                rowname1 <- as.numeric(rownames(i_model)[1])
                if (ncol(i_model) == 0 || is.na(rowname1) == 
                    FALSE) {
                    rowname_flag <- FALSE
                }
                if (ncol(i_model) == 0 || (colnames(i_model))[1] == 
                    "V1") {
                    colname_flag <- FALSE
                }
                rowname_line <- paste("rowname::", rowname_flag, 
                                      sep = "")
                colname_line <- paste("colname::", colname_flag, 
                                      sep = "")
                rowcolnames_lines <- rbind(rowname_line, colname_line)
                write.table(rowcolnames_lines, file, sep = "\t", 
                            quote = FALSE, row.names = FALSE, col.names = FALSE, 
                            append = TRUE)
                if (rowname_flag == FALSE && colname_flag == 
                    FALSE) {
                    write.table(i_model, file, sep = "\t", quote = FALSE, 
                                row.names = FALSE, col.names = FALSE, append = TRUE)
                }
                else if (rowname_flag == TRUE && colname_flag == 
                         TRUE) {
                    new_col_names <- t(as.matrix(c("ID", colnames(i_model))))
                    write.table(new_col_names, file = file, sep = "\t", 
                                quote = FALSE, row.names = FALSE, col.names = FALSE, 
                                append = TRUE)
                    write.table(i_model, file = file, sep = "\t", 
                                quote = FALSE, row.names = TRUE, col.names = FALSE, 
                                append = TRUE)
                }
                else if (rowname_flag == FALSE && colname_flag == 
                         TRUE) {
                    write.table(i_model, file, sep = "\t", quote = FALSE, 
                                row.names = FALSE, col.names = TRUE, append = TRUE)
                }
                end_line <- paste("<end::", names(model)[i_name], 
                                  ">", sep = "")
                write.table(end_line, file, sep = "\t", quote = FALSE, 
                            row.names = FALSE, col.names = FALSE, append = TRUE)
            }
        }
        options(ops)
    }
    invisible(file)
}