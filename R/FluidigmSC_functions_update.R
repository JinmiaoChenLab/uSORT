#' Modify some fluidigmSC's functions
fluidigmSC_functions_update <- function() {
    assignInNamespace(
        "itn_isSubSetIgnoreCase",
        itn_isSubSetIgnoreCase,
        ns = "fluidigmSC",
        envir = as.environment("package:fluidigmSC")
    )
    assignInNamespace("itn_getMatchingOriginalName",
                      itn_getMatchingOriginalName,
                      ns = "fluidigmSC")
    assignInNamespace("itn_anova_by_factor", itn_anova_by_factor, ns = "fluidigmSC")
    assignInNamespace("identifyOutliers", identifyOutliers, ns = "fluidigmSC")
    assignInNamespace("itn_displayOutlierGUI", itn_displayOutlierGUI, ns =
                          "fluidigmSC")
    assignInNamespace("itn_ScatterPlot", itn_ScatterPlot, ns = "fluidigmSC")
    assignInNamespace("itn_saveData", itn_saveData, ns = "fluidigmSC")
}



#' copied from fluidigmSC
updateGeneListFromList <- function (exp, gene_list)
{
    if (missing(exp))
        stop("Provide exp data")
    else if (missing(gene_list))
        stop("No Gene List File or Data provided.")
    if (itn_isCorrectInputs(exp = exp, gene_list = gene_list) ==
        FALSE)
        stop(fldm_last_error_msg)
    gene_list[, 1] <- itn_getMatchingOriginalName(as.vector(gene_list[, 1]),
                                                  rownames(exp$org_data))
    if (class(gene_list) != "data.frame")
        gene_list <- as.data.frame(gene_list)
    gene_list <- itn_removeEmptyDataLine(gene_list)
    if (ncol(gene_list) == 1) {
        GroupID <- rep("Untitled", nrow(gene_list))
        gene_list <- cbind(gene_list, GroupID)
    }
    else if (ncol(gene_list) >= 2 && colnames(gene_list)[2] ==
             "GroupID") {
        gene_list <- gene_list[, 1:2]
        colnames(gene_list) <- c("GeneID", "GroupID")
    }
    else if (ncol(gene_list) >= 2 && colnames(gene_list)[2] !=
             "GroupID") {
        GroupID <- rep("Untitled", nrow(gene_list))
        gene_list <- as.data.frame(cbind(as.character(gene_list[,
                                                                1]), GroupID))
        colnames(gene_list) <- c("GeneID", "GroupID")
    }
    gene_ids <- rownames(exp$org_data)
    if (itn_isSubSetIgnoreCase(gene_list[, 1], gene_ids)) {
        exp$gene_list <- as.data.frame(unique(gene_list))
    }
    else if (!itn_isSubSetIgnoreCase(gene_list[, 1], gene_ids) &&
             itn_isSubSetIgnoreCase(gene_list[-1, 1], gene_ids)) {
        exp$gene_list <- as.data.frame(unique(gene_list[-1, ]))
    }
    else {
        diff_genes <- itn_setDiffIgnoreCase(gene_list[-1, 1],
                                            gene_ids)
        stop(paste("The following genes cannot be found in EXP data: ",
                   diff_genes, sep = ""))
    }
    updated_exp <- itn_generateLog2EXDataFromOrgData(exp)
    return(updated_exp)
}




#' Modified from fluidigmSC
itn_anova_by_factor <-
    function (exp,
              aov_factor = "",
              group_difference = TRUE,
              p_FDR = TRUE)
    {
        expression_data <- exp$log2ex_data
        pValues <-
            apply(
                expression_data,
                1,
                itn_anova_x,
                aovFactor = aov_factor,
                Tukey_HSD = TRUE,
                colName = FALSE
            )
        if (sum(is.na(pValues)) > 0)
            pValues[is.na(pValues)] <- 1
        pValues <- t(pValues)
        if (group_difference) {
            diff <-
                apply(
                    expression_data,
                    1,
                    itn_anova_x,
                    aovFactor = aov_factor,
                    Tukey_HSD = TRUE,
                    colName = FALSE,
                    Difference = TRUE
                )
            diff <- t(diff)
        }
        x1 <- as.numeric(expression_data[1,])
        colNames <-
            itn_anova_x(x1,
                        aovFactor = aov_factor,
                        Tukey_HSD = TRUE,
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


itn_displayOutlierGUI <-
    function (exp_trimmed_all,
              exp_file,
              exp_outliers_all,
              exp_threshold_all,
              current_group_no = 1)
    {
        pca_plot <- TRUE
        exp <- exp_trimmed_all[[current_group_no]]
        outliers <- exp_outliers_all[[current_group_no]]
        threshold <- exp_threshold_all[[current_group_no]]
        sample_groups <- names(exp_outliers_all)
        title <-
            paste("Outlier Identification : ", basename(exp_file),
                  sep = "")
        #     itn_newPlotWindow(width = 11, height = 4.5, title = title,
        #                       xpos = 0)
        pdf(
            'outlier.boxplot.pdf',
            width = 11,
            height = 4.5,
            title = title
        )
        display_window <- dev.cur()
        expression_data_all <- exp$log2ex_data
        expression_data_all <-
            apply(expression_data_all, 2, as.numeric)
        data_median <- apply(expression_data_all, 2, median)
        data_median_sorted <- sort(data_median, decreasing = T)
        data_median_rownames <-
            rownames(as.matrix(data_median_sorted))
        expression_data_all <-
            expression_data_all[, data_median_rownames]
        data_median <- data_median_sorted
        d <- data.frame(data_median)
        d$QC = (d >= threshold)
        colnames(d) <- c('median', 'passQC')
        write.table(
            d,
            row.names = T,
            sep = '\t',
            col.names = NA,
            file = paste0('outliers.data.median.', basename(exp_file))
        )
        sample_ids <- colnames(expression_data_all)
        outliers <- as.matrix(outliers)[, 1]
        outlier_index <-
            c(1:ncol(expression_data_all))[is.element(sample_ids,
                                                      outliers)]
        outlier_num <- length(outlier_index)
        outlier_color <- "red"
        non_outlier_color <- "blue"
        sample_outlier_group_all <-
            rep(non_outlier_color, ncol(expression_data_all))
        sample_outlier_group_all[outlier_index] <- outlier_color
        boxplot_median_color <- "white"
        sample_legend_top <- 17
        sample_legend_grid <- 0.05
        sample_sep_scale <- 0.8
        xyN_all <- c()
        xyN <- c()
        ret <- as.data.frame(c())
        selected_outliers <- ""
        myplot <- dev.cur()
        draw_plot <- TRUE
        done_id <- FALSE
        hide_save_exit <- FALSE
        init_flag <- TRUE
        ret_cancel <- FALSE
        while (done_id == FALSE) {
            expression_data <- expression_data_all
            sample_outlier_group <- sample_outlier_group_all
            if (length(xyN_all) > 0) {
                sample_outlier_group[xyN_all] <- non_outlier_color
                sample_outlier_group[setdiff(xyN_all, outlier_index)] <-
                    outlier_color
            }
            if (draw_plot == TRUE || hide_save_exit == TRUE) {
                layout.matrix <- matrix(c(2, 1), ncol = 2, byrow = TRUE)
                layout.widths <- c(4, 0.4)
                layout.heights <- 1
                layout(
                    layout.matrix,
                    widths = layout.widths,
                    heights = layout.heights,
                    respect = FALSE
                )
                par(mar = c(6.5, 0, 3, 1))
                frame()
                value_range <- par("usr")
                min_x <- value_range[1]
                max_x <- value_range[2]
                min_y <- value_range[3]
                max_y <- value_range[4]
                sample_top <- seq(min_y, max_y, sample_legend_grid *
                                      (max_y - min_y))[sample_legend_top]
                sample_sep <-
                    sample_sep_scale * sample_legend_grid *
                    (max_y - min_y)
                sample_pos <-
                    seq(sample_top - (length(sample_groups) -
                                          1) * sample_sep,
                        sample_top,
                        sample_sep)
                sample_pos <- rev(sample_pos)
                sample_col <- "black"
                if (length(sample_groups) > 1) {
                    for (i in 1:length(sample_groups)) {
                        sample_pch <- 1
                        sample_cex <- 0.7
                        sample_font <- 1
                        text_pos <- 0.1
                        if (current_group_no == i) {
                            sample_pch <- 16
                            sample_cex <- 0.75
                            sample_font <- 2
                        }
                        points(text_pos,
                               sample_pos[i],
                               pch = sample_pch,
                               col = sample_col)
                        text(
                            text_pos,
                            sample_pos[i],
                            sample_groups[i],
                            cex = sample_cex,
                            font = sample_font,
                            pos = 4
                        )
                    }
                    points(
                        text_pos,
                        sample_pos[1] + 4.5 * (sample_pos[1] -
                                                   sample_pos[2]),
                        pch = 15,
                        col = non_outlier_color
                    )
                    text(
                        text_pos,
                        sample_pos[1] + 4.5 * (sample_pos[1] -
                                                   sample_pos[2]),
                        "Normal",
                        cex = 1,
                        font = 1,
                        pos = 4
                    )
                    points(
                        text_pos,
                        sample_pos[1] + 3 * (sample_pos[1] -
                                                 sample_pos[2]),
                        pch = 15,
                        col = outlier_color
                    )
                    text(
                        text_pos,
                        sample_pos[1] + 3 * (sample_pos[1] -
                                                 sample_pos[2]),
                        "Outlier",
                        cex = 1,
                        font = 1,
                        pos = 4
                    )
                }
                bh <- 0.1
                bw <- 0.85
                bspace <- 0.2 * bh
                if (hide_save_exit == FALSE) {
                    rect(0, 0, bw, bh)
                    text(0.5 * bw, 0.5 * bh, "Cancel", cex = 1.1)
                    rect(0, bspace + bh, bw, 2 * bh + bspace)
                    text(0.5 * bw, 0.5 * bh + bh + bspace, "OK",
                         cex = 1.1)
                    pca_button_adj <- 0.05
                    rect(0,
                         bspace + 3 * bh - pca_button_adj,
                         bw,
                         4 * bh + bspace - pca_button_adj)
                    text(0.5 * bw,
                         0.5 * bh + 3 * bh + bspace - pca_button_adj,
                         "PCA",
                         cex = 1.1)
                }
                par(mar = c(7, 5, 3, 0))
                title <-
                    paste(
                        "Outlier Identification : ",
                        basename(exp_file),
                        "  (",
                        sample_groups[current_group_no],
                        ":",
                        length(exp$gene_list[, 1]),
                        " Genes)",
                        sep = ""
                    )
                boxplot(
                    expression_data,
                    axes = FALSE,
                    frame = TRUE,
                    main = title,
                    col = "white",
                    medcol = boxplot_median_color,
                    border = sample_outlier_group,
                    ylab = "Expression (Log2)"
                )
                points(
                    1:length(data_median),
                    data_median,
                    pch = 16,
                    cex = 1.2,
                    col = sample_outlier_group
                )
                axis(
                    1,
                    1:length(sample_ids),
                    sample_ids,
                    cex.axis = 0.75,
                    las = 2
                )
                axis(2)
                abline(h = threshold, lty = 2)

                outlier_idx <- sample_outlier_group == outlier_color
                if (length(outlier_idx) > 0)
                    exp_outliers_all[[current_group_no]] <-
                    as.data.frame(colnames(expression_data)[outlier_idx])
                else
                    exp_outliers_all[[current_group_no]] <- as.data.frame(c(NA))
                outlier_list <- c()
                for (i in 1:length(sample_groups)) {
                    if (sum(is.na(exp_outliers_all[[i]])) == 0)
                        outlier_list <-
                            rbind(outlier_list, cbind(
                                as.character(exp_outliers_all[[i]][,
                                                                   1]),
                                rep(sample_groups[i], nrow(exp_outliers_all[[i]]))
                            ))
                }
                outlier_list <- as.data.frame(outlier_list)
                if (nrow(outlier_list) > 0) {
                    colnames(outlier_list) <- c("SampleID", "GroupID")
                    rownames(outlier_list) <-
                        1:length(outlier_list[,
                                              1])
                }
                outlier_list <- list(outlier_list)
                value_range <- par("usr")
                min_x <- value_range[1]
                max_x <- value_range[2]
                min_y <- value_range[3]
                max_y <- value_range[4]
                sample_top <- seq(min_y, max_y, sample_legend_grid *
                                      (max_y - min_y))[sample_legend_top]
                sample_sep <-
                    sample_sep_scale * sample_legend_grid *
                    (max_y - min_y)
                sample_pos <-
                    seq(sample_top - (length(sample_groups) -
                                          1) * sample_sep,
                        sample_top,
                        sample_sep)
                sample_pos <- rev(sample_pos)
            }
            draw_plot <- TRUE
            if (hide_save_exit == TRUE) {
                if (ret_cancel)
                    ret <- FALSE
                else
                    ret <- outlier_list
                dev.off()
                return(ret)
            }
            #xy <- locator(n = 1, type = "n")
            xy <- data.frame('x' = 262.5875, 'y' = 2.202172)
            if (pca_plot == TRUE)
                xy <- data.frame('x' = 263.7925, 'y' = 4.733315)
            #xy$x <-262.5875
            #xy$y <-2.202172
            xP <- xy$x
            yP <- xy$y
            xyN <- round(xP)
            value_range <- par("usr")
            min_x <- value_range[1]
            max_x <- value_range[2]
            min_y <- value_range[3]
            max_y <- value_range[4]
            if (init_flag == TRUE || (xP < max_x && xP > min_x &&
                                      yP > min_y && yP < max_y)) {
                init_flag <- FALSE
                removed_n <- match(xyN, xyN_all, nomatch = 0)
                if (removed_n > 0 &&
                    xyN > 0 && xyN < length(sample_ids) +
                    1) {
                    xyN_all <- xyN_all[-removed_n]
                }
                else if (xyN > 0 && xyN < length(sample_ids) + 1) {
                    xyN_all <- c(xyN_all, xyN)
                }
            }
            else {
                draw_plot <- FALSE
            }
            yP_no <- abs(sample_pos - yP) < 0.5 * (sample_pos[1] -
                                                       sample_pos[2])
            if (!is.na(yP_no[1])) {
                if (xyN > max_x && sum(yP_no) > 0) {
                    current_groupno <- c(1:length(yP_no))[yP_no]
                    if (current_group_no != current_groupno) {
                        current_group_no <- current_groupno
                        dev.off(display_window)
                        current_group_no <- list(current_group_no,
                                                 exp_outliers_all)
                        return(current_group_no)
                    }
                }
            }
            #         if (xyN > max_x && xyN < max_x + 0.1 * (max_x - min_x) &&
            #             yP > min_y + (0.3 - pca_button_adj) * (max_y - min_y) &&
            #             yP < min_y + (0.415 - pca_button_adj) * (max_y -
            #                                                      min_y))
            if (pca_plot == TRUE) {
                selected_outliers <-
                    as.data.frame(exp_outliers_all[[current_group_no]])
                colnames(selected_outliers) <- "Outlier"
                if (is.na(selected_outliers[1, 1]))
                    selected_outliers <- FALSE
                itn_displayPCAWithOutlierList(exp, outliers = selected_outliers)
                hide_save_exit <- TRUE
                dev.set(which = myplot)
            }
            #         if (xyN > max_x && xyN < max_x + 0.1 * (max_x - min_x) &&
            #             yP > min_y + 0.115 * (max_y - min_y) && yP < min_y +
            #             0.23 * (max_y - min_y))
            if (pca_plot == FALSE) {
                hide_save_exit <- TRUE
            }
            if (xyN > max_x &&
                xyN < max_x + 0.1 * (max_x - min_x) &&
                yP > min_y &&
                yP < min_y + 0.115 * (max_y - min_y)) {
                hide_save_exit <- TRUE
                ret_cancel <- TRUE
            }

        }

        return(ret)
    }


itn_getMatchingOriginalName <-
    function (new, original, ignore_case = FALSE)
    {
        if (ignore_case == FALSE)
            return(new)
        org_sub <- original[itn_isElementIgnoreCase(original, new,
                                                    ignore_case)]
        if (length(org_sub) != length(new))
            stop("Error in itn_getOriginalName.")
        updated_new <- original[toupper(original) %in% toupper(new)]
        updated_new_mtx <- as.matrix(updated_new)
        rownames(updated_new_mtx) <- toupper(updated_new)
        updated_new <- as.character(updated_new_mtx[toupper(new),])
        return(updated_new)
    }


itn_isSubSetIgnoreCase <- function (subset, set, ignore_case = FALSE)
{
    subset <- as.vector(subset)
    set <- as.vector(set)
    if (length(set) == 0 || length(subset) == 0) {
        return(FALSE)
    }
    if (ignore_case == TRUE) {
        subset <- tolower(subset)
        set <- tolower(set)
    }
    if (length((intersect(subset, set))) == length(subset)) {
        return(TRUE)
    }
    else {
        return(FALSE)
    }
}



itn_saveData <-
    function (data,
              file = TRUE,
              save_title = "",
              init_filename = "",
              back_to_r_console = FALSE)
    {
        model <- data
        #require("tcltk", quiet = TRUE)
        is_sc_object <- FALSE
        if (save_title == "") {
            save_title = "Save Object..."
        }
        if (is.list(data) &&
            sum(match(names(data), "obj_type", nomatch = 0)) >
            0) {
            is_sc_object <- TRUE
        }
        if (file == TRUE)
            file <- paste0(save_title, '.fso')
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
                if (ncol(data) == 0 ||
                    is.na(as.numeric(rownames(data)[1])) ==
                    FALSE) {
                    rowname_flag <- FALSE
                }
                if (ncol(data) == 0 ||
                    (colnames(data))[1] == "V1") {
                    colname_flag <- FALSE
                }
                if (rowname_flag == TRUE && colname_flag == TRUE) {
                    new_col_names <- t(as.matrix(c("ID", colnames(
                        data
                    ))))
                    write.table(
                        new_col_names,
                        file = file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        append = FALSE
                    )
                    write.table(
                        data,
                        file = file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = TRUE,
                        col.names = FALSE,
                        append = TRUE
                    )
                }
                else {
                    write.table(
                        data,
                        file = file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = rowname_flag,
                        col.names = colname_flag,
                        append = FALSE
                    )
                }
            }
            else {
                write.table(
                    "fluidigmSC Object File Version 1.0",
                    file,
                    sep = "\t",
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    append = FALSE
                )
                for (i_name in 1:length(model)) {
                    write.table(
                        "",
                        file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        append = TRUE
                    )
                    begin_line <-
                        paste("<start::", names(model)[i_name],
                              ">", sep = "")
                    write.table(
                        begin_line,
                        file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        append = TRUE
                    )
                    i_model <- model[[i_name]]
                    rowname_flag <- TRUE
                    colname_flag <- TRUE
                    rowname1 <- as.numeric(rownames(i_model)[1])
                    if (ncol(i_model) == 0 || is.na(rowname1) ==
                        FALSE) {
                        rowname_flag <- FALSE
                    }
                    if (ncol(i_model) == 0 ||
                        (colnames(i_model))[1] ==
                        "V1") {
                        colname_flag <- FALSE
                    }
                    rowname_line <- paste("rowname::", rowname_flag,
                                          sep = "")
                    colname_line <- paste("colname::", colname_flag,
                                          sep = "")
                    rowcolnames_lines <-
                        rbind(rowname_line, colname_line)
                    write.table(
                        rowcolnames_lines,
                        file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        append = TRUE
                    )
                    if (rowname_flag == FALSE && colname_flag ==
                        FALSE) {
                        write.table(
                            i_model,
                            file,
                            sep = "\t",
                            quote = FALSE,
                            row.names = FALSE,
                            col.names = FALSE,
                            append = TRUE
                        )
                    }
                    else if (rowname_flag == TRUE && colname_flag ==
                             TRUE) {
                        new_col_names <- t(as.matrix(c(
                            "ID", colnames(i_model)
                        )))
                        write.table(
                            new_col_names,
                            file = file,
                            sep = "\t",
                            quote = FALSE,
                            row.names = FALSE,
                            col.names = FALSE,
                            append = TRUE
                        )
                        write.table(
                            i_model,
                            file = file,
                            sep = "\t",
                            quote = FALSE,
                            row.names = TRUE,
                            col.names = FALSE,
                            append = TRUE
                        )
                    }
                    else if (rowname_flag == FALSE &&
                             colname_flag ==
                             TRUE) {
                        write.table(
                            i_model,
                            file,
                            sep = "\t",
                            quote = FALSE,
                            row.names = FALSE,
                            col.names = TRUE,
                            append = TRUE
                        )
                    }
                    end_line <-
                        paste("<end::", names(model)[i_name],
                              ">", sep = "")
                    write.table(
                        end_line,
                        file,
                        sep = "\t",
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        append = TRUE
                    )
                }
            }
            options(ops)
        }
        invisible(file)
    }


itn_ScatterPlot <-
    function (x,
              y,
              xy_annotation,
              id_clr_symb = NULL,
              xlab,
              ylab,
              locate = TRUE,
              plot_title = "",
              x_greater_thre = NA,
              x_less_thre = NA,
              y_greater_thre = NA,
              y_less_thre = NA,
              xlim = NA,
              ylim = NA,
              show_group = FALSE,
              group_color = FALSE,
              max_group_length = 15)
    {
        if (!require(SDMTools, quietly = TRUE))
            itn_LoadPackageMsg()
        if (class(x) != "numeric" || class(y) != "numeric")
            stop("Numeric values are needed.")
        if (length(x) != length(y))
            stop("The two numeric vectors must have the same length.")
        if (class(xy_annotation) != "data.frame")
            stop("Annotation in data frame is needed.")
        if (length(x) != nrow(xy_annotation))
            stop(
                "The number of rows in annotation data frame must be the same as the lengths ofthe two numeric vectors."
            )
        if (ncol(xy_annotation) < 2) {
            if (ncol(xy_annotation) == 0)
                xy_annotation <- as.data.frame(1:length(x))
            xy_annotation <- cbind(xy_annotation, rep("UNTITLE",
                                                      nrow(xy_annotation)))
            colnames(xy_annotation)[2] <- "GroupID"
            xy_annotation <- as.data.frame(xy_annotation)
        }
        if (!is.na(xlim[1])) {
            x[x < xlim[1]] <- xlim[1]
            x[x > xlim[2]] <- xlim[2]
        }
        if (!is.na(ylim[1])) {
            y[y < ylim[1]] <- ylim[1]
            y[y > ylim[2]] <- ylim[2]
        }
        x_y <- cbind(x, y)
        xy_ids <- as.character(xy_annotation[, 1])
        group_id <- as.character(xy_annotation[, 2])
        xy_group_num <- length(unique(group_id))
        xy_group_symbol_type <-
            itn_getGraphSymbols(unique(group_id),
                                id_clr_symb)
        if (class(group_color) == "logical")
            group_color <-
            itn_getColors(unique(group_id), id_clr_symb)
        unique_group_id <- unique(group_id)
        if (length(group_color) != length(unique_group_id))
            stop("group_color number should be the same as group number!")
        xy_color <- rep("", length(group_id))
        xy_symbol <- rep(1, length(group_id))
        for (k in 1:length(unique_group_id)) {
            xy_color[group_id == unique_group_id[k]] <- group_color[k]
            xy_symbol[group_id == unique_group_id[k]] <-
                xy_group_symbol_type[k]
        }
        #itn_newPlotWindow(title = plot_title)
        pdf(paste0(plot_title, '.pdf'))
        located_xyIDs = NULL
        assign("fldm_located_xyIDs", located_xyIDs, env = .GlobalEnv)
        locate_choice <- 1
        polygon_point_num <- 20
        button_name <- c("Done", "Save", "Clear", "Circle", "Point")
        button_position_no <- c(1, 2, 4, 6, 7)
        button_height <- 0.08
        button_width <- 0.8
        pixel_dist <- 8
        label_cex <- 0.6
        button_cex <- 0.9
        legend_cex <- 0.7
        max_legend_nchar <-
            max(nchar(as.character(unique_group_id)))
        if (max_legend_nchar > 20)
            legend_cex <- 0.4
        else if (max_legend_nchar > 15)
            legend_cex <- 0.5
        else if (max_legend_nchar > 10)
            legend_cex <- 0.6
        group_unique <- unique(group_id)
        group_id_trim <- substr(group_unique, 1, max_group_length)
        group_id_omit <- rep("...", length(group_unique))
        group_id_omit[nchar(as.character(group_unique)) <= max_group_length] <-
            ""
        group_id_legend <-
            paste(group_id_trim, group_id_omit, sep = "")
        xP_all <- c()
        yP_all <- c()
        xyN_all <- c()
        out_range <- 0
        current_button <- button_name[5]
        click_status <- button_name[5]
        hide_save_exit <- FALSE
        done_locate <- FALSE
        draw_plot = TRUE
        while (done_locate == FALSE) {
            if (draw_plot) {
                layout.matrix <- matrix(c(2, 1), ncol = 2)
                layout.widths <- c(6, 1)
                layout.heights <- 5
                layout(
                    layout.matrix,
                    widths = layout.widths,
                    heights = layout.heights,
                    respect = FALSE
                )
                par(mar = c(5, 0, 3, 1))
                frame()
                value_range <- par("usr")
                min_x <- value_range[1]
                max_x <- value_range[2]
                min_y <- value_range[3]
                max_y <- value_range[4]
                button_down <-
                    seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no]
                button_up <-
                    seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no] +
                    button_height * (max_y - min_y)
                if (hide_save_exit == FALSE && locate == TRUE) {
                    for (i in 1:length(button_name)) {
                        bname <- button_name[i]
                        rect(0, button_down[i], button_width, button_up[i])
                        if (click_status == button_name[i])
                            button_font <- 2
                        else
                            button_font <- 1
                        text(
                            0.5 * button_width,
                            0.5 * (button_down[i] +
                                       button_up[i]),
                            bname,
                            cex = 1.1,
                            font = button_font
                        )
                    }
                }
                if (xy_group_num > 1)
                    legend(
                        "topleft",
                        legend = group_id_legend,
                        col = group_color,
                        pt.lwd = 1,
                        pch = xy_group_symbol_type,
                        cex = legend_cex,
                        box.lwd = 0,
                        border = "white",
                        box.col = NA,
                        title = ""
                    )
                par(mar = c(5, 5, 3, 0))
                plot(
                    x,
                    y,
                    main = plot_title,
                    xlab = xlab,
                    ylab = ylab,
                    type = "n"
                )
                value_range <- par("usr")
                min_x <- value_range[1]
                max_x <- value_range[2]
                min_y <- value_range[3]
                max_y <- value_range[4]
                index_black <- c()
                x_defined <-
                    !is.na(x_less_thre) || !is.na(x_greater_thre)
                y_defined <-
                    !is.na(y_less_thre) || !is.na(y_greater_thre)
                if (x_defined || y_defined) {
                    x_index_black <- c(1:length(xy_color))
                    if (!is.na(x_less_thre))
                        x_index_black <-
                            intersect(x_index_black, c(1:length(xy_color))[x >
                                                                               x_less_thre])
                    if (!is.na(x_greater_thre))
                        x_index_black <-
                            intersect(x_index_black, c(1:length(xy_color))[x <
                                                                               x_greater_thre])
                    y_index_black <- c(1:length(xy_color))
                    if (!is.na(y_less_thre))
                        y_index_black <-
                        intersect(y_index_black, c(1:length(xy_color))[y >
                                                                           y_less_thre])
                    if (!is.na(y_greater_thre))
                        y_index_black <-
                        intersect(y_index_black, c(1:length(xy_color))[y <
                                                                           y_greater_thre])
                    if (x_defined && y_defined)
                        index_black <-
                        union(x_index_black, y_index_black)
                    else if (x_defined && !y_defined)
                        index_black <- x_index_black
                    else
                        index_black <- y_index_black
                }
                if (length(index_black) > 0)
                    xy_color[index_black] <- "black"
                rect(min_x, min_y, max_x, max_y, col = itn_getBackgroundColor())
                abline(h = 0,
                       v = 0,
                       col = "white")
                points(
                    x,
                    y,
                    lwd = 1,
                    cex = 0.8,
                    col = xy_color,
                    pch = xy_symbol
                )
                if (!is.na(x_less_thre))
                    abline(v = x_less_thre, lty = 2)
                if (!is.na(x_greater_thre))
                    abline(v = x_greater_thre, lty = 2)
                if (!is.na(y_less_thre))
                    abline(h = y_less_thre, lty = 2)
                if (!is.na(y_greater_thre))
                    abline(h = y_greater_thre, lty = 2)
                label_position_all <- rep(4, length(x))
                label_position_all[x - min_x > 0.75 * (max_x - min_x)] <-
                    2
                plot_range <- par("plt")
                pixel_size <- dev.size("px")
                pixel_size[1] <- pixel_size[1] * (plot_range[2] -
                                                      plot_range[1])
                pixel_size[2] <- pixel_size[2] * (plot_range[4] -
                                                      plot_range[3])
                if (show_group == "") {
                    text(x, y, xy_ids[1:length(xy_ids)], pos = label_position_all,
                         cex = label_cex)
                }
                if (show_group != "" && show_group != FALSE) {
                    if (show_group == "all") {
                        displayed_no <- c(1:length(xy_ids))
                    }
                    else {
                        displayed_no <- c(1:length(xy_ids))[as.vector(group_id) ==
                                                                show_group]
                        if (length(displayed_no) == 0)
                            stop(paste("No sample for group ", show_group,
                                       sep = ""))
                    }
                    text(x[displayed_no],
                         y[displayed_no],
                         xy_ids[displayed_no],
                         pos = label_position_all[displayed_no],
                         cex = label_cex)
                }
                button_down <-
                    seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no]
                button_up <-
                    seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no] +
                    button_height * (max_y - min_y)
            }
            draw_plot = TRUE
            circling_cancelled = FALSE
            if (locate == FALSE) {
                done_locate = TRUE
            }
            else if (hide_save_exit == TRUE) {
                if (length(xyN_all) > 0) {
                    text(
                        xP_all,
                        yP_all,
                        xy_ids[xyN_all],
                        pos = label_position_all[xyN_all],
                        cex = label_cex,
                        font = 3
                    )
                }
                done_locate = TRUE
            }
            else {
                if (length(xyN_all) > 0) {
                    text(
                        xP_all,
                        yP_all,
                        xy_ids[xyN_all],
                        pos = label_position_all[xyN_all],
                        cex = label_cex,
                        font = 3
                    )
                    cex = 1.2
                    points(x[xyN_all],
                           y[xyN_all],
                           pch = 16,
                           cex = cex,
                           col = xy_color[xyN_all])
                    points(x[xyN_all],
                           y[xyN_all],
                           pch = xy_symbol[xyN_all],
                           cex = cex,
                           col = xy_color[xyN_all])
                    located_xyIDs = xy_annotation[xyN_all, 1:2]
                    assign("fldm_located_xyIDs", located_xyIDs, env = .GlobalEnv)
                }
                if (click_status == button_name[4]) {
                    locate_choice <- polygon_point_num
                }
                if (locate_choice > 1) {
                    xp <- c()
                    yp <- c()
                    for (j in 1:locate_choice) {
                        xyp <- locator(n = 1, type = "o")
                        xp <- c(xp, xyp$x)
                        yp <- c(yp, xyp$y)
                        if (length(xp) > 1) {
                            points(xp[(length(xp) - 1):length(xp)], yp[(length(xp) -
                                                                            1):length(xp)], type = "l")
                            xs <-
                                (xp[1] - min_x) * pixel_size[1] / (max_x -
                                                                       min_x)
                            ys <-
                                (yp[1] - min_y) * pixel_size[2] / (max_y -
                                                                       min_y)
                            xe <-
                                (xp[length(xp)] - min_x) * pixel_size[1] / (max_x -
                                                                                min_x)
                            ye <-
                                (yp[length(xp)] - min_y) * pixel_size[2] / (max_y -
                                                                                min_y)
                            start_point <- c(xs, ys)
                            end_point <- c(xe, ye)
                            done_circling <- dist(rbind(start_point,
                                                        end_point)) <= pixel_dist
                            if (done_circling == TRUE) {
                                Sys.sleep(0.5)
                                break
                            }
                        }
                        if (xyp$x > max_x ||
                            xyp$x < min_x || xyp$y >
                            max_y || xyp$y < min_y) {
                            locate_choice <- 1
                            if (j > 1)
                                circling_cancelled = TRUE
                            break
                        }
                    }
                    xy <- list(x = xp, y = yp)
                    if (locate_choice == 1)
                        xy <- xyp
                }
                else {
                    xy <- locator(n = 1, type = "n")
                }
                xP <- xy$x
                yP <- xy$y
                xP_yP <- cbind(xP, yP)
                button_clicked = FALSE
                last_x <- xP[length(xP)]
                last_y <- yP[length(yP)]
                if (last_x > max_x &&
                    last_x < max_x + 0.125 * (max_x -
                                              min_x)) {
                    button_no <- c(1:length(button_down))[button_down <
                                                              last_y &
                                                              last_y < button_up]
                    if (length(button_no) > 0) {
                        button_clicked = TRUE
                        current_button <- button_name[button_no]
                        if (current_button == button_name[5]) {
                            locate_choice <- 1
                            click_status <- button_name[5]
                        }
                        if (current_button == button_name[4]) {
                            locate_choice <- polygon_point_num
                            click_status <- button_name[4]
                        }
                        if (current_button == button_name[3]) {
                            xP_all <- c()
                            yP_all <- c()
                            xyN_all <- c()
                            located_xyIDs = NULL
                            assign("fldm_located_xyIDs",
                                   located_xyIDs,
                                   env = .GlobalEnv)
                        }
                        if (current_button == button_name[2]) {
                            if (located_xyIDs == "" || nrow(located_xyIDs) ==
                                0) {
                                tkmessageBox(
                                    title = "fluidigmSC Analysis",
                                    message = "Nothing selected to save.",
                                    icon = "info",
                                    type = "ok"
                                )
                            }
                            else {
                                file <- tclvalue(
                                    tkgetSaveFile(
                                        title = "Save Selected IDs to File",
                                        filetypes = "{{Delimited File} {.txt}}"
                                    )
                                )
                                if (file != "") {
                                    file_ext = tools::file_ext(file)
                                    if (file_ext == "")
                                        file = paste(file, "txt", sep = ".")
                                    write.table(
                                        located_xyIDs,
                                        file = file,
                                        sep = "\t",
                                        quote = FALSE,
                                        row.names = FALSE,
                                        col.names = TRUE
                                    )
                                }
                            }
                        }
                        if (current_button == button_name[1]) {
                            hide_save_exit <- TRUE
                        }
                    }
                }
                if (button_clicked == FALSE && circling_cancelled ==
                    FALSE) {
                    if (locate_choice > 1) {
                        xy_in_polygon <- pnt.in.poly(x_y, xP_yP)
                        xyN <- c(1:length(x))[xy_in_polygon[, 3] ==
                                                  1]
                    }
                    else {
                        xyP <- rbind(c(xP, yP), x_y)
                        x_yP1 <-
                            (xyP[, 1] - min_x) * pixel_size[1] / (max_x -
                                                                      min_x)
                        x_yP2 <-
                            (xyP[, 2] - min_y) * pixel_size[2] / (max_y -
                                                                      min_y)
                        xyP <- cbind(x_yP1, x_yP2)
                        xyP_xy_dist <-
                            apply(xyP[-1,], 1, itn_twoVectorDistance,
                                  y = xyP[1,])
                        xyN <-
                            c(1:length(x))[xyP_xy_dist == min(xyP_xy_dist)]
                        if (min(xyP_xy_dist) >= pixel_dist) {
                            draw_plot <- FALSE
                        }
                    }
                    if (locate_choice > 1 ||
                        min(xyP_xy_dist) < pixel_dist) {
                        removed_n <- match(xyN, xyN_all, nomatch = 0)
                        removed_n <- removed_n[removed_n > 0]
                        if (length(removed_n) > 0) {
                            xP_all <- xP_all[-removed_n]
                            yP_all <- yP_all[-removed_n]
                            xyN_all <- xyN_all[-removed_n]
                        }
                        else {
                            xP_all <- c(xP_all, x[xyN])
                            yP_all <- c(yP_all, y[xyN])
                            xyN_all <- c(xyN_all, xyN)
                        }
                    }
                }
            }
        }
        dev.off()
        return(located_xyIDs)
    }
