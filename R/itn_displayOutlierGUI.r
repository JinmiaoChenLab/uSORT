itn_displayOutlierGUI <- 
function (exp_trimmed_all, exp_file, exp_outliers_all, exp_threshold_all, 
          current_group_no = 1) 
{
    pca_plot <- TRUE
    exp <- exp_trimmed_all[[current_group_no]]
    outliers <- exp_outliers_all[[current_group_no]]
    threshold <- exp_threshold_all[[current_group_no]]
    sample_groups <- names(exp_outliers_all)
    title <- paste("Outlier Identification : ", basename(exp_file), 
                   sep = "")
#     itn_newPlotWindow(width = 11, height = 4.5, title = title, 
#                       xpos = 0)
    pdf('outlier.boxplot.pdf',width = 11, height = 4.5, title = title)
    display_window <- dev.cur()
    expression_data_all <- exp$log2ex_data
    expression_data_all <- apply(expression_data_all, 2, as.numeric)
    data_median <- apply(expression_data_all, 2, median)
    data_median_sorted <- sort(data_median, decreasing = T)
    data_median_rownames <- rownames(as.matrix(data_median_sorted))
    expression_data_all <- expression_data_all[, data_median_rownames]
    data_median <- data_median_sorted
    d<-data.frame(data_median)
    d$QC=(d>=threshold)
    colnames(d)<-c('median','passQC')
    write.table(d, row.names = T,sep='\t',col.names = NA,
                file=paste0('outliers.data.median.',basename(exp_file)))
    sample_ids <- colnames(expression_data_all)
    outliers <- as.matrix(outliers)[, 1]
    outlier_index <- c(1:ncol(expression_data_all))[is.element(sample_ids, 
                                                               outliers)]
    outlier_num <- length(outlier_index)
    outlier_color <- "red"
    non_outlier_color <- "blue"
    sample_outlier_group_all <- rep(non_outlier_color, ncol(expression_data_all))
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
            sample_outlier_group[setdiff(xyN_all, outlier_index)] <- outlier_color
        }
        if (draw_plot == TRUE || hide_save_exit == TRUE) {
            layout.matrix <- matrix(c(2, 1), ncol = 2, byrow = TRUE)
            layout.widths <- c(4, 0.4)
            layout.heights <- 1
            layout(layout.matrix, widths = layout.widths, heights = layout.heights, 
                   respect = FALSE)
            par(mar = c(6.5, 0, 3, 1))
            frame()
            value_range <- par("usr")
            min_x <- value_range[1]
            max_x <- value_range[2]
            min_y <- value_range[3]
            max_y <- value_range[4]
            sample_top <- seq(min_y, max_y, sample_legend_grid * 
                                  (max_y - min_y))[sample_legend_top]
            sample_sep <- sample_sep_scale * sample_legend_grid * 
                (max_y - min_y)
            sample_pos <- seq(sample_top - (length(sample_groups) - 
                                                1) * sample_sep, sample_top, sample_sep)
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
                    points(text_pos, sample_pos[i], pch = sample_pch, 
                           col = sample_col)
                    text(text_pos, sample_pos[i], sample_groups[i], 
                         cex = sample_cex, font = sample_font, pos = 4)
                }
                points(text_pos, sample_pos[1] + 4.5 * (sample_pos[1] - 
                                                            sample_pos[2]), pch = 15, col = non_outlier_color)
                text(text_pos, sample_pos[1] + 4.5 * (sample_pos[1] - 
                                                          sample_pos[2]), "Normal", cex = 1, font = 1, 
                     pos = 4)
                points(text_pos, sample_pos[1] + 3 * (sample_pos[1] - 
                                                          sample_pos[2]), pch = 15, col = outlier_color)
                text(text_pos, sample_pos[1] + 3 * (sample_pos[1] - 
                                                        sample_pos[2]), "Outlier", cex = 1, font = 1, 
                     pos = 4)
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
                rect(0, bspace + 3 * bh - pca_button_adj, bw, 
                     4 * bh + bspace - pca_button_adj)
                text(0.5 * bw, 0.5 * bh + 3 * bh + bspace - pca_button_adj, 
                     "PCA", cex = 1.1)
            }
            par(mar = c(7, 5, 3, 0))
            title <- paste("Outlier Identification : ", basename(exp_file), 
                           "  (", sample_groups[current_group_no], ":", 
                           length(exp$gene_list[, 1]), " Genes)", sep = "")
            boxplot(expression_data, axes = FALSE, frame = TRUE, 
                    main = title, col = "white", medcol = boxplot_median_color, 
                    border = sample_outlier_group, ylab = "Expression (Log2)")
            points(1:length(data_median), data_median, pch = 16, 
                   cex = 1.2, col = sample_outlier_group)
            axis(1, 1:length(sample_ids), sample_ids, cex.axis = 0.75, 
                 las = 2)
            axis(2)
            abline(h = threshold, lty = 2)
            
            outlier_idx <- sample_outlier_group == outlier_color
            if (length(outlier_idx) > 0) 
                exp_outliers_all[[current_group_no]] <- as.data.frame(colnames(expression_data)[outlier_idx])
            else exp_outliers_all[[current_group_no]] <- as.data.frame(c(NA))
            outlier_list <- c()
            for (i in 1:length(sample_groups)) {
                if (sum(is.na(exp_outliers_all[[i]])) == 0) 
                    outlier_list <- rbind(outlier_list, cbind(as.character(exp_outliers_all[[i]][, 
                                                                                                 1]), rep(sample_groups[i], nrow(exp_outliers_all[[i]]))))
            }
            outlier_list <- as.data.frame(outlier_list)
            if (nrow(outlier_list) > 0) {
                colnames(outlier_list) <- c("SampleID", "GroupID")
                rownames(outlier_list) <- 1:length(outlier_list[, 
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
            sample_sep <- sample_sep_scale * sample_legend_grid * 
                (max_y - min_y)
            sample_pos <- seq(sample_top - (length(sample_groups) - 
                                                1) * sample_sep, sample_top, sample_sep)
            sample_pos <- rev(sample_pos)
        }
        draw_plot <- TRUE
        if (hide_save_exit == TRUE) {
            if (ret_cancel) 
                ret <- FALSE
            else ret <- outlier_list
            dev.off() 
            return(ret)
        }
        #xy <- locator(n = 1, type = "n")
        xy <- data.frame('x'=262.5875, 'y'=2.202172)
        if(pca_plot == TRUE)xy <- data.frame('x'=263.7925, 'y'=4.733315)
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
            if (removed_n > 0 && xyN > 0 && xyN < length(sample_ids) + 
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
        if(pca_plot==TRUE){
            selected_outliers <- as.data.frame(exp_outliers_all[[current_group_no]])
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
        if(pca_plot==FALSE){
            hide_save_exit <- TRUE
        }
        if (xyN > max_x && xyN < max_x + 0.1 * (max_x - min_x) && 
            yP > min_y && yP < min_y + 0.115 * (max_y - min_y)) {
            hide_save_exit <- TRUE
            ret_cancel <- TRUE
        }
        
    }
    
    return(ret)
}