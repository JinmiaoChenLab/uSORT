itn_ScatterPlot <- 
function (x, y, xy_annotation, id_clr_symb = NULL, xlab, ylab, 
          locate = TRUE, plot_title = "", x_greater_thre = NA, x_less_thre = NA, 
          y_greater_thre = NA, y_less_thre = NA, xlim = NA, ylim = NA, 
          show_group = FALSE, group_color = FALSE, max_group_length = 15) 
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
        stop("The number of rows in annotation data frame must be the same as the lengths ofthe two numeric vectors.")
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
    xy_group_symbol_type <- itn_getGraphSymbols(unique(group_id), 
                                                id_clr_symb)
    if (class(group_color) == "logical") 
        group_color <- itn_getColors(unique(group_id), id_clr_symb)
    unique_group_id <- unique(group_id)
    if (length(group_color) != length(unique_group_id)) 
        stop("group_color number should be the same as group number!")
    xy_color <- rep("", length(group_id))
    xy_symbol <- rep(1, length(group_id))
    for (k in 1:length(unique_group_id)) {
        xy_color[group_id == unique_group_id[k]] <- group_color[k]
        xy_symbol[group_id == unique_group_id[k]] <- xy_group_symbol_type[k]
    }
    #itn_newPlotWindow(title = plot_title)
    pdf(paste0(plot_title,'.pdf'))
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
    max_legend_nchar <- max(nchar(as.character(unique_group_id)))
    if (max_legend_nchar > 20) 
        legend_cex <- 0.4
    else if (max_legend_nchar > 15) 
        legend_cex <- 0.5
    else if (max_legend_nchar > 10) 
        legend_cex <- 0.6
    group_unique <- unique(group_id)
    group_id_trim <- substr(group_unique, 1, max_group_length)
    group_id_omit <- rep("...", length(group_unique))
    group_id_omit[nchar(as.character(group_unique)) <= max_group_length] <- ""
    group_id_legend <- paste(group_id_trim, group_id_omit, sep = "")
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
            layout(layout.matrix, widths = layout.widths, heights = layout.heights, 
                   respect = FALSE)
            par(mar = c(5, 0, 3, 1))
            frame()
            value_range <- par("usr")
            min_x <- value_range[1]
            max_x <- value_range[2]
            min_y <- value_range[3]
            max_y <- value_range[4]
            button_down <- seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no]
            button_up <- seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no] + 
                button_height * (max_y - min_y)
            if (hide_save_exit == FALSE && locate == TRUE) {
                for (i in 1:length(button_name)) {
                    bname <- button_name[i]
                    rect(0, button_down[i], button_width, button_up[i])
                    if (click_status == button_name[i]) 
                        button_font <- 2
                    else button_font <- 1
                    text(0.5 * button_width, 0.5 * (button_down[i] + 
                                                        button_up[i]), bname, cex = 1.1, font = button_font)
                }
            }
            if (xy_group_num > 1) 
                legend("topleft", legend = group_id_legend, col = group_color, 
                       pt.lwd = 1, pch = xy_group_symbol_type, cex = legend_cex, 
                       box.lwd = 0, border = "white", box.col = NA, 
                       title = "")
            par(mar = c(5, 5, 3, 0))
            plot(x, y, main = plot_title, xlab = xlab, ylab = ylab, 
                 type = "n")
            value_range <- par("usr")
            min_x <- value_range[1]
            max_x <- value_range[2]
            min_y <- value_range[3]
            max_y <- value_range[4]
            index_black <- c()
            x_defined <- !is.na(x_less_thre) || !is.na(x_greater_thre)
            y_defined <- !is.na(y_less_thre) || !is.na(y_greater_thre)
            if (x_defined || y_defined) {
                x_index_black <- c(1:length(xy_color))
                if (!is.na(x_less_thre)) 
                    x_index_black <- intersect(x_index_black, c(1:length(xy_color))[x > 
                                                                                        x_less_thre])
                if (!is.na(x_greater_thre)) 
                    x_index_black <- intersect(x_index_black, c(1:length(xy_color))[x < 
                                                                                        x_greater_thre])
                y_index_black <- c(1:length(xy_color))
                if (!is.na(y_less_thre)) 
                    y_index_black <- intersect(y_index_black, c(1:length(xy_color))[y > 
                                                                                        y_less_thre])
                if (!is.na(y_greater_thre)) 
                    y_index_black <- intersect(y_index_black, c(1:length(xy_color))[y < 
                                                                                        y_greater_thre])
                if (x_defined && y_defined) 
                    index_black <- union(x_index_black, y_index_black)
                else if (x_defined && !y_defined) 
                    index_black <- x_index_black
                else index_black <- y_index_black
            }
            if (length(index_black) > 0) 
                xy_color[index_black] <- "black"
            rect(min_x, min_y, max_x, max_y, col = itn_getBackgroundColor())
            abline(h = 0, v = 0, col = "white")
            points(x, y, lwd = 1, cex = 0.8, col = xy_color, 
                   pch = xy_symbol)
            if (!is.na(x_less_thre)) 
                abline(v = x_less_thre, lty = 2)
            if (!is.na(x_greater_thre)) 
                abline(v = x_greater_thre, lty = 2)
            if (!is.na(y_less_thre)) 
                abline(h = y_less_thre, lty = 2)
            if (!is.na(y_greater_thre)) 
                abline(h = y_greater_thre, lty = 2)
            label_position_all <- rep(4, length(x))
            label_position_all[x - min_x > 0.75 * (max_x - min_x)] <- 2
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
                text(x[displayed_no], y[displayed_no], xy_ids[displayed_no], 
                     pos = label_position_all[displayed_no], cex = label_cex)
            }
            button_down <- seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no]
            button_up <- seq(min_y, max_y, 0.1 * (max_y - min_y))[button_position_no] + 
                button_height * (max_y - min_y)
        }
        draw_plot = TRUE
        circling_cancelled = FALSE
        if (locate == FALSE) {
            done_locate = TRUE
        }
        else if (hide_save_exit == TRUE) {
            if (length(xyN_all) > 0) {
                text(xP_all, yP_all, xy_ids[xyN_all], pos = label_position_all[xyN_all], 
                     cex = label_cex, font = 3)
            }
            done_locate = TRUE
        }
        else {
            if (length(xyN_all) > 0) {
                text(xP_all, yP_all, xy_ids[xyN_all], pos = label_position_all[xyN_all], 
                     cex = label_cex, font = 3)
                cex = 1.2
                points(x[xyN_all], y[xyN_all], pch = 16, cex = cex, 
                       col = xy_color[xyN_all])
                points(x[xyN_all], y[xyN_all], pch = xy_symbol[xyN_all], 
                       cex = cex, col = xy_color[xyN_all])
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
                        xs <- (xp[1] - min_x) * pixel_size[1]/(max_x - 
                                                                   min_x)
                        ys <- (yp[1] - min_y) * pixel_size[2]/(max_y - 
                                                                   min_y)
                        xe <- (xp[length(xp)] - min_x) * pixel_size[1]/(max_x - 
                                                                            min_x)
                        ye <- (yp[length(xp)] - min_y) * pixel_size[2]/(max_y - 
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
                    if (xyp$x > max_x || xyp$x < min_x || xyp$y > 
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
            if (last_x > max_x && last_x < max_x + 0.125 * (max_x - 
                                                            min_x)) {
                button_no <- c(1:length(button_down))[button_down < 
                                                          last_y & last_y < button_up]
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
                        assign("fldm_located_xyIDs", located_xyIDs, 
                               env = .GlobalEnv)
                    }
                    if (current_button == button_name[2]) {
                        if (located_xyIDs == "" || nrow(located_xyIDs) == 
                            0) {
                            tkmessageBox(title = "fluidigmSC Analysis", 
                                         message = "Nothing selected to save.", 
                                         icon = "info", type = "ok")
                        }
                        else {
                            file <- tclvalue(tkgetSaveFile(title = "Save Selected IDs to File", 
                                                           filetypes = "{{Delimited File} {.txt}}"))
                            if (file != "") {
                                file_ext = tools::file_ext(file)
                                if (file_ext == "") 
                                    file = paste(file, "txt", sep = ".")
                                write.table(located_xyIDs, file = file, 
                                            sep = "\t", quote = FALSE, row.names = FALSE, 
                                            col.names = TRUE)
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
                    x_yP1 <- (xyP[, 1] - min_x) * pixel_size[1]/(max_x - 
                                                                     min_x)
                    x_yP2 <- (xyP[, 2] - min_y) * pixel_size[2]/(max_y - 
                                                                     min_y)
                    xyP <- cbind(x_yP1, x_yP2)
                    xyP_xy_dist <- apply(xyP[-1, ], 1, itn_twoVectorDistance, 
                                         y = xyP[1, ])
                    xyN <- c(1:length(x))[xyP_xy_dist == min(xyP_xy_dist)]
                    if (min(xyP_xy_dist) >= pixel_dist) {
                        draw_plot <- FALSE
                    }
                }
                if (locate_choice > 1 || min(xyP_xy_dist) < pixel_dist) {
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