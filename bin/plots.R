#!/usr/bin/env Rscript

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(scales))

# Function to generate plots
generate_plots <- function(table, benchmark, type, filter, stats) {
    if (type != "None" && filter != "None") {
        table <- table[table$Type == type & table$Filter == filter, ]
        name1 <- paste(type, "_", filter, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(type, "_", filter, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(type, "_", filter, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else if (stats != "None") {
        table <- table[table$StatsType == stats, ]
        name1 <- paste(stats, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(stats, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(stats, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else if (type != "None") {
        table <- table[table$Type == type, ]
        name1 <- paste(type, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(type, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(type, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else {
        name1 <- paste("f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste("variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste("pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c( "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("F1"), ]
    metric_data$value <- as.numeric(as.character(metric_data$value))
    tp_data$value <- as.numeric(as.character(tp_data$value))

    print(tp_data)
    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("F1"))

    axis_text_size <- 12
    axis_title_size <- 14
    point_size <- 3
    line_size <- 1.2
    title_size <- 16
    facet_text_size <- 14
    legend_text_size <- 12
    legend_title_size <- 14

    # Visualize TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool, group = interaction(variable, Tool))) +
        geom_line(size = line_size) +
        geom_point(size = point_size) +
        labs(x = "Tool", y = "Value", color = "Tool", title = "Variant Comparison Metrics") +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            strip.text = element_text(size = facet_text_size, face = "bold")) +
        scale_y_continuous(labels = scales::label_number(), limits = c(0, NA))

    # Visualize f1
    f1_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = Tool)) +
        geom_point(size = point_size) +
        labs(x = "Tool", y = "F1 Score", title = "F1 Score") +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
            ) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0, 1)) +
            scale_color_discrete()

    # Visualize Precision vs Recall
    pr_plot <- ggplot(table) +
        geom_point(aes(x = Recall, y = Precision, color = Tool), size = point_size) +
        labs(x = "Recall", y = "Precision", title = "Precision vs Recall") +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size)) +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1))

    # Save the plots
    tryCatch({
        if (!is.null(f1_plot)) {
            ggsave(name1, f1_plot, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving metric plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(tp_plot)) {
            ggsave(name2, tp_plot, width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving TP plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(pr_plot)) {
            ggsave(name3, pr_plot, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving precision recall plot: ", conditionMessage(e))
    })
}

# Main script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please input the .summary file and benchmark", call. = FALSE)
}

table <- read.csv(args[1])
benchmark <- args[2]

if (benchmark == "happy") {
    generate_plots(table, benchmark, "SNP", "PASS", "None")
    generate_plots(table, benchmark, "SNP", "ALL", "None")
    generate_plots(table, benchmark, "INDEL", "PASS", "None")
    generate_plots(table, benchmark, "INDEL", "ALL", "None")
}else if (benchmark == "wittyer") {
    generate_plots(table, benchmark, "None", "None", "Base")
    generate_plots(table, benchmark, "None", "None", "Event")
}else if (benchmark == "concordance") {
    generate_plots(table, benchmark, "SNP", "None", "None")
    generate_plots(table, benchmark, "INDEL", "None", "None")
}else {
    if (benchmark == "rtgtools") {
        table <- table[table$Threshold == "None", ]
    }
    generate_plots(table, benchmark, "None", "None", "None")
}
