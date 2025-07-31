#!/usr/bin/env Rscript

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

# Function to generate plots
generate_plots <- function(table, benchmark, type, filter, stats) {
    # Melt the data for easier plotting
    ## where type and filter are both none

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
    } else {
        name1 <- paste("f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste("variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste("pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c("TP_base", "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("F1"), ]
    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_base", "TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("F1"))

    # Visualize TP_base, TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool, group = interaction(variable, Tool))) +
        geom_line() +
        geom_point() +
        labs(x = "Tool", y = "Value", color = "Tool") +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5))

    # Visualize f1
    f1_plot <- ggplot(metric_data, aes(x = Tool, y = value)) +
        geom_point() +
        labs(x = "Tool", y = "f1") +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5))

    # Visualize Precision vs Recall
    pr_plot <- ggplot(table) +
        geom_point(aes(x = Recall, y = Precision, color = Tool)) +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"))

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

}else {
    if (benchmark == "rtgtools") {
        table <- table[table$Threshold == "None", ]
    }
    generate_plots(table, benchmark, "None", "None", "None")
}
