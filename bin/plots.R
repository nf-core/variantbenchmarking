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

    if (type != "None" && filter != "None" ){
        table = table[table$Type == type & table$Filter == filter, ]
        title1 = paste("Type=",type," Filter=",filter, " | TP/FP/FN by tool", sep="")
        title2 = paste("Type=",type," Filter=",filter, " | Precision, Recall, and F1 by Tool", sep="")
        name1 = paste(type, "_", filter, "_metric_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 = paste(type, "_", filter, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    else if (stats != "None" ){
        table = table[table$StatsType == stats, ]
        title1 = paste("StatsType=",stats, " | TP/FP/FN by tool", sep="")
        title2 = paste("StatsType=",stats, " | Precision, Recall, and F1 by Tool", sep="")
        name1 = paste(stats, "_metric_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 = paste(stats, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    else{
        title1 = paste("TP/FP/FN by tool", sep="")
        title2 = paste("Precision, Recall, and F1 by Tool", sep="")
        name1 = paste("metric_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 = paste("variants_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c("TP_base", "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("Precision", "Recall", "F1"), ]
    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_base", "TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("Precision", "Recall", "F1"))

    # Visualize TP_base, TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool, group = interaction(variable, Tool))) +
        geom_line() +
        geom_point() +
        labs(title = title1, x = "Tool", y = "Value", color = "Tool") +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5))

    # Visualize Precision, Recall, and F1 in separate plots with white background
    metric_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = Tool, shape = variable,  linetype = variable, group = interaction(variable, Tool))) +
        geom_point() +
        labs(title = title2, x = "Tool", y = "Value", color = "Metric", linetype = "Metric") +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5))

    # Save the plots
    tryCatch({
        if (!is.null(metric_plot)) {
            ggsave(name1, metric_plot, width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        } else {
            ggsave(name1, metric_plot, width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving metric plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(tp_plot)) {
            ggsave(name2, tp_plot, width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        } else {
            ggsave(name2,tp_plot, width = 10, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving TP plot: ", conditionMessage(e))
    })
}

# Main script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please input the .summary file and benchmark", call. = FALSE)
}

table <- read.csv(args[1])
benchmark <- args[2]

if (benchmark == "happy"){
    generate_plots(table, benchmark, "SNP", "PASS", "None")
    generate_plots(table, benchmark, "SNP", "ALL", "None")
    generate_plots(table, benchmark, "INDEL", "PASS", "None")
    generate_plots(table, benchmark, "INDEL", "ALL", "None")
}else if (benchmark == "wittyer"){
    generate_plots(table, benchmark, "None", "None", "Base")
    generate_plots(table, benchmark, "None", "None", "Event")

}else{
    generate_plots(table, benchmark, "None", "None", "None")
}
