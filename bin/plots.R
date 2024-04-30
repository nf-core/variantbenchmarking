#!/usr/bin/env Rscript

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

# Function to generate plots
generate_plots <- function(table, benchmark) {
    # Melt the data for easier plotting
    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c("TP_base", "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("Precision", "Recall", "F1"), ]
    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_base", "TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("Precision", "Recall", "F1"))

    # Visualize TP_base, TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = variable, group = interaction(variable, Tool))) +
        geom_line() +
        geom_point() +
        labs(title = "TP_base, TP_comp, FP, and FN by Tool", x = "Tool", y = "Value", color = "Metric") +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(legend.position = "top", panel.background = element_rect(fill = "white"))

    # Visualize Precision, Recall, and F1 in separate plots with white background
    metric_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = variable, linetype = variable, group = interaction(variable, Tool))) +
        geom_line() +
        geom_point() +
        labs(title = "Precision, Recall, and F1 by Tool", x = "Tool", y = "Value", color = "Metric", linetype = "Metric") +
        theme_minimal() +
        theme(legend.position = "top", panel.background = element_rect(fill = "white"))

    # Save the plots
    tryCatch({
        if (!is.null(metric_plot)) {
            png(paste("metric_by_tool_", benchmark, ".png", sep = ""), width = 10, height = 6, units = "in", res = 300, type = "cairo")
            print(metric_plot)
            dev.off()
        } else {
            png(paste("metric_by_tool_", benchmark, ".png", sep = ""), width = 10, height = 6, units = "in", res = 300, type = "cairo")
            dev.off()
        }
    }, error = function(e) {
        message("Error occurred while saving metric plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(tp_plot)) {
            png(paste("variants_by_tool_", benchmark, ".png", sep = ""), width = 16, height = 16, units = "in", res = 300, type = "cairo")
            print(tp_plot)
            dev.off()
        } else {
            png(paste("variants_by_tool_", benchmark, ".png", sep = ""), width = 16, height = 16, units = "in", res = 300, type = "cairo")
            dev.off()
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

generate_plots(table, benchmark)
