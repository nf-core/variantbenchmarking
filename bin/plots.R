#!/usr/bin/env Rscript

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

# Function to generate plots
generate_plots <- function(table, benchmark, type, filter) {
    # Melt the data for easier plotting

    if (type != "None" | filter != "None" ){
        table = table[table$Type == type & table$Filter == filter, ]
    }
    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c("TP_base", "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("Precision", "Recall", "F1"), ]
    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_base", "TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("Precision", "Recall", "F1"))

    if (type != "None" | filter != "None" ){
        title = paste("Type=",type," Filter=",filter, " | TP/FP/FN by tool", sep="")
    }else{
        title = paste("TP/FP/FN by tool", sep="")
    }

    # Visualize TP_base, TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool, group = interaction(variable, Tool))) +
        geom_line() +
        geom_point() +
        labs(title = title, x = "Tool", y = "Value", color = "Tool") +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(legend.position = "right", panel.background = element_rect(fill = "white"))

    if (type != "None" | filter != "None" ){
        title = paste("Type=",type," Filter=",filter, " | Precision, Recall, and F1 by Tool", sep="")
    }else{
        title = paste("Precision, Recall, and F1 by Tool", sep="")
    }
    # Visualize Precision, Recall, and F1 in separate plots with white background
    metric_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = Tool, shape = variable,  linetype = variable, group = interaction(variable, Tool))) +
        geom_point() +
        labs(title = title, x = "Tool", y = "Value", color = "Metric", linetype = "Metric") +
        theme_minimal() +
        theme(legend.position = "right", panel.background = element_rect(fill = "white"))

    # Save the plots
    if (type != "None" | filter != "None" ){
        name = paste(type, "_", filter, "_metric_by_tool_", benchmark, ".png", sep = "")
    }else{
        name = paste("metric_by_tool_", benchmark, ".png", sep = "")
    }
    tryCatch({
        if (!is.null(metric_plot)) {
            png(name, width = 10, height = 6, units = "in", res = 300, type = "cairo")
            print(metric_plot)
            dev.off()
        } else {
            png(name, width = 10, height = 6, units = "in", res = 300, type = "cairo")
            dev.off()
        }
    }, error = function(e) {
        message("Error occurred while saving metric plot: ", conditionMessage(e))
    })
    if (type != "None" | filter != "None" ){
        name = paste(type, "_", filter, "_variants_by_tool_", benchmark, ".png", sep = "")
    }else{
        name = paste("variants_by_tool_", benchmark, ".png", sep = "")
    }
    tryCatch({
        if (!is.null(tp_plot)) {
            png(name, width = 10, height = 6, units = "in", res = 300, type = "cairo")
            print(tp_plot)
            dev.off()
        } else {
            png(name, width = 10, height = 6, units = "in", res = 300, type = "cairo")
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

if (benchmark == "happy"){
    generate_plots(table, benchmark, "SNP", "PASS")
    generate_plots(table, benchmark, "SNP", "ALL")
    generate_plots(table, benchmark, "INDEL", "PASS")
    generate_plots(table, benchmark, "INDEL", "ALL")
}
else{
    generate_plots(table, benchmark, "None", "None")
}
