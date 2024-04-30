#!/usr/bin/env Rscript


# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please input the .summary file and benchmark", call. = FALSE)
}

table     <- read.csv(args[1])
benchmark <- args[2]

# Load necessary libraries
library(ggplot2)

# Melt the data for easier plotting
input_data_melted <- melt(table, id.vars = "Tool")

tp_data <- input_data_melted[input_data_melted$variable %in% c("TP_base", "TP_comp", "FP", "FN"), ]
metric_data <- input_data_melted[input_data_melted$variable %in% c("Precision", "Recall", "F1"), ]
# Specify the order of levels for the variable aesthetic
tp_data$variable <- factor(tp_data$variable, levels = c("TP_base", "TP_comp", "FP", "FN"))
metric_data$variable <- factor(metric_data$variable, levels = c("Precision", "Recall", "F1"))

# Visualize TP_base, TP_comp, FP, and FN in separate plots
tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = variable)) +
    geom_line() +
    geom_point() +  # Add points for clarity
    labs(title = "TP_base, TP_comp, FP, and FN by Tool", x = "Tool", y = "Value", color = "Metric") +
    theme_minimal() +
    theme(legend.position = "top", panel.background = element_rect(fill = "white"))


# Visualize Precision, Recall, and F1 in separate plots with white background
metric_plot <- ggplot(metric_data, aes(x = Tool, y = value, fill = variable)) +
geom_bar(stat = "identity") +
    labs(title = "Precision, Recall, and F1 by Tool", x = "Tool", y = "Value", fill = "Metric") +
    facet_wrap(~variable, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none", panel.background = element_rect(fill = "white"),
        # Adjust plot dimensions
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12))
# Adjust the dimensions of the plot
# Save the plots (adjust file names as needed)
ggsave(paste("metric_by_tool_",benchmark,".png", sep=""), metric_plot, width = 10, height = 6, units = "in", dpi = 300)
ggsave(paste("variants_by_tool_",benchmark,".png", sep=""), tp_plot, width = 16, height = 16, units = "in",dpi = 300)
