#!/usr/bin/env Rscript

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(scales))

seaborn_cb_palette <- c(
    "#0173b2", "#d55e00", "#029e73", "#de8f05", "#cc78bc",
    "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"
)

# Dictionary to translate raw pipeline tags into presentation-ready titles
term_mapping <- c(
    "GT"         = "Genotype Matching",
    "BASEPAIR"   = "Sequence Matching",
    "Snv"        = "SNVs",
    "SNP"        = "SNVs",
    "INDEL"      = "Indels",
    "Insertion"  = "Insertions",
    "Deletion"   = "Deletions",
    "JointIndel" = "Joint Indels",
    "PASS"       = "PASS Filter Only"
)

# Function to generate plots 
generate_plots <- function(table, benchmark, file_prefix, title_suffix, show_labels) {
    
    # Construct filenames based on the cleaned prefix
    if (file_prefix != "" && file_prefix != "overall") {
        name1 <- paste(file_prefix, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(file_prefix, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(file_prefix, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else {
        name1 <- paste("f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste("variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste("pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    }
    
    # Construct clean titles
    title_tp <- ifelse(title_suffix == "Overall", "Variant Comparison Metrics", paste("Variant Comparison Metrics -", title_suffix))
    title_f1 <- ifelse(title_suffix == "Overall", "F1 Score", paste("F1 Score -", title_suffix))
    title_pr <- ifelse(title_suffix == "Overall", "Precision vs Recall", paste("Precision vs Recall -", title_suffix))

    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c( "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("F1"), ]
    metric_data$value <- as.numeric(as.character(metric_data$value))
    tp_data$value <- as.numeric(as.character(tp_data$value))

    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("F1"))

    # Sizes
    axis_text_size <- 16
    axis_title_size <- 16
    point_size <- 4
    line_size <- 1.2
    title_size <- 18
    facet_text_size <- 16
    legend_text_size <- 14
    legend_title_size <- 16
    label_text_size <- 5

    # Visualize TP_comp, FP, and FN in separate plots
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool, group = interaction(variable, Tool))) +
        geom_line(size = line_size) +
        geom_point(size = point_size) +
        labs(x = "Tool", y = "Value", color = "Tool", title = title_tp) +
        facet_wrap(~variable, scales = "free_y") +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5, size = axis_text_size, face = "bold"),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size, face = "bold"),
            axis.title.y = element_text(size = axis_title_size, face = "bold"),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            strip.text = element_text(size = facet_text_size, face = "bold")) +
        scale_y_continuous(labels = scales::label_number(), limits = c(0, NA), expand = expansion(mult = c(0.05, 0.15))) +
        scale_color_manual(values = seaborn_cb_palette)

    if (show_labels) {
        tp_plot <- tp_plot + geom_text(aes(label = value), vjust = -1.0, size = label_text_size, color = "black")
    }

    # Visualize f1
    f1_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = Tool)) +
        geom_point(size = point_size) +
        labs(x = "Tool", y = "F1 Score", title = title_f1) +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 30, hjust = 0.5, size = axis_text_size, face = "bold"),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size, face = "bold"),
            axis.title.y = element_text(size = axis_title_size, face = "bold"),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
            ) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0, 1)) +
            scale_color_manual(values = seaborn_cb_palette)

    if (show_labels) {
        f1_plot <- f1_plot + geom_text(aes(label = round(value, 3)), vjust = -1.0, size = label_text_size, color = "black")
    }

    # Visualize Precision vs Recall
    pr_plot <- ggplot(table) +
        geom_point(aes(x = Recall, y = Precision, color = Tool), size = point_size) +
        labs(x = "Recall", y = "Precision", title = title_pr) +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size, face = "bold"),
            axis.title.y = element_text(size = axis_title_size, face = "bold"),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size)) +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_color_manual(values = seaborn_cb_palette)

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
show_point_labels <- "--labels" %in% args
args <- args[args != "--labels"]

if (length(args) < 2) {
    stop("Please input the .summary file and benchmark", call. = FALSE)
}

table <- read.csv(args[1])
benchmark <- args[2]

# Apply any tool specific pre filtering
if (benchmark == "rtgtools" && "Threshold" %in% colnames(table)) {
    table <- table[table$Threshold == "None", ]
}

# Define the columns we want to split our plots by
possible_groupings <- c("Type", "Filter", "StatsType", "Comparison", "Region")

# Find which of those columns actually exist in the current tool's table
split_cols <- intersect(colnames(table), possible_groupings)

# Automatically split and plot
if (length(split_cols) > 0) {
    
    groups <- interaction(table[, split_cols, drop = FALSE], drop = TRUE, sep = "_")
    split_data <- split(table, groups)
    
    for (raw_prefix in names(split_data)) {
        sub_table <- split_data[[raw_prefix]]
        
        # Split the string to filter out filler words
        parts <- unlist(strsplit(raw_prefix, "_"))
        
        # Remove case insensitive matches for all and none
        valid_parts <- parts[!tolower(parts) %in% c("all", "none")]
        
        if (length(valid_parts) > 0) {
            # File prefix remains raw and filesystem-safe
            clean_prefix <- paste(valid_parts, collapse = "_")
            
            # Create a pretty title by checking the dictionary mapping
            pretty_parts <- sapply(valid_parts, function(x) {
                if (x %in% names(term_mapping)) {
                    return(term_mapping[[x]])
                } else {
                    return(x) # Fallback to original text if not in dictionary
                }
            })
            
            # Join the pretty parts with a pipe or dash for clean readability
            title_suffix <- paste(pretty_parts, collapse = " | ")
        } else {
            clean_prefix <- "overall"
            title_suffix <- "Overall"
        }
        
        # Remove illegal filename characters just in case
        clean_prefix <- gsub("[^A-Za-z0-9_]", "_", clean_prefix)
        
        generate_plots(sub_table, benchmark, clean_prefix, title_suffix, show_point_labels)
    }
} else {
    generate_plots(table, benchmark, "overall", "Overall", show_point_labels)
}