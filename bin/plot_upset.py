#!/usr/bin/env python3

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

"""
Variant Upset Plot Generator (Interactive Plotly HTML)
Command-line tool to create interactive upset plots from FP/FN/TP results using pure Plotly.
"""

import argparse
import pandas as pd
import sys
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict, Counter

def parse_file(file_path, category_name):
    """
    Parse a VCF-like CSV file, extract tool detection information,
    and return a dictionary of sets.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return {}

    try:
        # Read only the header first to identify columns
        header = pd.read_csv(file_path, sep=',', nrows=0)
        gt_cols = [col for col in header.columns if '_GT' in col]
        
        if not gt_cols:
            print(f"Warning: No GT columns detected in {file_path}.")
            return {}

        required_cols = ['CHROM', 'POS', 'REF', 'ALT']
        if not all(c in header.columns for c in required_cols):
            print(f"Warning: Missing required columns in {file_path}")
            return {}

        # Only load the columns we actually need to save RAM
        cols_to_use = required_cols + gt_cols
        
        # Use memory efficient data types
        dtypes = {
            'CHROM': 'category',
            'POS': 'int32',
            'REF': 'category',
            'ALT': 'category'
        }
        
        df = pd.read_csv(file_path, sep=',', usecols=cols_to_use, dtype=dtypes)
    except pd.errors.EmptyDataError:
        print(f"Warning: File {file_path} is empty.")
        return {}
    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        return {}

    contents = {}
    non_variant_genotypes = {'0/0', '0|0', './.', '.|.', '0'}

    # Vectorized ID generation: Reuses the same strings for all tools instead of creating millions of new tuples
    variant_ids_series = df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + ':' + df['REF'].astype(str) + ':' + df['ALT'].astype(str)

    for col in gt_cols:
        name_part = col.split('_GT')[0]
        tool_name = name_part.split('.')[0]

        upset_key = f"{tool_name}_{category_name}"

        # Vectorized masking directly on the series
        gt_mask = ~df[col].isin(non_variant_genotypes)
        variant_ids = set(variant_ids_series[gt_mask])

        if variant_ids:
            contents[upset_key] = variant_ids

    print(f"File: {file_path}, Detected tools: {[k.split('_')[0] for k in contents.keys()]}")
    total_variants = sum(len(s) for s in contents.values())
    print(f"Found {total_variants} total variants in {file_path} after filtering.")

    return {k: v for k, v in contents.items() if v}

def create_interactive_upset(plot_data, title, output_file):
    """
    Creates an interactive Plotly-based UpSet plot and saves it as a self-contained HTML file.
    """
    if not plot_data or not any(plot_data.values()):
        print(f"Warning: No data for {title} plot. Skipping.")
        return

    if len(plot_data) < 2:
        print(f"Warning: Only one tool found in {title} group. Skipping Upset plot.")
        return

    categories = list(plot_data.keys())

    # Tally intersections natively and clear sets from memory iteratively.
    variant_masks = defaultdict(int)
    for i, cat in enumerate(categories):
        for var in plot_data[cat]:
            variant_masks[var] |= (1 << i)
        # Explicitly free memory for this tool's variant set once it is tallied
        plot_data[cat].clear()
        
    mask_counts = Counter(variant_masks.values())
    del variant_masks # Free the tally dictionary

    # Build the aggregated dataframe manually
    rows = []
    for mask, count in mask_counts.items():
        row = {'intersection_size': count}
        for i, cat in enumerate(categories):
            row[cat] = bool(mask & (1 << i))
        rows.append(row)

    if not rows:
        print(f"Warning: Invalid or empty data returned for {title}. Skipping.")
        return

    df_upset = pd.DataFrame(rows)
    df_upset = df_upset[df_upset['intersection_size'] > 0].sort_values('intersection_size', ascending=False).reset_index(drop=True)

    num_intersections = len(df_upset)

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.12,
        row_heights=[0.55, 0.45]
    )

    x_indices = list(range(num_intersections))
    sizes = df_upset['intersection_size'].tolist()

    # Top bar chart for intersection sizes
    fig.add_trace(
        go.Bar(
            x=x_indices,
            y=sizes,
            marker_color='black',
            text=[f"<b>{s}</b>" for s in sizes],
            textposition='outside',
            showlegend=False,
            hovertemplate="Intersection Size: %{y}<extra></extra>"
        ),
        row=1, col=1
    )

    # Connecting lines for multi-category intersections
    for idx, row in df_upset.iterrows():
        x_idx = idx
        active_cats = [i for i, cat in enumerate(categories) if row[cat]]
        if len(active_cats) > 1:
            fig.add_trace(
                go.Scatter(
                    x=[x_idx, x_idx],
                    y=[min(active_cats), max(active_cats)],
                    mode='lines',
                    line=dict(color='black', width=3),
                    showlegend=False,
                    hoverinfo='skip'
                ),
                row=2, col=1
            )

    # Bottom dot matrix for category presence/absence
    for i, cat in enumerate(categories):
        y_vals = [i] * num_intersections
        marker_colors = []
        marker_symbols = []
        hover_texts = []

        for idx, row in df_upset.iterrows():
            present = row[cat]
            if present:
                marker_colors.append('black')
                marker_symbols.append('circle')
                hover_texts.append('Present')
            else:
                marker_colors.append('#dcdcdc')
                marker_symbols.append('circle')
                hover_texts.append('Absent')

        fig.add_trace(
            go.Scatter(
                x=x_indices,
                y=y_vals,
                mode='markers',
                marker=dict(size=12, color=marker_colors, symbol=marker_symbols),
                name=cat,
                text=hover_texts,
                hovertemplate=f"<b>{cat}</b>: %{{text}}<extra></extra>",
                showlegend=False
            ),
            row=2, col=1
        )

    fig.update_layout(
        template="plotly_white",
        font=dict(size=16, color="black"),
        height=750, # Explicitly taller figure to prevent row label crowding
        margin=dict(l=220, r=40, t=40, b=40),
        hovermode="x unified"
    )

    fig.update_xaxes(showticklabels=False, row=1, col=1)

    max_size = max(sizes) if sizes else 10
    fig.update_yaxes(
        title_text="<b>Intersection size</b>",
        tickfont=dict(size=16),
        title_font=dict(size=18),
        showgrid=True,
        range=[0, max_size * 1.3],
        row=1, col=1
    )

    fig.update_xaxes(showticklabels=False, row=2, col=1)
    fig.update_yaxes(
        tickvals=list(range(len(categories))),
        ticktext=[f"<b>{cat}</b>" for cat in categories],
        tickfont=dict(size=16),
        title_font=dict(size=18),
        autorange="reversed",
        showgrid=True,
        row=2, col=1
    )

    plot_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{title}</title>
        <style>
            body {{ font-family: -apple-system, sans-serif; background-color: #f8f9fa; margin: 20px; }}
            .plot-container {{ width: 100%; margin-bottom: 20px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); box-sizing: border-box; }}
            .html-title {{ font-size: 26px; font-weight: bold; color: #2c3e50; margin-bottom: 15px; text-align: center; word-wrap: break-word; }}
        </style>
    </head>
    <body>
        <div class="plot-container">
            <div class="html-title">{title}</div>
            {plot_html}
        </div>
    </body>
    </html>
    """

    with open(output_file, "w") as f:
        f.write(html_template)
    print(f"Plot saved to {output_file}")

def create_grouped_upset_plots(all_data, prefix=None, title="Variant Detection Upset Plots"):
    """
    Create two separate interactive upset plots and save them to HTML files.
    """
    if not all_data:
        print("Error: No data to plot.")
        return

    master_data = {}
    for data in all_data:
        master_data.update(data)

    if not master_data or not any(master_data.values()):
        print("Warning: Master data dictionary is empty. No plots generated.")
        return

    base_prefix = prefix if prefix else "upset_plot"

    # Group 1: TP_Base and FN
    group1_data = {k: v for k, v in master_data.items() if 'TP_Base' in k or 'FN' in k}
    if group1_data and any(group1_data.values()):
        out_html = f"{base_prefix}_tp_fn_mqc.html"
        create_interactive_upset(group1_data, f"{title} - TP_Base + FN", out_html)
    else:
        print("Warning: No TP_Base and/or FN variants found. Skipping TP_Base + FN plot.")

    # Group 2: TP_Comp and FP
    group2_data = {k: v for k, v in master_data.items() if 'TP_Comp' in k or 'FP' in k}
    if group2_data and any(group2_data.values()):
        out_html = f"{base_prefix}_tp_fp_mqc.html"
        create_interactive_upset(group2_data, f"{title} - TP_Comp + FP", out_html)
    else:
        print("Warning: No TP_Comp and/or FP variants found. Skipping TP_Comp + FP plot.")

def main():
    parser = argparse.ArgumentParser(description='Create interactive HTML upset plots from variant benchmarking results')
    parser.add_argument('--fp', help='False positives VCF file')
    parser.add_argument('--fn', help='False negatives VCF file')
    parser.add_argument('--tp-base', help='True positives base VCF file')
    parser.add_argument('--tp-comp', help='True positives comparison VCF file')
    parser.add_argument('--output', '-o', help='Output file prefix (optional)')
    parser.add_argument('--title', help='Plot title', default='Variant Detection Upset Plots')

    args = parser.parse_args()

    if not any([args.fp, args.fn, args.tp_base, args.tp_comp]):
        parser.print_help()
        sys.exit(1)

    all_data = []

    file_map = {
        args.fp: 'FP',
        args.fn: 'FN',
        args.tp_base: 'TP_Base',
        args.tp_comp: 'TP_Comp'
    }

    for file_path, category in file_map.items():
        if file_path:
            print(f"Processing {category} file: {file_path}")
            data = parse_file(file_path, category)
            if data and any(data.values()):
                all_data.append(data)

    if not all_data:
        print("Error: No valid data found in any input files")
        sys.exit(1)

    create_grouped_upset_plots(all_data, args.output, args.title)

if __name__ == "__main__":
    main()