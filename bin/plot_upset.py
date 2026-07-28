#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

"""
Variant Upset Plot Generator
Command-line tool to create upset plots from FP/FN/TP results
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
import sys
import os


def parse_file(file_path, category_name):
    """
    Parse a VCF-like CSV file, extract tool detection information,
    and return a dictionary of sets.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return {}

    try:
        df = pd.read_csv(file_path, sep=',')
    except pd.errors.EmptyDataError:
        print(f"Warning: File {file_path} is empty.")
        return {}
    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        return {}

    gt_cols = [col for col in df.columns if '_GT' in col]

    if not gt_cols:
        print(f"Warning: No GT columns detected in {file_path}.")
        return {}

    contents = {}
    non_variant_genotypes = {'0/0', '0|0', './.', '.|.', '0'}

    for col in gt_cols:
        name_part = col.split('_GT')[0]
        tool_name = name_part.split('.')[0]

        upset_key = f"{tool_name}_{category_name}"

        required_cols = ['CHROM', 'POS', 'REF', 'ALT']
        if not all(c in df.columns for c in required_cols):
            print(f"Warning: Missing required columns in {file_path}")
            continue

        gt_mask = ~df[col].isin(non_variant_genotypes)

        variant_ids = set(tuple(row) for row in df.loc[gt_mask, required_cols].astype(str).values)

        if variant_ids:
            contents[upset_key] = variant_ids

    print(f"File: {file_path}, Detected tools: {[k.split('_')[0] for k in contents.keys()]}")
    total_variants = sum(len(s) for s in contents.values())
    print(f"Found {total_variants} total variants in {file_path} after filtering.")

    return {k: v for k, v in contents.items() if v}

def create_grouped_upset_plots(all_data, prefix=None, title="Variant Detection Upset Plots"):
    """
    Create two separate upset plots and save them to different files.
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

    def plot_group(plot_data, suffix, group_title):
        if not plot_data or not any(plot_data.values()):
            print(f"Warning: No data for {group_title} plot. Skipping.")
            return

        if len(plot_data) < 2:
            print(f"Warning: Only one tool found in {group_title} group. Skipping Upset plot.")
            return

        upset_data = from_contents(plot_data)

        if not isinstance(upset_data, pd.DataFrame) or upset_data.empty:
            print(f"Warning: Invalid or empty data returned for {group_title}. Skipping.")
            return

        num_categories = len(upset_data.index.get_level_values(0).unique())
        fig_width = max(10, num_categories * 0.8)
        fig_height = 8

        fig = plt.figure(figsize=(fig_width, fig_height))

        upset = UpSet(upset_data,
                      subset_size='count',
                      show_counts=True,
                      sort_by='cardinality',
                      min_subset_size=1)

        upset.plot(fig=fig)
        plt.tight_layout()

        plt.suptitle(f"{title}\n{group_title}", fontsize=16, fontweight='bold', y=1.02)

        plt.gca().tick_params(axis='y', labelsize=10)

        if prefix:
            output_file = f"{prefix}_{suffix}_mqc.png"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        plt.close(fig)

    # Group 1: TP_Base and FN
    group1_data = {k: v for k, v in master_data.items() if 'TP_Base' in k or 'FN' in k}
    if group1_data and any(group1_data.values()):
        plot_group(group1_data, 'tp_fn', 'TP_Base + FN')
    else:
        print("Warning: No TP_Base and/or FN variants found. Skipping TP_Base + FN plot.")

    # Group 2: TP_Comp and FP
    group2_data = {k: v for k, v in master_data.items() if 'TP_Comp' in k or 'FP' in k}
    if group2_data and any(group2_data.values()):
        plot_group(group2_data, 'tp_fp', 'TP_Comp + FP')
    else:
        print("Warning: No TP_Comp and/or FP variants found. Skipping TP_Comp + FP plot.")

def main():
    parser = argparse.ArgumentParser(description='Create upset plots from variant benchmarking results')
    parser.add_argument('--fp', help='False positives VCF file')
    parser.add_argument('--fn', help='False negatives VCF file')
    parser.add_argument('--tp-base', help='True positives base VCF file')
    parser.add_argument('--tp-comp', help='True positives comparison VCF file')
    parser.add_argument('--output', '-o', help='Output plot file (optional)')
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
