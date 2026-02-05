#!/usr/bin/env python
# A script to parse, analyze, and plot indel or structural variant lenght from VCF and CSV files.

# Author: Victor Perez

import matplotlib.pyplot as plt
import os
import argparse
import numpy as np
import gzip
import csv
import pandas as pd
import matplotlib.transforms as mtrans

def get_sample_name(vcf_file):
    base_name = os.path.basename(vcf_file)
    return base_name.split('.')[0]

def parse_vcf_data(vcf_files):
    """
    Parses multiple VCF files, extracts SVLEN values, and determines the maximum SV length.

    Returns:
        tuple: A tuple containing a dictionary of SV data and the maximum SV length found.
    """
    svlen_data = {}
    max_svlen = 0

    for vcf_file in vcf_files:
        print(f"Processing VCF file: {vcf_file}")

        file_name = get_sample_name(vcf_file)
        svlen_data[file_name] = {'positive': [], 'negative': []}

        try:
            if vcf_file.endswith('.gz'):
                file_handle = gzip.open(vcf_file, 'rt')
            else:
                file_handle = open(vcf_file, 'r')

            with file_handle as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    info_column = parts[7]
                    svlen_val = None
                    sv_type = None

                    if 'SVLEN=' in info_column:
                        info_fields = info_column.split(';')
                        for field in info_fields:
                            if field.startswith('SVLEN='):
                                svlen_str = field.split('=')[1]
                                try:
                                    if ',' in svlen_str:
                                        svlen_str = svlen_str.split(',')[0]
                                    svlen_val = int(svlen_str)
                                except ValueError:
                                    print(f"Could not parse SVLEN value: {svlen_str} in {vcf_file}")

                    if svlen_val is None and 'SVTYPE=' in info_column:
                        info_fields = info_column.split(';')
                        for field in info_fields:
                            if field.startswith('SVTYPE='):
                                sv_type = field.split('=')[1]

                        if sv_type == 'DEL':
                            if 'END=' in info_column:
                                for field in info_fields:
                                    if field.startswith('END='):
                                        end_str = field.split('=')[1]
                                        try:
                                            end_pos = int(end_str)
                                            svlen_val = end_pos - int(parts[1])
                                        except ValueError:
                                            print(f"Could not parse END value: {end_str} in {vcf_file}")
                                            svlen_val = None

                            if svlen_val is None:
                                svlen_val = len(parts[4]) - len(parts[3])

                        elif sv_type == 'INS':
                            svlen_val = len(parts[4]) - len(parts[3])

                    if svlen_val is None:
                        ref_len = len(parts[3])
                        alt_len = len(parts[4])
                        if ref_len > alt_len:
                            svlen_val = -(ref_len - alt_len)
                        elif alt_len > ref_len:
                            svlen_val = alt_len - ref_len

                    if svlen_val is not None and svlen_val != 0:
                        abs_svlen = abs(svlen_val)
                        if abs_svlen > max_svlen:
                            max_svlen = abs_svlen
                        if svlen_val > 0:
                            svlen_data[file_name]['positive'].append(svlen_val)
                        else:
                            svlen_data[file_name]['negative'].append(svlen_val)

        except FileNotFoundError:
            print(f"Error: File not found at {vcf_file}")
            continue

    return svlen_data, max_svlen

def parse_csv_data(csv_files):
    """
    Parses multiple CSV files, extracts variant lengths from REF/ALT, and determines the maximum length.

    Returns:
        tuple: A tuple containing a dictionary of variant data and the maximum variant length found.
    """
    svlen_data = {}
    max_svlen = 0

    for csv_file in csv_files:
        print(f"Processing CSV file: {csv_file}")

        file_name = get_sample_name(csv_file)
        svlen_data[file_name] = {'positive': [], 'negative': []}

        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    svlen_val = None

                    ref = row.get('REF', '')
                    alt = row.get('ALT', '')

                    if ref and alt:
                        if len(alt) > len(ref):
                            svlen_val = len(alt) - len(ref)
                        else:
                            svlen_val = -(len(ref) - len(alt))

                    if svlen_val is None:
                        ref_truth = row.get('REF.truth', '')
                        alt_truth = row.get('ALT.truth', '')
                        if ref_truth and alt_truth:
                            if len(alt_truth) > len(ref_truth):
                                svlen_val = len(alt_truth) - len(ref_truth)
                            else:
                                svlen_val = -(len(ref_truth) - len(alt_truth))

                    if svlen_val is not None and svlen_val != 0:
                        abs_svlen = abs(svlen_val)
                        if abs_svlen > max_svlen:
                            max_svlen = abs_svlen
                        if svlen_val > 0:
                            svlen_data[file_name]['positive'].append(svlen_val)
                        else:
                            svlen_data[file_name]['negative'].append(svlen_val)

        except FileNotFoundError:
            print(f"Error: File not found at {csv_file}")
            continue

    return svlen_data, max_svlen

def human_format(num):
    """Converts a number into a readable string with units (bp, kb, Mb)."""
    val = abs(float(num))
    sign = "-" if num < 0 else ""
    if val == 0: return "0"
    if val < 1000: return f"{sign}{int(val)}bp"
    if val < 1000000: return f"{sign}{val/1000:.1f}kb"
    return f"{sign}{val/1000000:.1f}Mb"

def format_bp_label(bin_edges, half_open=True, format="sci", decimals=1):
    left = human_format(bin_edges[0])
    right = human_format(bin_edges[1])
    symbol = ")" if half_open else "]"
    return f"[{left}, {right}{symbol}"

def default_bins():
    """
    Creates a set of bins based on typical orders of magnitue(kilo,Mega,Giga,Tera)
    Returns:
        A list of bins as described above.

    """
    exps=[1,2,3,4,5,6,7,8]#bp,kbp,Mbp,Gbp,Tbp
    bins_pos = [10**i for i in exps]
    bins_neg = [-val for val in bins_pos]
    bins = sorted(bins_pos + bins_neg + [0])
    return bins

def filter_frame(df):
    """
    Filters a data frame based on the bin_label categories, if a category is empty for all
    samples it will be filtered out.
    Returns:
        A filtered pandas data frame.

    """
    group = df.groupby("bin_label", sort=False)["counts"]
    categ = group.sum().index
    idx_non_empty = group.sum().values > 0
    categ_filt = list(categ[idx_non_empty])
    df_upt = df.loc[df["bin_label"].isin(categ_filt)]
    return df_upt

def data2frame(sv_data, bin_edges="default"):
    """
    Takes the sv_data and calculates the histogram of the counts of insertions and deletions
    using the list of bin_edges.
    Returns:
        A pandas dataframe with the information organized for casting a bar plot.

    """
    data_table = {
        "sample": [],
        "bin_left_edge": [],
        "bin_right_edge": [],
        "bin_label": [],
        "type": [],
        "counts": []
    }

    if bin_edges == "default":
        bins = default_bins()
    elif isinstance(bin_edges, (list, np.ndarray)):
        bins = sorted(bin_edges)

    bin_left_edge = bins[:-1]
    bin_right_edge = bins[1::]
    bin_intervals = list(zip(bin_left_edge, bin_right_edge))
    """
    Note: In the following 2 lines the label of the intervals
    is formatted to comply with the np.histogram output, i.e.
    all bins are considered half-open [) except the right-most bin which is closed
    [].  Chek np.histogram documentation.
    """
    bin_labels = list(map(format_bp_label, bin_intervals[0:-1]))
    bin_labels.append(format_bp_label(bin_intervals[-1], half_open=False))
    mod_type = ["insertion" if l >= 0 else "deletion" for l, r in bin_intervals]

    for file_name, lengths_dict in sv_data.items():
        all_lengths = lengths_dict['positive'] + lengths_dict['negative']
        if all_lengths:
            counts, _ = np.histogram(all_lengths, bins=bins)
            data_table["sample"].extend([file_name] * len(counts))
            data_table["bin_left_edge"].extend(bin_left_edge)
            data_table["bin_right_edge"].extend(bin_right_edge)
            data_table["bin_label"].extend(bin_labels)
            data_table["type"].extend(mod_type)
            data_table["counts"].extend(counts)

    df_table = pd.DataFrame(data_table)
    return df_table

def plot_svlen_distributions(sv_data, bins, output_file, plot_title, show_labels=True):
    """
    CreateS a bar plot and writes it in the output_file path
    Returns:
        A .png figure saved in the output_file path.

    """
    df_table = data2frame(sv_data, bin_edges=bins)
    df_upt = filter_frame(df_table)

    if df_upt.empty:
        print("Warning: No data available to plot after filtering (counts are 0). Generating empty placeholder plot.")
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No Data Available", ha='center', va='center', fontsize=20)
        plt.title(plot_title)
        plt.savefig(output_file)
        return

    category_label = df_upt["bin_label"].unique()
    files = sorted(df_upt["sample"].unique())
    bar_height = {sample: df_upt.loc[df_upt["sample"] == sample, "counts"].values for sample in files}

    types_list = df_upt["type"].unique()

    width = 0.15
    group_spacing = 1.3
    x = np.arange(len(category_label)) * group_spacing

    xlabel_text = "Variant Length Range" 
    xlabel_future_action = None
    xlabel_position = 0

    if len(types_list) == 2:
        xlabel_text = "Deletions | Insertions"
        xlabel_future_action = "reposition"
        ref_file = df_upt["sample"].unique()[0]
        aux = df_upt.loc[df_upt["sample"] == ref_file]
        transition_index = (aux["type"].values != "deletion").argmax()
        xlabel_position = x[transition_index]

    elif len(types_list) == 1:
        xlabel_future_action = None
        if types_list[0] == "insertion":
            xlabel_text = "Insertions"
        elif types_list[0] == "deletion":
            xlabel_text = "Deletions"

    plt.style.use('seaborn-v0_8-colorblind')
    fig, ax = plt.subplots(figsize=(26, 16))
    
    multiplier = 0
    for attribute, measurement in bar_height.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute, alpha=1)
        
        if show_labels:
            ax.bar_label(rects, 
                         padding=10, 
                         rotation=45, 
                         fontsize=24, 
                         fontweight='bold',
                         clip_on=False)
        multiplier += 1

    ax.set_title(plot_title, fontsize=42, fontweight='bold', pad=60)
    ax.grid(True, which="major", linestyle='--', alpha=0.6)
    
    tick_pos = x + (width * (len(files)-1)/2)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(category_label, rotation=40, fontsize=32, ha='right')
    
    ax.set_xlabel(xlabel_text, fontsize=36, fontweight='bold')
    if xlabel_future_action == "reposition":
        trans = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
        ax.xaxis.set_label_coords(xlabel_position, -0.3, transform=trans)

    ax.set_yscale('log')
    headroom_factor = 20 if show_labels else 5
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * headroom_factor)
    
    ax.tick_params(axis='y', which='major', labelsize=32)
    ax.set_ylabel("Count (Log10)", fontsize=36, fontweight='bold')

    ax.legend(title="Tool",
          loc='upper left',
          title_fontsize=34,
          fancybox=True,
          shadow=True,
          bbox_to_anchor=(1, 1),
          prop={'size': 32}
          )
    ax.set_facecolor((0.95, 0.95, 0.95))
    
    fig.savefig(output_file, dpi=100, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot SVLEN distributions from one or more VCF or CSV files.")
    parser.add_argument('input_files', nargs='+', help="One or more VCF or CSV files to process.")
    parser.add_argument('--output', '-o', dest='output_file', default='svlen_distributions.png')
    parser.add_argument('--title', '-t', dest='plot_title', default='Structural Variant Length Distributions by Type')
    parser.add_argument('--bins', '-b', nargs='+', type=str, default=None, 
                        help="""
                        List of integer numbers representing the bin_edges for the plot.
                        If no list is given the default edges are bp,kbp,Mbp,Gbp,Tbp in both
                        negative(Deletions) and positive(Insertions) directions. The counts of
                        the histogram are carried out by numpy.histogram.
                        """) 

    parser.add_argument('--no-labels', action='store_true', help="Hide the count labels on top of the bars.")

    args = parser.parse_args()

    final_bins = None
    if args.bins:
        final_bins = []
        for b in args.bins:
            clean_bins = b.replace(',', ' ').split()
            final_bins.extend([int(x) for x in clean_bins])
        final_bins = sorted(list(set(final_bins)))

    vcf_files_to_parse = [f for f in args.input_files if f.endswith(('.vcf', '.vcf.gz'))]
    csv_files_to_parse = [f for f in args.input_files if f.endswith('.csv')]

    sv_data = {}
    max_len = 0

    if vcf_files_to_parse:
        vcf_data, vcf_max_len = parse_vcf_data(vcf_files_to_parse)
        sv_data.update(vcf_data)
        if vcf_max_len > max_len:
            max_len = vcf_max_len

    if csv_files_to_parse:
        csv_data, csv_max_len = parse_csv_data(csv_files_to_parse)
        sv_data.update(csv_data)
        if csv_max_len > max_len:
            max_len = csv_max_len

    if sv_data:
        bin_edges = final_bins if final_bins else "default"
        plot_svlen_distributions(sv_data, bin_edges, args.output_file, args.plot_title, show_labels=not args.no_labels)
    else:
        print("No valid input files found to plot.")