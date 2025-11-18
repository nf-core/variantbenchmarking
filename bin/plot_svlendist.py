#!/usr/bin/env python
# A script to parse, analyze, and plot indel or structural variant lenght from VCF and CSV files.

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

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


def format_bp_label(bin_edges,half_open=True,format="sci",decimals=1):
    """
    Creates a string using the edges of the bins.
    Returns:
        A string of the form [lef_edge,right_edge]

    """
    
    if format=="sci":
        format_type="E"
    elif format=="raw":
        format_type="f"
    
    if half_open:
        symbol=")"
    else:
        symbol="]"
    
    label="[{left},{right}{S}".format(left= f"{bin_edges[0]:.{decimals}{format_type}}",
                                    right=f"{bin_edges[1]:.{decimals}{format_type}}",
                                    S=symbol
                                    )
    return label

def default_bins():
    """
    Creates a set of bins based on typical orders of magnitue(kilo,Mega,Giga,Tera)
    Returns:
        A list of bins as described above.

    """
    exps=[1,2,3,6,9,12]#bp,kbp,Mbp,Gbp,Tbp
    bins_pos=[10**i for i in exps]
    bins_neg=[ -val for val in bins_pos ]
    bins=sorted(bins_pos+bins_neg+[0])
    return bins

def filter_frame(df):
    """
    Filters a data frame based on the bin_label categories, if a category is empty for all
    samples it will be filtered out.
    Returns:
        A filtered pandas data frame.

    """
    group=df.groupby("bin_label",sort=False)["counts"]
    categ=group.sum().index
    idx_non_empty=group.sum().values>0
    categ_filt=list(categ[idx_non_empty])
    df_upt=df.loc[df["bin_label"].isin(categ_filt)]
    return df_upt


def data2frame(sv_data,bin_edges="default"):
    """
    Takes the sv_data and calculates the histogram of the counts of insertions and deletions
    using the list of bin_edges. 
    Returns:
        A pandas dataframe with the information organized for casting a bar plot.

    """

    data_table={
                "sample":[],
                "bin_left_edge":[],
                "bin_right_edge":[],
                "bin_label":[],
                "type":[],
                "counts":[]
                }


    if bin_edges=="default":
        bins=default_bins()
    elif isinstance(bin_edges, (list, np.ndarray)):
        bins=sorted(bin_edges)
        
    bin_left_edge=bins[:-1]
    bin_right_edge=bins[1::]
    bin_intervals=list(zip(bin_left_edge,bin_right_edge))

    """
    Note: In the following 2 lines the label of the intervals
    is formatted to comply with the np.histogram output, i.e.
    all bins are considered half-open [) except the right-most bin which is closed
    [].  Chek np.histogram documentation.
    """
    bin_labels=list ( map( format_bp_label, bin_intervals[0:-1] ) )
    bin_labels.append( format_bp_label(bin_intervals[-1],half_open=False) )
    mod_type=["insertion" if l>=0 else "deletion" for l,r in bin_intervals ]


    for file_name, lengths_dict in sv_data.items():

        all_lengths = lengths_dict['positive'] + lengths_dict['negative']

        if all_lengths:
            counts, _ = np.histogram(all_lengths, bins=bins)
            data_table["sample"].extend([file_name]*len(counts))
            data_table["bin_left_edge"].extend(bin_left_edge)
            data_table["bin_right_edge"].extend(bin_right_edge)
            data_table["bin_label"].extend(bin_labels)
            data_table["type"].extend(mod_type)
            data_table["counts"].extend(counts)
            
    df_table=pd.DataFrame(data_table) 
    return df_table

def plot_svlen_distributions(sv_data, output_file, plot_title):
    """
    CreateS a bar plot and writes it in the output_file path
    Returns:
        A .png figure saved in the output_file path.

    """
    df_table=data2frame(sv_data,bin_edges="default")
    df_upt=filter_frame(df_table)
    
    category_label=df_upt["bin_label"].unique()
    files=df_upt["sample"].unique()
    bar_height={sample: df_upt.loc[ df_upt["sample"]==sample, "counts"].values for sample in files }

    types_list=df_upt["type"].unique()

    width = 0.11  # the width of the bars
    interval_capacity=1/width
    interval_load=len(category_label)

    if interval_load<interval_capacity:
        x = np.arange(len(category_label))  # the label locations
    else:
        x=(width*interval_load)*np.array(list(range(interval_load)))
    

    if len(types_list)==2:
        xlabel_text="Deletions | Insertions"
        xlabel_future_action="reposition"
        ref_file=df_upt["sample"].unique()[0]
        aux=df_upt.loc[df_upt["sample"]==ref_file]
        transition_index=(aux["type"].values!="deletion").argmax()
        xlabel_position=x[transition_index]
    
    elif len(types_list)==1:
        xlabel_future_action=None
        if types_list[0]=="insertion":
            xlabel_text="Insertions"
            
        elif types_list[0]=="deletion":
            xlabel_text="Deletions"
            
    
    plt.style.use('seaborn-v0_8-colorblind')
    fig, ax = plt.subplots(figsize=(14, 10),layout='constrained')
    multiplier = 0
    pads=[-1,2]
    for attribute, measurement in bar_height.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute,alpha=1)
        ax.bar_label(rects, padding=pads[multiplier%len(pads)])
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_title(plot_title,fontsize=18, fontweight='bold', pad=20)
    ax.grid(True, which="major", linestyle='--', alpha=0.6)
    ax.set_xticks(x + width, category_label,rotation=35,fontsize=14)
    ax.set_xlabel(xlabel_text, fontsize=16, fontweight='bold')
    if xlabel_future_action=="reposition":
        trans = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
        ax.xaxis.set_label_coords(xlabel_position, -0.3,transform=trans)

    ax.set_yscale('log')
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_ylabel("Count (Log10)", fontsize=16, fontweight='bold')

    ax.legend(title="Tool",
          loc='upper right',
          title_fontsize=16,
          fancybox=True,
          shadow=True,
          bbox_to_anchor=(1.4, 1),
          prop={'size': 14}
          )
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot SVLEN distributions from one or more VCF or CSV files.")
    parser.add_argument('input_files', nargs='+', help="One or more VCF or CSV files to process.")
    parser.add_argument('--output', '-o', dest='output_file', default='svlen_distributions.png',
                        help="The name of the output image file (default: svlen_distributions.png).")
    parser.add_argument('--title', '-t', dest='plot_title', default='Structural Variant Length Distributions by Type',
                        help="The title for the plot.")

    args = parser.parse_args()

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
        plot_svlen_distributions(sv_data, args.output_file, args.plot_title)
    else:
        print("No valid input files found to plot.")
