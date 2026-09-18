#!/usr/bin/env python3
# A script to parse, analyze, and plot indel or structural variant length from VCF and CSV files.

# Author: Victor Perez
# Adapted for Interactive HTML/Plotly by Kuebra Narci

import os
import argparse
import numpy as np
import gzip
import csv
import pandas as pd
import plotly.express as px

SEABORN_CB_PALETTE = [
    "#0173b2", "#d55e00", "#029e73", "#de8f05", "#cc78bc",
    "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"
]

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

def format_bp_label(bin_edges, half_open=True):
    """Formats the bin edges into readable bracket labels."""
    left = human_format(bin_edges[0])
    right = human_format(bin_edges[1])
    symbol = ")" if half_open else "]"
    return f"[{left}, {right}{symbol}"

def default_bins():
    """
    Creates a set of bins based on typical orders of magnitude (bp, kbp, Mbp, etc.).
    Returns:
        A list of bins as described above.
    """
    exps = [1, 2, 3, 4, 5, 6, 7, 8]
    bins_pos = [10**i for i in exps]
    bins_neg = [-val for val in bins_pos]
    bins = sorted(bins_pos + bins_neg + [0])
    return bins

def filter_frame(df):
    """
    Filters a data frame based on the bin_label categories. If a category is empty for all
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
    Creates an interactive Plotly bar plot and writes it to the output file path as an HTML report.
    """
    df_table = data2frame(sv_data, bin_edges=bins)
    df_upt = filter_frame(df_table).copy()

    if df_upt.empty:
        print("Warning: No data available to plot after filtering (counts are 0). Generating empty placeholder plot.")
        html_template = f"""
        <!DOCTYPE html>
        <html><body>
        <h1 style="text-align: center; font-family: sans-serif; color: #2c3e50;">No Data Available for {plot_title}</h1>
        </body></html>
        """
        with open(output_file, "w") as f:
            f.write(html_template)
        return

    def format_label(val):
        if val == 0:
            return ""
        if val >= 1000:
            return f"<b>{val/1000:g}k</b>"
        return f"<b>{val}</b>"

    if show_labels:
        df_upt["formatted_counts"] = df_upt["counts"].apply(format_label)
    else:
        df_upt["formatted_counts"] = ""

    unique_tools = sorted(df_upt["sample"].unique())
    global_color_map = {tool: SEABORN_CB_PALETTE[i % len(SEABORN_CB_PALETTE)] for i, tool in enumerate(unique_tools)}

    types_list = df_upt["type"].unique()
    xlabel_text = "Variant Length Range"
    if len(types_list) == 2:
        xlabel_text = "Deletions | Insertions"
    elif len(types_list) == 1:
        if types_list[0] == "insertion":
            xlabel_text = "Insertions"
        elif types_list[0] == "deletion":
            xlabel_text = "Deletions"

    fig = px.bar(
        df_upt,
        x="bin_label",
        y="counts",
        color="sample",
        barmode="group",
        text="formatted_counts",
        color_discrete_map=global_color_map,
        labels={"counts": "Variant Count", "bin_label": "Length Range", "sample": "Tool"}
    )

    fig.update_layout(
        template="plotly_white",
        font=dict(size=18, color="black"),
        legend=dict(
            title_text="<b>Tool</b>",
            font=dict(size=18),
            title_font=dict(size=20)
        ),
        margin=dict(l=80, r=40, t=40, b=80),
        hovermode="x unified"
    )

    tick_labels = df_upt["bin_label"].unique()
    fig.update_xaxes(
        title_text=f"<b>{xlabel_text}</b>",
        tickfont=dict(size=18),
        showline=True, linewidth=1.5, linecolor='black', mirror=True,
        tickvals=tick_labels,
        ticktext=[f"<b>{t}</b>" for t in tick_labels]
    )

    max_val = df_upt["counts"].max()
    upper_lim = np.log10(max_val) + 1.5 if max_val > 0 else 1

    fig.update_yaxes(
        title_text="<b>Variant Count (Log Scale)</b>",
        tickfont=dict(size=18),
        showline=True, linewidth=1.5, linecolor='black', mirror=True,
        type="log",
        range=[np.log10(0.5), upper_lim],
        tickformat="~s"
    )

    if show_labels:
        fig.update_traces(
            textposition='outside',
            textfont=dict(size=16, color="black")
        )

    plot_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{plot_title}</title>
        <style>
            body {{ font-family: -apple-system, sans-serif; background-color: #f8f9fa; margin: 20px; }}
            .plot-container {{ width: 100%; margin-bottom: 20px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); box-sizing: border-box; }}
            .html-title {{ font-size: 26px; font-weight: bold; color: #2c3e50; margin-bottom: 15px; text-align: center; word-wrap: break-word; }}
        </style>
    </head>
    <body>
        <div class="plot-container">
            <div class="html-title">{plot_title}</div>
            {plot_html}
        </div>
    </body>
    </html>
    """

    with open(output_file, "w") as f:
        f.write(html_template)
    print(f"Plot saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot SVLEN distributions from one or more VCF or CSV files.")
    parser.add_argument('input_files', nargs='+', help="One or more VCF or CSV files to process.")
    parser.add_argument('--output', '-o', dest='output_file', default='svlen_distributions_mqc.html')
    parser.add_argument('--title', '-t', dest='plot_title', default=None, help="Custom plot title")
    parser.add_argument('--vartype', dest='vartype', choices=['small', 'structural'], default='structural', help="Analysis type: 'small' for Indels, 'structural' for SVs.")
    parser.add_argument('--tag', dest='tag', default=None, help="Tag representing the data subset (e.g., FN, FP, TP_comp).")
    parser.add_argument('--bins', '-b', nargs='+', type=str, default=None,
                        help="""
                        List of integer numbers representing the bin_edges for the plot.
                        If no list is given the default edges are bp,kbp,Mbp,Gbp,Tbp in both
                        negative(Deletions) and positive(Insertions) directions. The counts of
                        the histogram are carried out by numpy.histogram.
                        """)
    parser.add_argument('--no-labels', action='store_true', help="Hide the count labels on top of the bars.")

    args = parser.parse_args()

    if args.plot_title:
        final_title = args.plot_title
    else:
        base_title = "Indel Length Distribution" if args.vartype == 'small' else "SV Length Distribution"
        if args.tag:
            final_title = f"{base_title} of {args.tag} Variants"
        else:
            final_title = base_title

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
        plot_svlen_distributions(sv_data, bin_edges, args.output_file, final_title, show_labels=not args.no_labels)
    else:
        print("No valid input files found to plot.")
