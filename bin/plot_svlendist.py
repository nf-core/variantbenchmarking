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
import plotly.graph_objects as go
from plotly.subplots import make_subplots

SEABORN_CB_PALETTE = [
    "#0173b2", "#d55e00", "#029e73", "#de8f05", "#cc78bc",
    "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"
]

def get_sample_name(vcf_file):
    base_name = os.path.basename(vcf_file)
    return base_name.split('.')[0]

# SV types with a defined direction: deletions are plotted as negative lengths,
# insertions and duplications as positive lengths. Sub-types such as DEL:ME, INS:ME:ALU
# or DUP:TANDEM are matched by their first component. All other types (BND, TRA, INV,
# CNV, CPX, ...) have no insertion/deletion direction and are counted as separate categories.
DIRECTION = {"DEL": -1, "INS": 1, "DUP": 1}

def _info_dict(info_column):
    info = {}
    for field in info_column.split(';'):
        if '=' in field:
            key, value = field.split('=', 1)
            info[key] = value
        else:
            info[field] = True
    return info

def _is_breakend(alt):
    """Breakend notation: N[chr:pos[, ]chr:pos]N, or single breakends .N / N. (not a lone '.')"""
    return '[' in alt or ']' in alt or (len(alt) > 1 and (alt.startswith('.') or alt.endswith('.')))

def _sv_type(info, alt):
    """
    SV type of a record: the SVTYPE INFO field if present, otherwise the symbolic ALT
    allele (<DEL>, <DUP:TANDEM>, ...) or BND for breakend notation. Returns None for
    sequence-resolved records without SVTYPE (e.g. small variants).
    """
    sv_type = info.get('SVTYPE')
    if sv_type and sv_type is not True:
        return str(sv_type).upper()
    if alt.startswith('<') and alt.endswith('>'):
        return alt[1:-1].upper()
    if _is_breakend(alt):
        return 'BND'
    return None

def vcf_record_length(parts):
    """
    Returns (signed_length, category) for one VCF record; exactly one of the two is None.

    - Records with a directional SV type (DEL, INS, DUP and their sub-types) get a length from
      SVLEN, END-POS or REF/ALT (in this order) and the sign from the type, because callers
      differ in the SVLEN sign convention (VCF 4.4 defines SVLEN as positive).
    - Records with any other SV type (BND, TRA, INV, CNV, CPX, ...) are returned as a category.
    - Records without SV type (sequence-resolved, e.g. small variants) keep the REF/ALT length,
      whose sign already gives the direction.
    """
    ref, alt = parts[3], parts[4].split(',')[0]
    info = _info_dict(parts[7])
    sv_type = _sv_type(info, alt)

    if sv_type is None:
        length = len(alt) - len(ref)
        return (length, None) if length != 0 else (None, None)

    base_type = sv_type.split(':')[0]
    if base_type not in DIRECTION:
        return None, base_type

    length = None
    if 'SVLEN' in info and info['SVLEN'] is not True:
        try:
            length = int(str(info['SVLEN']).split(',')[0])
        except ValueError:
            length = None
    if length is None and 'END' in info and info['END'] is not True and base_type != 'INS':
        try:
            length = int(info['END']) - int(parts[1])
        except ValueError:
            length = None
    if length is None and not alt.startswith('<') and not _is_breakend(alt):
        length = len(alt) - len(ref)
    if not length:
        return None, f"{base_type} (no length)"

    return DIRECTION[base_type] * abs(length), None

def parse_vcf_data(vcf_files):
    """
    Parses multiple VCF files, extracts signed variant lengths, and determines the maximum SV length.

    Returns:
        tuple: A tuple containing a dictionary of SV data and the maximum SV length found.
    """
    svlen_data = {}
    max_svlen = 0

    for vcf_file in vcf_files:
        print(f"Processing VCF file: {vcf_file}")

        file_name = get_sample_name(vcf_file)
        svlen_data[file_name] = {'positive': [], 'negative': [], 'other': {}}

        try:
            if vcf_file.endswith('.gz'):
                file_handle = gzip.open(vcf_file, 'rt')
            else:
                file_handle = open(vcf_file, 'r')

            with file_handle as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 8:
                        continue
                    svlen_val, category = vcf_record_length(parts)

                    if svlen_val is None:
                        if category:
                            other = svlen_data[file_name]['other']
                            other[category] = other.get(category, 0) + 1
                        continue

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

    # Non-directional SV types (BND, TRA, INV) are shown as separate categories after the length bins
    other_categories = sorted({c for d in sv_data.values() for c in d.get('other', {})})
    for file_name, lengths_dict in sv_data.items():
        for category in other_categories:
            data_table["sample"].append(file_name)
            data_table["bin_left_edge"].append(np.nan)
            data_table["bin_right_edge"].append(np.nan)
            data_table["bin_label"].append(category)
            data_table["type"].append("other")
            data_table["counts"].append(lengths_dict.get('other', {}).get(category, 0))

    df_table = pd.DataFrame(data_table)
    return df_table

def plot_svlen_distributions(sv_data, bins, output_file, plot_title, show_labels=True, is_structural=True):
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
            return f"<b>{val/1000:.1f}k</b>"
        return f"<b>{val}</b>"

    def size_label(row):
        """Unsigned size range, e.g. '100 bp - 1 kb'; the panel gives the direction."""
        if row["type"] == "other":
            return row["bin_label"]
        def fmt(v):
            v = abs(v)
            for unit, f in (("Mb", 1e6), ("kb", 1e3)):
                if v >= f:
                    return f"{v/f:g} {unit}"
            return f"{int(v)} bp"
        lo, hi = sorted([abs(row["bin_left_edge"]), abs(row["bin_right_edge"])])
        if lo == 0:
            return f"<{fmt(hi)}"
        return f"{fmt(lo)}\u2013{fmt(hi)}"

    df_upt["formatted_counts"] = df_upt["counts"].apply(format_label) if show_labels else ""
    df_upt["size_label"] = df_upt.apply(size_label, axis=1)
    # sort key: increasing absolute size within each panel, so deletions and insertions read the same way
    df_upt["order"] = df_upt.apply(
        lambda r: min(abs(r["bin_left_edge"]), abs(r["bin_right_edge"])) if r["type"] != "other" else 0, axis=1)

    positive_name = "Insertions/Duplications" if is_structural else "Insertions"
    panels = [(t, name) for t, name in (("deletion", "Deletions"), ("insertion", positive_name),
                                         ("other", "Other SV types")) if t in set(df_upt["type"])]
    panel_cats = []
    for t, _ in panels:
        sub = df_upt[df_upt["type"] == t].sort_values(["order", "size_label"])
        panel_cats.append(list(dict.fromkeys(sub["size_label"])))

    unique_tools = sorted(df_upt["sample"].unique())
    global_color_map = {tool: SEABORN_CB_PALETTE[i % len(SEABORN_CB_PALETTE)] for i, tool in enumerate(unique_tools)}

    fig = make_subplots(rows=1, cols=len(panels), shared_yaxes=True, horizontal_spacing=0.03,
                        column_widths=[max(len(c), 1.5) for c in panel_cats],
                        subplot_titles=[f"<b>{name}</b>" for _, name in panels])
    for col, ((t, _), cats) in enumerate(zip(panels, panel_cats), start=1):
        sub = df_upt[df_upt["type"] == t]
        for tool in unique_tools:
            d = sub[sub["sample"] == tool].set_index("size_label").reindex(cats)
            fig.add_trace(go.Bar(
                x=cats, y=d["counts"].fillna(0), name=tool, legendgroup=tool, showlegend=(col == 1),
                marker_color=global_color_map[tool], text=d["formatted_counts"].fillna(""),
                textposition="outside", cliponaxis=False, textfont=dict(size=13, color="black"),
                hovertemplate=f"{tool}<br>%{{x}}<br>Count: %{{y}}<extra></extra>"), row=1, col=col)

    max_val = df_upt["counts"].max()
    upper_lim = np.log10(max_val) + 0.6 if max_val > 0 else 1
    fig.update_layout(
        template="plotly_white", barmode="group", bargap=0.2,
        font=dict(size=16, color="black"),
        legend=dict(title_text="<b>Tool</b>", font=dict(size=16), title_font=dict(size=18)),
        margin=dict(l=80, r=40, t=60, b=80), hovermode="closest",
        uniformtext=dict(minsize=12, mode="show"))
    fig.update_xaxes(tickfont=dict(size=14), tickangle=-45,
                     showline=True, linewidth=1.5, linecolor="black", mirror=True)
    fig.update_yaxes(type="log", range=[-0.05, upper_lim], dtick=1, tickformat="~s",
                     showline=True, linewidth=1.5, linecolor="black", mirror=True, tickfont=dict(size=14))
    fig.update_yaxes(title_text="<b>Variant count (log scale)</b>", row=1, col=1)

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
    parser.add_argument('--vartype', dest='vartype', choices=['small', 'indel', 'structural' ,'copynumber'], default='structural', help="Analysis type: 'small' or 'indel' for indels, 'structural' or 'copynumber' for SVs.")
    parser.add_argument('--tag', dest='tag', default=None, help="Tag representing the data subset (e.g., FN, FP, TP_comp, TP_base).")
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
        plot_svlen_distributions(sv_data, bin_edges, args.output_file, final_title, show_labels=not args.no_labels,
                                 is_structural=args.vartype in ("structural", "copynumber"))
    else:
        print("No valid input files found to plot.")
