#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

import pandas as pd
from pybedtools import BedTool

def fix_chrom_prefix(chrom_series, reference):
    if reference == "GRCh38":
        return chrom_series.apply(lambda x: x if x.startswith("chr") else f"chr{x}")
    elif reference == "GRCh37":
        return chrom_series.apply(lambda x: x.replace("chr", "") if x.startswith("chr") else x)
    return chrom_series

def write_to_csv(bed_df, output_bed, reference):
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df["chrom"] = fix_chrom_prefix(bed_df["chrom"].astype(str), reference)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def convert_caveman_to_bed(caveman_file):
    """Convert CAVEMAN output format to BED format."""
    df = pd.read_csv(caveman_file, sep=",", header=None)
    df = df.iloc[:, 1: 8]
    bed_df = df.rename(columns={1: "chrom", 2: "start", 3: "end", 4: "major_cn_n", 6: "major_cn_t"})
    return bed_df

def convert_ascat_to_bed(ascat_file):
    """Convert ASCAT output format to BED format."""
    df = pd.read_csv(ascat_file, sep="\t")
    bed_df = df[["chr", "startpos", "endpos", "nMajor", "nMinor"]].copy()
    bed_df = bed_df.rename(columns={"chr": "chrom", "startpos": "start", "endpos": "end"})
    return bed_df

def convert_cnvkit_to_bed(cnvkit_file):
    """Convert CNVKIT output format to BED format."""
    df = pd.read_csv(cnvkit_file, sep="\t")
    bed_df = df[["chromosome", "start", "end", "cn", "gene"]].copy()
    bed_df["chrom"] = bed_df["chromosome"].astype(str)
    return bed_df

def convert_facest_to_bed(facest_file):
    """Convert FACETS output format to BED format."""
    df = pd.read_csv(facest_file, sep="\t")
    bed_df = df[["chrom", "start", "end", "tcn.em", "nhet"]].copy()
    return bed_df

def convert_controlfreec_to_bed(controlfreec_file):
    """Convert CONTROLFREEC output format to BED format."""
    df = pd.read_csv(controlfreec_file, sep="\t", header=None)
    bed_df = df.iloc[:, : 5]
    bed_df = bed_df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "cn", 4: "effect"})
    return bed_df

def compute_statistics(truth_file, test_file):
    """Compute intersection, missed regions, precision, recall, and F1-score."""
    bed1 = BedTool(truth_file)  # Ground truth
    bed2 = BedTool(test_file)   # Test file


    # Convert to sets of string lines for comparison
    tp_from_truth = bed1.intersect(bed2, u=True, f=args.min_overlap, r=True)
    bed1_lines = set(str(feature) for feature in bed1)
    bed2_lines = set(str(feature) for feature in bed2)
    tp_lines = set(str(feature) for feature in tp_from_truth)

    # False Negatives = in truth but not in TP
    fn_lines = bed1_lines - tp_lines

    # False Positives = in test but not in TP (from test's perspective)
    # So we need to find TP regions from bed2's point of view
    tp_from_test = bed2.intersect(bed1, u=True, f=args.min_overlap, r=True)
    tp_test_lines = set(str(feature) for feature in tp_from_test)
    fp_lines = bed2_lines - tp_test_lines

    # Compute false negatives (FN) and false positives (FP)
    TP_base = len(tp_lines)
    TP_comp = len(tp_test_lines)
    FN = len(fn_lines)
    FP = len(fp_lines)


    # Precision, Recall, F1-score
    precision = TP_comp / (TP_comp + FP) if (TP_comp + FP) > 0 else 0
    recall = TP_base / (TP_base + FN) if (TP_base + FN) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "TP_base": TP_base,
        "TP_comp": TP_comp,
        "FN": FN,
        "FP": FP,
        "Precision": precision,
        "Recall": recall,
        "F1": f1_score
    }, tp_from_truth, tp_from_test, fn_lines, fp_lines

def save_statistics(stats, output_prefix):
    """Save statistics to a CSV file."""

    df = pd.DataFrame([stats])
    df.to_csv(f"{output_prefix}_stats.csv", index=False)

def write_bed(features, output_path):
    features = list(features)  # Force materialization once

    # Try sorting based on whether it's Interval or str
    if features and hasattr(features[0], 'chrom'):
        features = sorted(features, key=lambda x: (x.chrom, x.start, x.end))
    else:
        features = sorted(features)

    with open(output_path, 'w') as out:
        for item in features:
            if isinstance(item, str):
                out.write(item if item.endswith('\n') else item + '\n')
            else:
                out.write(str(item))

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Compare a ground truth BED file with a test BED file and compute intersection statistics.")
    parser.add_argument("truth_file", help="Ground truth BED file")
    parser.add_argument("test_file", help="Test BED file (or FACETS file)")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("test_tool", choices=["bed", "facets", "controlfreec", "cnvkit", "caveman", "ascat"], default="bed", help="Specify the tool used for test file (default: bed)")
    parser.add_argument("genome", choices=["GRCh37", "GRCh38"], help="Reference genome build: GRCh37 (no 'chr' prefix) or GRCh38 (with 'chr' prefix)")
    parser.add_argument("min_overlap", default=0.000000001, help="Minimum overlap required as a fraction between the files")

    args = parser.parse_args()

    test_bed_file = args.test_file

    if args.test_tool != "bed":
        if args.test_tool == "facets":
            test_bed_file = f"{args.output_prefix}_converted.bed"
            df = convert_facest_to_bed(args.test_file)

        if args.test_tool == "controlfreec":
            test_bed_file = f"{args.output_prefix}_converted.bed"
            df = convert_controlfreec_to_bed(args.test_file)

        if args.test_tool == "cnvkit":
            test_bed_file = f"{args.output_prefix}_converted.bed"
            df = convert_cnvkit_to_bed(args.test_file)

        if args.test_tool == "caveman":
            test_bed_file = f"{args.output_prefix}_converted.bed"
            df = convert_caveman_to_bed(args.test_file)

        if args.test_tool == "ascat":
            test_bed_file = f"{args.output_prefix}_converted.bed"
            df = convert_ascat_to_bed(args.test_file)

        write_to_csv(df, test_bed_file, args.genome)

    stats, tp_from_truth, tp_from_test, fn_lines, fp_lines = compute_statistics(args.truth_file, test_bed_file)
    save_statistics(stats, args.output_prefix)
    write_bed(tp_from_truth, f"{args.output_prefix}_TP_base.bed")
    write_bed(tp_from_test, f"{args.output_prefix}_TP_comp.bed")
    write_bed(fn_lines, f"{args.output_prefix}_FN.bed")
    write_bed(fp_lines, f"{args.output_prefix}_FP.bed")
