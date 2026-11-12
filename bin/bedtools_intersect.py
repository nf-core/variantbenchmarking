#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

import argparse
import os
import pandas as pd
from pybedtools import BedTool

TOOL_CONVERTERS = {
    "facets": lambda file: pd.read_csv(file, sep="\t")[["chrom", "start", "end", "tcn.em", "nhet"]].copy(),
    "controlfreec": lambda file: pd.read_csv(file, sep="\t", header=None).iloc[:, :5].rename(columns={0: "chrom", 1: "start", 2: "end", 3: "cn", 4: "effect"}),
    "cnvkit": lambda file: pd.read_csv(file, sep="\t")[["chromosome", "start", "end", "cn", "gene"]].copy().rename(columns={"chromosome": "chrom"}),
    "caveman": lambda file: pd.read_csv(file, sep=",", header=None).iloc[:, 1:8].rename(columns={1: "chrom", 2: "start", 3: "end", 4: "major_cn_n", 6: "major_cn_t"}),
    "ascat": lambda file: pd.read_csv(file, sep="\t")[["chr", "startpos", "endpos", "nMajor", "nMinor"]].copy().rename(columns={"chr": "chrom", "startpos": "start", "endpos": "end"})
}

def fix_chrom_prefix(chrom_series, reference):
    if reference == "GRCh38":
        return chrom_series.apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")
    elif reference == "GRCh37":
        return chrom_series.apply(lambda x: str(x).replace("chr", "") if str(x).startswith("chr") else x)
    return chrom_series

def write_to_csv(bed_df, output_bed, reference):
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df["chrom"] = fix_chrom_prefix(bed_df["chrom"], reference)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def convert_to_bed(tool_name, input_file):
    converter = TOOL_CONVERTERS.get(tool_name)
    if not converter:
        raise ValueError(f"Unknown tool: {tool_name}")
    return converter(input_file)

def compute_statistics(truth_file, test_file, min_overlap):
    """
    Compute intersection, missed regions, precision, recall, and F1-score
    using BedTools.
    """
    bed1 = BedTool(truth_file)  # Ground truth
    bed2 = BedTool(test_file)   # Test file

    # True Positives from the perspective of the ground truth file
    tp_from_truth = bed1.intersect(bed2, u=True, f=min_overlap, r=True)
    # True Positives from the perspective of the test file
    tp_from_test = bed2.intersect(bed1, u=True, f=min_overlap, r=True)

    # Convert to sets of string lines for comparison
    bed1_lines = set(str(feature) for feature in bed1)
    bed2_lines = set(str(feature) for feature in bed2)
    tp_lines = set(str(feature) for feature in tp_from_truth)
    tp_test_lines = set(str(feature) for feature in tp_from_test)

    # Calculate FN and FP based on set differences
    FN_lines = bed1_lines - tp_lines
    FP_lines = bed2_lines - tp_test_lines

    TP_base = len(tp_lines)
    TP_comp = len(tp_test_lines)
    FN = len(FN_lines)
    FP = len(FP_lines)

    # Precision, Recall, F1-score
    precision = TP_comp / (TP_comp + FP) if (TP_comp + FP) > 0 else 0
    recall = TP_base / (TP_base + FN) if (TP_base + FN) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    stats = {
        "Tool": args.test_tool,
        "Analysis": "Bedtools_Intersect",
        "TP_base": TP_base,
        "TP_comp": TP_comp,
        "FN": FN,
        "FP": FP,
        "Precision": precision,
        "Recall": recall,
        "F1": f1_score
    }
    return stats, tp_from_truth, tp_from_test, FN_lines, FP_lines

def save_statistics(stats, output_prefix):
    df = pd.DataFrame([stats])
    df.to_csv(f"{output_prefix}_stats.csv", index=False)

def write_bed(features, output_path):
    features = list(features)
    # Sort based on chrom, start, end if possible
    if features and hasattr(next(iter(features)), 'chrom'):
        features.sort(key=lambda x: (x.chrom, x.start, x.end))
    else:
        features.sort()

    with open(output_path, 'w') as out:
        for item in features:
            if isinstance(item, str):
                out.write(item if item.endswith('\n') else item + '\n')
            else:
                out.write(str(item))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare a ground truth BED file with a test BED file and compute intersection statistics.")
    parser.add_argument("truth_file", help="Ground truth BED file")
    parser.add_argument("test_file", help="Test file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("test_tool", choices=["bed", "facets", "controlfreec", "cnvkit", "caveman", "ascat"], default="bed", help="Specify the tool used for test file (default: bed)")
    parser.add_argument("genome", choices=["GRCh37", "GRCh38"], help="Reference genome build: GRCh37 (no 'chr' prefix) or GRCh38 (with 'chr' prefix)")
    parser.add_argument("--min_overlap", type=float, default=1e-9, help="Minimum overlap required as a fraction between the files")

    args = parser.parse_args()

    test_bed_file = args.test_file
    if args.test_tool != "bed":
        test_bed_file = f"{args.output_prefix}_converted.bed"
        df = convert_to_bed(args.test_tool, args.test_file)
        write_to_csv(df, test_bed_file, args.genome)

    stats, tp_from_truth, tp_from_test, fn_lines, fp_lines = compute_statistics(args.truth_file, test_bed_file, args.min_overlap)

    # Save outputs
    save_statistics(stats, args.output_prefix)
    write_bed(tp_from_truth, f"{args.output_prefix}_TP_base.bed")
    write_bed(tp_from_test, f"{args.output_prefix}_TP_comp.bed")
    write_bed(fn_lines, f"{args.output_prefix}_FN.bed")
    write_bed(fp_lines, f"{args.output_prefix}_FP.bed")
