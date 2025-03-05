#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

import pandas as pd
from pybedtools import BedTool

def read_bed(file):
    """Read a BED file into a BedTool object."""
    return BedTool(file)

def convert_caveman_to_bed(caveman_file, output_bed):
    """Convert CAVEMAN output format to BED format."""
    df = pd.read_csv(caveman_file, sep=",", header=None)
    bed_df = df.rename(columns={1: "chrom", 2: "start", 3: "end", 4: "major_cn_n", 5: "minor_cn_n", 6: "major_cn_t", 7: "minor_cn_t"})
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def convert_cnvkit_to_bed(cnvkit_file, output_bed):
    """Convert CNVKIT output format to BED format."""
    df = pd.read_csv(cnvkit_file, sep="\t")
    bed_df = df[["chromosome", "start", "end", "cn", "gene"]].copy()
    bed_df["chrom"] = bed_df["chromosome"].astype(str)
    bed_df["gene"] = bed_df["gene"].astype(str)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def convert_facest_to_bed(facest_file, output_bed):
    """Convert FACETS output format to BED format."""
    df = pd.read_csv(facest_file, sep="\t")
    bed_df = df[["chrom", "start", "end", "tcn.em", "nhet"]].copy()
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def convert_controlfreec_to_bed(controlfreec_file, output_bed):
    """Convert CONTROLFREEC output format to BED format."""
    df = pd.read_csv(controlfreec_file, sep="\t", header=None)
    bed_df = df.iloc[:, : 5]
    bed_df = bed_df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "cn", 4: "effect"})
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    return output_bed

def compute_statistics(truth_file, test_file):
    """Compute intersection, missed regions, precision, recall, and F1-score."""
    bed1 = read_bed(truth_file)  # Ground truth
    bed2 = read_bed(test_file)   # Test file

    # Compute intersections
    intersect = bed1.intersect(bed2, u=True)
    TP = len(intersect)

    # Compute false negatives (FN) and false positives (FP)
    FN_regions = bed1.intersect(bed2, v=True)  # Regions in truth file missed by test file
    FP_regions = bed2.intersect(bed1, v=True)  # Extra regions in test file not in truth file
    FN = len(FN_regions)
    FP = len(FP_regions)

    # Precision, Recall, F1-score
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "TP": TP,
        "FN": FN,
        "FP": FP,
        "precision": precision,
        "recall": recall,
        "f1_score": f1_score
    }, intersect, FN_regions, FP_regions

def save_statistics(stats, output_prefix):
    """Save statistics to a text file."""
    with open(f"{output_prefix}_stats.txt", "w") as f:
        f.write("Intersection statistics:\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")

def save_bedfile(bedtool_obj, filename):
    """Save the BedTool object to a BED file."""
    bedtool_obj.saveas(filename)

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Compare a ground truth BED file with a test BED file and compute intersection statistics.")
    parser.add_argument("truth_file", help="Ground truth BED file")
    parser.add_argument("test_file", help="Test BED file (or FACETS file)")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("test_tool", choices=["bed", "facets", "controlfreec", "cnvkit", "caveman"], default="bed", help="Specify the tool used for test file (default: bed)")
    args = parser.parse_args()

    test_bed_file = args.test_file
    if args.test_tool == "facets":
        test_bed_file = f"{args.output_prefix}_converted.bed"
        convert_facest_to_bed(args.test_file, test_bed_file)

    if args.test_tool == "controlfreec":
        test_bed_file = f"{args.output_prefix}_converted.bed"
        convert_controlfreec_to_bed(args.test_file, test_bed_file)

    if args.test_tool == "cnvkit":
        test_bed_file = f"{args.output_prefix}_converted.bed"
        convert_cnvkit_to_bed(args.test_file, test_bed_file)

    if args.test_tool == "caveman":
        test_bed_file = f"{args.output_prefix}_converted.bed"
        convert_caveman_to_bed(args.test_file, test_bed_file)

    stats, intersect, missed, extra = compute_statistics(args.truth_file, test_bed_file)
    save_statistics(stats, args.output_prefix)
    save_bedfile(intersect, f"{args.output_prefix}_TP.bed")
    save_bedfile(missed, f"{args.output_prefix}_FN.bed")
    save_bedfile(extra, f"{args.output_prefix}_FP.bed")

    print("Statistics saved to", f"{args.output_prefix}_stats.txt")
    print("Intersected regions saved to", f"{args.output_prefix}_TP.bed")
    print("Missed regions saved to", f"{args.output_prefix}_FN.bed")
    print("Extra regions saved to", f"{args.output_prefix}_FP.bed")
    print("Converted BED file", f"{args.output_prefix}_converted.bed")

