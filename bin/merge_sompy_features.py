#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci
'''
Generates a CSV file from a VCF
Expected usage:
    $ python merge_sompy_features.py  $csvs --output ${prefix}.${meta.tag}.csv
Use --help for more information.
'''

import csv
import argparse
from collections import defaultdict

KEY_COLUMNS = ["CHROM", "POS"]

def load_csv_by_key(filepath):
    """
    Read a CSV file, find the dynamic GT column, and rename
    it with the extracted sample name and a '_GT' suffix.
    """
    with open(filepath, newline='') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames

        gt_column = next((col for col in header if col.endswith('.GT')), None)
        if not gt_column:
            raise ValueError(f"No column ending with '.GT' found in {filepath}")

        sample_suffix = gt_column.replace('.GT', '')
        data = {}
        for row in reader:
            key = tuple(row[k] for k in KEY_COLUMNS)
            processed_row = {
                "CHROM": row.get("CHROM"),
                "POS": row.get("POS"),
                "REF": row.get("REF", ""),
                "ALT": row.get("ALT", ""),
                "FILTER": row.get("FILTER", "")
            }
            processed_row[f"{sample_suffix}_GT"] = row.get(gt_column, "")

            if key not in data:
                data[key] = processed_row
            else:
                data[key].update(processed_row)

        return data, sample_suffix

def merge_dicts_by_key(dicts, sample_names):
    merged = defaultdict(dict)

    for d in dicts:
        for key, row in d.items():
            merged[key].update(row)

    for key in merged.keys():
        for sample in sample_names:
            gt_field = f"{sample}_GT"
            if gt_field not in merged[key]:
                merged[key][gt_field] = "./."

    return merged

def get_sample_names(files):
    """Extract sample names from a list of files."""
    sample_names = []
    for file in files:
        with open(file, newline='') as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            gt_column = next((col for col in header if col.endswith('.GT')), None)
            if gt_column:
                sample_names.append(gt_column.replace('.GT', ''))
    return sample_names

def write_merged_file(merged_data, output_file, sample_names, delimiter):
    """Write merged dictionary to file according to delimiter."""
    sorted_keys = sorted(merged_data.keys(), key=lambda x: (x[0], int(x[1])))

    fixed_fields = ["CHROM", "POS", "REF", "ALT", "FILTER"]
    dynamic_fields = sorted([f"{sample}_GT" for sample in sample_names])
    fieldnames = fixed_fields + dynamic_fields

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval='./.', delimiter = delimiter)
        writer.writeheader()
        for key in sorted_keys:
            writer.writerow({field: merged_data[key].get(field, "./.") for field in fieldnames})

def main():
    parser = argparse.ArgumentParser(
        description="Merge CSVs by CHROM, POS, and handle dynamic GT columns."
    )
    parser.add_argument("files", nargs='+', help="Input CSV files (e.g. *_TP.csv)")
    parser.add_argument("--output", required=True, help="Output merged TSV file")
    args = parser.parse_args()

    all_dicts = []
    all_sample_names = get_sample_names(args.files)

    for file in args.files:
        print(f"Processing {file}")
        sample_dict, _ = load_csv_by_key(file)
        all_dicts.append(sample_dict)

    merged = merge_dicts_by_key(all_dicts, all_sample_names)
    write_merged_file(merged, args.output, all_sample_names, delimiter = "\t")
    print(f"Merged CSV written to {args.output}")

if __name__ == "__main__":
    main()
