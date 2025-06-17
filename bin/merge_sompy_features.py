#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci
'''
Generates a CSV file from a VCF
Expected usage:
    $ python split_sompy_features.py <vcf_file> <prefix>
Use --help for more information.
'''

import csv
import argparse
from collections import defaultdict
import os

KEY_COLUMNS = ["CHROM", "POS", "tag"]
FIELDS_TO_EXTRACT = ["CHROM", "POS", "tag", "REF", "REF.truth", "ALT", "ALT.truth", "QUAL", "FILTER"]
FIELDS_TO_SUFFIX = ["REF", "ALT"]

def extract_sample_suffix(filename):
    """Extract sample suffix from filename (without extension)."""
    return os.path.splitext(os.path.basename(filename))[0]

def load_csv_by_key(filepath, suffix):
    """Read a CSV file, filter relevant fields, and suffix REF/ALT."""
    with open(filepath, newline='') as f:
        reader = csv.DictReader(f)
        data = {}
        for row in reader:
            key = tuple(row[k] for k in KEY_COLUMNS)
            filtered = {}

            for field in FIELDS_TO_EXTRACT:
                if field in FIELDS_TO_SUFFIX:
                    filtered[f"{field}_{suffix}"] = row.get(field, "")
                elif field in KEY_COLUMNS:
                    filtered[field] = row.get(field, "")
                else:
                    if field not in data.get(key, {}):
                        filtered[field] = row.get(field, "")

            if key not in data:
                data[key] = filtered
            else:
                data[key].update(filtered)

        return data

def merge_dicts_by_key(dicts):
    """Merge all dicts on shared key."""
    merged = defaultdict(dict)
    for d in dicts:
        for key, row in d.items():
            merged[key].update(row)
    return merged

def write_merged_csv(merged_data, output_file):
    """Write merged dictionary to CSV."""
    sorted_keys = sorted(merged_data.keys(), key=lambda x: (x[0], int(x[1])))

    # Determine full set of columns
    all_fields = set()
    for row in merged_data.values():
        all_fields.update(row.keys())

    # Reorder fields: key columns first, then others
    fieldnames = KEY_COLUMNS + sorted(all_fields - set(KEY_COLUMNS))

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted_keys:
            writer.writerow(merged_data[key])

def main():
    parser = argparse.ArgumentParser(
        description="Merge TP/FP/FN CSVs by CHROM,POS,tag, keep selected fields, and suffix REF/ALT from filename."
    )
    parser.add_argument("files", nargs='+', help="Input CSV files (e.g. *_TP.csv)")
    parser.add_argument("--output", required=True, help="Output merged CSV file")
    args = parser.parse_args()

    all_dicts = []
    for file in args.files:
        suffix = extract_sample_suffix(file)
        print(f"Processing {file} (sample: {suffix})")
        sample_dict = load_csv_by_key(file, suffix)
        all_dicts.append(sample_dict)

    merged = merge_dicts_by_key(all_dicts)
    write_merged_csv(merged, args.output)
    print(f"Merged CSV written to {args.output}")

if __name__ == "__main__":
    main()
