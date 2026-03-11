#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

import argparse
import sys
import gzip
import bisect

def open_file(file_path, mode='r'):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, f'{mode}t')
    else:
        return open(file_path, mode)

def parse_vcf(file_path):
    variants = set()
    try:
        with open_file(file_path) as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    variant_id = (parts[0], parts[1], parts[3], parts[4])
                    variants.add(variant_id)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    return variants

def parse_bed_file(file_path):
    regions = {}
    starts = []
    try:
        with open_file(file_path) as bed_file:
            for line in bed_file:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    starts.append(start)
                    if chrom not in regions:
                        regions[chrom] = []
                    regions[chrom].append((start, end))
        
        min_start = min(starts)
        for chrom in regions:
            regions[chrom].sort(key=lambda x: x[0])
        
        if min_start == 0:
            print("BED file appears to be 0-based. Converting to 1-based coordinates to comparing with VCF downstream.")
            for chrom in regions:
                regions[chrom] = [(start + 1, end) for (start, end) in regions[chrom]]
        elif min_start == 1:
            print("BED file appears to be 1-based.")
        else:
            print(f"ERROR: Unexpected minimum start in BED file: {min_start}. BED must be 0- or 1-based.", file=sys.stderr)
            sys.exit(1)
    
    except FileNotFoundError:
        print(f"Error: BED file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return regions

def check_1_based(position, chrom, source="VCF/BED"):
    if position < 1:
        print(f"ERROR: {source} has position < 1 at {chrom}:{position}. Must be 1-based.", file=sys.stderr)
        sys.exit(1)

def is_in_region(chrom, pos, bed_regions):
    if chrom not in bed_regions:
        return False

    chrom_regions = bed_regions[chrom]
    idx = bisect.bisect_left(chrom_regions, (pos, float('inf')))

    if idx < len(chrom_regions):
        start, end = chrom_regions[idx]
        if start <= pos < end: 
            return True
    if idx > 0:
        start, end = chrom_regions[idx - 1]
        if start <= pos < end:
            return True

    return False

def update_sample_name_in_header(line, new_sample_name):
    if line.startswith('#CHROM') and len(line.strip().split('\t')) >= 10:
        parts = line.strip().split('\t')
        parts[-1] = new_sample_name
        return '\t'.join(parts) + '\n'
    return line

def subtract_vcf_files(primary_vcf, to_subtract_vcf, output_vcf, bed_file=None, sample_name=None):
    """
    Subtracts variants from the primary VCF that exist in the second VCF,
    with an optional BED file for region filtering.
    """
    print(f"Parsing variants from {to_subtract_vcf}...")
    variants_to_subtract = parse_vcf(to_subtract_vcf)
    print(f"Found {len(variants_to_subtract)} variants to subtract.")

    bed_regions = None
    if bed_file:
        print(f"Parsing regions from BED file {bed_file}...")
        bed_regions = parse_bed_file(bed_file)
        print(f"Found regions on {len(bed_regions)} chromosomes in BED file.")

    print(f"Filtering variants from {primary_vcf}...")

    try:
        with open_file(primary_vcf) as primary_file, open_file(output_vcf, 'w') as out_file:
            for line in primary_file:
                if line.startswith('#'):
                    if sample_name:
                        line = update_sample_name_in_header(line, sample_name)
                    out_file.write(line)
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]

                variant_id = (chrom, str(pos), ref, alt)

                if variant_id in variants_to_subtract:
                    continue

                if bed_regions and not is_in_region(chrom, pos, bed_regions):
                    continue

                out_file.write(line)

        print(f"Successfully created filtered VCF at {output_vcf}.")

    except FileNotFoundError:
        print(f"Error: File not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during filtering: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Subtracts variants from one VCF file if they exist in another, with an optional BED file for region filtering.")
    parser.add_argument("primary_vcf", help="The VCF file to filter.")
    parser.add_argument("to_subtract_vcf", help="The VCF file containing variants to remove.")
    parser.add_argument("output_vcf", help="The name of the output VCF file.")
    parser.add_argument("--bed-file", help="Optional BED file to restrict the output to specific regions.", default=None)
    parser.add_argument("--sample-name", help="Optional new sample name to use in the output VCF header.", default=None)

    args = parser.parse_args()

    subtract_vcf_files(args.primary_vcf, args.to_subtract_vcf, args.output_vcf, args.bed_file, args.sample_name)

if __name__ == "__main__":
    main()
