#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci
'''
Generates a CSV file from a VCF
Expected usage:
    $ python vcf_to_csv.py <vcf_file> <output>
Use --help for more information.
'''
import sys
from argparse import ArgumentParser
import re

import csv

def parse_info_field(info):
    """Parse the INFO field of a VCF line into a dictionary."""
    info_dict = {}
    for entry in info.split(';'):
        key_value = entry.split('=')
        if len(key_value) == 2:
            info_dict[key_value[0]] = key_value[1]
        else:
            info_dict[key_value[0]] = True
    return info_dict

def extract_gt_from_sample(sample, format_field):
    """Extract GT value from the sample field."""
    format_fields = format_field.split(':')
    sample_values = sample.split(':')
    if 'GT' in format_fields:
        gt_index = format_fields.index('GT')
        return sample_values[gt_index]
    return './.'  # Default GT value if not present

def vcf_to_csv(vcf_file, csv_file):
    """Convert a VCF file to a CSV file."""
    with open(vcf_file, 'r') as vcf:
        headers = []
        sample_headers = []
        records = []
        include_supp_vec = False
        include_supp = False
        include_type_inferred = False
        include_svtype = False
        include_svlen = False

        for line in vcf:
            if line.startswith('##'):
                continue  # Skip meta-information lines
            elif line.startswith('#'):
                headers = line[1:].strip().split('\t')
                sample_headers = headers[9:]  # The sample headers start from the 10th column
            else:
                row = line.strip().split('\t')
                info_dict = parse_info_field(row[7])

                # Check for SUPP_VEC, SUPP, type_inferred, SVTYPE, and SVLEN in the INFO field
                if 'SUPP_VEC' in info_dict:
                    include_supp_vec = True
                if 'SUPP' in info_dict:
                    include_supp = True
                if 'type_inferred' in info_dict:
                    include_type_inferred = True
                if 'SVTYPE' in info_dict:
                    include_svtype = True
                if 'SVLEN' in info_dict:
                    include_svlen = True

                records.append((row, info_dict))

        # Write the header with optional fields
        headers_to_write = headers[:7]  # Only keep CHROM, POS, ID, REF, ALT, QUAL, FILTER
        if include_supp_vec:
            headers_to_write.append("SUPP_VEC")
        if include_supp:
            headers_to_write.append("SUPP")
        if include_type_inferred:
            headers_to_write.append("type_inferred")
        if include_svtype:
            headers_to_write.append("SVTYPE")
        if include_svlen:
            headers_to_write.append("SVLEN")
        headers_to_write.extend([f'{sample}_GT' for sample in sample_headers])

        with open(csv_file, 'w', newline='') as csvf:
            csv_writer = csv.writer(csvf)
            csv_writer.writerow(headers_to_write)

            for row, info_dict in records:
                row_to_write = row[:7]  # Only keep CHROM, POS, ID, REF, ALT, QUAL, FILTER
                if include_supp_vec:
                    row_to_write.append(info_dict.get('SUPP_VEC', ''))
                if include_supp:
                    row_to_write.append(info_dict.get('SUPP', ''))
                if include_type_inferred:
                    row_to_write.append(info_dict.get('type_inferred', ''))
                if include_svtype:
                    row_to_write.append(info_dict.get('SVTYPE', ''))
                if include_svlen:
                    row_to_write.append(info_dict.get('SVLEN', ''))
                format_field = row[8]
                gt_values = [extract_gt_from_sample(sample, format_field) for sample in row[9:]]
                row_to_write.extend(gt_values)

                csv_writer.writerow(row_to_write)

if __name__ == '__main__':

    # Parse arguments
    parser = ArgumentParser(description='Generates a CSV file from a VCF')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output', help='Output CSV file')
    args = parser.parse_args()

    vcf_to_csv(args.vcf_file, args.output)
