#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci
'''
Generates a CSV file from a VCF.
Expected usage:
    $ python split_sompy_features.py <vcf_file> <prefix>
Use --help for more information.
'''
import csv
import argparse

def split_csv_by_tag(input_file, prefix):
    """
    Splits a CSV file into TP, FP, and FN files based on the 'tag' column.
    The output files will contain only CHROM, POS, REF, ALT, FILTER, and a new prefix.GT column.
    """
    output_files = {
        'TP': f'{prefix}.TP_comp.csv',
        'FP': f'{prefix}.FP.csv',
        'FN': f'{prefix}.FN.csv'
    }

    try:
        with open(input_file, newline='') as infile:
            reader = csv.reader(infile)
            new_header = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', f'{prefix}.GT']

            writers = {}
            files = {}
            for tag, filename in output_files.items():
                f = open(filename, 'w', newline='')
                writer = csv.writer(f)
                writer.writerow(new_header)
                writers[tag] = writer
                files[tag] = f

            for row in reader:
                if len(row) > 3:
                    tag = row[3]
                    if tag in writers:
                        chrom = row[1]
                        pos = row[2]

                        if tag == 'TP' or tag == 'FP':
                            ref = row[4]   # REF
                            alt = row[6]   # ALT
                            filt = row[8]  # FILTER
                        elif tag == 'FN':
                            ref = row[5]   # REF.truth
                            alt = row[7]   # ALT.truth
                            filt = row[9]  # FILTER.truth

                        # Handle empty filter values
                        if not filt:
                            filt = '.'

                        gt = '1/1'
                        new_row = [chrom, pos, ref, alt, filt, gt]
                        writers[tag].writerow(new_row)

            for f in files.values():
                f.close()

        print("Done. Files created:")
        for filename in output_files.values():
            print(f"  - {filename}")

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    parser = argparse.ArgumentParser(description="Split a CSV file into TP, FP, and FN files based on the 'tag' column.")
    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("prefix", help="Prefix for the output CSV files (e.g., 'sample_data')")
    args = parser.parse_args()

    split_csv_by_tag(args.input_csv, args.prefix)

if __name__ == "__main__":
    main()
