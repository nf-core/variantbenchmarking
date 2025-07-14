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
import os

def split_csv_by_tag(input_file, prefix):
    output_files = {
        'TP': f'{prefix}_TP.csv',
        'FP': f'{prefix}_FP.csv',
        'FN': f'{prefix}_FN.csv'
    }

    try:
        with open(input_file, newline='') as infile:
            reader = csv.reader(infile)
            header = next(reader)

            # Prepare output writers
            writers = {}
            files = {}
            for tag, filename in output_files.items():
                f = open(filename, 'w', newline='')
                writer = csv.writer(f)
                writer.writerow(header)
                writers[tag] = writer
                files[tag] = f

            # Write rows to correct files
            for row in reader:
                if len(row) > 3:
                    tag = row[3]
                    if tag in writers:
                        writers[tag].writerow(row)

            # Close all output files
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
    parser.add_argument("prefix", help="Path to the input CSV file")
    args = parser.parse_args()

    split_csv_by_tag(args.input_csv, args.prefix)

if __name__ == "__main__":
    main()
