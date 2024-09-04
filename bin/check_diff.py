#!/usr/bin/env python

import pysam
from argparse import ArgumentParser

def compare_vcf(file1, file2):
    vcf1 = pysam.VariantFile(file1)
    vcf2 = pysam.VariantFile(file2)

    variants1 = {rec.id for rec in vcf1.fetch()}
    variants2 = {rec.id for rec in vcf2.fetch()}

    unique_to_file1 = variants1 - variants2
    unique_to_file2 = variants2 - variants1

    print(f"Unique to {file1}: {unique_to_file1}")
    print(f"Unique to {file2}: {unique_to_file2}")

if __name__ == '__main__':
    import os
    import sys

    # Parse arguments
    parser = ArgumentParser(description='Check differences')
    parser.add_argument('file1', help='VCF file')
    parser.add_argument('file2', help='Output VCF file')

    args = parser.parse_args()
    compare_vcf(args.file1, args.file2)
