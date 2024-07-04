#!/usr/bin/env python
# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License
'''
Generates a CSV file from an input VCF file
Expected usage:
    $ python vcf_to_csv.py <vcf_file> <output_file>
Use --help for more information.
'''
from argparse import ArgumentParser

if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor

    # Parse arguments
    parser = ArgumentParser(description='Generate CSV file from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('-f', '--fasta-ref', help='FASTA reference file')
    args = parser.parse_args()

    variants = []

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor(args.vcf_file,ensure_pairs=False,fasta_ref=args.fasta_ref)
    df = extractor.to_dataframe()
    # Insert id column in the first position
    df.insert(0, 'id', '')
    df['id'] = df['variant_record_obj'].apply(lambda x: x.id)
    df.drop(['variant_record_obj'], axis=1, inplace=True)
    df.to_csv(f'{args.output_file}', index=False)
