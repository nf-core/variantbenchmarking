#!/usr/bin/env python

# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License
'''
Generates a normalized VCF file from a VCF
Expected usage:
    $ python normalize.py <vcf_file> <output_vcf_file>
Use --help for more information.
'''
import sys
from argparse import ArgumentParser

import pysam


def _contig_to_int(contig):
    contig = contig.lower().replace('chr', '')
    if contig.isdigit():
        return int(contig)
    else:
        return 22 + ord(contig[0])


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor
    from variant_extractor.variants import VariantType

    # Parse arguments
    parser = ArgumentParser(description='Generate normalized VCF file from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_vcf_file', help='Output VCF file')
    parser.add_argument('fasta', help='Reference fasta file coupled with fai')

    args = parser.parse_args()

    # Extract header from original file
    with open(args.vcf_file, 'r') as f:
        with pysam.VariantFile(f) as input_file:
            input_file.header.add_meta('cmdline', ' '.join(sys.argv))
            header_str = str(input_file.header)

    # Open output file, write as stream
    with open(args.output_vcf_file, 'w') as output_vcf:
        # Write header
        output_vcf.write(header_str)
        print(f'Reading {args.vcf_file}...')
        # Open input file, read with variant_extractor
        extractor = VariantExtractor(args.vcf_file,ensure_pairs=False,fasta_ref=args.fasta)
        records = list(extractor)
        # Sort record by chromosome and position
        records.sort(key=lambda x: (_contig_to_int(x.contig), x.pos))
        for variant_record in records:
            output_vcf.write(str(variant_record)+'\n')
