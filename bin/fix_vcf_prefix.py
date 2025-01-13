#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci
'''
Generates a CSV file from a VCF
Expected usage:
    $ python fix_vcf_prefix.py <vcf_file> <output>
Use --help for more information.
'''
import os
import subprocess
import argparse
import shutil


def determine_genome_version(vcf_file):
    """
    Determine the genome version by inspecting chromosome naming in the VCF file.
    """
    with subprocess.Popen(
        f"bcftools view -h {vcf_file} | grep -m 1 '^##contig=<ID='",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    ) as proc:
        output, _ = proc.communicate()

    if not output:
        raise ValueError("Unable to determine genome version from contigs in the VCF header.")

    if "chr" in output:
        return "GRCh38"
    else:
        return "GRCh37"


def fix_vcf_prefix(input_vcf, output_vcf, rename_file, target_version):
    """
    Check the prefix of chromosome names in the VCF and fix it if necessary using the provided rename file.
    """
    # Determine genome version from input VCF
    current_version = determine_genome_version(input_vcf)
    print(f"Detected genome version: {current_version}")

    # If the current genome version matches the target, simply copy the input VCF
    if current_version == target_version:
        print(f"Genome version matches the target ({target_version}). Copying input VCF to output.")
        shutil.copy(input_vcf, output_vcf)
        shutil.copy(input_vcf + ".tbi", output_vcf + ".tbi")  # Copy index as well
        return

    # Verify the rename file exists
    if not os.path.isfile(rename_file):
        raise FileNotFoundError(f"Rename file '{rename_file}' not found.")

    # Use bcftools to rename chromosomes without modifying the header
    subprocess.check_call(
        f"bcftools annotate --rename-chrs {rename_file} --no-version {input_vcf} -Oz -o {output_vcf}",
        shell=True,
    )
    subprocess.check_call(f"bcftools index {output_vcf}", shell=True)
    print(f"Chromosome names updated and output written to {output_vcf}")


def main():
    parser = argparse.ArgumentParser(description="Check and fix VCF chromosome naming prefix.")
    parser.add_argument("input_vcf", help="Input VCF file (compressed and indexed).")
    parser.add_argument("output_vcf", help="Output VCF file.")
    parser.add_argument(
        "--rename-chr",
        required=False,
        help="Path to a file with chromosome rename mappings (for bcftools --rename-chrs).",
    )
    parser.add_argument(
        "--target-version",
        required=True,
        choices=["GRCh37", "GRCh38"],
        help="Target genome version (GRCh37 or GRCh38).",
    )
    args = parser.parse_args()

    fix_vcf_prefix(args.input_vcf, args.output_vcf, args.rename_chr, args.target_version)


if __name__ == "__main__":
    main()
