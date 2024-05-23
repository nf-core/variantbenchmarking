import argparse
import pysam
import logging
import os
import subprocess as sp

logging.basicConfig(
    format='%(levelname)-7s | %(asctime)s | %(message)s',
    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', required=True, help='Input VCF')
args = parser.parse_args()

in_vcf = pysam.VariantFile(args.vcf)
out_name = os.path.basename(args.vcf)
if out_name.endswith('.gz'):
    out_name = out_name[:-3]
if out_name.endswith('.vcf'):
    out_name = out_name[:-4]

out_header = in_vcf.header
out_vcf = pysam.VariantFile('{}.resolved.vcf.gz'.format(out_name), 'w', header=out_header)

logger.info('Only <DEL> will be sequence resolved, BND, <INS> and <DUP:TANDEM> are filtered out!')
for v in in_vcf:
    if v.info['SVTYPE'] == 'BND':
        continue
    else:
        if v.alts[0] == '<INS>':
            continue

        elif v.alts[0] == '<DUP:TANDEM>':
            continue

        elif v.alts[0] == '<DEL>':
            v.alleles = (fasta_f.fetch(v.chrom, v.pos - 1, v.stop), v.ref)
            out_vcf.write(v)
        else:
            out_vcf.write(v)
fasta_f.close()
in_vcf.close()
out_vcf.close()

sp.run(['tabix', out_vcf.filename])
logger.info('VCF indexed')

logger.info('DONE')
