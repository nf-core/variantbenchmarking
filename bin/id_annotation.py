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

anno_header = in_vcf.header

anno_vcf = pysam.VariantFile('{}.annotated.vcf.gz'.format(out_name), 'w', header=anno_header)

logger.info('Adding CHROM_POS_TYPE to ID')
for v in in_vcf:
    svlen=(len(v.alts[0])-len(v.ref))
    if (svlen < 0):
        string=str(v.chrom) +'_'+ str(v.pos) + '_DEL'
        v.id =string
    elif (svlen > 0):
        string=str(v.chrom) +'_'+ str(v.pos) + '_INS'
        v.id =string
    else:
        string=str(v.chrom) +'_'+ str(v.pos) + '_SNP'
        v.id =string


in_vcf.close()
anno_vcf.close()
logger.info('Finished')

sp.run(['tabix', anno_vcf.filename])
logger.info('VCF indexed')

logger.info('DONE')
