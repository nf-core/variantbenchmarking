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
parser.add_argument('--input', required=True, help='Input VCF')
parser.add_argument('--svlength', type=int, required=False, help='SV Lenght', default=50)
args = parser.parse_args()



in_vcf = pysam.VariantFile(args.input)
out_name = os.path.basename(args.input)
if out_name.endswith('.gz'):
    out_name = out_name[:-3]
if out_name.endswith('.vcf'):
    out_name = out_name[:-4]

anno_header = in_vcf.header
logger.info('Writing Header')
anno_header.info.add("SVTYPE","1","String","Type of structural variant")
anno_header.info.add("SVLEN",".","Integer","Difference in length between the REF and ALT alleles")

anno_vcf = pysam.VariantFile('{}.annotated.vcf.gz'.format(out_name), 'w', header=anno_header)

logger.info('Adding SVTYE/SVLEN to VCF')
counter=0
for v in in_vcf:
    svlen=(len(v.alts[0])-len(v.ref))
    limit=args.svlength

    if (svlen < -1*limit):
        v.info.update({'SVTYPE':"DEL"})
        v.info.update({'SVLEN':svlen})
        anno_vcf.write(v)
        counter=counter+1
    elif (svlen > limit):
        v.info.update({'SVTYPE':"INS"})
        v.info.update({'SVLEN':svlen})
        anno_vcf.write(v)
        counter=counter+1
    else:
        anno_vcf.write(v)

in_vcf.close()
anno_vcf.close()
logger.info('Finished adding SVTYE/SVLEN to VCF: {}'.format(counter))

sp.run(['tabix', anno_vcf.filename])
logger.info('VCF indexed')

logger.info('DONE')
