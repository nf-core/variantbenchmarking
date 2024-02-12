import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', required=True, help='Input VCF')
parser.add_argument('--repeat', required=True, help='Simple Repeat regions')
parser.add_argument('--repmask', required=True, help='Repeat Mask regions')
parser.add_argument('--segdup', required=True, help='Segmental Duplication regions')

args = parser.parse_args

def print_AF(in_vcf,samplerep,repeatmask,segdup):
    in_file = pysam.VariantFile(in_vcf)
    out_name = os.path.basename(in_vcf)
    if out_name.endswith('.gz'):
        out_name = out_name[:-3]
    if out_name.endswith('.vcf'):
        out_name = out_name[:-4]

    with open(out_name + ".txt","w") as out_file:

        print('SV_Type', 'SV', 'Length' , 'AF', 'Simple_Repeat','Repeat_Name','Repeat_Class','Segmental_Dup', file=out_file)
        for v in in_vcf:
            temp1=samplerep[(samplerep['chrom']== v.chrom) &( (( v.start > samplerep['chromStart']-1 ) & ( v.start < samplerep['chromEnd'] +1 ) )| (( v.stop > samplerep['chromStart'] -1) & ( v.stop < samplerep['chromEnd'] +1 ))) ]

            if len(temp1) > 0:
                simple_repeat=temp1.iloc[0].period
            else:
                simple_repeat='None'
            
            temp2=repeatmask[(repeatmask['genoName']== v.chrom) &( (( v.start > repeatmask['genoStart']-1 ) & ( v.start < repeatmask['genoEnd'] +1 ) )| (( v.stop > repeatmask['genoStart'] -1) & ( v.stop < repeatmask['genoEnd'] +1 ))) ]

            if len(temp2) > 0:
                repname=temp2.iloc[0].repName
                repclass=temp2.iloc[0].repClass
            else:
                repname='None'
                repclass='None'
                
            temp3=segdup[(segdup['chrom']== v.chrom) &( (( v.start > segdup['chromStart']-1 ) & ( v.start < segdup['chromEnd'] +1 ) )| (( v.stop > segdup['chromStart'] -1) & ( v.stop < segdup['chromEnd'] +1 ))) ]
            if len(temp3) > 0:
                segdups=temp3.iloc[0].fracMatch
            else:
                segdups='None'


            if len(v.ref) -  len(v.alts[0]) > 49:
                if len(v.ref) -  len(v.alts[0]) < 1001:
                    print('DEL', v.ref, len(v.alts[0])- len(v.ref) , v.info['AF'][0],simple_repeat,repname, repclass,segdups, file=out_file)
                else:
                    print('DEL', '...', len(v.alts[0])- len(v.ref) , v.info['AF'][0],simple_repeat,repname, repclass,segdups, file=out_file)
                    
            elif len(v.alts[0]) -  len(v.ref) > 49:
                if len(v.ref) -  len(v.alts[0]) < 1001:
                    print('INS', v.alts[0],len(v.alts[0]) -  len(v.ref), v.info['AF'][0],simple_repeat,repname, repclass,segdups, file=out_file)
                else:
                    print('INS', '...',len(v.alts[0]) -  len(v.ref), v.info['AF'][0],simple_repeat,repname, repclass,segdups, file=out_file)
                    
            else:
                continue


repeat_file = pd.read_csv(args.repeat, header=0, sep='\t')
repmask_file = pd.read_csv(args.repmask, header=0, sep='\t')
segdup_file = pd.read_csv(args.segdups, header=0, sep='\t')

print_AF(args.vcf, repeat_file,repmask_file,segdup_file )
