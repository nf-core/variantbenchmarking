//
// LIFTOVER_VCFS: SUBWORKFLOW TO LIFTOVER VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_LIFTOVERVCF           } from '../../../modules/nf-core/picard/liftovervcf'
include { GAWK as REFORMAT_HEADER      } from '../../../modules/nf-core/gawk'
include { BCFTOOLS_ANNOTATE            } from '../../../modules/nf-core/bcftools/annotate'
include { UCSC_LIFTOVER                } from '../../../modules/nf-core/ucsc/liftover'
include { GNU_SORT                     } from '../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE               } from '../../../modules/nf-core/bedtools/merge'
include { TABIX_BGZIPTABIX             } from '../../../modules/nf-core/tabix/bgziptabix'


workflow LIFTOVER_VCFS {
    take:
    ch_vcf          // channel: [val(meta), vcf]
    ch_bed          // channel: [bed]
    fasta           // reference channel [val(meta), ref.fa]
    chain           // chain channel [val(meta), chain.gz]
    rename_chr      // reference channel [val(meta), chrlist.txt]
    dictionary      // reference channel [val(meta), genome.dict]

    main:

    // Use picard liftovervcf tool to convert vcfs
    PICARD_LIFTOVERVCF(
        ch_vcf,
        dictionary,
        fasta,
        chain
    )

    // reformat header, convert PS TYPE integer to string after liftover
    REFORMAT_HEADER(
        PICARD_LIFTOVERVCF.out.vcf_lifted,
        [],
        false
    )

    TABIX_BGZIPTABIX(
        REFORMAT_HEADER.out.output
    )

    // rename chr after liftover
    BCFTOOLS_ANNOTATE(
        TABIX_BGZIPTABIX.out.gz_index.map{meta, vcf, tbi -> tuple(meta, vcf, tbi, [], [])},
        [],
        [],
        rename_chr.map{_meta, file -> file}
    )
    vcf_ch = BCFTOOLS_ANNOTATE.out.vcf

    // liftover high confidence bed file if given
    UCSC_LIFTOVER(
        ch_bed.map{file -> tuple([id: params.truth_id], file)},
        chain.map{_meta, file -> file}
    )

    // sort bed file
    GNU_SORT(
        UCSC_LIFTOVER.out.lifted
    )

    // merge the intersected regions
    BEDTOOLS_MERGE(
        GNU_SORT.out.sorted
    )
    bed_ch = BEDTOOLS_MERGE.out.bed

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    bed_ch      // channel: [val(meta), bed]
}
