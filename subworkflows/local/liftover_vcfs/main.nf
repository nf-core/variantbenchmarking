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
    ch_targets_bed  // channel: [bed]
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
        rename_chr.map{_meta, vcf -> vcf}
    )
    vcf_ch = BCFTOOLS_ANNOTATE.out.vcf

    // liftover bed files if given
    ch_targets_bed.map{bed -> tuple([id: "targets"], bed)}
        .mix(ch_bed.map{bed -> tuple([id: "regions"], bed)})
        .set{bed_ch}

    UCSC_LIFTOVER(
        bed_ch,
        chain.map{_meta, bed -> bed}
    )

    // sort bed file
    GNU_SORT(
        UCSC_LIFTOVER.out.lifted
    )

    // merge the intersected regions
    BEDTOOLS_MERGE(
        GNU_SORT.out.sorted
    )

    BEDTOOLS_MERGE.out.bed
        .filter{ meta, _bed -> meta.id == "targets" }
        .set{targets_ch}

    BEDTOOLS_MERGE.out.bed
        .filter{ meta, _bed -> meta.id == "regions" }
        .set{bed_ch}


    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    bed_ch      // channel: [val(meta), bed]
    targets_ch  // channel: [val(meta), bed]
}
