//
// LIFTOVER_VCFS_TRUTH: SUBWORKFLOW TO LIFTOVER TRUTH VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary'
include { PICARD_LIFTOVERVCF              } from '../../modules/nf-core/picard/liftovervcf'
include { REFORMAT_HEADER                 } from '../../modules/local/reformat_header.nf'
include { BCFTOOLS_RENAME_CHR             } from '../../modules/local/bcftools_rename_chr.nf'
include { UCSC_LIFTOVER                   } from '../../modules/nf-core/ucsc/liftover'
include { MODIFY_CHR_NOTATION             } from '../../modules/local/modify_chr_notation.nf'
include { BEDTOOLS_MERGE                  } from '../../modules/nf-core/bedtools/merge'


workflow LIFTOVER_VCFS_TRUTH {
    take:
    truth_ch        // channel: [val(meta), vcf]
    high_conf_ch    // channel: [val(meta), bed]
    liftover_genome // reference channel [val(meta), ref.fa]
    chain           // chain channel [val(meta), chain.gz]
    rename_chr      // reference channel [val(meta), chrlist.txt]

    main:

    versions = Channel.empty()

    //prepare dict file for liftover of vcf files
    PICARD_CREATESEQUENCEDICTIONARY(
        liftover_genome
    )
    versions = versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    // Use picard liftovervcf tool to convert vcfs
    PICARD_LIFTOVERVCF(
        truth_ch,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        liftover_genome,
        chain
    )
    versions = versions.mix(PICARD_LIFTOVERVCF.out.versions.first())
    vcf_ch   = PICARD_LIFTOVERVCF.out.vcf_lifted

    // reformat header, convert PS TYPE integer to string after liftover
    REFORMAT_HEADER(
        vcf_ch.map{meta, vcf -> tuple(meta, vcf, [])}
    )
    versions = versions.mix(REFORMAT_HEADER.out.versions.first())

    // rename chr after liftover
    BCFTOOLS_RENAME_CHR(
        REFORMAT_HEADER.out.gz_tbi,
        rename_chr
    )
    vcf_ch = BCFTOOLS_RENAME_CHR.out.vcf

    // liftover high confidence file if given
    UCSC_LIFTOVER(
        high_conf_ch,
        chain.map{meta, file -> file}
    )
    versions = versions.mix(UCSC_LIFTOVER.out.versions.first())

    // modify chr notation and sort file
    MODIFY_CHR_NOTATION(
        UCSC_LIFTOVER.out.lifted
    )
    versions = versions.mix(MODIFY_CHR_NOTATION.out.versions.first())

    // merge the intersected regions
    BEDTOOLS_MERGE(
        MODIFY_CHR_NOTATION.out.bed
    )
    versions = versions.mix(BEDTOOLS_MERGE.out.versions.first())

    bed_ch = BEDTOOLS_MERGE.out.bed

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    bed_ch      // channel: [val(meta), bed]
    versions
}
