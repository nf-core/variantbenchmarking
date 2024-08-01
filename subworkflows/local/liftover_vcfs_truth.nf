//
// LIFTOVER_VCFS_TRUTH: SUBWORKFLOW TO LIFTOVER TRUTH VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary'
include { PICARD_LIFTOVERVCF              } from '../../modules/nf-core/picard/liftovervcf'
include { REFORMAT_HEADER                 } from '../../modules/local/reformat_header.nf'
include { BCFTOOLS_RENAME_CHR             } from '../../modules/local/bcftools_rename_chr.nf'
include { UCSC_LIFTOVER                   } from '../../modules/nf-core/ucsc/liftover'


workflow LIFTOVER_VCFS_TRUTH {
    take:
    truth_ch     // channel: [val(meta), vcf]
    high_conf_ch // channel: [val(meta), bed]
    fasta        // reference channel [val(meta), ref.fa]
    chain        // chain channel [val(meta), chain.gz]
    rename_chr   // reference channel [val(meta), chrlist.txt]

    main:

    versions = Channel.empty()

    //prepare dict file for liftover
    PICARD_CREATESEQUENCEDICTIONARY(
        fasta
    )
    versions = versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    // Use picard liftovervcf tool to convert vcfs
    PICARD_LIFTOVERVCF(
        truth_ch,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        fasta,
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

    // correct high confidence file if given
    UCSC_LIFTOVER(
        high_conf_ch,
        chain.map{meta, file -> file}
    )
    versions = versions.mix(UCSC_LIFTOVER.out.versions.first())
    bed_ch = UCSC_LIFTOVER.out.lifted

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    bed_ch      // channel: [val(meta), bed]
    versions
}
