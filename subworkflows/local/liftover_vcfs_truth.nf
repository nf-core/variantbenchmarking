//
// LIFTOVER_VCFS_TRUTH: SUBWORKFLOW TO LIFTOVER TRUTH VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary'
include { PICARD_LIFTOVERVCF              } from '../../modules/nf-core/picard/liftovervcf'
include { REFORMAT_HEADER                 } from '../../modules/local/reformat_header.nf'
include { BCFTOOLS_RENAME_CHR             } from '../../modules/local/bcftools_rename_chr.nf'


workflow LIFTOVER_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    chain       // chain channel [val(meta), chain.gz]
    rename_chr  // reference channel [val(meta), chrlist.txt]

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

    // reformat header, convert PS TYPE integer to string.
    REFORMAT_HEADER(
        vcf_ch.map{meta, vcf -> tuple(meta, vcf, [])}
    )

    // rename chr
    BCFTOOLS_RENAME_CHR(
        REFORMAT_HEADER.out.gz_tbi,
        rename_chr
    )
    vcf_ch = BCFTOOLS_RENAME_CHR.out.gz_tbi.map{meta, vcf, index -> tuple(meta, vcf)}

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    versions
}
