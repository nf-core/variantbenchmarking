//
// LIFTOVER_VCFS_TRUTH: SUBWORKFLOW TO LIFTOVER TRUTH VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary'
include { PICARD_LIFTOVERVCF              } from '../../modules/nf-core/picard/liftovervcf'


workflow LIFTOVER_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    chain       // chain channel [val(meta), chain.gz]

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

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    versions
}
