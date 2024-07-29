//
// SUBSAMPLE_SOMATIC_VCFS_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES IN SOMATIC CASES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE   } from '../../modules/nf-core/bcftools/view'
include { TABIX_TABIX                                } from '../../modules/nf-core/tabix/tabix'


workflow SUBSAMPLE_SOMATIC_VCFS_TEST {
    take:
    input_ch    // channel: [val(meta),vcf,index]

    main:

    versions = Channel.empty()

    BCFTOOLS_VIEW_SUBSAMPLE(
        input_ch,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_VIEW_SUBSAMPLE.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_VIEW_SUBSAMPLE.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set{vcf_ch}

    emit:
    vcf_ch
    versions
}
