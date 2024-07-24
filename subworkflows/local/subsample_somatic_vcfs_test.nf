//
// SUBSAMPLE_SOMATIC_VCFS_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES IN SOMATIC CASES
//

params.options = [:]
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE   } from '../../modules/nf-core/bcftools/view'
include { TABIX_TABIX  } from '../../modules/nf-core/tabix/tabix'


workflow SUBSAMPLE_SOMATIC_VCFS_TEST {
    take:
    input_ch    // channel: [val(meta),vcf,index]

    main:

    versions=Channel.empty()

    BCFTOOLS_VIEW_SUBSAMPLE(
        input_ch,
        [],[],[]
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions)

    TABIX_TABIX(
        BCFTOOLS_VIEW_SUBSAMPLE.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_VIEW_SUBSAMPLE.out.vcf.join(TABIX_TABIX.out.tbi, by:0)
                            .set{vcf_ch}

    emit:
    vcf_ch
    versions
}
