//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE   } from '../../modules/nf-core/bcftools/view'
include { BGZIP_TABIX  } from '../../modules/local/bgzip_tabix.nf'


workflow SUBSAMPLE_VCF_TEST {
    take:
    input_ch    // channel: [val(meta), vcf]

    main:

    versions = Channel.empty()

    //
    // MODULE: BGZIP_TABIX
    //
    // zip and index input test files
    BGZIP_TABIX(
        input_ch
    )
    versions = versions.mix(BGZIP_TABIX.out.versions)

    BCFTOOLS_VIEW_SUBSAMPLE(
        BGZIP_TABIX.out.gz_tbi,
        [],[],[]
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions)
    vcf_ch   = BCFTOOLS_VIEW_SUBSAMPLE.out.vcf

    emit:
    vcf_ch
    versions
}
