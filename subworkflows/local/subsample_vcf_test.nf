//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE     } from '../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTERMISSING } from '../../modules/nf-core/bcftools/view'
include { BCFTOOLS_SORT  } from '../../modules/nf-core/bcftools/sort'


workflow SUBSAMPLE_VCF_TEST {
    take:
    input_ch    // channel: [val(meta), vcf]

    main:

    versions = Channel.empty()

    BCFTOOLS_SORT(
        input_ch
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    BCFTOOLS_VIEW_SUBSAMPLE(
        BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions.first())

    BCFTOOLS_VIEW_FILTERMISSING(
        BCFTOOLS_VIEW_SUBSAMPLE.out.vcf.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []

    )
    versions = versions.mix(BCFTOOLS_VIEW_FILTERMISSING.out.versions.first())
    vcf_ch   = BCFTOOLS_VIEW_FILTERMISSING.out.vcf

    emit:
    vcf_ch      // channel: [val(meta), vcf]
    versions
}
