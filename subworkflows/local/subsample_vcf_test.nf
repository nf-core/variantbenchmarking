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

    // sorts multisample vcf
    BCFTOOLS_SORT(
        input_ch
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    // Subsamples sample name for multisample vcfs
    BCFTOOLS_VIEW_SUBSAMPLE(
        BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions.first())

    // filters out ./. genotypes (remainings from multisample vcf)
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
