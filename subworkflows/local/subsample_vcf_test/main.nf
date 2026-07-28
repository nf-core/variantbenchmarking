//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort'

workflow SUBSAMPLE_VCF_TEST {
    take:
    input_ch    // channel: [val(meta), vcf]

    main:

    // sorts multisample vcf
    BCFTOOLS_SORT(
        input_ch
    )

    // Subsample sample name for multisample vcfs
    BCFTOOLS_VIEW_SUBSAMPLE(
        BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )
    vcf_ch = BCFTOOLS_VIEW_SUBSAMPLE.out.vcf

    emit:
    vcf_ch      // channel: [val(meta), vcf]
}
