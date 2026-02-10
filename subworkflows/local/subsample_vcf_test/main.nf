//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTERMISSING } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_PLUGINSETGT  } from '../../../modules/nf-core/bcftools/pluginsetgt/main'

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

    // Add the tools known not to have GT field here (only strelka for now)
    ch_branched_vcf = BCFTOOLS_VIEW_SUBSAMPLE.out.vcf
        .branch { meta, vcf ->
            needs_gt: meta.caller.toLowerCase().contains("strelka")
            ok:       true
        }

    BCFTOOLS_PLUGINSETGT(
        ch_branched_vcf.needs_gt.map{ meta, vcf -> tuple(meta, vcf, []) },
        Channel.value('q'),
        Channel.value('0p --include "1"'),
        [], // regions
        []  // targets
    )

    ch_for_filtering = BCFTOOLS_PLUGINSETGT.out.vcf.mix(ch_branched_vcf.ok)

    // filters out ./. genotypes (remaining from multisample vcf)
    BCFTOOLS_VIEW_FILTERMISSING(
        ch_for_filtering.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []

    )
    vcf_ch   = BCFTOOLS_VIEW_FILTERMISSING.out.vcf

    emit:
    vcf_ch      // channel: [val(meta), vcf]
}
