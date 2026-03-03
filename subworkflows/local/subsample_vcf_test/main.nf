//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort'
include { GAWK as ADD_GT_STRELKA} from '../../../modules/nf-core/gawk'

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

    // Add GT field using 
    ADD_GT_STRELKA(
        ch_branched_vcf.needs_gt,
        [],
        false
    )
    vcf_ch = ADD_GT_STRELKA.out.output.mix(ch_branched_vcf.ok)

    emit:
    vcf_ch      // channel: [val(meta), vcf]
}
