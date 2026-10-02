//
// SPLIT_SMALL_VARIANTS_TEST: SUBWORKFLOW TO SPLIT SMALL SOMATIC VARIANTS INTO SNV AND INDEL
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SNV          } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INDEL        } from '../../../modules/nf-core/bcftools/view'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_SNV   } from '../../../modules/nf-core/htslib/bgziptabix'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_INDEL } from '../../../modules/nf-core/htslib/bgziptabix'

workflow SPLIT_SMALL_VARIANTS_TEST {
    take:
    input_ch // channel: [val(meta), vcf, index]

    main:

    out_vcf_ch = channel.empty()
    // split small into snv and indel if somatic
    BCFTOOLS_VIEW_SNV(
        input_ch,
        [],
        [],
        [],
    )

    TABIX_BGZIPTABIX_SNV(
        BCFTOOLS_VIEW_SNV.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    TABIX_BGZIPTABIX_SNV.out.output
        .join(TABIX_BGZIPTABIX_SNV.out.index, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, index -> tuple(meta + [vartype: "snv"], vcf, index) }
        .set { split_snv_vcf }
    out_vcf_ch = out_vcf_ch.mix(split_snv_vcf)

    BCFTOOLS_VIEW_INDEL(
        input_ch,
        [],
        [],
        [],
    )

    TABIX_BGZIPTABIX_INDEL(
        BCFTOOLS_VIEW_INDEL.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    TABIX_BGZIPTABIX_INDEL.out.output
        .join(TABIX_BGZIPTABIX_INDEL.out.index, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, index -> tuple(meta + [vartype: "indel"], vcf, index) }
        .set { split_indel_vcf }
    out_vcf_ch = out_vcf_ch.mix(split_indel_vcf)

    emit:
    out_vcf_ch // channel: [val(meta), vcf, index]
}
