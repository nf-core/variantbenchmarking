//
// SPLIT_SMALL_VARIANTS_TEST: SUBWORKFLOW TO SPLIT SMALL SOMATIC VARIANTS INTO SNV AND INDEL
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SNV     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INDEL   } from '../../../modules/nf-core/bcftools/view'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SNV   } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INDEL } from '../../../modules/nf-core/tabix/bgziptabix'

workflow SPLIT_SMALL_VARIANTS_TEST {
    take:
    input_ch    // channel: [val(meta), vcf, index]

    main:

    versions   = Channel.empty()
    out_vcf_ch = Channel.empty()

    // split small into snv and indel if somatic
    BCFTOOLS_VIEW_SNV(
        input_ch,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_SNV.out.versions.first())

    TABIX_BGZIPTABIX_SNV(
        BCFTOOLS_VIEW_SNV.out.vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX_SNV.out.versions.first())

    TABIX_BGZIPTABIX_SNV.out.gz_tbi
        .map { meta, file, index -> tuple(meta + [vartype: "snv"], file, index) }
        .set{split_snv_vcf}
    out_vcf_ch = out_vcf_ch.mix(split_snv_vcf)

    BCFTOOLS_VIEW_INDEL(
        input_ch,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_INDEL.out.versions.first())

    TABIX_BGZIPTABIX_INDEL(
        BCFTOOLS_VIEW_INDEL.out.vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX_INDEL.out.versions.first())
    TABIX_BGZIPTABIX_INDEL.out.gz_tbi
        .map { meta, file, index -> tuple(meta + [vartype: "indel"], file, index) }
        .set{split_indel_vcf}
    out_vcf_ch = out_vcf_ch.mix(split_indel_vcf)

    emit:
    out_vcf_ch     // channel: [val(meta), vcf, index]
    versions       // channel: [versions.yml]
}
