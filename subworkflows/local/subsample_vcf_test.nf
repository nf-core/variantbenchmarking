//
// SUBSAMPLE_VCF_TEST: SUBWORKFLOW TO PREPARE SUBSET SAMPLES
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSAMPLE     } from '../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTERMISSING } from '../../modules/nf-core/bcftools/view'
include { TABIX_TABIX as TABIX_TABIX_1                 } from '../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_TABIX_2                 } from '../../modules/nf-core/tabix/tabix'
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

    TABIX_TABIX_1(
        BCFTOOLS_SORT.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_1.out.versions)

    BCFTOOLS_SORT.out.vcf.join(TABIX_TABIX_1.out.tbi, by:0)
                            .set{vcf_ch}

    BCFTOOLS_VIEW_SUBSAMPLE(
        vcf_ch,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_SUBSAMPLE.out.versions)

    TABIX_TABIX_2(
        BCFTOOLS_VIEW_SUBSAMPLE.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_2.out.versions)

    BCFTOOLS_VIEW_SUBSAMPLE.out.vcf.join(TABIX_TABIX_2.out.tbi, by:0)
                            .set{vcf_ch}

    BCFTOOLS_VIEW_FILTERMISSING(
        vcf_ch,
        [],
        [],
        []

    )
    versions = versions.mix(BCFTOOLS_VIEW_FILTERMISSING.out.versions)
    vcf_ch   = BCFTOOLS_VIEW_FILTERMISSING.out.vcf

    emit:
    vcf_ch
    versions
}
