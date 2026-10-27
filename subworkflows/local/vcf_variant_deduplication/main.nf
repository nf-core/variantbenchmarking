//
// VCF_VARIANT_DEDUPLICATION: deduplicate, sort and index vcf files
//

include { BCFTOOLS_SORT                     } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_NORM as BCFTOOLS_DEDUP   } from '../../../modules/nf-core/bcftools/norm'

workflow VCF_VARIANT_DEDUPLICATION {
    take:
    vcf_ch    // channel: [val(meta), vcf]
    fasta     // reference channel [val(meta), ref.fa]

    main:

    // Deduplicates variants at the same position test
    BCFTOOLS_DEDUP(
        vcf_ch,
        fasta
    )

    // sort vcf
    BCFTOOLS_SORT(
        BCFTOOLS_DEDUP.out.vcf
    )

    BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi)
        .set{ch_vcf}

    emit:
    ch_vcf      // channel: [ val(meta), vcf, index ]

}
