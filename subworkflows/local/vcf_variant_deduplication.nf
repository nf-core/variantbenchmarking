//
// VCF_VARIANT_DEDUPLICATION: deduplicate, sort and index vcf files
//

include { BCFTOOLS_SORT                     } from '../../modules/nf-core/bcftools/sort'
include { TABIX_TABIX                       } from '../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_NORM as BCFTOOLS_DEDUP   } from '../../modules/nf-core/bcftools/norm'

workflow VCF_VARIANT_DEDUPLICATION {
    take:
    vcf_ch    // channel: [val(meta), vcf]
    fasta     // reference channel [val(meta), ref.fa]

    main:

    versions = Channel.empty()
    versions = Channel.empty()

    // Deduplicates variants at the same position test
    BCFTOOLS_DEDUP(
        vcf_ch,
        fasta
    )
    versions = versions.mix(BCFTOOLS_DEDUP.out.versions.first())

    // sort vcf
    BCFTOOLS_SORT(
        BCFTOOLS_DEDUP.out.vcf
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set{ch_vcf}

    emit:
    ch_vcf      // channel: [ val(meta), vcf,index ]
    versions    // channel: [ versions.yml ]

}
