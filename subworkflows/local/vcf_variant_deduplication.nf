//
// VCF_VARIANT_DEDUPLICATION: deduplicate, sort and index vcf files
//

params.options = [:]

include { BCFTOOLS_SORT       } from '../../modules/nf-core/bcftools/sort'    addParams( options: params.options )
include { TABIX_TABIX         } from '../../modules/nf-core/tabix/tabix'       addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_DEDUP } from '../../modules/nf-core/bcftools/norm'     addParams( options: params.options )

workflow VCF_VARIANT_DEDUPLICATION {
    take:
    vcf_ch    // channel: [val(meta),vcf]
    fasta     // reference channel [val(meta), ref.fa]

    main:

    versions=Channel.empty()

    //
    // BCFTOOLS_DEDUP
    //
    // Deduplicates variants at the same position test
    BCFTOOLS_DEDUP(
        vcf_ch,
        fasta
    )
    versions = versions.mix(BCFTOOLS_DEDUP.out.versions)

    // sort vcf
    BCFTOOLS_SORT(
        BCFTOOLS_DEDUP.out.vcf
    )

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )

    BCFTOOLS_SORT.out.vcf.join(TABIX_TABIX.out.tbi, by:0)
                        .set{ch_vcf}


    emit:
    ch_vcf      // channel: [ val(meta), [ vcf ], [index] ]
    versions    // channel: [ versions.yml ]

}
