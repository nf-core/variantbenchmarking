//
// VCF_REHEADER_SAMPLENAME: reheader sample names when needed
//

include { TABIX_TABIX         } from '../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_REHEADER   } from '../../modules/nf-core/bcftools/reheader'

workflow VCF_REHEADER_SAMPLENAME {
    take:
    vcf_ch    // channel: [val(meta), vcf]
    fai       // reference channel [val(meta), ref.fai]

    main:

    versions=Channel.empty()

    // rename sample name
    BCFTOOLS_REHEADER(
        vcf_ch.map{meta, vcf -> tuple( meta, vcf, [], [] )},
        fai
        )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions)

    TABIX_TABIX(
        BCFTOOLS_REHEADER.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_REHEADER.out.vcf
                            .join(TABIX_TABIX.out.tbi)
                            .set{ch_vcf}

    emit:
    ch_vcf      // channel: [ val(meta), vcf, index ]
    versions    // channel: [ versions.yml ]

}
