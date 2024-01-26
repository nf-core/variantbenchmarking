//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BGZIP_TABIX      } from '../../modules/local/bgzip_tabix.nf'       addParams( options: params.options )
include { BCFTOOLS_NORM    } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { TABIX_TABIX      } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )

workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    ref         // reference channel [ref.fa, ref.fa.fai]
    main_chroms // channel" path(chrom sizes)

    main:

    versions=Channel.empty()

    //
    // PREPARE_VCFS
    //
    truth_ch.map { it -> tuple([id: params.sample],[caller:"truth"], it[0]) }
            .set{truth}

    // BGZIP if needed and index truth
    BGZIP_TABIX(
        truth
    )
    versions = versions.mix(BGZIP_TABIX.out.versions)

    // Normalize test
    BCFTOOLS_NORM(
        BGZIP_TABIX.out.gz_tbi,
        ref,
        [[],[]]
    )
    versions = versions.mix(BCFTOOLS_NORM.out.versions)

    TABIX_TABIX(
        BCFTOOLS_NORM.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_NORM.out.vcf.join(TABIX_TABIX.out.tbi)
                        .map{it -> tuple(it[0], it[2], it[4])}
                        .set{vcf_ch}

    emit:
    vcf_ch
    versions
}
