//
// SV_VARIANT_FILTERING: Filter SV using survivor and bcftools
//

params.options = [:]

include { SURVIVOR_FILTER     } from '../../modules/nf-core/survivor/filter'  addParams( options: params.options )
include { TABIX_BGZIP         } from '../../modules/nf-core/tabix/bgzip'      addParams( options: params.options )
include { TABIX_BGZIPTABIX    } from '../../modules/nf-core/tabix/bgziptabix'  addParams( options: params.options )


workflow SV_VARIANT_FILTERING {
    take:
    vcf_ch    // channel: [val(meta),vcf]

    main:

    versions=Channel.empty()

    TABIX_BGZIP(
        vcf_ch.map{it -> tuple( it[0], it[1])}
    )
    versions = versions.mix(TABIX_BGZIP.out.versions)

    //
    // MODULE: SURVIVOR_FILTER
    //
    // filters out smaller SVs than min_sv_size
    SURVIVOR_FILTER(
        TABIX_BGZIP.out.output.map{it -> tuple( it[0], it[1],[])},
        params.min_sv_size,
        params.max_sv_size,
        params.min_allele_freq,
        params.min_num_reads
    )
    versions = versions.mix(SURVIVOR_FILTER.out.versions)

    TABIX_BGZIPTABIX(
        SURVIVOR_FILTER.out.vcf
    )
    vcf_ch = TABIX_BGZIPTABIX.out.gz_tbi

    emit:
    vcf_ch
    versions
}
