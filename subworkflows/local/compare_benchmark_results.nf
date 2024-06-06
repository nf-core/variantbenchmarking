
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge tp/fp/fn results from different tools.
//

params.options = [:]

include { BCFTOOLS_QUERY    } from '../../modules/nf-core/bcftools/query'   addParams( options: params.options )
include { SURVIVOR_MERGE    } from '../../modules/nf-core/survivor/merge'   addParams( options: params.options )
include { TABIX_BGZIPTABIX  } from '../../modules/nf-core/tabix/bgziptabix' addParams( options: params.options )
include { TABIX_BGZIP as TABIX_BGZIP_BENCH } from '../../modules/nf-core/tabix/bgzip'     addParams( options: params.options )

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    input_ch    // channel: [val(meta), vcf.gz]
    fai

    main:
    versions   = Channel.empty()

    //
    // MODULE: TABIX_BGZIP
    //
    // unzip vcfs

    TABIX_BGZIP_BENCH(
        input_ch
    )
    versions = versions.mix(TABIX_BGZIP_BENCH.out.versions)

    TABIX_BGZIP_BENCH.out.output
                .groupTuple()
                .set{vcf_ch}

    //
    // MODULE: SURVIVOR_MERGE
    //
    // Merge Benchmark SVs from different tools
    SURVIVOR_MERGE(
        vcf_ch,
        1000,
        1,
        1,
        0,
        0,
        30
    )
    versions = versions.mix(SURVIVOR_MERGE.out.versions)

    TABIX_BGZIPTABIX(
        SURVIVOR_MERGE.out.vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX.out.versions)

    // add a plot to see better supported variants
    // https://github.com/fritzsedlazeck/SURVIVOR/wiki

    BCFTOOLS_QUERY(
        TABIX_BGZIPTABIX.out.gz_tbi,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_QUERY.out.versions)

    emit:
    versions
}
