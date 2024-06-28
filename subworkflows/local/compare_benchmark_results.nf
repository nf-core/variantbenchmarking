
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

params.options = [:]

include { SURVIVOR_MERGE    } from '../../modules/nf-core/survivor/merge'   addParams( options: params.options )
include { TABIX_BGZIPTABIX  } from '../../modules/nf-core/tabix/bgziptabix' addParams( options: params.options )
include { TABIX_TABIX       } from '../../modules/nf-core/tabix/tabix'      addParams( options: params.options )
include { BCFTOOLS_MERGE    } from '../../modules/nf-core/bcftools/merge'   addParams( options: params.options )
include { TABIX_BGZIP as TABIX_BGZIP_BENCH } from '../../modules/nf-core/tabix/bgzip'            addParams( options: params.options )
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SV    } from '../../modules/nf-core/bcftools/query'   addParams( options: params.options )
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SMALL } from '../../modules/nf-core/bcftools/query'   addParams( options: params.options )

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    small_ch    // channel: [val(meta), vcf.gz, index]
    sv_ch       // channel: [val(meta), vcf.gz]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions   = Channel.empty()

    // Small Variants
    //
    // MODULE: BCFTOOLS_MERGE
    //
    BCFTOOLS_MERGE(
        small_ch.groupTuple(),
        fasta,
        fai,
        []
    )
    versions = versions.mix(BCFTOOLS_MERGE.out.versions)

    TABIX_TABIX(
        BCFTOOLS_MERGE.out.merged_variants
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_QUERY_SMALL(
        BCFTOOLS_MERGE.out.merged_variants.join(TABIX_TABIX.out.tbi),
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_QUERY_SMALL.out.versions)

    // SV part
    //
    // MODULE: TABIX_BGZIP
    //
    // unzip vcfs

    TABIX_BGZIP_BENCH(
        sv_ch
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

    BCFTOOLS_QUERY_SV(
        TABIX_BGZIPTABIX.out.gz_tbi,
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_QUERY_SV.out.versions)

    emit:
    versions
}
