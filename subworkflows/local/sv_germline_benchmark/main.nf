//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

include { TRUVARI_BENCHMARK        } from '../../../subworkflows/local/truvari_benchmark'
include { SVANALYZER_BENCHMARK     } from '../../../subworkflows/local/svanalyzer_benchmark'
include { WITTYER                  } from '../../../modules/nf-core/wittyer'
include { TABIX_BGZIP as TABIX_BGZIP_QUERY  } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIP as TABIX_BGZIP_TRUTH  } from '../../../modules/nf-core/tabix/bgzip'

workflow SV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()
    tagged_variants = Channel.empty()

    // SV benchmarking
    if (params.method.contains('truvari')){
        // use truvari benchmarking
        TRUVARI_BENCHMARK(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(TRUVARI_BENCHMARK.out.versions)
        summary_reports = summary_reports.mix(TRUVARI_BENCHMARK.out.report)
        tagged_variants = tagged_variants.mix(TRUVARI_BENCHMARK.out.tagged_variants)
    }

    if (params.method.contains('svanalyzer') && params.variant_type != "copynumber"){
        // WARN: svbenchmark cannot be run with copynumber analysis
        // apply svanalyzer to benchmark SVs
        SVANALYZER_BENCHMARK(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(SVANALYZER_BENCHMARK.out.versions)
        summary_reports = summary_reports.mix(SVANALYZER_BENCHMARK.out.report)
        tagged_variants = tagged_variants.mix(SVANALYZER_BENCHMARK.out.tagged_variants)
    }

    if (params.method.contains('wittyer')){

        // unzip vcf.gz files
        TABIX_BGZIP_QUERY(
            input_ch.map { meta, vcf, _tbi, _truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                [ meta, vcf ]
            }
        )
        versions = versions.mix(TABIX_BGZIP_QUERY.out.versions.first())

        TABIX_BGZIP_TRUTH(
            input_ch.map { meta, _vcf, _tbi, truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                [ meta, truth_vcf ]
            }
        )
        versions = versions.mix(TABIX_BGZIP_TRUTH.out.versions.first())

        input_ch.map { meta, _vcf, _tbi, _truth_vcf, _truth_tbi, bed, _targets_bed  ->
                [ meta, bed ]
            }
            .set { bed }

        //
        // MODULE: WITTYER
        //
        WITTYER(
            TABIX_BGZIP_QUERY.out.output
                .join(TABIX_BGZIP_TRUTH.out.output, failOnDuplicate:true, failOnMismatch:true)
                .join(bed, failOnDuplicate:true, failOnMismatch:true)
        )
        versions = versions.mix(WITTYER.out.versions.first())

        WITTYER.out.report
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "wittyer"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)
    }

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    summary_reports // channel: [val(meta), reports]
    versions        // channel: [val(meta), versions.yml]
}
