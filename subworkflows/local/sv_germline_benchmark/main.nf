//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

include { TRUVARI_BENCHMARK      } from '../../../subworkflows/local/truvari_benchmark'
include { SVANALYZER_BENCHMARK   } from '../../../subworkflows/local/svanalyzer_benchmark'
include { WITTYER_BENCHMARK      } from '../../../subworkflows/local/wittyer_benchmark'

workflow SV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    versions        = channel.empty()
    summary_reports = channel.empty()
    tagged_variants = channel.empty()
    logs            = channel.empty()

    // SV benchmarking
    if (params.method.contains('truvari')){
        // use truvari benchmarking
        TRUVARI_BENCHMARK(
            input_ch,
            fasta,
            fai
        )
        summary_reports = summary_reports.mix(TRUVARI_BENCHMARK.out.report)
        tagged_variants = tagged_variants.mix(TRUVARI_BENCHMARK.out.tagged_variants)
        logs            = logs.mix(TRUVARI_BENCHMARK.out.logs)
    }

    if (params.method.contains('svanalyzer') && params.variant_type != "copynumber"){
        // WARN: svbenchmark cannot be run with copynumber analysis
        // apply svanalyzer to benchmark SVs
        SVANALYZER_BENCHMARK(
            input_ch,
            fasta,
            fai
        )
        summary_reports = summary_reports.mix(SVANALYZER_BENCHMARK.out.report)
        tagged_variants = tagged_variants.mix(SVANALYZER_BENCHMARK.out.tagged_variants)
    }

    if (params.method.contains('wittyer')){
        WITTYER_BENCHMARK(
            input_ch
        )
        summary_reports = summary_reports.mix(WITTYER_BENCHMARK.out.report)

    }

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    summary_reports // channel: [val(meta), reports]
    logs            // channel: [log.txt]
    versions        // channel: [val(meta), versions.yml]
}
