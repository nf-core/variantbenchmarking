//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

include { TRUVARI_BENCH          } from '../../../modules/nf-core/truvari/bench'
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
        TRUVARI_BENCH(
            input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                [ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed ]
            },
            fasta,
            fai
        )
        summary_reports = summary_reports.mix(TRUVARI_BENCH.out.summary
                            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "truvari"], file) }
                            .groupTuple()
                            .map { meta, files -> tuple(meta, files.flatten()) })
        tagged_variants = tagged_variants.mix(TRUVARI_BENCH.out.fn_vcf,
                                            TRUVARI_BENCH.out.fp_vcf,
                                            TRUVARI_BENCH.out.tp_base_vcf,
                                            TRUVARI_BENCH.out.tp_comp_vcf)
                                        .map { _meta, file ->
                                            def mapping = [
                                                'fn': 'FN',
                                                'fp': 'FP',
                                                'tp-base': 'TP_base',
                                                'tp-comp': 'TP_comp'
                                            ]
                                            def tag = file.getName().tokenize('.').find { token -> token in ['fn', 'fp', 'tp-base', 'tp-comp'] }
                                            def transformedTag = mapping[tag] ?: tag
                                            tuple([vartype: params.variant_type, benchmark_tool: "truvari", tag: transformedTag], file)
                                        }
        logs            = logs.mix(TRUVARI_BENCH.out.log)
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
