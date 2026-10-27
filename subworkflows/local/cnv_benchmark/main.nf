//
// CNV_BENCHMARK: SUBWORKFLOW FOR COPY NUMBER VARIATIONS
//

include { TRUVARI_BENCH         } from '../../../modules/nf-core/truvari/bench'
include { WITTYER_BENCHMARK     } from '../../../subworkflows/local/wittyer_benchmark'
include { RTGTOOLS_CNVEVAL      } from '../../../modules/nf-core/rtgtools/cnveval/main'

workflow CNV_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    summary_reports = channel.empty()
    tagged_variants = channel.empty()
    logs            = channel.empty()

    if (params.method.contains('truvari')){
        TRUVARI_BENCH(
            input_ch.map{ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed, _targets_bed  ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed ]
            },
            fasta,
            fai
        )
        summary_reports = summary_reports.mix(TRUVARI_BENCH.out.summary
                            .map { _meta, summary -> tuple([vartype: params.variant_type] + [benchmark_tool: "truvari"], summary) }
                            .groupTuple()
                            .map { meta, files -> tuple(meta, files.flatten()) })
        tagged_variants = tagged_variants.mix(TRUVARI_BENCH.out.fn_vcf,
                                            TRUVARI_BENCH.out.fp_vcf,
                                            TRUVARI_BENCH.out.tp_base_vcf,
                                            TRUVARI_BENCH.out.tp_comp_vcf)
                                        .map { _meta, vcfs ->
                                            def mapping = [
                                                'fn': 'FN',
                                                'fp': 'FP',
                                                'tp-base': 'TP_base',
                                                'tp-comp': 'TP_comp'
                                            ]
                                            def tag = vcfs.getName().tokenize('.').find { token -> token in ['fn', 'fp', 'tp-base', 'tp-comp'] }
                                            def transformedTag = mapping[tag] ?: tag
                                            tuple([vartype: params.variant_type, id: "truvari", tag: transformedTag], vcfs)
                                        }
        logs            = logs.mix(TRUVARI_BENCH.out.log)
    }

    if (params.method.contains('rtgtools')){
        RTGTOOLS_CNVEVAL(
            input_ch.map{ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed, _targets_bed  ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed ]
            }
        )
        summary_reports = summary_reports.mix(RTGTOOLS_CNVEVAL.out.summary
                            .map { _meta, summary -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], summary) }
                            .groupTuple()
                            .map { meta, files -> tuple(meta, files.flatten()) })
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
}
