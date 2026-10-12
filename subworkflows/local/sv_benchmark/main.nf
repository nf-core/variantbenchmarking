//
// SV_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

include { TRUVARI_BENCH          } from '../../../modules/nf-core/truvari/bench'
include { SVANALYZER_BENCHMARK   } from '../../../subworkflows/local/svanalyzer_benchmark'
include { WITTYER_BENCHMARK      } from '../../../subworkflows/local/wittyer_benchmark'
include { RTGTOOLS_BNDEVAL       } from '../../../modules/nf-core/rtgtools/bndeval'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_FN      } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TP_BASE } from '../../../modules/nf-core/bcftools/reheader'

workflow SV_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    summary_reports = channel.empty()
    tagged_variants = channel.empty()
    logs            = channel.empty()

    // SV benchmarking
    if (params.method.contains('truvari')){
        TRUVARI_BENCH(
            input_ch.map{ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed, _targets_bed  ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed ]
            },
            fasta,
            fai
        )
        summary_reports = summary_reports.mix(TRUVARI_BENCH.out.summary
                            .map { _meta, report -> tuple([vartype: params.variant_type] + [benchmark_tool: "truvari"], report) }
                            .groupTuple()
                            .map { meta, files -> tuple(meta, files.flatten()) })
        tagged_variants = tagged_variants.mix(TRUVARI_BENCH.out.fn_vcf,
                                            TRUVARI_BENCH.out.fp_vcf,
                                            TRUVARI_BENCH.out.tp_base_vcf,
                                            TRUVARI_BENCH.out.tp_comp_vcf)
                                        .map { _meta, report ->
                                            def mapping = [
                                                'fn': 'FN',
                                                'fp': 'FP',
                                                'tp-base': 'TP_base',
                                                'tp-comp': 'TP_comp'
                                            ]
                                            def tag = report.getName().tokenize('.').find { token -> token in ['fn', 'fp', 'tp-base', 'tp-comp'] }
                                            def transformedTag = mapping[tag] ?: tag
                                            tuple([vartype: params.variant_type, id: "truvari", tag: transformedTag], report)
                                        }
        logs = logs.mix(TRUVARI_BENCH.out.log)
    }

    if (params.method.contains('svanalyzer') ){
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

    if (params.method.contains("rtgtools")){
            
        // running bndeval might require svdecompose
        RTGTOOLS_BNDEVAL(
            input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                [ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed ]
            }
        )
        summary_reports  = summary_reports.mix(RTGTOOLS_BNDEVAL.out.summary
                                                    .map { _meta, summary -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], summary) }
                                                    .groupTuple())

        // Reheader sample name in base outputs
        BCFTOOLS_REHEADER_FN(
            RTGTOOLS_BNDEVAL.out.fn_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )

        BCFTOOLS_REHEADER_TP_BASE(
            RTGTOOLS_BNDEVAL.out.baseline_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )

        tagged_variants = tagged_variants.mix(
                BCFTOOLS_REHEADER_FN.out.vcf.map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "rtgtools"], vcf) },
                RTGTOOLS_BNDEVAL.out.fp_vcf.map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "rtgtools"], vcf) },
                BCFTOOLS_REHEADER_TP_BASE.out.vcf.map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "rtgtools"], vcf) },
                RTGTOOLS_BNDEVAL.out.tp_vcf.map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "rtgtools"], vcf) }
            )
    }

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    summary_reports // channel: [val(meta), reports]
    logs            // channel: [log.txt]
}
