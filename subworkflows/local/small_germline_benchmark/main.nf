//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

include { RTGTOOLS_VCFEVAL    } from '../../../modules/nf-core/rtgtools/vcfeval'
include { HAPPY_BENCHMARK     } from '../../../subworkflows/local/happy_benchmark'

workflow SMALL_GERMLINE_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]
    falsepositive_bed  // reference channel [val(meta), bed]
    stratification_bed // reference channel [val(meta), bed files]
    stratification_tsv // reference channel [val(meta), tsv]

    main:

    summary_reports = channel.empty()
    tagged_variants = channel.empty()

    if (params.method.contains('rtgtools')){

        RTGTOOLS_VCFEVAL(
            input_ch,
            sdf
        )

        summary_reports = summary_reports.mix(RTGTOOLS_VCFEVAL.out.summary
                                         .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
                                         .groupTuple())
        tagged_variants = tagged_variants.mix(RTGTOOLS_VCFEVAL.out.fn_vcf.join(RTGTOOLS_VCFEVAL.out.fn_tbi),
                                            RTGTOOLS_VCFEVAL.out.fp_vcf.join(RTGTOOLS_VCFEVAL.out.fp_tbi),
                                            RTGTOOLS_VCFEVAL.out.baseline_vcf.join(RTGTOOLS_VCFEVAL.out.baseline_tbi),
                                            RTGTOOLS_VCFEVAL.out.tp_vcf.join(RTGTOOLS_VCFEVAL.out.tp_tbi))
                                        .map { _meta, file , index ->
                                            def mapping = [
                                                'fn': 'FN',
                                                'fp': 'FP',
                                                'tp-baseline': 'TP_base',
                                                'tp': 'TP_comp'
                                            ]
                                            def tag = file.getName().tokenize('.').find { token -> token in ['fn', 'fp', 'tp-baseline', 'tp'] }
                                            def transformedTag = mapping[tag] ?: tag
                                            tuple([vartype: params.variant_type, id: "rtgtools", tag: transformedTag], file , index)
                                        }

    }

    if (params.method.contains('happy')){
        HAPPY_BENCHMARK(
            input_ch,
            fasta,
            fai,
            falsepositive_bed,
            stratification_bed,
            stratification_tsv
        )
        summary_reports = summary_reports.mix(HAPPY_BENCHMARK.out.summary_reports)
        tagged_variants = tagged_variants.mix(HAPPY_BENCHMARK.out.tagged_variants)
    }
    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]
}
